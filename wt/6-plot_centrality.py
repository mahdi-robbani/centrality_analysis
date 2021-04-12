import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import MDAnalysis as mda
import argparse
import pickle

WT_DICT = {"Degree": 0.030864,
              "Betweenness": 0.144122,
              "Closeness": 0.067009,
              "Eigenvector": 0.519175}

RANGE_DICT = {"Degree": 0.035,
              "Betweenness": 0.15,
              "Closeness": 0.10,
              "Eigenvector": 0.55,
              "CF_betweenness_1" : 0.75,
              "CF_betweenness_2" : 0.75,
              "CF_betweenness_3" : 0.75,
              "CF_closeness_1" : 0.135,
              "CF_closeness_2" : 0.135,
              "CF_closeness_3" : 0.135,
              }

RES_DICT = {
    'ALA':  'A','ARG':	'R','ASN':	'N','ASP':	'D',
    'ASX':	'B','CYS':	'C','GLU':	'E','GLN':	'Q',
    'GLX':	'Z','GLY':	'G','HIS':	'H','ILE':	'I',
    'LEU':	'L','LYS':	'K','MET':	'M','PHE':	'F',
    'PRO':	'P','SER':	'S','THR':	'T','TRP':	'W',
    'TYR':	'Y','VAL':	'V'
    }

# data preprocessing functions

def get_res_column(data):
    return [RES_DICT[data['name'][i]] + n[1:] for i, n in enumerate(data['node'])]

def get_new_col_names(data):
    return [f"CF_{'_'.join(c.split('_')[2:])}" if "current_flow" in c else c.capitalize() for c in data.columns]

def load_centrality(file):
    # load file
    data = pd.read_csv(file, sep ="\t")
    # add residue column
    data['residue'] = get_res_column(data)
    # get new names
    new_colnames = get_new_col_names(data)
    # repalce col names
    data = data.rename(columns = dict(zip(data.columns, new_colnames)))
    return data

def load_sasa(file):
    data = pd.read_csv(file, sep = "\t")
    data['SASA'] = data['Mean']
    data['Node'] = data['ID'].apply(lambda x: 'A' + str(x))
    return data

def combine_cent_sasa(cent_file, sasa_file):
    cent = load_centrality(cent_file)
    sasa = load_sasa(sasa_file)
    sasa = sasa.drop(columns = ['Name'])
    data = cent.merge(sasa, on = 'Node')
    return data

# plotting functions

def plot_centrality_vs_residues(data, columns, sds, fname, out_dir, share_y, max_range = None, size = (9, 7)):
    nrows = len(columns)//2 if len(columns) % 2 == 0 else len(columns)//2 + 1
    for sd in sds:
        fig, axs = plt.subplots(nrows, 2, figsize = size, sharex= True, sharey= share_y)
        for i, ax in enumerate(axs.flat):
            if i < len(columns):
                ax.plot(data[columns[i]], label = columns[i])
                ax.set_title(columns[i])
                if max_range == "max":
                    ax.set_ylim(top = 1)
                elif max_range == "range_dict":
                    ax.set_ylim(top = RANGE_DICT[columns[i]])
                mean_std = [data[columns[i]].mean(), data[columns[i]].std()]
                cutoff = mean_std[0] + sd*mean_std[1]
                for j, val in enumerate(data[columns[i]]):
                    if val > cutoff:
                        ax.annotate(data['Residue'][j], (j, val))
                ax.hlines(y = cutoff, xmin = 0, xmax = len(data['Node']), linestyles = "dashed")
                if i % 2 == 0:
                    ax.set(ylabel='Centrality Value')
                if i == (nrows*2)-2 or i == (nrows*2)-1:
                    ax.set(xlabel='Residues')
                #ax.set(xlabel = 'residues', ylabel='Centrality Value')
                ax.grid()
            else:
                ax.set_visible(False)
        fig.tight_layout()
        plt.savefig(f"{out_dir}/centrality_{sd}_{fname}.pdf")
        plt.clf()

def heatmap(data, colnames, fname, oudir):
    cor = data[colnames].corr().round(2)
    plt.figure(figsize=(10,8))
    sns.heatmap(cor, cmap="RdBu_r", annot = True)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/correlation_{fname}.pdf")
    plt.clf()


# set directory and files
cent_file = "centrality/wt.txt"
sasa_file = "sasa.txt"
psn_file = "hc_graph_filtered.dat"
pdb_file = "model0_A.pdb"
out_dir = "plots/"

# load and fix df
data = combine_cent_sasa(cent_file, sasa_file)
# get colnames
basic_cols = ['Degree', 'Betweenness', 'Closeness', 'Eigenvector']
cf_cols = [c for c in data.columns if "CF_" in c]
cf_cols = sorted(cf_cols, key = lambda c : c.split('_')[2:])

#plot toggles
plot = False
CENT_B = plot
CENT_CF = plot
CORR = plot
NETWORK = True


# plot centrality vs residues

sds = [3]

# plot basic
if CENT_B:
    print("plotting centrality vs residues for basic cenralities")
    # plot custom plots (defined by range dict)
    plot_centrality_vs_residues(data, basic_cols, sds, "basic", out_dir, share_y = False, max_range = "range_dict")
    # Plot absolute plots (max range [0-1])
    plot_centrality_vs_residues(data, basic_cols, sds, "basic_m", out_dir, share_y = False, max_range = "max")
# plot only flow centralities
if CENT_CF:
    print("plotting centrality vs residues for CF cenralities")
    # Plot relative plots (default range)
    plot_centrality_vs_residues(data, cf_cols, sds, "cf", out_dir, share_y = False, max_range = "range_dict")
    # Plot absolute plots (max range [0-1])
    plot_centrality_vs_residues(data, cf_cols, sds, "cf_m", out_dir, share_y = False, max_range = "max")

# plot heatmap
if CORR:
    #plot heatmap
    print("plotting heatmap for all cenralities")
    sasa_cols = basic_cols + cf_cols + ['SASA']
    heatmap(data, sasa_cols, "sasa", out_dir)


def build_graph(fname, pdb = None):
    """Build a graph from the provided matrix"""

    try:
        adj_matrix = np.loadtxt(fname)
    except:
        errstr = f"Could not load file {fname} or wrong file format."
        raise ValueError(errstr)
    # if the user provided a reference structure
    if pdb is not None:
        try:
            # generate a Universe object from the PDB file
            u = mda.Universe(pdb)
        except FileNotFoundError:
            raise FileNotFoundError(f"PDB not found: {pdb}")
        except:
            raise Exception(f"Could not parse pdb file: {pdb}")
        # generate identifiers for the nodes of the graph
        identifiers = [f"{r.segment.segid}{r.resnum}" for r in u.residues]
    # if the user did not provide a reference structure
    else:
        # generate automatic identifiers going from 1 to the
        # total number of residues considered
        identifiers = [f"_{i}" for i in range(1, adj_matrix.shape[0]+1)]
    
    # generate a graph from the data loaded
    G = nx.Graph(adj_matrix)
    # set the names of the graph nodes (in place)
    node_names = dict(zip(range(adj_matrix.shape[0]), identifiers))
    nx.relabel_nodes(G, mapping = node_names, copy = False)
    # return the idenfiers and the graph
    return identifiers, G

# save pos
def save_pos(pos):
    with open('pos.bin', 'wb') as f:
        pickle.dump(pos, f)

# load pos
def load_pos(file):
    with open(file, 'rb') as f:
        data = pickle.load(f)
    return data


def get_graph_pos(psn, pdb):
    nodes, G = build_graph(psn, pdb)
    weighted_edges = [(u, v, d["weight"]) for u, v, d in G.edges(data=True)]
    # H = nx.Graph()
    # H.add_weighted_edges_from(weighted_edges)
    H = G
    # k = 0.5, iterations = 100
    pos = nx.spring_layout(H, k = 0.5, seed = 1, iterations = 100)
    save_pos(pos)
    return H, pos

def plot_graph(G, pos, df, measure, out_dir):
    nodes = G.nodes()
    lab = {node:node[1:] for node in nodes}
    #lab = {node:COVERT_DICT[node] for node in nodes}
    #weights = df[measure]
    weights = df[df['Node'].isin(nodes)][measure]
    fig, ax = plt.subplots(figsize = (14,10))
    ec = nx.draw_networkx_edges(G, pos)
    nc = nx.draw_networkx_nodes(G, pos, node_size = 300, nodelist = nodes, node_color = weights, cmap = plt.cm.RdBu_r, edgecolors='black', vmax = WT_DICT[measure])
    plt.colorbar(nc)
    nx.draw_networkx_labels(G, pos, labels = lab, font_size=9, font_color = 'black')
    plt.axis('off')
    #plt.tight_layout()
    plt.title(f"{measure} centrality")
    plt.savefig(f'{out_dir}{measure}_graph.pdf')
    plt.clf()

G, pos = get_graph_pos(psn_file, pdb_file)
if NETWORK:
    for measure in basic_cols:
        print(f"Plotting network: {measure}")
        plot_graph(G, pos, data, measure, out_dir)
        if measure == "Betweenness":
            break


# if args.graph:
#     data_largest_comp, pos, C = get_plot_requirements(matrix = args.inp_dat, 
#                                                       pdb = args.inp_pdb, 
#                                                       cent_df = data)
#     for col in basic:
#         print(f"Plotting network for all basic components: {col}")
#         plot_graph(G = C, pos = pos, df = data_largest_comp, 
#                    measure = col, out_dir = out_dir)


# if args.size:
#     nodes, G = build_graph(args.inp_dat, args.inp_pdb)
#     comps = nx.algorithms.components.connected_components(G)
#     comps_size = sorted([len(c) for c in comps], reverse= True)
#     print("Plotting size of components")
#     plt.plot(comps_size)
#     for i, length in enumerate(comps_size):
#         if length > 1:
#             plt.annotate(length, (i, length))
#     plt.savefig(f"{out_dir}comp_size.png")
    


# def get_df_basic(fname):
#     df = pd.read_csv(fname, sep ="\t")
#     colnames = ["degree", "betweenness", "closeness", "eigenvector"]
#     df = df[colnames]
#     return df

# def get_diff_df(wt, mut, abs = True):
#     node_wt = pd.read_csv(wt, sep = "\t")['node']
#     node_mut = pd.read_csv(mut, sep = "\t")['node']
#     assert node_wt.shape == node_mut.shape
#     wt_df = get_df_basic(wt)
#     mt_df = get_df_basic(mut)
#     diff_df = (wt_df - mt_df).abs() if abs else wt_df - mt_df
#     diff_df['node'] = node_wt
#     return diff_df

# if args.d_graph:
#     if args.wt and args.mut:
#         diff_df = get_diff_df(args.wt, args.mut)
#         data_largest_comp, pos, C = get_plot_requirements(matrix = args.inp_dat, 
#                                                           pdb = args.inp_pdb, 
#                                                           cent_df = diff_df)
#         for col in basic:
#             print(f"Plotting network for all basic components: {col}")
#             plot_graph(G = C, pos = pos, df = data_largest_comp, 
#                     measure = col, out_dir = f"{out_dir}diff_")
        
#     else:
#         print("ERROR: Specify WT and Mut files. Exiting...")
#         exit(1)





# def node_convert_dict(data):
#     new_nodes = [RES_DICT[data['name'][i]] + n[1:] for i, n in enumerate(data['node'])]
#     node_convert_dict = dict(zip(data['node'], new_nodes))
#     return node_convert_dict

# # network plots

# def remove_isolates(G):
#     """ Get graph without useless nodes
#     """

#     # create duplicate graph so the original is unaffected
#     H = G.copy()
#     # isolates are nodes with zero edges
#     isolates = nx.isolates(G)
#     # remove from duplicate graph
#     H.remove_nodes_from(list(isolates))
#     return H

# def keep_largest_comonent(G):
#     H = G.copy()
#     components = nx.algorithms.components.connected_components(G)
#     remove = sorted(components, key = lambda x: len(x))[:-1]
#     for nodes in remove:
#         H.remove_nodes_from(nodes)
#     return H

# def get_reduced_df(df, H):
#     pos = nx.spring_layout(H, k = 0.5)
#     node_list = H.nodes()
#     reduced_df = df.loc[df['node'].isin(node_list)]
#     return reduced_df, pos



# def get_plot_requirements(matrix, pdb, cent_df):
#     nodes, G = build_graph(matrix, pdb)
#     H = remove_isolates(G)
#     C = keep_largest_comonent(G)
#     data_largest_comp, pos = get_reduced_df(cent_df, C)
#     return data_largest_comp, pos, C

# Helper functions for network

