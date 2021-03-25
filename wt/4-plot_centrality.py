import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import MDAnalysis as mda
import argparse

# Helper functions
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

def replace_dict(names):
    new_names = []
    for name in names:
        if "_" in name:
            name = name.split("_")
            new_name = '_'.join((name[0], name[-2], name[-1]))
            new_names.append(new_name)
        else:
            new_names.append(name)
    return new_names

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--centrality_file", dest="inp_cent", default = "centrality/flow.txt")
parser.add_argument("-d", "--dat_file", dest="inp_dat", default = "hc_graph_filtered.dat")
parser.add_argument("-p", "--pdb_file", dest="inp_pdb", default = "model0_A.pdb")
parser.add_argument("-b", "--basic", dest = "basic", default = False, action="store_true")
parser.add_argument("-f", "--flow", dest = "flow", default = False, action="store_true")
parser.add_argument("-m", "--heatmap", dest = "heatmap", default = False, action="store_true")
parser.add_argument("-g", "--graph", dest = "graph", default = False, action="store_true")
parser.add_argument("-s", "--comp_size", dest = "size", default = False, action="store_true")
parser.add_argument("-o", "--output", dest = "output", default="plots")
args = parser.parse_args()


# set directory
file_path = args.inp_cent
out_dir = f"{args.output}/"


# plotting functions

def plot_centrality_vs_residues(data, columns, sds, fname, out_dir, size = (9, 7)):
    nrows = len(columns)//2 if len(columns) % 2 == 0 else len(columns)//2 + 1
    for sd in sds:
        fig, axs = plt.subplots(nrows, 2, figsize = size, sharex= True, sharey= False)
        for i, ax in enumerate(axs.flat):
            if i < len(columns):
                ax.plot(data[columns[i]], label = columns[i])
                ax.set_title(columns[i])
                mean_std = [data[columns[i]].mean(), data[columns[i]].std()]
                cutoff = mean_std[0] + sd*mean_std[1]
                for j, val in enumerate(data[columns[i]]):
                    if val > cutoff:
                        ax.annotate(data['node'][j], (j, val))
                ax.hlines(y = cutoff, xmin = 0, xmax = len(data['node']), linestyles = "dashed")
                ax.set(xlabel = 'residues', ylabel='Centrality Value')
            else:
                ax.set_visible(False)
        # for ax in axs.flat:
        #     ax.set(ylabel='Centrality Value')
        #plt.xlabel("Residues")
        #plt.ylabel("Centrality value")
        #plt.legend()
        plt.savefig(f"{out_dir}/centrality_{sd}_{fname}.pdf")
        plt.clf()

def heatmap(data, colnames, fname, oudir):
    cor = data[colnames].corr().round(2)
    plt.figure(figsize=(8,8))
    ax = sns.heatmap(cor, square = True, annot = True, cmap="viridis", vmin=0, vmax=1)
    ax.set_yticks(np.arange(0, len(colnames)+1))
    # print(np.arange(0.5, len(cor.index), 1))
    # plt.xticks(np.arange(0, len(cor.index), 1), cor.index)
    # plt.yticks(np.arange(-0.5, len(cor.index), 1), cor.index)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/correlation_{fname}.pdf")

    plt.clf()

#load file
data = pd.read_csv(file_path, sep ="\t")
# remove node names, res name, and hubs
colnames = list(data.columns)[2:]
new_colnames = replace_dict(colnames)
data = data.rename(columns = dict(zip(colnames, new_colnames)))

#plot centrality vs residues
sds = [3]
# plot basic
basic = ['degree', 'betweenness', 'closeness', 'eigenvector']
if args.basic:
    plot_centrality_vs_residues(data, basic, sds, "basic", out_dir)
# plot only flow centralities
flow = ['current_betweenness_c0', 'current_closeness_c0', 'current_betweenness_c1', 'current_closeness_c1', 'current_betweenness_c2', 'current_closeness_c2']
if args.flow:
    plot_centrality_vs_residues(data, flow, sds, "flow", out_dir)
if args.heatmap:
    #plot heatmap
    heatmap(data, new_colnames, "flow", out_dir) #xticklabels

# network plots

def remove_isolates(G):
    """Takes in a graph and returns a new graph where all the nodes with 
    zero edges have been removed.
    """

    # create duplicate graph so the original is unaffected
    H = G.copy()
    # isolates are nodes with zero edges
    isolates = nx.isolates(G)
    # remove from duplicate graph
    H.remove_nodes_from(list(isolates))
    return H

def keep_largest_comonent(G):
    H = G.copy()
    components = nx.algorithms.components.connected_components(G)
    remove = sorted(components, key = lambda x: len(x))[:-1]
    for nodes in remove:
        H.remove_nodes_from(nodes)
    return H

def get_reduced_df(df, H):
    pos = nx.spring_layout(H, k = 0.5)
    node_list = H.nodes()
    reduced_df = df.loc[df['node'].isin(node_list)]
    return reduced_df, pos


def plot_graph(G, df, measure, out_dir):
    nodes = df['node']
    lab = {node:node[1:] for node in nodes}
    weights = df[measure]
    ec = nx.draw_networkx_edges(G, pos)
    nc = nx.draw_networkx_nodes(G, pos, nodelist = nodes, node_color = weights, cmap = plt.cm.hot, edgecolors='black')
    plt.colorbar(nc)
    nx.draw_networkx_labels(G, pos, labels = lab, font_size=9, font_color ='midnightblue')
    plt.axis('off')
    plt.title(f"{measure} centrality")
    plt.savefig(f'{out_dir}{measure}_graph.png')
    plt.clf()



nodes, G = build_graph(args.inp_dat, args.inp_pdb)
H = remove_isolates(G)
C = keep_largest_comonent(G)
data_largest_comp, pos = get_reduced_df(data, C)

if args.graph:
    for col in basic:
        plot_graph(C, data_largest_comp, col, out_dir)

comps = nx.algorithms.components.connected_components(G)
comps_size = sorted([len(c) for c in comps], reverse= True)
if args.size:
    plt.plot(comps_size)
    for i, length in enumerate(comps_size):
        if length > 1:
            plt.annotate(length, (i, length))
    plt.savefig(f"{out_dir}comp_size.png")
    
print("plot compelte")
