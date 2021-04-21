import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import MDAnalysis as mda
import argparse
import pickle
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from collections import Counter
from pyinteraph import path_analysis as pa

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

def load_ddg(file):
    data = pd.read_csv(file, sep = "\t")
    data['Node'] = data['Node'].apply(lambda x: 'A' + str(x))
    data.columns = ["Node", "Res", "DDG"]
    data = data[["Node", "DDG"]]
    return data

def combine_cent_sasa(cent_file, sasa_file, ddg_file):
    cent = load_centrality(cent_file)
    sasa = load_sasa(sasa_file)
    ddg = load_ddg(ddg_file)
    sasa = sasa.drop(columns = ['Name'])
    data = cent.merge(sasa, on = 'Node')
    data = data.merge(ddg, on = 'Node')
    print(data)
    return data

# plotting functions

def plot_centrality_vs_residues(data, columns, sds, fname, out_dir, share_y, max_range = None, size = (9, 7)):
    node_dict = {}
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
                node_list = []
                for j, val in enumerate(data[columns[i]]):
                    if val > cutoff:
                        ax.annotate(data['Residue'][j], (j, val))
                        node_list.append(data['Residue'][j])
                node_dict.update({columns[i] : node_list})
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
        plt.savefig(f"{out_dir}/centrality_{sd}_{fname}{w}.pdf")
        plt.clf()
    return node_dict

def heatmap(data, colnames, fname, oudir):
    cor = data[colnames].corr().round(2)
    plt.figure(figsize=(10,8))
    sns.heatmap(cor, cmap="RdBu_r", annot = True)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/correlation_{fname}{w}.pdf")
    plt.clf()


# set directory and files
w = False
if w:
    cent_file = "centrality/wt_w.txt"
    w = "_w"
else:
    cent_file = "centrality/wt.txt"
    w = ""
#cent_file = "unfiltered/wt_unfilterd.txt"
sasa_file = "sasa.txt"
ddg_file = "../../mutatex/mut_data/mean_per_pos_df.txt"
psn_file = "hc_graph_filtered.dat"
#psn_file = "hc_graph.dat"
pdb_file = "model0_A.pdb"
out_dir = "plots/"
#out_dir = "unfiltered/"


# load and fix df
data = combine_cent_sasa(cent_file, sasa_file, ddg_file)
# get colnames
basic_cols = ['Degree', 'Betweenness', 'Closeness', 'Eigenvector']
cf_cols = [c for c in data.columns if "CF_" in c]
cf_cols = sorted(cf_cols, key = lambda c : c.split('_')[2:])

#plot toggles
plot = False
CENT_B = plot
CENT_CF = plot
CORR = True
NETWORK = plot


def get_count(node_dict):
    print(node_dict)
    node_list = []
    for nodes in node_dict.values():
        node_list += nodes
    cnt = Counter(node_list)
    node_list = [n for n, c in cnt.items()]
    return node_list

# plot centrality vs residues

sds = [3]

# plot basic
if CENT_B:
    print("plotting centrality vs residues for basic cenralities")
    # plot custom plots (defined by range dict)
    wt_nodes = plot_centrality_vs_residues(data, basic_cols, sds, "basic", out_dir, share_y = False, max_range = "range_dict")
    print(get_count(wt_nodes))
    plot_centrality_vs_residues(data, basic_cols, sds, "basic_f", out_dir, share_y = False)
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
    sasa_cols = basic_cols + cf_cols + ['SASA', 'DDG']
    heatmap(data, sasa_cols, "ddg", out_dir)

# save pos
def save_pos(pos, name = f"pos{w}.bin"):
    with open(name, 'wb') as f:
        pickle.dump(pos, f)

# load pos
def load_pos(file):
    with open(file, 'rb') as f:
        data = pickle.load(f)
    return data


def get_graph_pos(psn, pdb):
    nodes, names, G = pa.build_graph(psn, pdb)
    weighted_edges = [(u, v, d["weight"]) for u, v, d in G.edges(data=True)]
    # H = nx.Graph()
    # H.add_weighted_edges_from(weighted_edges)
    H = G
    # k = 0.5, iterations = 100
    pos = nx.spring_layout(H, k = 0.55, seed = 1, iterations = 100)
    save_pos(pos)
    return H, pos

def plot_graph(G, pos, df, measure, out_dir, r_dict = True):
    vmax = WT_DICT[measure] if r_dict else None
    r = "" if r_dict else "_f"
    isolates = list(nx.isolates(G))
    connected = [node for node in G.nodes() if node not in isolates]
    isolates_lab = {node:node[1:] for node in isolates}
    connected_lab = {node:node[1:] for node in connected}
    active_site = [55,60,102,113,122,126]
    acitve_lab = {f"A{str(node)}" : str(node) for node in active_site}
    #lab = {node:COVERT_DICT[node] for node in nodes}
    #weights = df[measure]
    c_weight = df[df['Node'].isin(connected)][measure]
    i_weight = df[df['Node'].isin(isolates)][measure]
    fig, ax = plt.subplots(figsize = (14,10))
    ec = nx.draw_networkx_edges(G, pos)
    colors = ["blue", "purple", "red"]
    RdPuBu = LinearSegmentedColormap.from_list("RdPuBu", colors)
    nx.draw_networkx_nodes(G, pos, node_size = 250, nodelist = isolates, node_color = 'white', edgecolors='black', vmax = vmax)
    nc = nx.draw_networkx_nodes(G, pos, node_size = 450, nodelist = connected, node_color = c_weight, cmap = RdPuBu, edgecolors='black', vmax = vmax)
    plt.colorbar(nc)
    nx.draw_networkx_labels(G, pos, labels = isolates_lab, font_size=9, font_color = 'black')
    nx.draw_networkx_labels(G, pos, labels = connected_lab, font_size=9, font_color = 'white')
    nx.draw_networkx_labels(G, pos, labels = acitve_lab, font_size=9, font_color = 'cyan')
    plt.axis('off')
    #plt.tight_layout()
    plt.title(f"{measure} centrality")
    plt.savefig(f'{out_dir}{measure}_graph{r}{w}.pdf')
    plt.clf()

G, pos = get_graph_pos(psn_file, pdb_file)
if NETWORK:
    for measure in basic_cols:
        print(f"Plotting network: {measure}")
        plot_graph(G, pos, data, measure, out_dir, r_dict = False)
        # if measure == "Betweenness":
        #     break

