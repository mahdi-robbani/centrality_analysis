import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import MDAnalysis as mda
import pickle
from matplotlib.colors import LinearSegmentedColormap
from collections import Counter
from pyinteraph import path_analysis as pa
import adjustText

# WT_DICT = {"Degree": 0.030864,
#               "Betweenness": 0.144122,
#               "Closeness": 0.067009,
#               "Eigenvector": 0.519175}

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

CUTOFF_DICT = {"Degree": 5,
              "Betweenness": 5,
              "Closeness": 5,
              "Eigenvector": 5,
              "CF_betweenness_1" : 3,
              "CF_betweenness_2" : 3,
              "CF_betweenness_3" : 3,
              "CF_betweenness_4" : 3,
              "CF_betweenness_5" : 3,
              "CF_closeness_1" : 3,
              "CF_closeness_2" : 3,
              "CF_closeness_3" : 3,
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
    data = data.drop(["Res"], axis = 1)
    return data

def combine_cent_sasa(cent_file, sasa_file = False, ddg_file = False):
    data = load_centrality(cent_file)
    if sasa_file:
        sasa = load_sasa(sasa_file)
        sasa = sasa.drop(columns = ['Name'])
        data = data.merge(sasa, on = 'Node')
    if ddg_file:
        ddg = load_ddg(ddg_file)
        data = data.merge(ddg, on = 'Node')
    return data

# plotting functions

def plot_cent_vs_res_multiplot(data, columns, n, out_dir, ext):
    labeled_residues = ['A6','A29','A63','A66','A66','A99']
    # add index to df
    data['index'] = list(range(len(data['Node'])))
    #get subset of dataframe with only labeled residue
    data_lab = data[data['Node'].isin(labeled_residues)]
    data_other = data[~data['Node'].isin(labeled_residues)]
    # set figure size and grid
    fig, axs = plt.subplots(3, 1, figsize = (9, 7), sharex= True, sharey= False)
    for i, ax in enumerate(axs.flat):
        ax.grid(alpha = 0.5)
        # make scatter plot
        #plot normal points
        ax.scatter(x = data_other['index'], y = data_other[columns[i]], edgecolors = 'black', alpha = 0.7)
        # plot labeled points
        ax.scatter(x = data_lab['index'], y = data_lab[columns[i]], edgecolors = 'black', alpha = 0.7, color = "purple")
        # value of top n value
        top_n_cutoff = sorted(data[columns[i]], reverse = True)[n:n+1][0] - 1e-06
        max_val = data[columns[i]].max()
        # label top n
        for j, val in enumerate(data[columns[i]]):
            x_pos = j + (0.015 * max_val)
            y_pos = val  + (0.015 * max_val)
            # pick all top values, ignore value if zero
            if val > top_n_cutoff and val > 0:
                ax.annotate(data['Residue'][j], (x_pos, y_pos), alpha = 0.8)
            # label active site residues
            if data['Node'][j] in labeled_residues:
                ax.annotate(data['Residue'][j], (x_pos, y_pos), color = "purple", alpha = 0.8)
        # set limits and labels
        if i == 1:
            ax.set(ylabel='Centrality Value')
        ax.set_ylim(bottom = 0, top = (max_val * 1.1))
    plt.xlabel("Residue Index")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/{columns[0]}_centrality_top{n}{ext}.pdf")
    plt.clf()

def plot_cent_vs_res(data, column, n, out_dir, ext):
    labeled_residues = ['A6','A29','A63','A66','A66','A99']
    # add index to df
    data['index'] = list(range(len(data[column])))
    #get subset of dataframe with only labeled residue
    data_lab = data[data['Node'].isin(labeled_residues)]
    data_other = data[~data['Node'].isin(labeled_residues)]
    # set figure size and grid
    plt.figure(figsize=(9,7))
    plt.grid(alpha = 0.5)
    # make scatter plot
    #plot normal points
    plt.scatter(x = data_other['index'], y = data_other[column], edgecolors = 'black', alpha = 0.7)
    # plot labeled points
    plt.scatter(x = data_lab['index'], y = data_lab[column], edgecolors = 'black', alpha = 0.7, color = "purple")
    # value of top n value
    top_n_cutoff = sorted(data[column], reverse = True)[n:n+1][0] - 1e-06
    texts = []
    # label top n
    for i, val in enumerate(data[column]):
        # pick all top values, ignore value if zero
        if data['Node'][i] in labeled_residues:
            texts.append(plt.text(i, val, data['Residue'][i], color = "purple"))
        elif val > top_n_cutoff and val > 0:
            texts.append(plt.text(i, val, data['Residue'][i]))
        # label active site residues

    # set limits and labels
    plt.ylim(bottom = 0)
    plt.xlabel("Residue Index")
    plt.ylabel("Centrality Value")
    plt.tight_layout()
    # adjust text
    adjustText.adjust_text(texts)
    plt.savefig(f"{out_dir}/{column}_centrality_top{n}{ext}.pdf")
    plt.clf()


def heatmap(data, colnames, fname, oudir):
    #ignore glycines
    data = data.loc[data['Name'] != "GLY"]
    new_colnames = [f"Mutant Count > {col.split('_')[1]}" if col.split('_')[0] == "Count" else col for col in colnames]
    col_dict = {col : f"Mutant Count > {col.split('_')[1]}" for col in data.columns if col.split('_')[0] == "Count"}
    data = data.rename(columns = col_dict)
    cor = data[new_colnames].corr().round(2)
    plt.figure(figsize=(10,8))
    sns.heatmap(cor, cmap="RdBu_r", annot = True)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/correlation_{fname}.pdf")
    plt.clf()


# set directory and files
cent_file = "centrality/wt.txt"
sasa_file = "sasa.txt"
ddg_dict = {"mean" : "../mutatex/mut_data/mean_per_pos_df.txt",
            "median" : "../mutatex/mut_data/median_per_pos_df.txt",
            "count_3" : "../mutatex/mut_data/Count_3_df.txt"}
ddg_file = "../mutatex/mut_data/per_pos_df.txt"
psn_file = "hc_graph_filtered.dat"
pdb_file = "model0_A.pdb"
out_dir = "plots/"

# load and fix df
# replace sasa_file/dgg_file with None if not present
data = combine_cent_sasa(cent_file, sasa_file, ddg_file)
# #ignore glycines
# data = data.loc[data['Name'] != "GLY"]
# data = data.reset_index()
# get colnames
basic_cols = ['Degree', 'Betweenness', 'Closeness', 'Eigenvector']
cf_cols = [c for c in data.columns if "CF_" in c]
#cf_cols = sorted(cf_cols, key = lambda c : c.split('_')[2:])
#cf_cols = cf_cols[:6]
cf_bet_cols = [c for c in cf_cols if "betweenness" in c][:3]
cf_close_cols = [c for c in cf_cols if "closeness" in c][:3]

#plot toggles
plot = False
CENT_B = True
CENT_CF = True
CORR = plot
NETWORK = plot

def get_count(node_dict):
    #print(node_dict)
    node_list = []
    for nodes in node_dict.values():
        node_list += nodes
    cnt = Counter(node_list)
    node_list = [n for n, c in cnt.items()]
    return node_list

# plot centrality vs residues

# plot basic
if CENT_B:
    print("plotting centrality vs residues")
    for cent in basic_cols:
        print(f"plotting {cent}")
        plot_cent_vs_res(data, cent, 5, "plots", "")
    #plot_centrality_vs_residues(data, basic_cols, sds, "basic_f", out_dir, share_y = False)
    # Plot absolute plots (max range [0-1])
    #plot_centrality_vs_residues(data, basic_cols, sds, "basic_m", out_dir, share_y = False, max_range = "max")
    # #plot sasa and ddg values
    # ddg_cols = ["SASA", "Mean DDG", "Median DDG", "Count_3"]
    # plot_centrality_vs_residues(data, ddg_cols, sds, "ddg", out_dir, share_y = False)

# plot only flow centralities
if CENT_CF:
    print("plotting centrality vs residues for CF cenralities")
    # Plot relative plots (default range)
    print(f"plotting {cf_bet_cols}")
    plot_cent_vs_res_multiplot(data, cf_bet_cols, 5, "plots", "")
    print(f"plotting {cf_close_cols}")
    plot_cent_vs_res_multiplot(data, cf_close_cols, 5, "plots", "")
    #plot_centrality_vs_residues(data, cf_cols, sds, "cf", out_dir, share_y = False)
    # # Plot absolute plots (max range [0-1])
    # plot_centrality_vs_residues(data, cf_cols, sds, "cf_m", out_dir, share_y = False, max_range = "max")

# plot heatmap
if CORR:
    #plot heatmap
    print("plotting heatmap for all cenralities")
    cent_cols = basic_cols + cf_cols
    heatmap(data, cent_cols, "noG", out_dir)
    print("plotting heatmap for all sasa cenralities")
    sasa_cols = basic_cols + ['SASA', 'Mean DDG', 'Median DDG', 'Count_3']
    heatmap(data, sasa_cols, "ddg_noG", out_dir)
    # print("plotting heatmap for all count cenralities")
    # sasa_cols = basic_cols + ['SASA', 'Mean DDG', 'Median DDG'] + [f"Count_{i}" for i in range(0, 12)]
    # heatmap(data, sasa_cols, "ddg_count_noG", out_dir)

# save pos
def save_pos(pos, name = f"pos.bin"):
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
    #pos = load_pos("../ccmpsn/pos.bin") #CHANGE
    #print(pos)
    return H, pos

def plot_graph(G, pos, df, measure, out_dir, ext, r_dict = False):
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
    plt.savefig(f'{out_dir}{measure}_graph{r}_{ext}.pdf')
    plt.clf()

G, pos = get_graph_pos(psn_file, pdb_file)
if NETWORK:
    for measure in basic_cols + ['Mean DDG', 'Median DDG', 'Count_3']:
        print(f"Plotting network: {measure}")
        plot_graph(G, pos, data, measure, out_dir, "", r_dict = False)
        # if measure == "Betweenness":
        #     break

