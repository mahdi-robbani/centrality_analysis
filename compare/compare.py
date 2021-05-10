import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import MDAnalysis as mda
from Bio import PDB
import argparse
from collections import Counter
from pyinteraph import path_analysis as pa
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

RES_DICT = {
    'ALA':  'A','ARG':	'R','ASN':	'N','ASP':	'D',
    'ASX':	'B','CYS':	'C','GLU':	'E','GLN':	'Q',
    'GLX':	'Z','GLY':	'G','HIS':	'H','ILE':	'I',
    'LEU':	'L','LYS':	'K','MET':	'M','PHE':	'F',
    'PRO':	'P','SER':	'S','THR':	'T','TRP':	'W',
    'TYR':	'Y','VAL':	'V'
    }

wt = "../wt/centrality/wt.txt"
mt = "../s99t/centrality/s99t.txt"
pdb = "../wt/model0_A.pdb"
psn = "../wt/hc_graph_filtered.dat"


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

def get_diff_df(wt_df, mt_df, scale, abs = False):
    colnames = ["Degree", "Betweenness", "Closeness", "Eigenvector"]
    wt_df = wt_df[colnames]
    mt_df = mt_df[colnames]
    if scale:
        wt_df = wt_df / wt_df.max()
        mt_df = mt_df / mt_df.max()
    diff = mt_df - wt_df
    diff = diff.abs() if abs else diff
    return diff

def get_diff_rank_df(wt_df, mt_df, abs = False):
    colnames = ["Degree", "Betweenness", "Closeness", "Eigenvector"]
    wt_df = wt_df[colnames].rank(method = 'min', ascending = False)
    mt_df = mt_df[colnames].rank(method = 'min', ascending = False)
    diff = mt_df - wt_df
    diff = diff.abs() if abs else diff
    return diff

wt_df = load_centrality(wt)
mt_df = load_centrality(mt)

diff_df = get_diff_df(wt_df, mt_df, scale = False, abs = False)
diff_df['Residue'] = wt_df['Residue']

def plot_centrality_vs_residues(data, columns, sds, fname = "", size = (9, 7)):
    node_dict = {}
    # make sure equal number of rows and cols
    nrows = len(columns)//2 if len(columns) % 2 == 0 else len(columns)//2 + 1
    # make one plot for each standard deviation
    for sd in sds:
        fig, axs = plt.subplots(nrows, 2, figsize = size, sharex= True, sharey= False)
        for i, ax in enumerate(axs.flat):
            if i < len(columns):
                ax.plot(data[columns[i]], label = columns[i])
                ax.set_title(columns[i])
                mean_std = [data[columns[i]].mean(), data[columns[i]].std()]
                cutoff_u = mean_std[0] + sd*mean_std[1]
                cutoff_l = mean_std[0] - sd*mean_std[1]
                node_list = []
                for j, val in enumerate(data[columns[i]]):
                    if val > cutoff_u or val < cutoff_l:
                        resname = data['Residue'][j]
                        if resname != 'S99':
                            ax.annotate(data['Residue'][j], (j, val))
                        else:
                            ax.annotate('T99', (j, val))
                        node_list.append(data['Residue'][j])
                node_dict.update({columns[i] : node_list})
                ax.hlines(y = cutoff_u, xmin = 0, xmax = len(data['Residue']), linestyles = "dashed")
                ax.hlines(y = cutoff_l, xmin = 0, xmax = len(data['Residue']), linestyles = "dashed")
                if i % 2 == 0:
                    ax.set(ylabel='Centrality Value')
                if i == (nrows*2)-2 or i == (nrows*2)-1:
                    ax.set(xlabel='Residues')
                ax.grid()
            else:
                ax.set_visible(False)
        fig.tight_layout()
        plt.savefig(f"{fname}.pdf")
        plt.clf()
    return node_dict


def get_count(node_dict):
    print(node_dict)
    node_list = []
    for nodes in node_dict.values():
        node_list += nodes
    cnt = Counter(node_list)
    node_list = [n for n, c in cnt.items()]
    return node_list

# Plot
#diff_nodes = plot_centrality_vs_residues(diff_df, list(diff_df.columns)[:-1], [3], "diff")
#print(get_count(diff_nodes))

# create pdb

def replace_bfac_column(pdb, vals, pdb_out):
    """Replace the column containing B-factors in a PDB with
    custom values."""

    # create tthe PDB parser
    parser = PDB.PDBParser()
    # get the protein structure
    structure = parser.get_structure("protein", pdb)
    io = PDB.PDBIO()
    chain_offset = 0
    for model in structure:
        for chain in model:
            for i, residue in enumerate(chain):
                for atom in residue:
                    # set the custom value
                    atom.set_bfactor(float(vals[i+chain_offset]))
            chain_offset += len(chain)
    # set the structure for the output
    io.set_structure(structure)
    # save the structure to a new PDB file
    io.save(pdb_out)

def write_pdb_files(df, pdb, fname):
    """Save a pdb file for every centrality measure in the input 
    centrality dictionary.
    """

    for col in df.columns:
        # Ignore residue name column
        if col != "node":
            # Create input array
            cent_array = df[col]
            # Replace column and save PDB file
            replace_bfac_column(pdb, cent_array, f"{fname}_{col}.pdb")

#write_pdb_files(diff_df, pdb, "pdb/diff")

# scaled_diff_df = get_diff_df(wt_df, mt_df, scale = True, abs = False)
# scaled_diff_df['Residue'] = wt_df['Residue']
# # Plot
# scaled_diff_nodes = plot_centrality_vs_residues(scaled_diff_df, list(scaled_diff_df.columns)[:-1], [3], "scaled_diff")
# print(get_count(scaled_diff_nodes))
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
    #save_pos(pos)
    #pos = load_pos("../ccmpsn/pos.bin") #CHANGE
    #print(pos)
    return H, pos

def plot_graph(G, pos, df, measure, out_dir, ext):
    isolates = list(nx.isolates(G))
    connected = [node for node in G.nodes() if node not in isolates]
    isolates_lab = {node:node[1:] for node in isolates}
    connected_lab = {node:node[1:] for node in connected}
    active_site = [55,60,102,113,122,126]
    acitve_lab = {f"A{str(node)}" : str(node) for node in active_site}
    #lab = {node:COVERT_DICT[node] for node in nodes}
    #weights = df[measure]
    #df[measure] = df[measure]
    c_weight = df[df['Residue'].isin(connected)][measure]
    fig, ax = plt.subplots(figsize = (14,10))
    ec = nx.draw_networkx_edges(G, pos)
    colors = ["blue", "purple", "red"]
    RdPuBu = LinearSegmentedColormap.from_list("RdPuBu", colors)
    nx.draw_networkx_nodes(G, pos, node_size = 250, nodelist = isolates, node_color = 'white', edgecolors='black')
    nc = nx.draw_networkx_nodes(G, pos, node_size = 450, nodelist = connected, node_color = c_weight, cmap = RdPuBu, edgecolors='black')
    plt.colorbar(nc)
    nx.draw_networkx_labels(G, pos, labels = isolates_lab, font_size=9, font_color = 'black')
    nx.draw_networkx_labels(G, pos, labels = connected_lab, font_size=9, font_color = 'white')
    nx.draw_networkx_labels(G, pos, labels = acitve_lab, font_size=9, font_color = 'cyan')
    plt.axis('off')
    #plt.tight_layout()
    plt.title(f"{measure} centrality")
    plt.savefig(f'{out_dir}{measure}_graph_{ext}.pdf')
    plt.clf()

G, pos = get_graph_pos(psn, pdb)

rank_diff_df = get_diff_rank_df(wt_df, mt_df)
rank_diff_df = rank_diff_df#.abs()#/rank_diff_df.max()
rank_diff_df['Residue'] = wt_df['Residue'].apply(lambda x: f"A{x[1:]}")


for measure in ["Degree", "Betweenness", "Closeness", "Eigenvector"]:
    plot_graph(G, pos, rank_diff_df, measure, "", "rank_diff")


