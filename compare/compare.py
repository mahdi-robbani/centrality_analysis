import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import MDAnalysis as mda
from Bio import PDB
import argparse
from collections import Counter

RES_DICT = {
    'ALA':  'A','ARG':	'R','ASN':	'N','ASP':	'D',
    'ASX':	'B','CYS':	'C','GLU':	'E','GLN':	'Q',
    'GLX':	'Z','GLY':	'G','HIS':	'H','ILE':	'I',
    'LEU':	'L','LYS':	'K','MET':	'M','PHE':	'F',
    'PRO':	'P','SER':	'S','THR':	'T','TRP':	'W',
    'TYR':	'Y','VAL':	'V'
    }

wt = "../wt/centrality/wt.txt"
w = True
if not w:
    mt = "../s99t/centrality/s99t.txt"
    w = ""
else:
    mt = "../wt/centrality/wt_w.txt"
    w = "_w"
pdb = "../wt/model0_A.pdb"


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

def get_diff_rank_df(wt_df, mt_df)

#node = pd.read_csv(wt, sep = "\t")['node']
wt_df = load_centrality(wt)
mt_df = load_centrality(mt)

print(wt_df)

#diff_df = get_diff_df(wt_df, mt_df, scale = False, abs = False)
#diff_df['Residue'] = wt_df['Residue']

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
        plt.savefig(f"{fname}{w}.pdf")
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



