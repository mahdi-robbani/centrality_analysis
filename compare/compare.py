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
mt = "../s99t/centrality/s99t.txt"
pdb = "../wt/model0_A.pdb"

# def get_df_basic(fname):
#     df = pd.read_csv(fname, sep ="\t")
#     colnames = ["degree", "betweenness", "closeness", "eigenvector"]
#     df = df[colnames]
#     return df

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

#node = pd.read_csv(wt, sep = "\t")['node']
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
diff_nodes = plot_centrality_vs_residues(diff_df, list(diff_df.columns)[:-1], [3], "diff")
print(get_count(diff_nodes))

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

scaled_diff_df = get_diff_df(wt_df, mt_df, scale = True, abs = False)
scaled_diff_df['Residue'] = wt_df['Residue']
# Plot
scaled_diff_nodes = plot_centrality_vs_residues(scaled_diff_df, list(scaled_diff_df.columns)[:-1], [3], "scaled_diff")
print(get_count(scaled_diff_nodes))

# old code
#write_pdb_files(scaled_diff_df, pdb, "pdb/scaled_diff")


# ratio_df = wt_df / mt_df
# ratio_df['node'] = node
# ratio_df_nodes = plot_centrality_vs_residues(ratio_df, list(ratio_df.columns)[:-1], [3], "ratio")
# print(get_count(ratio_df_nodes))


wt = ['T32', 'N102', 'N108', 'H126', 'H92', 'S99', 'T119', 'L122', 'V128']
mut = ['F22', 'T32', 'L98', 'N102', 'F129', 'V132', 'M136', 'V139', 'M142', 'H92', 'T99', 'T119', 'L122', 'V128']
u_diff = ['A82', 'A102', 'A111', 'A100', 'A108', 'A84', 'A125', 'A32', 'A85', 'A86', 'A113', 'A127']
s_diff = ['A82', 'A102', 'A111', 'A99', 'A100', 'A108', 'A113', 'A127', 'A66', 'A84', 'A125', 'A32', 'A85', 'A86']

all_nodes = wt + mut + s_diff
all_nodes = [n[1:] for n in all_nodes]
print(Counter(all_nodes))
#print(Counter(s_diff))

# all_nodes = [wt, mut, u_diff, s_diff] 
# for lst in all_nodes:
#     lst = sorted([n[1:] for n in lst])
#     print(lst)

print(wt_df)
plt.plot(wt_df['CF_closeness_1'])
plt.savefig("tst.png")



