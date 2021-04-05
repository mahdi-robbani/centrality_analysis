import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import MDAnalysis as mda
from Bio import PDB
import argparse
from collections import Counter


wt = "../wt/centrality/flow.txt"
mt = "../s99t/centrality/s99t.txt"
pdb = "../wt/model0_A.pdb"

def get_df_basic(fname):
    df = pd.read_csv(fname, sep ="\t")
    colnames = ["degree", "betweenness", "closeness", "eigenvector"]
    df = df[colnames]
    return df

def get_diff_df(wt_df, mt_df, scale, abs = False):
    if scale:
        wt_df = wt_df / wt_df.max()
        mt_df = mt_df / mt_df.max()
    diff = mt_df - wt_df
    diff = diff.abs() if abs else diff
    return diff

node = pd.read_csv(wt, sep = "\t")['node']
wt_df = get_df_basic(wt)
mt_df = get_df_basic(mt)

diff_df = get_diff_df(wt_df, mt_df, scale = False, abs = False)
diff_df['node'] = node

def plot_centrality_vs_residues(data, columns, sds, fname = "", size = (9, 7)):
    node_dict = {}
    # make sure equal number of rows and cols
    nrows = len(columns)//2 if len(columns) % 2 == 0 else len(columns)//2 + 1
    # make one plot for each standard deviation
    for sd in sds:
        fig, axs = plt.subplots(nrows, 2, figsize = size, sharex= True, sharey= True)
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
                        ax.annotate(data['node'][j], (j, val))
                        node_list.append(data['node'][j])
                node_dict.update({columns[i] : node_list})
                ax.hlines(y = cutoff_u, xmin = 0, xmax = len(data['node']), linestyles = "dashed")
                ax.hlines(y = cutoff_l, xmin = 0, xmax = len(data['node']), linestyles = "dashed")
                ax.set(xlabel = 'residues', ylabel='Centrality Value')
            else:
                ax.set_visible(False)
        # for ax in axs.flat:
        #     ax.set(ylabel='Centrality Value')
        #plt.xlabel("Residues")
        #plt.ylabel("Centrality value")
        #plt.legend()
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
scaled_diff_df['node'] = node
scaled_diff_nodes = plot_centrality_vs_residues(scaled_diff_df, list(scaled_diff_df.columns)[:-1], [3], "scaled_diff")
print(get_count(scaled_diff_nodes))

#write_pdb_files(scaled_diff_df, pdb, "pdb/scaled_diff")


# ratio_df = wt_df / mt_df
# ratio_df['node'] = node
# ratio_df_nodes = plot_centrality_vs_residues(ratio_df, list(ratio_df.columns)[:-1], [3], "ratio")
# print(get_count(ratio_df_nodes))







