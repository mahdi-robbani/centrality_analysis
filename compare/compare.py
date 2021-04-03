import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import MDAnalysis as mda
from Bio import PDB
import argparse


wt = "../wt/centrality/flow.txt"
mt = "../s99t/centrality/s99t.txt"
pdb = "../wt/model0_A.pdb"

def get_df_basic(fname):
    df = pd.read_csv(fname, sep ="\t")
    colnames = ["degree", "betweenness", "closeness", "eigenvector"]
    df = df[colnames]
    return df

def get_diff_df(wt_df, mt_df, scale, abs = False):
    w_max = wt_df.max()
    m_max = mt_df.max()
    scale_val = w_max/m_max if scale else 1
    diff = (mt_df * scale) - wt_df
    diff = diff.abs() if abs else diff
    return diff

node = pd.read_csv(wt, sep = "\t")['node']
wt_df = get_df_basic(wt)
mt_df = get_df_basic(mt)


diff_df = get_diff_df(wt_df, mt_df, scale = False)
diff_df['node'] = node

def plot_centrality_vs_residues(data, columns, sds, fname = "", size = (9, 7)):
    nodes = []
    nrows = len(columns)//2 if len(columns) % 2 == 0 else len(columns)//2 + 1
    for sd in sds:
        fig, axs = plt.subplots(nrows, 2, figsize = size, sharex= True, sharey= True)
        for i, ax in enumerate(axs.flat):
            if i < len(columns):
                ax.plot(data[columns[i]], label = columns[i])
                ax.set_title(columns[i])
                mean_std = [data[columns[i]].mean(), data[columns[i]].std()]
                cutoff = mean_std[0] + sd*mean_std[1]
                for j, val in enumerate(data[columns[i]]):
                    if abs(val) > cutoff:
                        ax.annotate(data['node'][j], (j, val))
                        nodes.append(data['node'][j])
                ax.hlines(y = cutoff, xmin = 0, xmax = len(data['node']), linestyles = "dashed")
                ax.hlines(y = -cutoff, xmin = 0, xmax = len(data['node']), linestyles = "dashed")
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


plot_centrality_vs_residues(diff_df, list(diff_df.columns)[:-1], [3], "diff")

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

scaled_diff_df = get_diff_df(wt_df, mt_df, scale = True)
scaled_diff_df['node'] = node
plot_centrality_vs_residues(scaled_diff_df, list(scaled_diff_df.columns)[:-1], [3], "scaled_diff")






