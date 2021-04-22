from pyinteraph import path_analysis as pa
from edgeswap import EdgeSwapGraph
import networkx as nx
import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', dest = "psn", help='input psn')
parser.add_argument('-p', dest = "pdb", help='input pdb')
parser.add_argument('-n', dest = "n", type = int, default = 10, help='number of iterations')
parser.add_argument('-o', dest = "out", help='output random mat')
args = parser.parse_args()


psn = args.psn#"hc_graph_filtered.dat"
pdb = args.pdb#"model0_A.pdb"
n = args.n
out = args.out

_, _, G = pa.build_graph(psn, pdb = pdb)
randomG = EdgeSwapGraph(G).randomize_by_edge_swaps(n)
randomG_mat = nx.to_numpy_matrix(randomG)
np.savetxt(f"{args.out}.dat", randomG_mat)

