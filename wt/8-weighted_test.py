import numpy as np
import pandas as pd
import networkx as nx
#from pyinteraph import graph_analysis as ga
import networkx.algorithms.centrality as nxc
import MDAnalysis as mda
import matplotlib.pyplot as plt
import seaborn as sns

RES_DICT = {
    'ALA':  'A','ARG':	'R','ASN':	'N','ASP':	'D',
    'ASX':	'B','CYS':	'C','GLU':	'E','GLN':	'Q',
    'GLX':	'Z','GLY':	'G','HIS':	'H','ILE':	'I',
    'LEU':	'L','LYS':	'K','MET':	'M','PHE':	'F',
    'PRO':	'P','SER':	'S','THR':	'T','TRP':	'W',
    'TYR':	'Y','VAL':	'V'
    }

def get_identifiers(pdb):
    u = mda.Universe(pdb)
    identifiers = [f"{RES_DICT[r.resname]}{r.resnum}" for r in u.residues]
    return identifiers

def build_graph(data, pdb):
    identifiers = get_identifiers(pdb)
    G = nx.Graph(data)
    node_names = dict(zip(range(data.shape[0]), identifiers))
    nx.relabel_nodes(G, mapping = node_names, copy = False)
    return G, identifiers

unf = "hc_graph.dat"
fil = "hc_graph_filtered.dat"
pdb = "model0_A.pdb"
OUT_DIR = "test/"

U = np.loadtxt(unf)
U = U/np.max(U)
lg_U = -np.log(U + 10e-6)
G, ids = build_graph(U, pdb)
print(len(G.edges()))
lg_G, _ = build_graph(lg_U, pdb)
print(len(lg_G.edges()))

sns.heatmap(U, annot = True)
plt.savefig(f"{OUT_DIR}U")
plt.clf()

cent1 = nxc.betweenness_centrality(G)
print(cent1)
cent2 = nxc.betweenness_centrality(lg_G, weight = "weight")
print(cent2)


POS = nx.spring_layout(G, k = 0.55, seed = 1, iterations = 100)
#print(min(cent2.values()))

def plot_network(G, cent_dict, fname):
    weights = list(cent_dict.values())
    fig, ax = plt.subplots(figsize = (14,10))
    nx.draw_networkx_edges(G, POS)
    nc = nx.draw_networkx_nodes(G, POS, node_color = weights, cmap = plt.cm.viridis, edgecolors='black')
    plt.colorbar(nc)
    nx.draw_networkx_labels(G, POS)
    plt.axis('off')
    plt.savefig(f'{fname}.png')
    plt.clf()

max_edges = [(u,v) for u,v,d in lg_G.edges(data=True) if d["weight"] == 11.512925464970229]
lg_G.remove_edges_from(max_edges)
#print(len(lg_G.edges))

plot_network(G, cent1, f"{OUT_DIR}bet_G")
plot_network(lg_G, cent2, f"{OUT_DIR}bet_G_lg")

# edge_weights = [d["weight"] for u,v,d in lg_G.edges(data=True)]
# #print(max(edge_weights))
# plt.plot(max_weights)
# plt.savefig(f"{OUT_DIR}weights.png")

