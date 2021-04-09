# usage: graph_analysis.py [-h] [-r TOPOLOGY] [-a DAT] [-c] [-u]
#                          [-k HUBS_CUTOFF] [-cb COMPONENTS_PDB] [-ub HUBS_PDB]
#                          [-d]

# PyInteraph network analysis module.

# optional arguments:
#   -h, --help            show this help message and exit
#   -r TOPOLOGY, --reference TOPOLOGY
#                         Reference topology file
#   -a DAT, --adj-matrix DAT
#                         Input graph file
#   -c, --components      Calculate connected components
#   -u, --hubs            Calculate hubs
#   -k HUBS_CUTOFF, --hubs-cutoff HUBS_CUTOFF
#                         Minimum number of connections for hubs (default: 3)
#   -cb COMPONENTS_PDB, --components-pdb COMPONENTS_PDB
#                         Save connected components ID in PDB file
#   -ub HUBS_PDB, --hubs-pdb HUBS_PDB
#                         Save hub degrees in PDB file
#   -d, --write-paths     Write the paths found as matrices

python ../../pyinteraph2/pyinteraph/graph_analysis.py -r model0_A.pdb -a hc_graph_filtered.dat -c -cb cc.pdb


