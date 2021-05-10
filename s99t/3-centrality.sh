dir=../../pyinteraph2/pyinteraph
script=$dir/centrality_analysis.py
inp=hc_graph_filtered.dat
pdb_A=model0_A.pdb
# remember to activate virtual environment

#node centralities (use edge weights)
python $script -i $inp -r $pdb_A -e -c all -o centrality/s99t -p
#edge centralities
#python $script -i $inp -r $pdb_A -c edge_betweenness -o centrality/edge


