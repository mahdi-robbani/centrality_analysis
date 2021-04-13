dir=../../pyinteraph2/pyinteraph
script=$dir/centrality_analysis.py
inp=hc_graph_filtered.dat
unf=hc_graph.dat
pdb=model0.pdb 
pdb_A=model0_A.pdb
# remember to activate virtual environment

# get all centralities
#python $script -i $inp -r $pdb_A -c all -o centrality/wt #-p
python $script -i $unf -r $pdb_A -w -c all -o unfiltered/wt_unf_weight
python $script -i $unf -r $pdb_A -c all -o unfiltered/wt_unf


#python $script -i $inp -r $pdb_A -c edge_betweenness -o centrality/edge