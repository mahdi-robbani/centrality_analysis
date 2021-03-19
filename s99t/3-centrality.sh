dir=../../pyinteraph2/pyinteraph
script=$dir/centrality_analysis.py
inp=hc_graph_filtered.dat
pdb=model0.pdb 
pdb_A=model0_A.pdb
# remember to activate virtual environment

#Test all non group centralities
python $script -i $inp -r $pdb_A -c all -o centrality/s99t -p
#python $script -i $inp -r $pdb_A -c group -g A40:A60 -o centrality/group
