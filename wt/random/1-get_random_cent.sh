graph_script=../9-random_graph.py
cent_script=../../../pyinteraph2/pyinteraph/centrality_analysis.py
psn=../hc_graph_filtered.dat
pdb=../model0_A.pdb
n_swaps=20
out=endpoints #CHANGE FILENAME



for i in {1..200}
do
    echo run $i
    python $graph_script -i $psn -p $pdb -n $n_swaps -o tmp
    python $cent_script -i tmp.dat -r $pdb -c degree betweenness closeness eigenvector -e -o tmp
    tail -n +2 tmp.txt >> $out
done

