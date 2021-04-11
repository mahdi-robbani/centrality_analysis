b=betweenness
c=closeness
d=degree

python sort.py -i wt.txt -s $c $b $d name node
head -n 5 sort.tmp
python sort.py -i wt.txt -s $b $c $d name node
head -n 5 sort.tmp
python sort.py -i wt.txt -s $d $c $b name node
head -n 5 sort.tmp