tail -n +6 model0.pdb | head -n -2 | sed 's/./A/22' > mid.tmp
head -n 5 model0.pdb > top.tmp
tail -n 2 model0.pdb > bot.tmp
cat top.tmp mid.tmp bot.tmp > model0_A.pdb
rm top.tmp mid.tmp bot.tmp
