import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest = "input", type = str)
parser.add_argument("-s", dest = "sort", nargs = "+")
args = parser.parse_args()

data = pd.read_csv(args.input, delimiter = "\t")
data = data[args.sort].sort_values(args.sort, ascending = False)
data.to_csv("sort.tmp", sep = "\t")
