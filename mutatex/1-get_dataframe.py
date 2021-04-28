import pandas as pd
import os

data_dir = "mut_data/final_averages"

files = os.listdir(data_dir)
files = sorted(files, key = lambda x: int(x[2:]))
names = [f[0] + f[2:] for f in files]

mean_df = pd.DataFrame(columns = names)
for i, f in enumerate(files):
    path = f"{data_dir}/{f}"
    df = pd.read_csv(path, sep ="\s", header =None, skiprows= 1)
    mean_df[names[i]] = df[0]

mean_df.to_csv("mut_data/mean_df.txt", sep = "\t", index = False)



