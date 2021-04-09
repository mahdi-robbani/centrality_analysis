import os
import numpy as np
import pandas as pd

name = []
res_id = []
mean = []
sd = []
maxim = []
minim = []
for res in os.listdir(os.getcwd()+'/SOLV/XVG'):
    name.append(''.join(s for s in res[:-4] if s.isalpha()))
    res_id.append(int(''.join(s for s in res[:-4] if s.isdigit())))
    sasa = np.loadtxt("SOLV/XVG/"+res, comments=["@", "#"])[:, 1]
    mean.append(sasa.mean())
    sd.append(sasa.std())
    maxim.append(sasa.max())
    minim.append(sasa.min())
    print(f"{res} completed")
    # if len(name) == 10:
    #     break

dct = {"Name": name, "ID": res_id, "Mean" : mean, "SD": sd, "Max": maxim, "Min" : minim}
df = pd.DataFrame(dct).sort_values(by="ID")
df.to_csv("sasa.txt", sep="\t", index = False)

