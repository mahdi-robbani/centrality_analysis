import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sasa = pd.read_csv("sasa.txt", sep = "\t")
sasa['node'] = sasa['ID'].apply(lambda x: 'A' + str(x))
print(list(sasa['node']))
cent = pd.read_csv("centrality/wt.txt", sep = "\t")
print(list(cent['node']))

cent_cols = ["node", "degree", "betweenness", "closeness", "eigenvector"]
data = pd.concat([sasa, cent[cent_cols]], axis = 1)
cor_cols = ["Mean", "degree", "betweenness", "closeness", "eigenvector"]
corr_data = data[cor_cols]
corr_data.columns = ["SASA", "Degree", "Betweenness", "Closeness", "Eigenvector"]
# sasa_col = -1 * corr_data['SASA'] / corr_data['SASA'].max()
# corr_data['SASA'] = sasa_col

def heatmap(data, colnames, fname, out_dir):
    cor = data[colnames].corr().round(2)
    plt.figure(figsize=(10,8))
    sns.heatmap(cor, cmap="viridis", annot = True)#, vmin = 0, vmax = 1)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/correlation_{fname}.pdf")
    plt.clf()

#heatmap(corr_data, corr_data.columns, "sasa", "plots")
#print(corr_data)


