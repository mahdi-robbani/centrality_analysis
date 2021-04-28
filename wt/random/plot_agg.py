import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp

OUT_DIR = ""
cent_file = "one_edge_2.txt"
cols = ["degree", "betweenness", "closeness", "eigenvector"]
cent_data = pd.read_csv(cent_file, sep = "\t", header= None)#[cols]
cent_data.columns = ['Node', 'Res'] + cols
cent_data = cent_data[cols]


def density_plot(data, measure, ext):
    plot = sns.distplot(data[measure])
    plot.figure.savefig(f"{OUT_DIR}{measure}_dist{ext}.png")
    plt.clf()


for measure in cent_data.columns:
    print(measure)
    density_plot(cent_data, measure, "_one_edge_2")

cypa_file = "../centrality/wt.txt"
cypa_data = pd.read_csv(cypa_file, sep = "\t",)[cols]
n = cypa_data.shape[0] + cent_data.shape[0]
print("\n")
for measure in cols:
    print(measure)
    stat, pval = sp.stats.kruskal(cypa_data[measure], cent_data[measure])
    eta = (stat - 1)/(n - 2)
    print("Statistic:", stat)
    print("eta:", eta)
    print("Pval:", pval)


