import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp

OUT_DIR = ""
cent_file = "agg.txt"
cols = ["degree", "betweenness", "closeness", "eigenvector"]
cent_data = pd.read_csv(cent_file, sep = "\t", header= None)#[cols]
cent_data.columns = ['Node', 'Res'] + cols
cent_data = cent_data[cols]


def density_plot(data, measure, ext):
    plot = sns.distplot(data[measure])
    plot.figure.savefig(f"{OUT_DIR}{measure}_dist{ext}.png")
    plt.clf()


# for measure in cent_data.columns:
#     print(measure)
#     density_plot(cent_data, measure, "_agg")

cypa_data = pd.read_csv("../centrality/wt.txt", sep = "\t",)[cols]

for measure in cols:
    print(measure)
    stat, pval = sp.stats.kruskal(cypa_data[measure], cent_data[measure])
    print("Statistic:", stat)
    print("Pval:", pval)


