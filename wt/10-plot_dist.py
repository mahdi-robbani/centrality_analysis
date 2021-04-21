import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

OUT_DIR = "plots/"
cent_file = "centrality/wt.txt"
cols = ["degree", "betweenness", "closeness", "eigenvector"]
cent_data = pd.read_csv(cent_file, sep = "\t")[cols]


def density_plot(data, measure, ext):
    plot = sns.distplot(data[measure])
    plot.figure.savefig(f"{OUT_DIR}{measure}_dist{ext}.png")
    plt.clf()


for measure in cent_data.columns:
    print(measure)
    density_plot(cent_data, measure, "_original")




