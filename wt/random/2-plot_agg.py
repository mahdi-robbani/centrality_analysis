import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp


# def density_plot(data, measure, ext):
#     plot = sns.distplot(data[measure])
#     plot.figure.savefig(f"{OUT_DIR}{ext}.png")
#     plt.clf()

def violin_plot(data, ext):
    plot = sns.violinplot(data = data, x = "Measure", y  = "Value")
    plt.tight_layout()
    plot.figure.savefig(f"cent_dist_{ext}.pdf")
    plt.clf()

#set variables
OUT_DIR = ""
cent_file = "endpoints"
cols = ["degree", "betweenness", "closeness", "eigenvector"]
# load file
cent_data = pd.read_csv(cent_file, sep = "\t", header= None)
cent_data.columns = ['Node', 'Res'] + cols
cent_data = cent_data[cols]
#convert to long format
cent_data = pd.melt(cent_data, value_vars = cols, var_name = 'Measure', value_name = "Value")
# plot
violin_plot(cent_data, "random")

#plot original file
original_file = "../centrality/wt.txt"
original_data = pd.read_csv(original_file, sep = "\t")
original_data = original_data[cols]
original_long = pd.melt(original_data, value_vars = cols, var_name = 'Measure', value_name = "Value")
violin_plot(original_long, "original")

#combined
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 6), sharey=True)
ax1.set_title('Random network')
sns.violinplot(ax = ax1, data = cent_data, x = "Measure", y  = "Value")
ax2.set_title('Original network')
sns.violinplot(ax = ax2, data = original_long, x = "Measure", y  = "Value")
plt.tight_layout
plt.savefig("combined.pdf")

# for measure in cent_data.columns:
#     print(measure)
#     density_plot(cent_data, measure, "_one_edge_2")

# cypa_file = "../centrality/wt.txt"
# cypa_data = pd.read_csv(cypa_file, sep = "\t",)[cols]
# n = cypa_data.shape[0] + cent_data.shape[0]
# print("\n")
# for measure in cols:
#     print(measure)
#     stat, pval = sp.stats.kruskal(cypa_data[measure], cent_data[measure])
#     eta = (stat - 1)/(n - 2)
#     print("Statistic:", stat)
#     print("eta:", eta)
#     print("Pval:", pval)


