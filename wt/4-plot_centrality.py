import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr

# set directory
file_path = "centrality/flow.txt"
out_dir = "plots/"

def replace_dict(names):
    new_names = []
    for name in names:
        if "_" in name:
            name = name.split("_")
            new_name = '_'.join((name[0], name[-2], name[-1]))
            new_names.append(new_name)
        else:
            new_names.append(name)
    return new_names

def plot_centrality_vs_residues(data, columns, sds, fname = "", size = (9, 7)):
    nrows = len(columns)//2 if len(columns) % 2 == 0 else len(columns)//2 + 1
    for sd in sds:
        fig, axs = plt.subplots(nrows, 2, figsize = size, sharex= True, sharey= False)
        for i, ax in enumerate(axs.flat):
            if i < len(columns):
                ax.plot(data[columns[i]], label = columns[i])
                ax.set_title(columns[i])
                mean_std = [data[columns[i]].mean(), data[columns[i]].std()]
                cutoff = mean_std[0] + sd*mean_std[1]
                for j, val in enumerate(data[columns[i]]):
                    if val > cutoff:
                        ax.annotate(data['node'][j], (j, val))
                ax.hlines(y = cutoff, xmin = 0, xmax = len(data['node']), linestyles = "dashed")
                ax.set(xlabel = 'residues', ylabel='Centrality Value')
            else:
                ax.set_visible(False)
        # for ax in axs.flat:
        #     ax.set(ylabel='Centrality Value')
        #plt.xlabel("Residues")
        #plt.ylabel("Centrality value")
        #plt.legend()
        plt.savefig(f"{out_dir}/centrality_{sd}_{fname}.pdf")
        plt.clf()

def heatmap(data, colnames, fname=""):
    cor = data[colnames].corr().round(2)
    plt.figure(figsize=(8,8))
    ax = sns.heatmap(cor, square = True, annot = True, cmap="viridis", vmin=0, vmax=1)
    ax.set_yticks(np.arange(0, len(colnames)+1))
    # print(np.arange(0.5, len(cor.index), 1))
    # plt.xticks(np.arange(0, len(cor.index), 1), cor.index)
    # plt.yticks(np.arange(-0.5, len(cor.index), 1), cor.index)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/correlation_{fname}.pdf")

    plt.clf()

#load file
data = pd.read_csv(file_path, sep ="\t")
# remove node names, res name, and hubs
colnames = list(data.columns)[2:]
new_colnames = replace_dict(colnames)
data = data.rename(columns = dict(zip(colnames, new_colnames)))

#plot centrality vs residues
sds = [3]
# plot basic
basic = ['degree', 'betweenness', 'closeness', 'eigenvector']
plot_centrality_vs_residues(data, basic, sds, "basic")
# plot only flow centralities
flow = ['current_betweenness_c0', 'current_closeness_c0', 'current_betweenness_c1', 'current_closeness_c1', 'current_betweenness_c2', 'current_closeness_c2']
plot_centrality_vs_residues(data, flow, sds, "flow")

#plot heatmap
heatmap(data, new_colnames, "flow") #xticklabels

#combine connected component cols

def combine_cols(data, cc_cols):
    for col in cc_cols:
        individual_cols = [name for name in new_colnames if col in name]
        combined = data[individual_cols].sum(axis = 1)
        data[f"{col}_c"] = combined
        x = individual_cols + [f"{col}_c"]
        data = data.drop(columns = individual_cols)
    return data

del_cols = ["current_betweenness", "current_closeness"] #maybe add communicability later on
data_combined = combine_cols(data, del_cols)
combined_colnames = list(data_combined.columns)[2:]
#plot 2
heatmap(data_combined, combined_colnames, "flow_combined")
plot_centrality_vs_residues(data_combined, combined_colnames, sds, "flow_combined", (10, 20))





# def plot_centrality_vs_residues(data, columns, sds, fname = ""):
#     #BACKUP
#     for sd in sds:
#         fig, axs = plt.subplots(len(columns), figsize = (10, 20), sharex= True, sharey= False)
#         for i, colname in enumerate(columns):
#             axs[i].plot(data[colname], label = colname)
#             axs[i].set_title(colname)
#             mean_std = [data[colname].mean(), data[colname].std()]
#             cutoff = mean_std[0] + sd*mean_std[1]
#             for j, val in enumerate(data[colname]):
#                 if val > cutoff:
#                     axs[i].annotate(data['node'][j], (j, val))
#             axs[i].hlines(y = cutoff, xmin = 0, xmax = len(data['node']), linestyles = "dashed") 
#         # for ax in axs.flat:
#         #     ax.set(ylabel='Centrality Value')
#         plt.xlabel("Residues")
#         #plt.ylabel("Centrality value")
#         #plt.legend()
#         plt.savefig(f"{out_dir}/centrality_{sd}_{fname}.pdf")
#         plt.clf()