import pandas as pd
import matplotlib.pyplot as plt

wt = "../wt/centrality/flow.txt"
mt = "../s99t/centrality/s99t.txt"

def get_df_basic(fname):
    df = pd.read_csv(fname, sep ="\t")
    colnames = ["degree", "betweenness", "closeness", "eigenvector"]
    df = df[colnames]
    return df

node = pd.read_csv(wt, sep = "\t")['node']
wt_df = get_df_basic(wt)
mt_df = get_df_basic(mt)


diff_df = (wt_df - mt_df).abs()
diff_df['node'] = node

def plot_centrality_vs_residues(data, columns, sds, fname = "", size = (9, 7)):
    nrows = len(columns)//2 if len(columns) % 2 == 0 else len(columns)//2 + 1
    for sd in sds:
        fig, axs = plt.subplots(nrows, 2, figsize = size, sharex= True, sharey= True)
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
        plt.savefig(f"{fname}.pdf")
        plt.clf()


plot_centrality_vs_residues(diff_df, list(diff_df.columns)[:-1], [3], "diff")