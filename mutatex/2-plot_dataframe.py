import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def heatmap(data, fname):
    plt.figure(figsize=(10,8))
    sns.heatmap(data, cmap="viridis")
    plt.tight_layout()
    plt.savefig(f"plots/{fname}.png")
    plt.clf()

def split_heatmap(data, num_splits, ext):
    start = 0
    step = len(data.columns)//num_splits + 1
    end = step
    for i in range(num_splits):
        cols = data.columns[start:end]
        # PLOT HEATMAP
        heatmap(data[cols], f"heatmap_{i}_{ext}")
        print(start, end)
        start +=step
        end +=step
        if end > ncol:
            end = ncol

def series_to_df(series, colname):
    df = series.to_frame()
    df.columns = [colname]
    df['Res'] = df.index
    df['Node'] = [int(r[1:]) for r in df['Res']]
    df = df.sort_values(by = "Node")
    df = df[["Node","Res", colname]]
    return df

data = pd.read_csv("mut_data/mean_df.txt", sep = "\t")
nrow, ncol = data.shape
order = list("GAVLIMFWPSTCYNQDEKRHsyp")
idx_dict = dict(zip(range(nrow), order))
data = data.rename(index = idx_dict)
ind_list = list("GAVLIMFWPSTCYNQDEKRH")
data = data.loc[ind_list]


# get mean_mean df
mean_series = data.mean(axis = 0).sort_values(ascending = False)
mean_df = series_to_df(mean_series, "Mean DDG")
mean_df.to_csv("mut_data/mean_per_pos_df.txt", sep = "\t", index = False)

# get median_mean df
median_series = data.median(axis = 0).sort_values(ascending = False)
median_df = series_to_df(median_series, "Median DDG")
median_df.to_csv("mut_data/median_per_pos_df.txt", sep = "\t", index = False)

# get count df
dgg_cutoff = 3
count_name = f"Count_{dgg_cutoff}"
count_series = (data > dgg_cutoff).sum(axis = 0).sort_values(ascending = False)
count_df = series_to_df(count_series, count_name)
count_df.to_csv(f"{count_name}_df.txt", sep = "\t", index = False)

#PLOT
# 5 splits = 32.6
num_splits = 5
# unchanged df
#split_heatmap(data, 5, "original")


#modify df values
cutoff = 3
data[data < cutoff] = 0
data[data > 0] = 1
# plot new

#split_heatmap(data, 5, f"bool_{cutoff}")

sum_df = data.sum(axis = 0).sort_values(ascending = False)
#sum_df.to_csv(f"mut_data/sum_cutoff_df.txt", sep = "\t", index = True)






