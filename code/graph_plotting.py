# Dependencies
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'Lucida Grande' 
plt.rcParams['font.size'] = 12               
plt.rcParams['axes.titlesize'] = 14          
plt.rcParams['axes.labelsize'] = 12          
plt.rcParams['xtick.labelsize'] = 10         
plt.rcParams['ytick.labelsize'] = 10         
plt.rcParams['legend.fontsize'] = 10 

plt.rcParams['axes.grid'] = True
plt.rcParams['axes.grid.axis'] = 'y'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.linewidth'] = 0.5
plt.rcParams['grid.alpha'] = 0.6

orange = "#FFA500"
teal = "#008080"

# Scatterplot of observed mutations
def scatterplot_observed (data: pd.DataFrame, sorted: bool = False):
    # Scatterplot of observed mutations
    if (sorted):
        data = data.sort_values(by="M. tuber", ascending=True).reset_index(drop=True)
    
    plt.figure(figsize=(8, 5))
    plt.scatter(x = data[data["M. tuber"]!= 0].index, 
                y = data[data["M. tuber"]!= 0]['M. tuber'], 
                color = orange,
                label = "Observed Mutations",)
    
    plt.scatter(x = data[data["M. tuber"] == 0].index,
                y = data[data["M. tuber"] == 0]['M. tuber'],
                color = teal, 
                alpha = 0.5,
                s = 10,
                label = "Zero observations")

    #plt.ylim(0, 800)   
    plt.ylabel("Number of Observed Mutations in MTB")
    if (sorted):
        plt.xlabel("All Records of Amino Acid Substitutions Sorted by Number of Observed Mutations")
    else:
        plt.xlabel("All Records of Amino Acid Substitutions Sorted by Codon Name")

    plt.xticks([])

    plt.tight_layout()
    plt.legend()
    plt.show()

# Histogram of all observed mutations
def histogram_all (data: pd.DataFrame, log_scale : bool):
    plt.figure(figsize=(7,5))
    
    counts, bins, patches = plt.hist(data['M. tuber'], bins=30, color=orange, edgecolor='white')

    # Najdeme bin, který obsahuje 0, a přebarvíme ho
    for i in range(len(bins) - 1):
        # Zkontrolujeme, zda je 0 uvnitř intervalu tohoto binu
        if bins[i] <= 0 < bins[i+1]:
            patches[i].set_facecolor(teal)
            #patches[i].set_edgecolor('teal')
            patches[i].set_alpha(0.8) 

    if (log_scale):
        plt.yscale('log')
        plt.yticks([1, 5, 10, 50, 100], [1, 5, 10, 50, 100])
    
    plt.xlabel('Number of Observed Mutations in MTB')
    plt.ylabel('Count')

    #plt.grid(False)
    #plt.grid(axis='y', linestyle='--', linewidth=0.5, alpha=0.6)

    plt.tight_layout()
    plt.show()

# Generate a boxplot
def boxplot(table : pd.DataFrame, log_scale : bool):
    plt.figure(figsize=(7,5))

    if log_scale:
        # log(1 + x) transformace
        y_data = np.log1p(table['M. tuber'])
        yticks = np.log1p([0, 1, 10, 100, 500, 1000])
        ylabels = [0, 1, 10, 100, 500, 1000]
    else:
        y_data = table['M. tuber']
        yticks = None
        ylabels = None

    sns.boxplot(x = table['Mutation_Pathways'],
                y = y_data,
                linewidth=1, 
                palette = [orange, "salmon", teal])

    if (log_scale):
        plt.yticks(yticks, ylabels)

    # Labels
    plt.xlabel("Mutation Pathways")
    plt.ylabel("Observed Mutations")

    plt.show()

# Boxplot divided by frequencies
# Generate a boxplot
def boxplot_freq(table : pd.DataFrame, log_scale : bool):
    plt.figure(figsize=(7,5))

    if log_scale:
        # log(1 + x) transformace
        y_data = np.log1p(table['LeftSide'])
        yticks = np.log1p([0, 1, 10, 100, 500, 1000])
        ylabels = [0, 1, 10, 100, 500, 1000]
    else:
        y_data = table['LeftSide']
        yticks = None
        ylabels = None

    sns.boxplot(x = table['Mutation_Pathways'],
                y = y_data,
                linewidth=1, 
                palette = [teal, "cornflowerblue", orange])

    if (log_scale):
        plt.yticks(yticks, ylabels)

    # Labels
    plt.xlabel("Genetic Code Pathways") # fontsize = 10
    plt.ylabel("Observed Mutations") # fontsize = 10

    plt.show()

# Generate a regression plot
def plot_regression(data: pd.DataFrame, x_vals: np.ndarray, y_pred: np.ndarray):
    plt.figure(figsize=(8, 5))
    
    y_data = np.log1p(data['LeftSide'])
    yticks = np.log1p([0, 1, 10, 100, 500, 800])
    ylabels = [0, 1, 10, 100, 500, 800]

    x_data = data['Mutation_Pathways'].copy()
    
    # Add small random noise (jitter) where x == 1
    jitter_strength = 0.03  # How much jitter
    jitter = np.random.uniform(-jitter_strength, jitter_strength, size=x_data.shape)
    x_data_jittered = x_data + np.where(x_data == 1, jitter, 0)

    x_ticks = [1, 2, 3]
    x_labels = [1, 2, 3]

    plt.scatter(x_data_jittered, y_data, alpha=0.6, label="Observed data", color = teal)

    plt.plot(x_vals, y_pred, color=orange, linewidth=2, label="Fitted regression line")

    plt.xlabel("Genetic Code Pathways")
    plt.ylabel("Observed Mutations")
    plt.legend()

    plt.yticks(yticks, ylabels)
    plt.xticks(x_ticks, x_labels)

    plt.tight_layout()
    plt.grid(True, axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.show()

# Generating barplots for the codons of interest
def histogram_codon (data: pd.DataFrame, codon : str):
    data_subset = data[data['Codon'] == codon]

    unique_aa = data_subset['Alt_aa'].unique()
    colors = [teal, "salmon", "cornflowerblue", "mediumseagreen", "hotpink", orange, "mediumpurple", "darkseagreen", "firebrick"] #plt.colormaps.get_cmap('Set2')  # Use a colormap with enough colors

    plt.bar(data_subset['Alt_aa'], data_subset['M. tuber'], color=[colors[i] for i in range(len(unique_aa))])
    plt.xlabel('Mutated Amino Acid')
    plt.ylabel('Number of Observed Point Mutations')

    plt.gca().yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.show()

def hist_all_alt (data: pd.DataFrame):
    # Histogram of observed mutations
    plt.figure(figsize=(8, 5))
    sns.histplot(data["M. tuber"], bins=20, kde=True, color="red")
    plt.xlabel("Observed Mutations (M. tuber)")
    plt.ylabel("Frequency")
    plt.title("Histogram of Observed Mutations")
    #plt.grid(True)
    plt.show()

def combine_distribution_graph(data: pd.DataFrame, sorted: bool = False, log_scale: bool = False):
    # Standalone function that directly plots both scatterplot and histogram without calling subfunctions
    fig, axs = plt.subplots(2, 1, figsize=(8, 10), sharex=False)

    # === Scatterplot ===
    if sorted:
        data = data.sort_values(by="M. tuber", ascending=True).reset_index(drop=True)

    axs[0].scatter(
        x=data[data["M. tuber"] != 0].index,
        y=data[data["M. tuber"] != 0]["M. tuber"],
        color=orange,
        label="Observed Mutations"
    )
    axs[0].scatter(
        x=data[data["M. tuber"] == 0].index,
        y=data[data["M. tuber"] == 0]["M. tuber"],
        color=teal,
        alpha=0.5,
        s=10,
        label="Zero observations"
    )
    axs[0].set_ylabel("Number of Observed Mutations in MTB")
    xlabel = "All Records of Amino Acid Substitutions Sorted by Number of Observed Mutations" if sorted else "All Records of Amino Acid Substitutions Sorted by Codon Name"
    axs[0].set_xlabel(xlabel)
    axs[0].set_xticks([])
    axs[0].legend()
    axs[0].grid(True, axis='y', linestyle='--', alpha=0.4)

    # === Histogram ===
    counts, bins, patches = axs[1].hist(data["M. tuber"], bins=30, color=orange, edgecolor='white')
    for i in range(len(bins) - 1):
        if bins[i] <= 0 < bins[i+1]:
            patches[i].set_facecolor(teal)
            patches[i].set_alpha(0.8)

    if log_scale:
        axs[1].set_yscale('log')
        axs[1].set_yticks([1, 5, 10, 50, 100])
        axs[1].set_yticklabels([1, 5, 10, 50, 100])

    axs[1].set_xlabel("Number of Observed Mutations in MTB")
    axs[1].set_ylabel("Count")
    axs[1].grid(True, axis='y', linestyle='--', alpha=0.4)

    plt.tight_layout()
    plt.show()

def transition_bias(data: pd.DataFrame):
    plt.figure(figsize=(5,4))
    sns.barplot(x=data.index, y=data.values, palette = [teal, orange])
    plt.ylabel('Total Events')
    plt.xlabel('Mutation Type')
    plt.tight_layout()
    plt.show()

def substitution_bias(data: pd.DataFrame):
    plt.figure(figsize=(10,6))
    sns.barplot(x=data.index, y=data.values, palette=[teal, "cornflowerblue", orange, "#4EA72E", "#EF857D", "#156082"])
    plt.ylabel('Total Events')
    plt.xlabel('Single-nucleotide Substitution Type')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()