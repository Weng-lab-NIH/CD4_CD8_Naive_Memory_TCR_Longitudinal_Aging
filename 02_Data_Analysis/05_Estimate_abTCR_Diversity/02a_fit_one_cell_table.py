"""
Info: Load in one cells table file and fit AB-B = x*A
    Reduces each cells table sample to # unique alpha, betas, TCRs. Combines samples with < 2000 cells. 

USE THIS COMMAND: python 02a_fit_one_cell_table.py 01_fit_input_data\cells_table_with_mixcr_cd4.tsv empty
Author: Jeffrey Cifello
"""
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import argparse, os
from scipy.optimize import minimize
from scipy.stats import pearsonr
import matplotlib.lines as mlines
import pickle

def get_subset_ua1(a1, a2):
    """Get a-prefactor for each sample, to average after"""
    if "" in a2:
        a2.remove("")
    if "None|None|None" in a2:
        a2.remove("None|None|None")
    all_a = a1 | a2
    a1_only = a1 - a2
    # fraction of alphas belonging only to primary
    pref = float(len(a1_only)) / len(all_a)
    
    # 
    # pref2 = float(len(a1)) / len(all_a)
    
    return pref

def get_sample_data(cells_csv):
    """Consolidate cell information to a smaller df with
    richness per chain per sample"""
    cells_df = pd.read_csv(cells_csv, sep="\t", low_memory=False)
    cells_df.fillna("", inplace=True)
    cells_df = cells_df.loc[cells_df["Source"] != "fontenot" ]
    samples = pd.unique(cells_df['Subject'])
    large_samples, small_samples = [], []
    combined_sources = set()
    small_samples_w_as = []

    print("original samples: {}".format(len(samples)))
    cells_df["small_w_2a"] = [ False for i in range(len(cells_df)) ]
    a_prefs = []
    small_a_prefs = []
    all_a_prefs = []
    for sample in samples:
        subset = cells_df.loc[cells_df['Subject']==sample]
        subset_idx = subset.index
        if len(subset) > 2000:
            large_samples.append(sample)
            if len(set(subset["Alpha 2"].tolist())) > 1:
                new_a1_fraction = get_subset_ua1(set(subset["Alpha 1"].tolist()), set(subset["Alpha 2"].tolist()))
                a_prefs.append(new_a1_fraction)
                # print(sample, new_a1_fraction)
                all_a_prefs.append(new_a1_fraction)

        else:
            small_samples.append(sample)
            cells_df.loc[subset_idx, 'Subject'] = 'combined'
            combined_sources.add(cells_df.loc[subset_idx, 'Source'].tolist()[0])

            if len(set(subset["Alpha 2"].tolist())) > 1:
                cells_df.loc[subset_idx, "small_w_2a"] = [ True for i in range(len(subset))]
                new_a1_fraction = get_subset_ua1(set(subset["Alpha 1"].tolist()), set(subset["Alpha 2"].tolist()))
                all_a_prefs.append(new_a1_fraction)
                # print(sample, new_a1_fraction, "small!")

    a_prefactor = np.mean(a_prefs)

    if len(small_samples) > 0:
        large_samples.append('combined')
        combined_sources = '&'.join(list(combined_sources))
        cells_df.loc[cells_df['Subject']=='combined', 'Source']=combined_sources

    chains_df = pd.DataFrame(index=large_samples, columns=['Source', 'Unique Alphas', 'Unique Betas', 'Unique TCRs'])

    for sample in large_samples:
        subset = cells_df.loc[cells_df['Subject']==sample]
        chains_df.loc[sample, 'Source'] = subset['Source'].tolist()[0]
        chains_df.loc[sample, 'Unique Alphas'] = len(pd.unique(subset['Alpha 1']))
        chains_df.loc[sample, 'Unique Betas'] = len(pd.unique(subset['Beta']))
        chains_df.loc[sample, 'Unique TCRs'] = len(pd.unique(subset['TCR']))

    return chains_df, all_a_prefs

def get_r_squared(y, predicted_y):
    """Manually calculate the R^2 value. This must be done manually because 
    the linregress function cannot be used. We are not allowing an intercept
    in this fitting. 
    """
    SStot = sum((y - np.mean(y))**2)
    SSreg = sum((predicted_y - np.mean(y))**2)
    Rsq = SSreg / SStot

    print("Rsq: {0:.3f}".format(Rsq))

def plot_fit(input_df, slope, p_val, output_dir, tag):
    """Create the plot for the fit, alongside the input data. 
    Also save the table used to generate this plot. 
    """
    max_chains = [ max(input_df.loc[idx, "Unique Alphas"], input_df.loc[idx, "Unique Betas"]) for idx in input_df.index ]
    max_chains = pd.Series(max_chains, index=input_df.index)
    x = input_df["Unique Alphas"]+input_df["Unique Betas"]
    y = input_df["Unique TCRs"] - max_chains
    fit_y = [min(slope*np.array(x)), max(slope*np.array(x))]
    fit_x = [min(x), max(x)]

    x = np.array(x.tolist())
    y = np.array(y.tolist())
    fit_x = np.array(fit_x)
    fit_y = np.array(fit_y)

    ## remove some outliers
    make_log=False
    if make_log:
        x = np.log10(x)
        y = np.log10(y)
        fit_x = np.log10(fit_x)
        fit_y = np.log10(fit_y)

    out_name=os.path.join(output_dir, "{}_tag.png".format(tag))
    plt.figure(figsize=(6.4, 6.4))
    plt.plot(fit_x, fit_y, "r")
    plt.plot(x, y, "ro")
    # plt.ylabel("\alpha\betaTCR - Max(TCR\alpha, TCR\beta)$")
    # plt.xlabel("TCR\alpha + TCR\beta$")
    plt.title(tag)
    plt.savefig(out_name, dpi=300)
    
    out_df = "{}\{}_fit.csv".format(output_dir, tag)
    fitted = pd.DataFrame(data={"X":x,
                            "Y":y,
                            "Fit_x":np.linspace(fit_x[0], fit_x[1], len(x)),
                            "Fit_y":np.linspace(fit_y[0], fit_y[1], len(y)) })
    fitted.to_csv(out_df, index=False)

def add_new_data_as_small(new_dat):
    """Take in the new data and return it as X,Y lists."""
    fraction_primary_list=[]
    print(new_dat)
    print(os.path.exists(new_dat))
    new_df = pd.read_csv(new_dat)
    
    large_samples_x=[]
    large_samples_y=[]
    large_samples_b=[]
    small_total_x, small_total_y, small_total_b = 0,0,0
    for i in range(len(new_df)):
        new_x=new_df.loc[i,"A1"]
        new_y=new_df.loc[i,"AB"]-new_df.loc[i,"A1"]
        new_b=new_df.loc[i,"B"]
        
        new_a1_frac = new_df.loc[i,"A1"] / (new_df.loc[i,"A1"] + new_df.loc[i,"A2"])
        fraction_primary_list.append(new_a1_frac)
        if new_df.loc[i, "cells"]>3000:
            large_samples_x.append(new_x)
            large_samples_y.append(new_y)
            large_samples_b.append(new_b)
        else:
            small_total_x += new_x
            small_total_y += new_y
            small_total_b += new_b
            if small_total_y >= 3000:
                large_samples_x.append(small_total_x)
                large_samples_y.append(small_total_y)
                large_samples_b.append(small_total_b)
                small_total_x=0
                small_total_y=0
                small_total_b=0
    
    return large_samples_x, large_samples_y, large_samples_b, fraction_primary_list

if __name__ == '__main__':
    for cell_type in ["cd4", "cd8"]:
        new_data="new_data/new_{}_data.csv".format(cell_type)
        new_x, new_y, new_b, a1_frac_list_new = add_new_data_as_small(new_data)
        
        cells_table="old_fit_data\cells_table_with_mixcr_{}.tsv".format(cell_type.lower())

        chains_df, a1_frac_list = get_sample_data(cells_table)
        
        chains_df["Subject"] = chains_df.index

        ## create a chains table for the new data
        #Index(['Source', 'Unique Alphas', 'Unique Betas', 'Unique TCRs', 'Subject']
        new_data_df = pd.DataFrame()
        new_data_df["Source"] = ["NewPub" for i in new_x]
        new_data_df["Unique Alphas"] = new_x
        new_data_df["Unique Betas"] = new_b
        new_data_df["Unique TCRs"] = np.array(new_y)+np.array(new_x)
        new_data_df["Subject"] = [str(i) for i in range(len(new_x))]

        chains_df = pd.concat([chains_df, new_data_df])

        max_chains = [ max(chains_df.loc[idx, "Unique Alphas"], chains_df.loc[idx, "Unique Betas"]) for idx in chains_df.index ]
        max_chains = pd.Series(max_chains, index=chains_df.index)
        Y = chains_df['Unique TCRs'] - max_chains
        X = chains_df['Unique Alphas']+chains_df['Unique Betas']

        X = np.array(X.tolist())
        Y = np.array(Y.tolist())

        # Y = np.concatenate([np.array(Y.tolist()),new_y])
        # X = np.concatenate([np.array(X.tolist()),new_x])
        
        ##
        all_a1_fractions= a1_frac_list + a1_frac_list_new
        average_a1_fraction = np.mean(all_a1_fractions)
        print("{} used to calculate prefactor!".format(len(all_a1_fractions)))
        ##
        num_before = len(X)
        if cell_type=="cd8":
            chains_df = chains_df[Y<2250]
            X = X[np.where(Y<2250)]
            Y = Y[np.where(Y<2250)]
            
        print("{} were removed!".format(num_before-len(X)))
        print("Fitting to {} samples!".format(len(X)))

        print("Fitting on {} samples.".format(len(chains_df)))
        print("Total cells: {}".format(sum(chains_df["Unique TCRs"])))

        fn = lambda a : sum((Y - (a*(X)))**2)
        slope = minimize(fn, 1.0).x[0]
        r, p_value = pearsonr(Y, X*slope)

        print("SLOPE: {0:.4f}".format(slope))
        print("PREFACTOR: {0:.3f}".format(average_a1_fraction))
        print("Pearson r: {0:.3f}".format(r))
        print("p-value: {0:.3e}".format(p_value))

        get_r_squared(Y, slope*X)

        plot_fit(chains_df, slope, p_value, "new_outs", cell_type)

        ## save alpha fractions
        
        with open("new_outs/{}_all_alpha_ratios.p".format(cell_type), "wb") as f:
            pickle.dump(all_a1_fractions, f)
