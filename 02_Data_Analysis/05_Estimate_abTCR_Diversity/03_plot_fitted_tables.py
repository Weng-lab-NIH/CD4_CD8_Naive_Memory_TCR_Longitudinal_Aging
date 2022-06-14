"""
Info: Load in tables

"""
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import argparse, os

def plot_one_fit(inp_df, my_color="blue", out_directory="./", out_name = "my_plot", log=False):
    if log:
        x = np.log10(inp_df["X"])
        y = np.log10(inp_df["Y"])
        fit_x = np.log10(inp_df["Fit_x"])
        fit_y = np.log10(inp_df["Fit_y"])
    else:
        x = inp_df["X"]
        y = inp_df["Y"]
        fit_x = inp_df["Fit_x"]
        fit_y = inp_df["Fit_y"]

    plt.clf()
    plt.figure(figsize=(6.4, 6.4))
    plt.plot(x, y, color=my_color, linestyle="None", marker="o")
    plt.plot(fit_x, fit_y, color=my_color)
    plt.savefig(os.path.join(out_directory, out_name+".png"))

def plot_both_fits(cd4, cd8, out_directory = "./", out_name = "my_log", log=False):
    if log:
        cd4_x = np.log10(cd4["X"])
        cd4_y = np.log10(cd4["Y"])
        cd4_fit_x = np.log10(cd4["Fit_x"])
        cd4_fit_y = np.log10(cd4["Fit_y"])
        cd8_x = np.log10(cd8["X"])
        cd8_y = np.log10(cd8["Y"])
        cd8_fit_x = np.log10(cd8["Fit_x"])
        cd8_fit_y = np.log10(cd8["Fit_y"])
    else:
        cd4_x = cd4["X"]
        cd4_y = cd4["Y"]
        cd4_fit_x = cd4["Fit_x"]
        cd4_fit_y = cd4["Fit_y"]
        cd8_x = cd8["X"]
        cd8_y = cd8["Y"]
        cd8_fit_x = cd8["Fit_x"]
        cd8_fit_y = cd8["Fit_y"]

    plt.clf()
    plt.figure(figsize=(6.4, 6.4))
    plt.plot(cd4_fit_x, cd4_fit_y, color="blue")
    plt.plot(cd8_fit_x, cd8_fit_y, color="red")
    plt.plot(cd4_x, cd4_y, color="blue", linestyle="None", marker="o")
    plt.plot(cd8_x, cd8_y, color="red", linestyle="None", marker="o")

    plt.savefig(os.path.join(out_directory, out_name+".png"), dpi=400)

if __name__ == "__main__":
    cd4_results = "C:/Users/Jeffrey/Documents/Programming/NIH_2021/sun paper/new_estimate_diversity/new_outs/cd4_fit.csv"
    cd8_results = "C:/Users/Jeffrey/Documents/Programming/NIH_2021/sun paper/new_estimate_diversity/new_outs/cd8_fit.csv"

    #C:/Users/Jeffrey/Documents/Programming/NIH_2021/sun paper/estimate_tcrs_diversity/02_model_fitting_results/fit_tables
    out_dir = "C:/Users/Jeffrey/Documents/Programming/NIH_2021/sun paper/new_estimate_diversity/new_outs/both"
    cd4_df = pd.read_csv(cd4_results)
    cd8_df = pd.read_csv(cd8_results)

    plot_one_fit(cd4_df, "blue", out_dir, "CD4", False)
    plot_one_fit(cd4_df, "blue", out_dir, "log_CD4", True)

    plot_one_fit(cd8_df, "red", out_dir, "CD8", False)
    plot_one_fit(cd8_df, "red", out_dir, "log_CD8", True)

    plot_both_fits(cd4_df, cd8_df, out_dir, "CD4CD8", False)
    plot_both_fits(cd4_df, cd8_df, out_dir, "log_CD4CD8", True)
