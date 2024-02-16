import sys
import argparse
import operator

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

def plotProbTable(table, ax, title, zoomed = False, xlim_ratio=3):
    if zoomed:
        scale_factor = 1e3
        ax.set_ylabel("Frequency (Kb)",fontsize=14)
    else:
        scale_factor = 1e6
        ax.set_ylabel("Frequency (Mb)", fontsize=14)
    freq = table["freq"] / scale_factor
    fit = table["fit"] / scale_factor
    ax.plot(table["#coverage"], freq, color="black", label="actual")
    ax.plot(table["#coverage"], fit, color="blue", label="fit")
    ax.plot(table["#coverage"], table["error"] * fit, color="red", label="error")
    ax.plot(table["#coverage"], table["duplicated"] * fit, color="orange", label="duplicated")
    ax.plot(table["#coverage"], table["haploid"] * fit, color="green", label="haploid")
    ax.plot(table["#coverage"], table["collapsed"] * fit, color="purple", label="collapsed")
    max_x = np.argmax(table["freq"]) * xlim_ratio
    if zoomed:
        max_y = 0.1 + 4 * max(np.max(table["collapsed"] * fit), \
                              np.max(table["duplicated"] * fit), \
                              np.max(table["error"] * fit))
    else:
        max_y = 1.25 * np.max(freq)
    ax.set_xlim(0, max_x)
    ax.set_ylim(0, max_y)
    ax.set_xlabel("Coverage", fontsize=14)
    ax.tick_params(axis='x', labelsize= 14)
    ax.tick_params(axis='y', labelsize= 14)
    ax.set_title(title, fontsize=16)
    ax.legend()

def plotPairDist(table, pdf, suptitle):
    fig = plt.figure(figsize=(12,6))
    ax1 = plt.axes([0.1, 0.2, 0.35, 0.5])
    ax2 = plt.axes([0.6, 0.2, 0.35, 0.5])
    ax3 = plt.axes([0.85, 0.8, 0.1, 0.1])
    plotProbTable(table, ax1, "Coverage Distribution")
    plotProbTable(table, ax2, "Coverage Distribution (Zoomed in)", zoomed = True)
    plotSummaryTable(ax3, table)
    fig.suptitle(suptitle, fontsize=16)
    pdf.savefig(fig)
    plt.close()

def plotTitlePage(title, pdf):
    fig = plt.figure(figsize=(12,6))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    ax.text(0.5, 0.5, title, fontsize=48, va='center', ha='center')
    ax.axis("off")
    pdf.savefig(fig)
    plt.close()

def plotSummaryTable(ax, table):
    ax.axis('off')
    ax.axis('tight')
    # Calculate total
    fit_total = sum(table["fit"])
    freq_total = int(sum(table["freq"]))
    error_total = sum(table["fit"] * table["error"])
    duplicated_total = sum(table["fit"] * table["duplicated"])
    haploid_total = sum(table["fit"] * table["haploid"])
    collapsed_total = sum(table["fit"] * table["collapsed"])
    # Calculate fractions
    error_per = round(error_total / fit_total *100, 3)
    duplicated_per = round(duplicated_total / fit_total * 100, 3)
    haploid_per = round(haploid_total / fit_total * 100, 3)
    collapsed_per = round(collapsed_total / fit_total * 100, 3)
    data = np.array([round(freq_total / 1e6, 3), error_per, duplicated_per, haploid_per, collapsed_per])
    df = pd.DataFrame(data)
    ax.table(cellText=df.values, rowLabels=["Size (Mb)",
                                            "Error (%)",
                                            "Duplicated (%)",
                                            "Haploid (%)",
                                            "Collapsed (%)"], loc='center', fontsize=16)


def main():
    parser = argparse.ArgumentParser(description='Generate a pdf showing the whole-genome and contig-specific coverage distributions')
    parser.add_argument('--table', type=str, help='path to the whole genome probability table file')
    parser.add_argument('--dir', type=str, help='directory that contains contig-specific probability table files')
    parser.add_argument('--pdf', type=str, help='path to output pdf')
    parser.add_argument('--diploid', action='store_true', help='set True if the reference is diploid')
    parser.add_argument('--hap1Pattern', type=str, help='pattern for finding hap1 contig names (like #1# or pat or h1tg)', default="h1tg")
    parser.add_argument('--hap2Pattern', type=str, help='pattern for finding hap2 contig names (like #2# or mat or h2tg)', default="h2tg")
    parser.set_defaults(diploid=False)
    args = parser.parse_args()
    tablePath = args.table
    tablesDir = args.dir
    outputPath = args.pdf
    isDiploid = args.diploid
    hap1Pattern = args.hap1Pattern
    hap2Pattern = args.hap2Pattern

    import matplotlib.backends.backend_pdf
    pdf = matplotlib.backends.backend_pdf.PdfPages(outputPath)
    prefix = os.path.basename(tablePath)[:-len(".table")]

    # Plot the whole-genome coverage distribution
    plotTitlePage("Whole Genome", pdf)
    table = pd.read_csv(tablePath, sep="\t")
    plotPairDist(table, pdf, "{}\n{}".format(prefix, "whole genome"))
    # Plot contig-specific coverage distributions
    from os import listdir
    from os.path import isfile, join
    tableFiles = [join(tablesDir, f) for f in listdir(tablesDir)]
    # x.split("_")[-2] is the start position of the window
    # "_".join(x.split("_")[0:-2]) is everything before that including the contig name
    #tableFiles= sorted(tableFiles, key = lambda x:("_".join(x.split("_")[0:-2]) , int(x.split("_")[-2])))
    if isDiploid:
        hap1TableFiles = [tableFile for tableFile in tableFiles if hap1Pattern in tableFile]
        hap2TableFiles = [tableFile for tableFile in tableFiles if hap2Pattern in tableFile]
        tableFilesPerHap = [[hap1TableFiles, "Hap1"], [hap2TableFiles, "Hap2"]]
    else:
        tableFilesPerHap = [[tableFiles, "All"]]
    windowname = "##"
    plotTitlePage("Window-Specific Coverage Distributions", pdf)
    for tableFilesOneHap, hapName in tableFilesPerHap:
        plotTitlePage(f"{hapName} Contigs", pdf)
        for p in tableFilesOneHap:
            table = pd.read_csv(p, sep="\t")
            filename = os.path.basename(p)
            # Extract contig name from file name
            windowname = filename[(len(prefix) + 1):-len(".table")]
            attrbs = windowname.strip().split("_")
            plotPairDist(table, pdf, "{}\n{}\n{}-{}".format(prefix, "_".join(attrbs[0:-2]), attrbs[-2], attrbs[-1]))
    pdf.close()

if __name__ == "__main__":
    main()
