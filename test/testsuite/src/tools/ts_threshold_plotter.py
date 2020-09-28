#!/usr/bin/env python

"""
COSMO Technical Testsuite

This script plots tolerance files.
"""

import numpy as np
import pandas as pd
import ts_thresholds
import matplotlib
import glob
import shutil
import os
import argparse
from itertools import islice, cycle

# information
__author__      = "Santiago Moreno"
__email__       = "cosmo-wg6@cosmo.org"
__maintainer__  = "xavier.lapillonne@meteoswiss.ch"

defaultInput = ["old/", "new/"]
defaultOutput = "plots/"

matplotlib.use("Agg")

def default_plot_creation():
    for filename in glob.iglob(args.inputFiles[0] + "/*"):
        filename = os.path.basename(filename)
        old_thresholds = ts_thresholds.Thresholds(args.inputFiles[0] + filename)
        new_thresholds = ts_thresholds.Thresholds(args.inputFiles[1] + filename)
        old_dict = old_thresholds.to_dict()
        new_dict = new_thresholds.to_dict()
        plot_data = {}

        for key, value in new_dict.items():
            if key not in ["minval", "steps"]:
               plot_data["new "  + key] = new_dict[key]
        for key, value in old_dict.items():
            if key not in ["minval", "steps"]:
               plot_data["old "  + key] = old_dict[key]
        old_dict_steps = old_dict["steps"]

        _plot(plot_data, filename, old_dict_steps)

def create_plot_from_files(files):
    plot_data = {}
    for file in files:
        filename = os.path.basename(file)
        threshold = ts_thresholds.Thresholds(file)
        threshold_dict = threshold.to_dict()
        for key, value in threshold_dict.items():
            if key not in ["minval", "steps", "CHKDAT"]:
                print("dbg: " + str(value) + " " + str(threshold_dict["steps"]))
                plot_data[key] = pd.Series(value, index=threshold_dict["steps"])
    _plot(plot_data, filename)

def _plot(plot_data,  title, index = ""):
    if index == "":
        df = pd.DataFrame(plot_data)
    else:
        df = pd.DataFrame(plot_data, index=index)
    my_colors = list(islice(cycle(['DarkGreen', 'DarkBlue', 'Yellow', 'Gray', 'DarkOrange', "Red", "Pink", "LightBlue", "Purple", "Pink"]), None, len(plot_data)))
    plot = df.plot(logy = True, marker = ".", title=title, figsize=(10,10), color=my_colors)
    plot.set_xlabel("steps")
    plot.set_ylabel("thresholds")
    plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig = plot.get_figure()
    fig.savefig(defaultOutput + title + ".png", bbox_inches='tight')
    matplotlib.pyplot.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A tool to create linegraphs of two \
             threshold files to visualize their differences")
    parser.add_argument("-i", "--input", dest="inputFiles", action="store", nargs="*", default=defaultInput,
             help="defines which files will be compared to each other\
             theres no limit to the amount of which can be defined. One can also define two folder names as arguments so \
             that the content of those two folders will get compared to each other. The files will only be compared if they have the same names")
    parser.add_argument("-o", "--output", dest="outputFolder", action="store",
             default=defaultOutput, help="define the destination folder for the generated pictures")
    args = parser.parse_args()
    os.makedirs(args.outputFolder, exist_ok=True)
    if (len(args.inputFiles) == 2) and (os.path.isdir(args.inputFiles[0])) and (os.path.isdir(args.inputFiles[1])):
        default_plot_creation()
    else:
        create_plot_from_files(args.inputFiles)


