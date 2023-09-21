'''
Note about the script:
It is a Python script to plot Levenstein distances between different lengths of barcodes to identify the right length of barcode and Levenstein distance to use for clustering. 

Output of this script:
    1. A pdf file containing the LV distance (of the full or partial barcode) histogram (log scale) for 3 samples

Command to run this script: python3 LvHistogramStarcode.py <experiment file> <Sample name> -length <30/40/50> --d <LV Distance>
'''
import regex as re
from argparse import ArgumentParser
import os
import numpy as np
import random
import itertools
import numpy as np
import Levenshtein
import matplotlib.pyplot as plt
import pandas as pd
import gc
from matplotlib.backends.backend_pdf import PdfPages

#Command line parser
parser = ArgumentParser()
parser.add_argument("experiment", help = "Specify the path to the experiment directory")
parser.add_argument("sampleName", help = "Specify the list of sample names")
parser.add_argument("--length", help = "Specify which starcode length has to be plotted", choices=["30","40","50"], default=50)
parser.add_argument("--d", help = "Specify the distance used for Levenshtein distance based clustering in starcode.", default="8")
args = parser.parse_args()

experimentDirectory = args.experiment

SampleList = args.sampleName.split(",")
for sample in SampleList:
    os.chdir('..')
    os.chdir(os.path.join(experimentDirectory, "analyzed", sample,"starcode")) #identifying the input file location and moving to it

    # Assigning the length of substring to determine LV distance 
    Barcode_raw = [] #to store all the Barcodes

    # Open required file and parse it to create an array of all barcodes
    readCountsFile = sample + "_Barcode{}_d{}.txt".format(args.length, args.d)

    with open(readCountsFile, "r") as reads:   
        for line in reads:
            Barcode_raw.append(line.split()[0])

    Barcode = np.unique(Barcode_raw)
    del(Barcode_raw)
    gc.collect()

    # Sampling the barcode to measure Levenshtein distance
    seed_value = 42
    random.seed(seed_value)
    Sample_1 = np.random.choice(Barcode, size=5000, replace=False)
    Sample_2 = np.random.choice(Barcode, size=5000, replace=False)
    Sample_3 = np.random.choice(Barcode, size=5000, replace=False)
    del(Barcode)
    gc.collect()

    # Calculate the string distance matrix using Levenshtein distance (edit distance)
    print("Beginning to find Levenshtein distance for the Barcode")

    # Sample 1
    pairs_1 = itertools.product(Sample_1, Sample_1)
    pairs_1 = [(pair[0], pair[1]) for pair in pairs_1 if pair[0] != pair[1]]
    matrix_1 = [Levenshtein.distance(pair[0], pair[1]) for pair in pairs_1]
    del(pairs_1)
    matrix_1 = np.array(matrix_1)
    del Sample_1
    gc.collect()
    print("Done for Sampling 1")

    # Sample 2
    pairs_2 = itertools.product(Sample_2, Sample_2)
    pairs_2 = [(pair[0], pair[1]) for pair in pairs_2 if pair[0] != pair[1]]
    matrix_2 = [Levenshtein.distance(pair[0], pair[1]) for pair in pairs_2]
    del(pairs_2)
    matrix_2 = np.array(matrix_2)
    del Sample_2
    gc.collect()
    print("Done for Sampling 2")

    # Sample 3
    pairs_3 = itertools.product(Sample_3, Sample_3)
    pairs_3 = [(pair[0], pair[1]) for pair in pairs_3 if pair[0] != pair[1]]
    matrix_3 = [Levenshtein.distance(pair[0], pair[1]) for pair in pairs_3]
    del(pairs_3)
    matrix_3 = np.array(matrix_3)
    del Sample_3
    gc.collect()
    print("Done for Sampling 3")

    # Save the merged files of plots together
    pdf_filename = sample + '_LV_distance_{}.pdf'.format(str(args.length))
    pdf = PdfPages(pdf_filename)

    x_tick = np.arange(0, max(np.max(matrix_1), np.max(matrix_2), np.max(matrix_3)), 4) 

    plt.figure(figsize=(6, 4))
    n_1, bins_1, patches_1 = plt.hist(matrix_1, np.max(matrix_1), density = True)
    plt.title('Sample 1')
    plt.yscale('log')
    plt.xlabel('LV Distance')
    plt.xticks(x_tick)
    plt.ylabel('Fraction')
    pdf.savefig()
    plt.close()

    plt.figure(figsize=(6, 4))
    n_2, bins_2, patches_2 = plt.hist(matrix_2, np.max(matrix_2), density = True)
    plt.title('Sample 2')
    plt.xticks(x_tick)
    plt.yscale('log')
    plt.xlabel('LV Distance')
    plt.ylabel('Fraction')
    pdf.savefig()
    plt.close()

    plt.figure(figsize=(6, 4))
    n_3, bins_3, patches_3 = plt.hist(matrix_3, np.max(matrix_3), density = True)
    plt.title('Sample 3')
    plt.xticks(x_tick)
    plt.yscale('log')
    plt.xlabel('LV Distance')
    plt.ylabel('Fraction')
    pdf.savefig()
    plt.close()
    pdf.close()



