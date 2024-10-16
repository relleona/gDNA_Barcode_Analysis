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
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages


def create_histogram(matrix, title, pdf, x_tick):
    # Create a figure for the normal scale histogram
    plt.figure(figsize=(6, 4))
    num_bins = np.max(matrix)  # Using the maximum value as the number of bins
    plt.hist(matrix, bins=num_bins, density=True, linewidth=10)  # Reduced linewidth for better visibility
    plt.title(title)
    plt.xlabel('LV Distance')
    plt.ylabel('Fraction')
    plt.xticks(x_tick)  # Assuming x_tick is defined elsewhere
    plt.tight_layout()
    pdf.savefig()  # Save the normal scale histogram
    plt.close()

    # Create a figure for the logscale histogram
    plt.figure(figsize=(6, 4))
    plt.hist(matrix, bins=num_bins, density=True, linewidth=10)  # Same histogram, but we'll change the y-axis scale
    plt.title(title + " (logscale)")
    plt.yscale('log')  # Set y-axis to logarithmic scale
    plt.xlabel('LV Distance')
    plt.ylabel('Fraction (log scale)')
    plt.xticks(x_tick)  # Assuming x_tick is defined elsewhere
    plt.tight_layout()
    pdf.savefig()  # Save the logscale histogram
    plt.close()


def analyze_LV(Sample, csvfilepath, sample_number):
    # Generate all pairs, including self-pairs
    pairs_1 = itertools.product(Sample, Sample)
    
    # Filter out self-pairs
    pairs_1 = [(pair[0], pair[1]) for pair in pairs_1 if pair[0] != pair[1]]
    
    # Calculate Levenshtein distances
    matrix_1 = [Levenshtein.distance(pair[0], pair[1]) for pair in pairs_1]
    
    # Free up memory
    del(pairs_1)
    
    # Convert to numpy array
    matrix_1 = np.array(matrix_1)
    
    # Saving the matrix to reproduce data
    pd.DataFrame(matrix_1).to_csv(csvfilepath, index=False, header=False)
    
    # Calculating mean
    sample1_mean = sum(matrix_1) / len(matrix_1)
    
    # Clean up to free memory
    del Sample
    gc.collect()
    
    print(f"Done for Sampling {sample_number}")

    return sample1_mean, matrix_1


if __name__ == "__main__":
    #Command line parser
    # Arguments for the script
    parser = ArgumentParser()
    parser.add_argument("pathExperiment", help = "Specify the path to the directory containing txt file needed to be analyzed.")
    parser.add_argument("pathtosample", help = "Specify the sample that needs to have its LV distance analyzed.")
    parser.add_argument("pathtosummary", help = "Specify the path to the summary file containing LV mean distance.")
    parser.add_argument("sampleName", help = "Specify the name of sample directory containing the results of step 1 and step 2." ,type=str)
    parser.add_argument("-inputFraction", help = "Specify if the histogram has to be constructed for the full barcode or for a partial part. Options: full, partial", choices=["full", "partial"], default="full")
    parser.add_argument("-inputLength", help= "If you have chosen partial in the inputFraction, enter the length of barcode you want to construct the histogram for.", type=int, default = 50)
    args = parser.parse_args()

    # Initialize the folders
    if not os.path.exists(os.path.join(args.pathExperiment,"LV_Analysis", f"{args.sampleName}_matrix")):
        os.makedirs(os.path.join(args.pathExperiment,"LV_Analysis", f"{args.sampleName}_matrix"))

    # Change directory to the pathExperiment 
    os.chdir(args.pathExperiment)

    # Assigning the length of substring to determine LV distance 
    Barcode_raw = [] #to store all the Barcodes

    with open(args.pathtosample, "r") as reads:   
        for line in reads:
            Barcode_raw.append(line.split()[0])

    Barcode = np.unique(Barcode_raw)

    del(Barcode_raw)
    gc.collect()
    
    # Make samples from the data 
    seed_value = 42
    random.seed(seed_value)
    Sample_1 = np.random.choice(Barcode, size=5000, replace=False)
    Sample_2 = np.random.choice(Barcode, size=5000, replace=False)
    Sample_3 = np.random.choice(Barcode, size=5000, replace=False)
    del(Barcode)
    gc.collect()

    # Combine Samples into a list 
    Samples = [Sample_1, Sample_2, Sample_3]

    # Calculate the string distance matrix using Levenshtein distance (edit distance)
    print("Beginning to find Levenshtein distance for the {} Barcode for sample {}".format(args.inputFraction, args.sampleName))

    # Initialize the matrix 
    sample_means = []
    matrices = []

    for count, samples in enumerate(Samples):
        csvpath = os.path.join("LV_Analysis", f"{str(args.sampleName)}_matrix", f"sample_{count}.csv")
        mean, matrix1 = analyze_LV(samples, csvpath, count+1)
        sample_means.append(mean)
        matrices.append(matrix1)
    
    if args.inputFraction == "partial": 
        pdf_filename = 'LV_Analysis/' + args.sampleName + '_LV_distance_{}.pdf'.format(args.inputLength)
    else: 
        pdf_filename = 'LV_Analysis/' + args.sampleName + '_LV_distance_full.pdf'

    pdf = PdfPages(pdf_filename)

    x_tick = np.arange(0, max(np.max(matrices[0]), np.max(matrices[1]), np.max(matrices[2])), 4) 

    # Create and save the plots
    with PdfPages(pdf_filename) as pdf:
        create_histogram(matrices[0], 'Sample 1', pdf, x_tick)
        create_histogram(matrices[1], 'Sample 2', pdf, x_tick)
        create_histogram(matrices[2], 'Sample 3', pdf, x_tick)

    print(f"PDF saved as {pdf_filename}")
    
    summary_file = args.pathtosummary
    #Keep all the mean data saved 
    with open(summary_file, "a") as summary:
        summary.write("Sample : {}\n".format(args.sampleName))
        summary.write("Sample 1 mean {}\n".format(sample_means[0]))
        summary.write("Sample 2 mean {}\n".format(sample_means[1]))
        summary.write("Sample 3 mean {}\n\n".format(sample_means[2]))

    print(summary)


