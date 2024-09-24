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
parser.add_argument("-d", help = "Specify the distance used for Levenshtein distance based clustering in starcode.", default="8")
args = parser.parse_args()


def create_histogram(matrix, title, pdf):
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

experimentDirectory = args.experiment
#Edited because there is an error 
if (args.sampleName == "Multiple_Samples"): 
    SampleList = [args.sampleName]
else: 
    SampleList = args.sampleName.split(",")

for sample in SampleList:
    # os.chdir('..')
    os.chdir(os.path.join(experimentDirectory, "analyzed", sample,"starcode")) #identifying the input file location and moving to it

    if not os.path.exists(os.path.join("matrix")):
        os.makedirs(os.path.join("matrix"))

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
    pd.DataFrame(matrix_1).to_csv(os.path.join("matrix", "sample_1.csv"), index=False, header=False)
    sample1_mean = sum(matrix_1) / len(matrix_1)
    del Sample_1
    gc.collect()
    print("Done for Sampling 1")

    # Sample 2
    pairs_2 = itertools.product(Sample_2, Sample_2)
    pairs_2 = [(pair[0], pair[1]) for pair in pairs_2 if pair[0] != pair[1]]
    matrix_2 = [Levenshtein.distance(pair[0], pair[1]) for pair in pairs_2]
    del(pairs_2)
    matrix_2 = np.array(matrix_2)
    pd.DataFrame(matrix_2).to_csv(os.path.join("matrix", "sample_2.csv"), index=False, header=False)
    sample2_mean = sum(matrix_2) / len(matrix_2)
    del Sample_2
    gc.collect()
    print("Done for Sampling 2")

    # Sample 3
    pairs_3 = itertools.product(Sample_3, Sample_3)
    pairs_3 = [(pair[0], pair[1]) for pair in pairs_3 if pair[0] != pair[1]]
    matrix_3 = [Levenshtein.distance(pair[0], pair[1]) for pair in pairs_3]
    del(pairs_3)
    matrix_3 = np.array(matrix_3)
    pd.DataFrame(matrix_3).to_csv(os.path.join("matrix", "sample_3.csv"), index=False, header=False)
    sample3_mean = sum(matrix_3) / len(matrix_3)
    del Sample_3
    gc.collect()
    print("Done for Sampling 3")

    # Save the merged files of plots together
    pdf_filename = sample + '_LV_distance_{}.pdf'.format(str(args.length))
    pdf = PdfPages(pdf_filename)

    x_tick = np.arange(0, max(np.max(matrix_1), np.max(matrix_2), np.max(matrix_3)), 4) 
    # Create and save the plots
    with PdfPages(pdf_filename) as pdf:
        create_histogram(matrix_1, 'Sample 1', pdf)
        create_histogram(matrix_2, 'Sample 2', pdf)
        create_histogram(matrix_3, 'Sample 3', pdf)

    print(f"PDF saved as {pdf_filename}")
    # plt.figure(figsize=(6, 4))
    # n_1, bins_1, patches_1 = plt.hist(matrix_1, np.max(matrix_1), density = True)
    # plt.title('Sample 1')
    # plt.yscale('log')
    # plt.xlabel('LV Distance')
    # plt.xticks(x_tick)
    # plt.ylabel('Fraction')
    # pdf.savefig()
    # plt.close()

    # plt.figure(figsize=(6, 4))
    # n_2, bins_2, patches_2 = plt.hist(matrix_2, np.max(matrix_2), density = True)
    # plt.title('Sample 2')
    # plt.xticks(x_tick)
    # plt.yscale('log')
    # plt.xlabel('LV Distance')
    # plt.ylabel('Fraction')
    # pdf.savefig()
    # plt.close()

    # plt.figure(figsize=(6, 4))
    # n_3, bins_3, patches_3 = plt.hist(matrix_3, np.max(matrix_3), density = True)
    # plt.title('Sample 3')
    # plt.xticks(x_tick)
    # plt.yscale('log')
    # plt.xlabel('LV Distance')
    # plt.ylabel('Fraction')
    # pdf.savefig()
    # plt.close()
    # pdf.close()

    #Keep all the mean data saved 
    summary_file = args.sampleName + "_meandistance.txt"
    with open(summary_file, "w") as summary:
        summary.write("Sample 1 mean {}\n".format(sample1_mean))
        summary.write("Sample 2 mean {}\n".format(sample2_mean))
        summary.write("Sample 3 mean {}\n".format(sample3_mean))

    print(summary)


#Now process Levensthein distances for the final data as well if its in mulriple samples 
if (args.sampleName == "Multiple_Samples"): 

    Separatedsamples = os.path.join(experimentDirectory, "analyzed",args.sampleName,"separated")
    os.chdir(Separatedsamples)

    # File setup - outside the loop
    summary_file = "final_meandistance.txt"

    # Check if the file exists, if not, create it with a header
    if not os.path.exists(summary_file):
        with open(summary_file, "w") as summary:
            summary.write("Summary of Mean Levenshtein Distances\n")
            summary.write("==========================\n\n")

    # Corrected pattern using raw string    
    pattern = re.compile(rf"_final{args.length}_{args.d}\.txt$")

    #Find the file with the pattern
    for root, dirs, files in os.walk(Separatedsamples):
        for file in files:
            if pattern.search(file):
                file_path = os.path.join(root, file)
                print(f"Processing file: {file_path}")
                # Assigning the length of substring to determine LV distance 
                Barcode_raw = [] #to store all the Barcodes

                # Open required file and parse it to create an array of all barcodes
                with open(file_path, 'r') as reads:
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

                base_name = file.split("_final")[0]

                # Make the matrix data and save it here 
                if not os.path.exists(os.path.join(f"{base_name}_matrix")):
                    os.makedirs(os.path.join(f"{base_name}_matrix"))

                # Sample 1
                pairs_1 = itertools.product(Sample_1, Sample_1)
                pairs_1 = [(pair[0], pair[1]) for pair in pairs_1 if pair[0] != pair[1]]
                matrix_1 = [Levenshtein.distance(pair[0], pair[1]) for pair in pairs_1]
                del(pairs_1)
                matrix_1 = np.array(matrix_1)
                pd.DataFrame(matrix_1).to_csv(os.path.join(f"{base_name}_matrix", "sample_1.csv"), index=False, header=False)
                sample1_mean = sum(matrix_1) / len(matrix_1)
                del Sample_1
                gc.collect()
                print(f"Done for Sampling 1 for {base_name}")

                # Sample 2
                pairs_2 = itertools.product(Sample_2, Sample_2)
                pairs_2 = [(pair[0], pair[1]) for pair in pairs_2 if pair[0] != pair[1]]
                matrix_2 = [Levenshtein.distance(pair[0], pair[1]) for pair in pairs_2]
                del(pairs_2)
                matrix_2 = np.array(matrix_2)
                pd.DataFrame(matrix_2).to_csv(os.path.join(f"{base_name}_matrix", "sample_2.csv"), index=False, header=False)
                sample2_mean = sum(matrix_2) / len(matrix_2)
                del Sample_2
                gc.collect()
                print(f"Done for Sampling 2 for {base_name}")

                # Sample 3
                pairs_3 = itertools.product(Sample_3, Sample_3)
                pairs_3 = [(pair[0], pair[1]) for pair in pairs_3 if pair[0] != pair[1]]
                matrix_3 = [Levenshtein.distance(pair[0], pair[1]) for pair in pairs_3]
                del(pairs_3)
                matrix_3 = np.array(matrix_3)
                pd.DataFrame(matrix_3).to_csv(os.path.join(f"{base_name}_matrix", "sample_3.csv"), index=False, header=False)
                sample3_mean = sum(matrix_3) / len(matrix_3)
                del Sample_3
                gc.collect()
                print(f"Done for Sampling 3 for {base_name}")

                # Save the merged files of plots together
                pdf_filename = file + '_LV_distance_{}.pdf'.format(str(args.length))
                pdf = PdfPages(pdf_filename)

                x_tick = np.arange(0, max(np.max(matrix_1), np.max(matrix_2), np.max(matrix_3)), 4) 

                # Create and save the plots
                with PdfPages(pdf_filename) as pdf:
                    create_histogram(matrix_1, 'Sample 1', pdf)
                    create_histogram(matrix_2, 'Sample 2', pdf)
                    create_histogram(matrix_3, 'Sample 3', pdf)

                print(f"PDF saved as {pdf_filename}")

                # plt.figure(figsize=(6, 4))
                # n_1, bins_1, patches_1 = plt.hist(matrix_1, np.max(matrix_1), density = True)
                # plt.title('Sample 1')
                # # plt.yscale('log')
                # plt.xlabel('LV Distance')
                # plt.xticks(x_tick)
                # plt.ylabel('Fraction')
                # pdf.savefig()
                # plt.close()

                # plt.figure(figsize=(6, 4))
                # n_2, bins_2, patches_2 = plt.hist(matrix_2, np.max(matrix_2), density = True)
                # plt.title('Sample 2')
                # plt.xticks(x_tick)
                # # plt.yscale('log')
                # plt.xlabel('LV Distance')
                # plt.ylabel('Fraction')
                # pdf.savefig()
                # plt.close()

                # plt.figure(figsize=(6, 4))
                # n_3, bins_3, patches_3 = plt.hist(matrix_3, np.max(matrix_3), density = True)
                # plt.title('Sample 3')
                # plt.xticks(x_tick)
                # # plt.yscale('log')
                # plt.xlabel('LV Distance')
                # plt.ylabel('Fraction')
                # pdf.savefig()
                # plt.close()
                # pdf.close()

                #Keep all the mean data saved 
                with open(summary_file, "a") as summary:
                    summary.write("Sample : {}\n".format(base_name))
                    summary.write("Sample 1 mean {}\n".format(sample1_mean))
                    summary.write("Sample 2 mean {}\n".format(sample2_mean))
                    summary.write("Sample 3 mean {}\n\n".format(sample3_mean))

                print(f"Finish recording {base_name}")