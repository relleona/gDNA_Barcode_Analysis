'''
Note about the script:
It is a Python script to plot Levenstein distances between different lengths of barcodes to identify the right length of barcode and Levenstein distance to use for clustering. 

Output of this script:
    1. 3 txt files containing 30, 40 and 50 starting bases of the barcode
    2. A CSV file containing the full barcode, 30, 40 and 50 starting bases of the barcode
    3. A pdf file containing the LV distance (of the full or partial barcode) histogram (log scale) for 3 samples

Command to run this script: python3 LvHistogram.py <experiment file> <Sample name> -inputFraction <full/partial> -inputLength <length>
'''
import gzip
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
parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory")
parser.add_argument("sampleName", help = "Specify the name of sample directory containing the results of step 1")
parser.add_argument("-inputFraction", help = "Specify if the histogram has to be constructed for the full barcode or for a partial part. Options: full, partial", choices=["full", "partial"], default="full")
parser.add_argument("-inputLength", help= "If you have chosen partial in the inputFraction, enter the length of barcode you want to construct the histogram for.", type=int)

args = parser.parse_args()

experimentDirectory = args.pathExperiment
os.chdir('..')
os.chdir(os.path.join(experimentDirectory, "analyzed", args.sampleName)) #identifying the input file location and moving to it
if not os.path.exists("LV_Analysis"):
    os.makedirs("LV_Analysis")

os.chdir(os.path.join("extractedBarcodeData")) 

# Assigning the length of substring to determine LV distance 
lengths = [30, 40, 50]

if(args.inputFraction == "partial"):
    length = int(args.inputLength)
    frac = 1
else: 
    frac = 0

Barcode_raw = [] #to store all the Barcodes

# Open required file and parse it to create an array of all barcodes

readCountsFile = args.sampleName + "_readCountsOnly.gz"
if not os.path.exists(readCountsFile):
    readCountsFile = args.sampleName + "_readCountsOnly_liberal.gz"

counts = []

with gzip.open(readCountsFile, 'rt') as reads:   
    for line in reads:
        try:
            Barcode_raw.append(line.split()[0])
            counts.append(int(line.split()[1]))
        except:
            print("error",line)


Barcode = np.array(Barcode_raw)

del(Barcode_raw)
gc.collect()
os.chdir("..")

# Save the barcodes of varying lengths separately
Barcode_1 = []
Barcode_2 = []
Barcode_3 = []
if(frac):
    Barcode_frac = []
for i in Barcode:
    Barcode_1.append(i[:lengths[0]])
    Barcode_2.append(i[:lengths[1]])
    Barcode_3.append(i[:lengths[2]])
    if(frac):
        Barcode_frac.append(i[:length])

maps = {"Sequence": np.array(Barcode), "Barcode_" + str(lengths[0]): np.array(Barcode_1), "Barcode_" + str(lengths[1]): np.array(Barcode_2), "Barcode_" + str(lengths[2]): np.array(Barcode_3), "Counts":np.array(counts)}
Barcode_new = pd.DataFrame.from_dict(maps)
Barcode_new.to_csv("LV_Analysis/" + args.sampleName + "_FullBarcode.csv", index = False, sep='\t')
m1 = {"Sequence": Barcode_1, "Counts": counts}
Barcode_1 = pd.DataFrame.from_dict(m1)
Barcode_1.to_csv("LV_Analysis/" + args.sampleName + "_Barcode_" + str(lengths[0]) + ".txt", index = False, header=False, sep='\t')
m2 = {"Sequence": Barcode_2, "Counts": counts}
Barcode_2 = pd.DataFrame.from_dict(m2)
Barcode_2.to_csv("LV_Analysis/" + args.sampleName + "_Barcode_" + str(lengths[1]) + ".txt", index = False, header=False, sep='\t')
m3 = {"Sequence": Barcode_3, "Counts": counts}
Barcode_3 = pd.DataFrame.from_dict(m3)
Barcode_3.to_csv("LV_Analysis/" + args.sampleName + "_Barcode_" + str(lengths[2]) + ".txt", index = False, header=False, sep='\t')

if(frac):
    Barcode = Barcode_frac
    del Barcode_frac

del Barcode_1, Barcode_2, Barcode_3
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
print("Beginning to find Levenshtein distance for the {} Barcode".format(args.inputFraction))

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
if(frac):
    pdf_filename = 'LV_Analysis/' + args.sampleName + '_LV_distance_{}.pdf'.format(args.inputLength)
else:
    pdf_filename = 'LV_Analysis/' + args.sampleName + '_LV_distance_full.pdf'
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



