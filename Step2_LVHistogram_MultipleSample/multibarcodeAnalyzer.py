import gzip
import regex as re
from argparse import ArgumentParser
import os
import numpy as np
import random
import itertools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gc
from matplotlib.backends.backend_pdf import PdfPages
# from Step2functions import analyze_LV, create_histogram

# Arguments for the script
parser = ArgumentParser()
parser.add_argument("experimentPath", help = "Specify the path to the experiment directory", type = str)
parser.add_argument("sampleNames", help = "Array of samples within the Experiment data.", type = str)
parser.add_argument("-inputFraction", help = "Specify if the histogram has to be constructed for the full barcode or for a partial part. Options: full, partial", choices=["full", "partial"], default="full")
parser.add_argument("-inputLength", help= "If you have chosen partial in the inputFraction, enter the length of barcode you want to construct the histogram for.", type=int, default = 50)

args = parser.parse_args()

# #Move into current directory
# os.chdir('..')

# Using the inputs to determine if we want partial or full and outputs txt and csv files for further processing. 
if(args.inputFraction == "partial"):
    inputlength = int(args.inputLength)
    frac = 1
else: 
    frac = 0

# Initiate an empty dataframe to combine all samples into Multiple_Samples 
barcodeMultipleSamples = pd.DataFrame()

# Assigning different lengths of substring to later determine LV distance 
lengths = [30, 40, 50]

# Separate the samples by a comma and make it into a list
SampleList = args.sampleNames.split(",")

for sample in SampleList:
    #Move into analyzed directory
    os.chdir(os.path.join(args.experimentPath, "analyzed", sample)) #identifying the input file location and moving to it

    #Make LV_Analysis and matrix folder
    if not os.path.exists("LV_Analysis"):
        os.makedirs("LV_Analysis")

    ## Extract the barcode from the extracted barcode data 
    #Move to extractdedBarcodeData directory
    os.chdir(os.path.join("extractedBarcodeData")) 

    Barcode_raw = [] #to store all the Barcodes
    counts = [] # to store the counts

    # Open required file and parse it to create an array of all barcodes
    readCountsFile = sample + "_readCountsOnly.gz" #both
    if not os.path.exists(readCountsFile):
        readCountsFile = sample + "_readCountsOnly_liberal.gz" #before

    #Open the file and append to the data 
    with gzip.open(readCountsFile, 'rt') as reads:   
        for line in reads:
            try:
                Barcode_raw.append(line.split()[0])
                counts.append(int(line.split()[1]))
            except:
                print("error",line)

    # Make it into a numpy array 
    Barcode = np.array(Barcode_raw)

    # Save some data 
    del(Barcode_raw)
    gc.collect()

    # Moves "up" one level in the directory structure.
    os.chdir("..")

    # Save the barcodes of varying lengths separately
    Barcode_1 = []
    Barcode_2 = []
    Barcode_3 = []
    if(frac):
        Barcode_frac = []

    # Read every line and make the lengths as mentioned 
    for i in Barcode:
        Barcode_1.append(i[:lengths[0]])
        Barcode_2.append(i[:lengths[1]])
        Barcode_3.append(i[:lengths[2]])
        if(frac):
            Barcode_frac.append(i[:inputlength])

    
    # Make an array of the Sample Names only that is the same length as the total barcodes in the Barcode variable
    sample_array = np.full(len(Barcode_1), sample)

    # if the inputLength is not within the list of length, then there should be another column that includes that fraction
    if frac == 1 and inputlength not in lengths:
        maps = {"Sequence": np.array(Barcode), "Barcode_" + str(lengths[0]): np.array(Barcode_1), "Barcode_" + str(lengths[1]): np.array(Barcode_2), "Barcode_" + str(lengths[2]): np.array(Barcode_3), "Barcode_" + str(inputlength): np.array(Barcode_frac), "Counts":np.array(counts), "Sample":sample_array}
                # Saving the barcode for starcode processing later on
        m5 = {"Sequence": Barcode_frac, "Counts": counts}
        Barcode_5 = pd.DataFrame.from_dict(m5)
        Barcode_5.to_csv("LV_Analysis/" + sample + "_Barcode_" + str(inputlength) + ".txt", index = False, header=False, sep='\t')
    else: 
        maps = {"Sequence": np.array(Barcode), "Barcode_" + str(lengths[0]): np.array(Barcode_1), "Barcode_" + str(lengths[1]): np.array(Barcode_2), "Barcode_" + str(lengths[2]): np.array(Barcode_3), "Counts":np.array(counts), "Sample":sample_array}
    
    # Saving all the data for later/further processing 
    m1 = {"Sequence": Barcode_1, "Counts": counts}
    Barcode_1 = pd.DataFrame.from_dict(m1)
    Barcode_1.to_csv("LV_Analysis/" + sample + "_Barcode_" + str(lengths[0]) + ".txt", index = False, header=False, sep='\t')
    m2 = {"Sequence": Barcode_2, "Counts": counts}
    Barcode_2 = pd.DataFrame.from_dict(m2)
    Barcode_2.to_csv("LV_Analysis/" + sample + "_Barcode_" + str(lengths[1]) + ".txt", index = False, header=False, sep='\t')
    m3 = {"Sequence": Barcode_3, "Counts": counts}
    Barcode_3 = pd.DataFrame.from_dict(m3)
    Barcode_3.to_csv("LV_Analysis/" + sample + "_Barcode_" + str(lengths[2]) + ".txt", index = False, header=False, sep='\t')
    m4 = {"Sequence": Barcode, "Counts": counts}
    Barcode_4 = pd.DataFrame.from_dict(m4)
    Barcode_4.to_csv("LV_Analysis/" + sample + "_Barcode_full.txt", index = False, header=False, sep='\t')

    # Make a full barcode csv including all the different lengths 
    Barcode_new = pd.DataFrame.from_dict(maps)
    Barcode_new.to_csv("LV_Analysis/" + sample + "_AllBarcode.csv", index = False, sep='\t')

    # Iterating through the samples to extract the barcodes and combine into a single dataframe
    if barcodeMultipleSamples.empty:
        barcodeMultipleSamples = Barcode_new
    else:
        barcodeMultipleSamples = pd.concat([barcodeMultipleSamples, Barcode_new], ignore_index=True)

# creating output directory
os.chdir(os.path.join(args.experimentPath, "analyzed"))
if(not os.path.exists(os.path.join("Multiple_Samples", "LV_Analysis"))):
    os.mkdir(os.path.join("Multiple_Samples", "LV_Analysis"))
os.chdir(os.path.join("Multiple_Samples", "LV_Analysis"))

# Outside of the loop 
# Saving the various combined barcode files
Barcode_30 = barcodeMultipleSamples[['Barcode_30','Counts']]
Barcode_30.to_csv('Multiple_Samples_Barcode_30.txt', index=False, header=False, sep="\t")
Barcode_40 = barcodeMultipleSamples[['Barcode_40','Counts']]
Barcode_40.to_csv('Multiple_Samples_Barcode_40.txt', index=False, header=False, sep="\t")
Barcode_50 = barcodeMultipleSamples[['Barcode_50','Counts']]
Barcode_50.to_csv('Multiple_Samples_Barcode_50.txt', index=False, header=False, sep="\t")

if frac == 1 and inputlength not in lengths:
    string = "Barcode_" + str(args.inputLength)
    Barcode_frac = barcodeMultipleSamples[[string,'Counts']]
    Barcode_frac.to_csv("Multiple_Samples_" + "Barcode_" + str(args.inputLength) + ".txt", index=False, header=False, sep="\t")

Barcode_100 = barcodeMultipleSamples[['Sequence','Counts']]
Barcode_100.to_csv("Multiple_Samples_Barcode_full.txt",index=False, header=False, sep="\t")

print(barcodeMultipleSamples.columns)

barcodeMultipleSamples.to_csv("Multiple_Samples_AllBarcode.csv",index=False, sep="\t") 