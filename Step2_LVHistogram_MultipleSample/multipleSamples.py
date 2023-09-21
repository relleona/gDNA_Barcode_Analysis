'''
Note about the script:
It is a Python script to combine multiple twin samples so that they can be combinedly analyzed in Starcode for accurate identification of redundant barcodes

Output of this script:
    1. 3 txt files containing 30, 40 and 50 starting bases of the barcode of multiple samples
    2. A CSV file containing the full barcode, 30, 40 and 50 starting bases of the barcode and sample name

Command to run this script: python3 <path to multipleSamples.py> <Path to experiment directory> <Sample name 1> <Sample name 2> <Sample name 3> 
'''
import os
import pandas as pd
from argparse import ArgumentParser

# Arguments for the script
parser = ArgumentParser()
parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory")
parser.add_argument("sampleNames", nargs='+', help = "Specify the name of sample directory containing the results of step 1 to combine")

args = parser.parse_args()

# Function to extract the barcodes for a given sample
def reader_barcode(sampleName):
    experimentDirectory = args.pathExperiment
    os.chdir(os.path.join(experimentDirectory, "analyzed", sampleName, "LV_Analysis")) #identifying the input file location and moving to it
    
    with open(sampleName + "_FullBarcode.csv", "r") as BarcodeCSV:
        BarcodeCSV = pd.read_csv(BarcodeCSV, sep="\t")
        BarcodeCSV['Sample'] = sampleName
    # os.chdir('../../../..')
    os.chdir(experimentDirectory)
    return BarcodeCSV

# Iterating through the samples to extract the barcodes and combine into a single dataframe
barcodeMultipleSamples = pd.DataFrame()
for sample in args.sampleNames:
    if barcodeMultipleSamples.empty:
        barcodeMultipleSamples = reader_barcode(sample)
    else:
        barcodeMultipleSamples = pd.concat([barcodeMultipleSamples, reader_barcode(sample)])

# creating output directory
os.chdir(os.path.join(args.pathExperiment, "analyzed"))
if(not os.path.exists("Multiple_Samples")):
    os.mkdir("Multiple_Samples")
os.chdir(os.path.join("Multiple_Samples"))

# Saving the various combined barcode files
Barcode_30 = barcodeMultipleSamples[['Barcode_30','Counts']]
Barcode_30.to_csv('Multiple_Samples_Barcode_30.txt', index=False, header=False, sep="\t")
Barcode_40 = barcodeMultipleSamples[['Barcode_40','Counts']]
Barcode_40.to_csv('Multiple_Samples_Barcode_40.txt', index=False, header=False, sep="\t")
Barcode_50 = barcodeMultipleSamples[['Barcode_50','Counts']]
Barcode_50.to_csv('Multiple_Samples_Barcode_50.txt', index=False, header=False, sep="\t")
print(barcodeMultipleSamples.columns)
barcodeMultipleSamples.to_csv("MultipleSamples_FullBarcode.csv",index=False, sep="\t")
