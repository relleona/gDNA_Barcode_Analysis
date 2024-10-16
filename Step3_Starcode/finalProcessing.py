'''
Note about the script:
This is a python script to update the barcode csv with starcoode merged barcodes
Output files:
	1. The barcode and its substrings are updated in the csv file and that is saved as a csv file
    2. Separate textfiles are created to store the specified length of barcode for each sample separately if they were combined

command to run this script: python3 finalProcessing.py <experiment name> <combined - yes/no> <sample name> -d <distance used for merging in Starcode> -length <Length of barcode for final analysis>
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

import os
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from utilityFunctionsv2 import Barcode_scanner, process_barcodes

parser = ArgumentParser()
parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory")
parser.add_argument("combined", help = "Specify if the starcode needs to run on combined file of multiple samples(yes); otherwise no", choices=["yes", "no"] )
parser.add_argument("sampleName", help = "Sample Name" )
parser.add_argument("-d", help = "LV Distance used for clustering using starcode.", default="8" )
parser.add_argument("-length", help = "Length of Barcode to be used for analysis", default=50 )
parser.add_argument("--f", "-Fraction",help = "Specify if its running starcode for the full barcode or for a partial part. Options: full, partial", choices=["full", "partial"], default="full")


args = parser.parse_args()

inputDirectory = args.pathExperiment
sampleName = args.sampleName

os.chdir("..")

lengths = ["30","40","50"]

if str(args.length) not in lengths: 
    if str(args.f) == "full": 
        # Will indicate that we are using the whole sequence. 
        # It will also to be inputted into the name of the post starcode txt to be read. 
        lengths.insert(0, "full")
    else: 
        # Will include the lengths that are not listed above as shown is AllBarcode.csv
	    lengths.append(str(args.length))

else: 
    # Will include first column that are not listed above as shown in AllBarcode.csv so that the column can be indexed properly
    lengths.insert(0, "sequence")

# Initialized all the paths needed to read the inputs
CombinedFilePath = os.path.join(args.pathExperiment, "analyzed",sampleName,"LV_Analysis",str(sampleName) + "_AllBarcode.csv")
Combined_Barcode = pd.read_csv(CombinedFilePath, sep="\t")
StarcodeInputPath = os.path.join(args.pathExperiment, "analyzed",sampleName,"starcode")
os.chdir(StarcodeInputPath)

if(args.combined == "yes"):

    # Make a new file that indicates all the new barcodes
    for ind, i in enumerate(lengths):
        print("Combining", i)
        Combined_Barcode = Barcode_scanner(Combined_Barcode, sampleName, i, str(args.d), ind)

    print("Combined")
    os.chdir("..")
    Combined_Barcode.to_csv("updatedAllBarcode.csv", index=False, sep="\t")

    # Make a file to store the separated samples, within Multipe_Samples 
    seperatesamplepath = os.path.join(args.pathExperiment, "analyzed",sampleName,"separated")
    if not os.path.exists(seperatesamplepath):
	    os.makedirs(seperatesamplepath)
    
    # Go into the newly made separated folders  
    os.chdir(seperatesamplepath)
    process_barcodes(Combined_Barcode,lengths.index(str(args.length)), sampleName,args.length, str(args.d))

else:
    # CombinedFilePath = os.path.join(args.pathExperiment, "analyzed",args.sampleName,"LV_Analysis",str(args.sampleName) + "_FullBarcode.csv")
    # Combined_Barcode = pd.read_csv(CombinedFilePath, sep="\t")
    # StarcodeInputPath = os.path.join(args.pathExperiment, "analyzed",sampleName,"starcode")
    # os.chdir(StarcodeInputPath)

    for ind, i in enumerate(lengths):
        Combined_Barcode = Barcode_scanner(Combined_Barcode, sampleName, i, str(args.d), ind)
    os.chdir("..")

    Combined_Barcode.to_csv("updatedAllBarcode.csv", index=False)

    # process_barcodes(Combined_Barcode,lengths.index(str(args.length)), sampleName,args.length, str(args.d))




