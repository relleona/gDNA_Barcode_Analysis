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
from utilityFunctions import Barcode_scanner, Separate_samples

parser = ArgumentParser()
parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory")
parser.add_argument("combined", help = "Specify if the starcode needs to run on combined file of multiple samples(yes); otherwise no", choices=["yes", "no"] )
parser.add_argument("sampleName", help = "Sample Name" )
parser.add_argument("-d", help = "LV Distance used for clustering using starcode.", default="8" )
parser.add_argument("-length", help = "Length of Barcode to be used for analysis", default=50 )

args = parser.parse_args()

inputDirectory = args.pathExperiment
sampleName = args.sampleName

os.chdir("..")

if(args.combined == "yes"):
    CombinedFilePath = os.path.join(args.pathExperiment, "analyzed",args.sampleName,"MultipleSamples_FullBarcode.csv")
    Combined_Barcode = pd.read_csv(CombinedFilePath, sep="\t")
    StarcodeInputPath = os.path.join(args.pathExperiment, "analyzed",args.sampleName,"starcode")
    os.chdir(StarcodeInputPath)
    length = ["30", "40", "50"]
    for ind, i in enumerate(length):
        print("Combining", i)
        Combined_Barcode = Barcode_scanner(Combined_Barcode, sampleName, i, str(args.d), ind)
    print("Combined")
    os.chdir("..")
    Combined_Barcode.to_csv("updatedAllBarcode.csv", index=False, sep="\t")
    Separate_samples(Combined_Barcode, 1, args.sampleName,args.length, str(args.d))
else:
    CombinedFilePath = os.path.join(args.pathExperiment, "analyzed",args.sampleName,"LV_Analysis",str(args.sampleName) + "_FullBarcode.csv")
    Combined_Barcode = pd.read_csv(CombinedFilePath, sep="\t")
    StarcodeInputPath = os.path.join(args.pathExperiment, "analyzed",sampleName,"starcode")
    os.chdir(StarcodeInputPath)
    length = ["30", "40","50"]
    for ind, i in enumerate(length):
        Combined_Barcode = Barcode_scanner(Combined_Barcode, args.sampleName, i, str(args.d), ind)
    os.chdir("..")
    Combined_Barcode.to_csv("updatedAllBarcode.csv", index=False)
    Separate_samples(Combined_Barcode, 0, args.sampleName,args.length)




