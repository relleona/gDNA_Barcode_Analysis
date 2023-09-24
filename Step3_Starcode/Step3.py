#!/home/keerthana/miniconda3/envs/Barcode_extraction/bin/python
import os, subprocess
from argparse import ArgumentParser
import numpy

#Command line parser
parser = ArgumentParser()

parser.add_argument("pathScript", help = "Specify the path containing the script files.", type = str)
parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory", type = str)
parser.add_argument("combinedSample", help = "Specify if the starcode needs to run on combined file of multiple samples(yes); otherwise no", choices=["yes", "no"] )
parser.add_argument("lengthStarcode", help = "Length of barcode to be used for analysis" )
parser.add_argument("distanceStarcode", help = "Specify the Levenshtein distance for clustering/merging sequences in Starcode" )
parser.add_argument("threadsStarcode", help = "Number of threads to use for Starcode" )
parser.add_argument("sampleArrayStarcode", help = "List of samples where barcodes have to be merged. For multiple samples, sample name is Multiple_Samples." )
args = parser.parse_args()

# Running Starcode on the samples
pathStarcodeEnvelope = os.path.join(args.pathScript, "Step3_Starcode","starcodeEnvelope.py")
commandStarcodeEnvelope = ["python3", pathStarcodeEnvelope, args.pathExperiment, args.pathScript, args.combinedSample,"-length", args.lengthStarcode, "-d",args.distanceStarcode,"-thread", args.threadsStarcode,"-sampleArray", args.sampleArrayStarcode]
print(commandStarcodeEnvelope)
subprocess.call(commandStarcodeEnvelope)

#Creating a LV Distance Histogram after merging to verify merging is accurate
pathLvHistogram = os.path.join(args.pathScript, "Step3_Starcode","LvHistogramStarcode.py")
commandLvHistogramFinal = ["python3", pathLvHistogram, args.pathExperiment, args.sampleArrayStarcode, "--length",  args.lengthStarcode, "-d",args.distanceStarcode]
subprocess.call(commandLvHistogramFinal)




