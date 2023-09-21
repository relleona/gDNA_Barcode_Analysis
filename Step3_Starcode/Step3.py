import os, subprocess
from argparse import ArgumentParser

#Command line parser
parser = ArgumentParser()

parser.add_argument("scriptPath", help = "Specify the path containing the script files.", type = str)
parser.add_argument("experimentPath", help = "Specify the path to the experiment directory", type = str)
parser.add_argument("combinedSample", help = "Specify if the starcode needs to run on combined file of multiple samples(yes); otherwise no", choices=["yes", "no"] )
parser.add_argument("lengthStarcode", help = "Length of barcode to be used for analysis" )
parser.add_argument("distanceStarcode", help = "Specify the Levenshtein distance for clustering/merging sequences in Starcode" )
parser.add_argument("threadsStarcode", help = "Number of threads to use for Starcode" )
parser.add_argument("sampleArrayStarcode", help = "List of samples wwhere barcodes have to be merged. For multiple samples, sample name is Multiple_Samples." )
args = parser.parse_args()

# Running Starcode on the samples
pathStarcodeEnvelope = os.path.join(args.scriptPath, "Step3_Starcode","starcodeEnvelope.py")
commandStarcodeEnvelope = ["python3", pathStarcodeEnvelope, args.experimentPath, args.scriptPath, args.combinedSample,"-length", args.lengthStarcode, "-d",args.distanceStarcode,"-thread", args.threadsStarcode,"-sampleArray", args.sampleArrayStarcode]
subprocess.call(commandStarcodeEnvelope)

#Creating a LV Distance Histogram after merging to verify merging is accurate
pathLvHistogram = os.path.join(args.scriptPath, "Step3_Starcode","LvHistogramStarcode.py")
commandLvHistogramFinal = ["python3", pathStarcodeEnvelope, args.experimentPath, args.sampleArrayStarcode, "--length",  args.lengthStarcode, "-d",args.distanceStarcode]
subprocess.call(commandLvHistogramFinal)




