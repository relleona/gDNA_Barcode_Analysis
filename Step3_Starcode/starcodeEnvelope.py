'''
Note about the script:
This is a python script to submit all the samples to Starcode to match for optimal Levenstein distance and merge similar barcodes to reduce redundancy
Output files:
	1. For each sample, 3 files are created where merged barcodes of varying lengths (50,40,30) are written
	2. The barcode and its substrings are updated in the csv file and that is saved as a new csv file

command to run this script: python3 <path to starcodeEnvelope.py> <Path to Experiment directory> <Path to script directory> <combined-yes/no> -length <Length of Barcode to use for final analysis> -d <LV_distance for Starcode> -thread <Number of threads to use>  -sampleName <sampleName1>   <sampleName2> 
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

from argparse import ArgumentParser
import os
import glob
import subprocess

parser = ArgumentParser()
parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory")
parser.add_argument("pathScript", help = "Specify the path to the script containing directory")
parser.add_argument("combined", help = "Specify if the starcode needs to run on combined file of multiple samples(yes); otherwise no", choices=["yes", "no"] )
parser.add_argument("-length", help = "Length of barcode to be used for analysis", default=50,type = str )
parser.add_argument("-d","--distance", help = "Specify the Levenshtein distance for clustering/merging sequences in Starcode", type = str, default = 8)
parser.add_argument("-thread", help = "Number of threads to use for Starcode", default=4,type = str )
parser.add_argument("-sampleArray", help = "Number of threads to use for Starcode", default=4,type = str )
args = parser.parse_args()

#  Location for the starcode script
starcodeScriptPath = os.path.join(args.pathScript,"Step3_Starcode","starcodeRun.py" )
# Location of the final processing script to update the csv file with new barcodes after merging by Stacode
processingScriptPath = os.path.join(args.pathScript,"Step3_Starcode","finalProcessing.py")

if(args.combined == "yes"):
	command = ["python3", starcodeScriptPath, args.pathExperiment, args.pathScript, "yes", "Multiple_Samples","-d", str(args.distance), "-thread", str(args.thread)]
	print(command)
	subprocess.call(command)
	processing = ["python3", processingScriptPath, args.pathExperiment, "yes", "Multiple_Samples","-d", str(args.distance), "-length",str(args.length) ]
	subprocess.call(processing)

else: 
		samples = args.sampleArray.split(",")
		#Run starcode on all sample files with same settings
		for index, sample in enumerate(samples):
			command = ["python3", starcodeScriptPath, args.pathExperiment, args.pathScript,"no", sample, "-d", str(args.distance), "-thread", str(args.thread)]
			print(command)
			subprocess.call(command)
			processing = ["python3", processingScriptPath, args.pathExperiment,"no", sample,"-d", str(args.distance), "-length",str(args.length)]
			subprocess.call(processing)


			
