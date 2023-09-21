'''
Note about the script:
This is a python script to submit a sample to Starcode to match for optimal Levenstein distance and merge similar barcodes to reduce redundancy
Output files:
	1. 3 files are created where merged barcodes of varying lengths (50,40,30) are written

command to run this script: python3 starcodeRun.py <path to experiment directory> <path to script directory> <combined - yes/no> <sample name> -d <distance> -thread <Number of threads to use>
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
from argparse import ArgumentParser
import os
import subprocess


#Command line parser
parser = ArgumentParser()
parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory")
parser.add_argument("pathScript", help = "Specify the path to the script containing directory")
parser.add_argument("combined", help = "If combines file of multiple samples has to be run, yes; otherwise, no", choices = ["yes", "no"])
parser.add_argument("sampleName", help = "Specify the name of sample directory containing the text files of barcodes")
parser.add_argument("-d", "--distance", help = "Specify the Levenshtein distance for clustering/merging sequences", type = str, default = 8) 
parser.add_argument("-thread", help = "Number of threads to use for Starcode", default='4',type = str )
args = parser.parse_args()

#The location where starcode is compiled and stored
starcodeLocation = os.path.join(args.pathScript,"Step3_Starcode","starcode","starcode")

threads = args.thread
os.chdir("..")

lengths = [30,40,50]
if(args.combined == "no"):
	outFileDirectory = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'starcode')		
	if not os.path.exists(outFileDirectory):
			os.makedirs(outFileDirectory)

	for i in lengths:
		samplePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'LV_Analysis/{}_Barcode_{}.txt'.format(args.sampleName, i))
		outfilePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'starcode', '{}_Barcode{}_d{}.txt'.format(args.sampleName, i, args.distance))
		starcodeCommand = [starcodeLocation, '-d', args.distance, '-t', threads, '-i', samplePath, "-s", "--seq-id", '-o', outfilePath]
		subprocess.run(starcodeCommand)
else:
	outFileDirectory = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'starcode')
	if not os.path.exists(outFileDirectory):
		os.makedirs(outFileDirectory)
	for i in lengths:
		samplePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, '{}_Barcode_{}.txt'.format(args.sampleName, i))
		outfilePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'starcode', '{}_Barcode{}_d{}.txt'.format(args.sampleName, i, args.distance))
		starcodeCommand = [starcodeLocation, '-d', args.distance, '-t', threads, '-i', samplePath, "-s", "--seq-id", '-o', outfilePath]
		subprocess.run(starcodeCommand)		



