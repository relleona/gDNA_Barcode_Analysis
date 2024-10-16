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
# parser.add_argument("combined", help = "If combines file of multiple samples has to be run, yes; otherwise, no", choices = ["yes", "no"])
parser.add_argument("sampleName", help = "Specify the name of sample directory containing the text files of barcodes")
parser.add_argument("-d", "--distance", help = "Specify the Levenshtein distance for clustering/merging sequences", type = str, default = "8") 
parser.add_argument("-thread", help = "Number of threads to use for Starcode", default='4',type = str )
parser.add_argument("-length", help = "Length of barcode to be used for analysis", default="50",type = str )
parser.add_argument("-Fraction",help = "Specify if its running starcode for the full barcode or for a partial part. Options: full, partial", choices=["full", "partial"], default="full")


args = parser.parse_args()

#The location where starcode is compiled and stored
starcodeLocation = os.path.join(args.pathScript,"Step3_Starcode","starcode","starcode")

threads = args.thread
os.chdir("..")

# These are the lengths that are normally used or initiated
lengths = ["30","40","50"]

# These are the input length  that is not within the measurements mentioned above
if int(args.length) not in lengths: 
	if args.Fraction =="full": 
		lengths.append("full")
	else:
		lengths.append(str(args.length))

# The outputfile of the data 
outFileDirectory = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'starcode')		
if not os.path.exists(outFileDirectory):
	os.makedirs(outFileDirectory)

# Starcode should be done on every one of the different lengths 
for i in lengths:
	# If the length that is used is the entire length of the barcode, it will be titled full instead of the length 
	if i =="full": 
		samplePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'LV_Analysis/{}_Barcode_full.txt'.format(args.sampleName))
		outfilePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'starcode', '{}_Barcodefull_d{}.txt'.format(args.sampleName, args.distance))

	else: 
		samplePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'LV_Analysis/{}_Barcode_{}.txt'.format(args.sampleName, i))
		outfilePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'starcode', '{}_Barcode{}_d{}.txt'.format(args.sampleName, i, args.distance))

	# Command includes location, distance to be collapsed, threads used, path to barcode, path to output 
	starcodeCommand = [starcodeLocation, '-d', args.distance, '-t', threads, '-i', samplePath, "-s", "--seq-id", '-o', outfilePath]
	subprocess.run(starcodeCommand)







# else:
# 	outFileDirectory = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'starcode')
# 	if not os.path.exists(outFileDirectory):
# 		os.makedirs(outFileDirectory)
	
# 	for i in lengths:
# 		samplePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, '{}_Barcode_{}.txt'.format(args.sampleName, i))
# 		outfilePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'starcode', '{}_Barcode{}_d{}.txt'.format(args.sampleName, i, args.distance))
# 		starcodeCommand = [starcodeLocation, '-d', args.distance, '-t', threads, '-i', samplePath, "-s", "--seq-id", '-o', outfilePath]
# 		subprocess.run(starcodeCommand)		



# for i in lengths:
#     samplePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'LV_Analysis/{}_Barcode_{}.txt'.format(args.sampleName, i))
#     outfilePath = os.path.join(args.pathExperiment, "analyzed", args.sampleName, 'starcode', '{}_Barcode{}_d{}.txt'.format(args.sampleName, i, args.distance))
    
#     print(f"Sample path: {samplePath}")
#     print(f"Output path: {outfilePath}")
    
#     if not os.path.exists(samplePath):
#         print(f"Error: Input file not found: {samplePath}")
#         continue
    
#     starcodeCommand = [starcodeLocation, '-d', args.distance, '-t', threads, '-i', samplePath, "--seq-id", '-o', outfilePath]
    
#     print(f"Running command: {' '.join(starcodeCommand)}")
    
#     try:
#         result = subprocess.run(starcodeCommand, check=True, capture_output=True, text=True)
#         print(f"Starcode output: {result.stdout}")
#     except subprocess.CalledProcessError as e:
#         print(f"Error running starcode: {e}")
#         print(f"Starcode stderr: {e.stderr}")
#     except Exception as e:
#         print(f"Unexpected error: {e}")

# print("Script completed.")