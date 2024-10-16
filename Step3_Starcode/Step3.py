#!/home/keerthana/miniconda3/envs/Barcode_extraction/bin/python
import os, subprocess
from argparse import ArgumentParser
import numpy
from datetime import datetime
import regex as re


#Command line parser
parser = ArgumentParser()

parser.add_argument("pathScript", help = "Specify the path containing the script files.", type = str)
parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory", type = str)
parser.add_argument("combinedSample", help = "Specify if the starcode needs to run on combined file of multiple samples(yes); otherwise no", choices=["yes", "no"] )
parser.add_argument("lengthStarcode", help = "Length of barcode to be used for analysis" )
parser.add_argument("distanceStarcode", help = "Specify the Levenshtein distance for clustering/merging sequences in Starcode" )
parser.add_argument("threadsStarcode", help = "Number of threads to use for Starcode" )
parser.add_argument("sampleArrayStarcode", help = "List of samples where barcodes have to be merged. For multiple samples, sample name is Multiple_Samples." )
parser.add_argument("Fraction",help = "Specify if its running starcode for the full barcode or for a partial part. Options: full, partial", choices=["full", "partial"], default="full")

args = parser.parse_args()

# # Running Starcode on the samples
# pathStarcodeEnvelope = os.path.join(args.pathScript, "Step3_Starcode","starcodeEnvelope.py")
# commandStarcodeEnvelope = ["python3", pathStarcodeEnvelope, args.pathExperiment, args.pathScript, args.combinedSample,"-length", args.lengthStarcode, "-d",args.distanceStarcode,"-thread", args.threadsStarcode,"-sampleArray", args.sampleArrayStarcode]
# print(commandStarcodeEnvelope)
# subprocess.call(commandStarcodeEnvelope)

# parser = ArgumentParser()
# parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory")
# parser.add_argument("pathScript", help = "Specify the path to the script containing directory")
# parser.add_argument("combined", help = "Specify if the starcode needs to run on combined file of multiple samples(yes); otherwise no", choices=["yes", "no"] )
# parser.add_argument("-length", help = "Length of barcode to be used for analysis", default=50,type = str )
# parser.add_argument("-d","--distance", help = "Specify the Levenshtein distance for clustering/merging sequences in Starcode", type = str, default = 8)
# parser.add_argument("-thread", help = "Number of threads to use for Starcode", default=4,type = str )
# parser.add_argument("-sampleArray", help = "Number of threads to use for Starcode", default=4,type = str )
# args = parser.parse_args()


# Location for the starcode script
starcodeScriptPath = os.path.join(args.pathScript,"Step3_Starcode","starcodeRun.py" )
# Location of the final processing script to update the csv file with new barcodes after merging by Stacode
processingScriptPath = os.path.join(args.pathScript,"Step3_Starcode","finalProcessing.py")
# Location to LV Distance Histogram script
pathLvHistogram = os.path.join(args.pathScript,"LVHistogram.py")

# Move to the Experiment folder 
os.chdir(args.pathExperiment)

# Create txt file in the Experiment folder 
summary_file = "LVmeandistance.txt"

# Generate a timestamp for the output filename
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

# Indicate that it is pre_starcode 
with open(summary_file, "a") as summary:
    summary.write(f"\nSummary of Mean Levenshtein Distances {timestamp}\n")
    summary.write("==========================\n\n")

    summary.write(f"Poststarcode {timestamp}\n")
    summary.write("------------------------------------\n\n")

# Save the path to this summary file 
pathtosummary = os.path.join(args.pathExperiment,summary_file)

if(args.combinedSample == "no"):
	samples = args.sampleArrayStarcode.split(",")

	# Run starcode and processing on all sample files with same settings
	# Iterate to the sample files that are separated
	for index, sample in enumerate(samples):

		# Run starcode 
		command = ["python3", starcodeScriptPath, args.pathExperiment, args.pathScript, sample, "-d", str(args.distanceStarcode), "-thread", str(args.threadsStarcode), "-length", str(args.lengthStarcode), "-Fraction", str(args.Fraction)]
		print(command)
		subprocess.call(command)

		# Run the finalProcessing command: include dataset with the new collapesed barcodes called "updatedAllBarcode.csv" 
		processing = ["python3", processingScriptPath, args.pathExperiment,"no", sample,"-d", str(args.distanceStarcode), "-length",str(args.lengthStarcode), "--f", str(args.Fraction)]
		print(processing)
		subprocess.call(processing)

		# Path to the sample that we want to do an LV analysis on 
		if args.Fraction == "partial" :
			pathtosample = os.path.join(args.pathExperiment, "analyzed", sample, "starcode", "{}_Barcode{}_d{}.txt".format(sample, args.lengthStarcode, args.distanceStarcode))
		else :
			pathtosample = os.path.join(args.pathExperiment, "analyzed", sample, "starcode", "{}_Barcodefull_d{}.txt".format(sample, args.distanceStarcode))
		
		# Path to the folder containing the sample
		pathtofolder = os.path.join(args.pathExperiment, "analyzed", sample, "starcode")

		# Run LV Histogram on the newly made files using starcode  
		commandLvHistogram = ["python3", pathLvHistogram, pathtofolder, pathtosample, pathtosummary, sample, "-inputFraction", args.Fraction, "-inputLength",args.lengthStarcode]
		print(commandLvHistogram)
		subprocess.call(commandLvHistogram)



    # #Creating a LV Distance Histogram after merging to verify merging is accurate
    # commandLvHistogramFinal = ["python3", pathLvHistogram, args.pathExperiment, args.sampleArrayStarcode, "--length",  args.lengthStarcode, "-d",args.distanceStarcode]
    # subprocess.call(commandLvHistogramFinal)

else: 
	# Run starcode 
	command = ["python3", starcodeScriptPath, args.pathExperiment, args.pathScript, "Multiple_Samples","-d", str(args.distanceStarcode), "-thread", str(args.threadsStarcode), "-length",str(args.lengthStarcode),  "-Fraction", args.Fraction ]
	print(command)
	subprocess.call(command)

	# Run the finalProcessing command : include dataset with the new collapesed barcodes called "updatedAllBarcode.csv", the separated samples folder and samples separated with new collapsed barcodes
	processing = ["python3", processingScriptPath, args.pathExperiment, "yes", "Multiple_Samples","-d", str(args.distanceStarcode), "-length",str(args.lengthStarcode), "--f", str(args.Fraction)]
	print(processing)
	subprocess.call(processing)

	# Get the LV distance of the Multiple_Samples data 
	if args.Fraction == "partial" :
		pathtosample = os.path.join(args.pathExperiment, "analyzed", "Multiple_Samples", "starcode", "Multiple_Samples_Barcode{}_d{}.txt".format(args.lengthStarcode, args.distanceStarcode))
	else :
		pathtosample = os.path.join(args.pathExperiment, "analyzed", "Multiple_Samples", "starcode", "Multiple_Samples_Barcodefull_d{}.txt".format(args.distanceStarcode))

	# Path to the folder containing the sample	
	pathtofolder = os.path.join(args.pathExperiment, "analyzed", "Multiple_Samples", "starcode")

	# Run LVHistogram 
	commandLvHistogram = ["python3", pathLvHistogram, pathtofolder, pathtosample, pathtosummary, "Multiple_Samples", "-inputFraction", args.Fraction, "-inputLength",args.lengthStarcode]
	print(commandLvHistogram)
	subprocess.call(commandLvHistogram)

	# Get the LV distance of each of the separated samples in the Multiple_Samples data 
	Separatedsamples = os.path.join(args.pathExperiment, "analyzed", "Multiple_Samples" ,"separated")
	os.chdir(Separatedsamples)
	for root, dirs, files in os.walk(Separatedsamples):
		for file in files:
			# Corrected pattern using raw string
			pattern = rf"_final{args.lengthStarcode}_{args.distanceStarcode}\.txt$"
			if re.search(pattern, file):
				file_path = os.path.join(root, file)
				print(f"Processing file: {file_path}")

				# Get the name of the file being processed 
				pattern = f"_final{args.lengthStarcode}_{args.distanceStarcode}.txt"
				# Using string split
				extracted = file.split(pattern)[0]

				# Run LVHistogram 	
				commandLvHistogram = ["python3", pathLvHistogram, Separatedsamples, file_path, pathtosummary, extracted , "-inputFraction", args.Fraction, "-inputLength",args.lengthStarcode]
				print(commandLvHistogram)
				subprocess.call(commandLvHistogram)

	# Creating a LV Distance Histogram after merging to verify merging is accurate
    # commandLvHistogramFinal = ["python3", pathLvHistogram, args.pathExperiment, "Multiple_Samples", "--length",  args.lengthStarcode, "-d",args.distanceStarcode]
    # subprocess.call(commandLvHistogramFinal)





	