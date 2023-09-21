''' 
Note about the script:
The python script to parse the sample FASTQ files, merge lanes and save it into gz files.

Required Packages: Biopython, Regex(not standard Python regex), numpy, pandas, os, glob 

The output of this script will be: 
	1.An 'analyzed' folder with a sub folder for each sample
	2.The folder for each sample contains an 'extractedBarcodes' folder which will contain 4 files:
		a. sampleName_counts or sampleName_counts_liberal - contains the barcode sequence, number of barcodes counted and number of UMIs associated
		b. sampleName_UMICounts or sampleName_UMICounts_liberal - contains the barcode sequence along with number of associated UMIs
		c. sampleName_UMIs or sampleName_UMIs_liberal - contains the barcode sequence, its quality score for first instance and the list of associated UMIs
		d. sampleName_readCountsOnly or sampleName_readCountsOnly_liberal - contains the barccode sequence and its count number
	3. If excluded reads is true, more files which contain details about bad samples are also outputted to the same folder. It will contain details of barcodes missing vector sequences, having bad quality score or having four repeated bases contiguously
	4. A summary text file containing details about length of barcode dictionary, number of UMIs and number of bad Barcodes classified into different issues. 

Command to run this file: python3 <path to Envelope.py> <path to Experiment directory> <path to Script> "--staggerFile" <path to Stagger File> "-r", "--checkVector", args.checkVector,"-l", args.barcodeLength, "-Q",args.minPhred ,"-e", args.excludeReads

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'''

import os
import subprocess
import csv
from argparse import ArgumentParser
import glob

#Command line parser 
parser = ArgumentParser()
parser.add_argument("experiment", help = "Specify the path to the experiment directory")
parser.add_argument("pathScript", help = "Specify the path containing the script files.", type = str)
parser.add_argument("--pathStaggerFile", help = "Specify the path to the file containing information on length of stagger for each sample")
parser.add_argument("-r", "--includeReads", help = "If specified, output additional tables with barcodes and only read counts or UMI counts. Otherwise outputs only one table with both counts", action = 'store_true')
parser.add_argument("-checkVector", help = "Option to check vector sequence before or on both sides of the barcode sequence.", default = "both", choices = ["both", "before"]) 
parser.add_argument("-barcodeLength", help = "If check_vector before specified, input here your desired barcode length. Default is 100", default = "100", type = str) 
parser.add_argument("-Q", "--minPhred", help = "Specify the minimum phredscore required to include a readout. Filters reads with more than 5 bases before the barcode with low phredscore.", default = "14", type = str) 
parser.add_argument("-e", "--excludedReads", help = "If specified, output txt.gz files containing reads excluded from the UMI and count files.")
args = parser.parse_args()

# pathExtractionScript is the location to the file that has the script to extract barcode
pathExtractionScript = os.path.join(args.pathScript,"Step1_extractBarcode","parseFastqMain.py" )


#Move to experiment directory
os.chdir(args.experiment)

# Including Stagger details for each sample
if args.pathStaggerFile is not None:
	samples = []
	staggers = []
	with open(args.pathStaggerFile, 'r') as file:
		tmp = csv.reader(file)
		for line in tmp:
			samples.append(line[0])
			staggers.append(str(line[1])) 
else:
	samples = [os.path.basename(i) for i in glob.glob("raw/*")]
	staggers = ['0'] * len(samples)

print(samples)

#Format additional arguments
additionalArguments = ["-checkVector", args.checkVector, "--minPhred", args.minPhred]
if args.includeReads:
	additionalArguments.extend(["--includeReads"])
if args.checkVector == "before":
	additionalArguments.extend(["-barcodeLength", args.barcodeLength])
if args.excludedReads == "True":
	additionalArguments.append(["--excludedReads"])

#Extract barcode, UMI and read counts for each sample file
#Extract barcode for each sample file
for index, i in enumerate(samples):
	command = ["python3", pathExtractionScript, args.experiment, i, "-s", staggers[index]] + additionalArguments
	print(command)
	subprocess.run(command)

