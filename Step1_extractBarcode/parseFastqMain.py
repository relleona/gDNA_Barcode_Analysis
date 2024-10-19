'''
Note about the script:
The python script to parse the barcodes from single end FASTQ files off NextSeq. 
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

Command to run this file: python3 <path to parseFastqMain.py> <path to Experiment> <Sample Name> -s <Stagger Length> -r --checkVector <both/before> --minPhred <min Quality>

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'''
# Importing necessary packages and functions

#!/opt/anaconda3/envs/Barcode_extraction/bin/python

# import sys
# print(sys.path)
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import regex as re
from collections import Counter
from argparse import ArgumentParser
import os, glob
import numpy as np
from extractionFunctions import parseBarcode_both,parseBarcode_before,\
													 writeOutFileBarcodeCounts,writeOutFileBadSeqRecord,\
													 writeOutFileBarcodeReadCounts,writeOutFileUMIs,\
													 count_read_UMI, writeOutFileBarcodeUMICounts


def check_file_created(filename):
    if os.path.exists(filename):
        print(f"File created successfully: {filename}")
    else:
        print(f"Failed to create file: {filename}")

#Command line parser
parser = ArgumentParser()
parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory")
parser.add_argument("sampleName", help = "Specify the name of sample directory containing the fastq.gz files")
parser.add_argument("-o", "--outFilePrefix", help = "Specify the output file prefix for table of barcodes and UMIs. If none specified, will use sampleDirectory") 
parser.add_argument("-s", "--stagger", help = "Specify the length of the stagger.", type=int, default = 0)
parser.add_argument("-r", "--includeReads", help = "If specified, output additional tables with barcodes and only read counts or UMI counts. Otherwise outputs only one table with both counts", action = 'store_true')
parser.add_argument("-checkVector", help = "Option to check vector sequence before or on both sides of the barcode sequence.", default = "both", choices = ["both", "before"]) 
parser.add_argument("-barcodeLength", help = "If checkVector before specified, input here your desired barcode length.", type = int) 
parser.add_argument("-Q", "--minPhred", help = "Specify the minimum phredscore required to include a readout. Filters reads with more than 5 bases before the barcode with low phredscore.", default = 14, type = int) 
parser.add_argument("-a", "--asciioffset", help = "If PhredScore has letters, ascii offset will be 33, otherwise it will be 64. Most recent version of Illumina uses Phred Score offset of 33. ", default = 33, type = int)
parser.add_argument("-e", "--excludeReads", help = "If specified, output txt.gz files containing reads excluded from the count files.", action = 'store_true')
args = parser.parse_args()

experimentDirectory = args.pathExperiment

#Making the directory for output files.
outFileDirectory = os.path.join(experimentDirectory, "analyzed", args.sampleName, 'extractedBarcodeData') 
if not os.path.exists(outFileDirectory):
    os.makedirs(outFileDirectory)

if args.outFilePrefix is not None:
    outFilePrefix = args.outFilePrefix
else:
    outFilePrefix = args.sampleName

if args.checkVector == 'both':
	outFileUMI = outFilePrefix + "_UMI.gz"
	outFileCounts = outFilePrefix + "_counts.gz"
	outFileUMICounts = outFilePrefix + "_UMICountsOnly.gz"
	outFileReadCounts = outFilePrefix + "_readCountsOnly.gz"
elif args.checkVector == 'before':
	outFileUMI = outFilePrefix + "_Index_liberal.gz"
	outFileCounts = outFilePrefix + "_counts_liberal.gz"
	outFileUMICounts = outFilePrefix + "_UMICountsOnly_liberal.gz"
	outFileReadCounts = outFilePrefix + "_readCountsOnly_liberal.gz"

outFileMissingBeforeBarcode = outFilePrefix + "_missingBeforeBarcode.gz"
outFileMissingAfterBarcode = outFilePrefix + "_missingAfterBarcode.gz"
outFileBadLength = outFilePrefix + "_badLength.gz"
outFileBadPhred = outFilePrefix + "_badPhred.gz"
staggerLength = args.stagger
minPhred = int(args.minPhred)
print(minPhred)
asciioffset = int(args.asciioffset)
print(asciioffset)
barcodeLength = int(args.barcodeLength)
print(barcodeLength)

# Moving to sample directory 
os.chdir(os.path.join(experimentDirectory, "raw", args.sampleName))

# Printing an update that the parsing for this sample has begun
print("Parsing sample {}".format(args.sampleName))

inFileNames = glob.glob("*fastq*")

#Filter the barcode
if args.checkVector == "both":
	barcode_dict, missingBeforeBarcode, missingAfterBarcode, badQscore, badLength, badBarcode, tot_reads = parseBarcode_both(inFileNames, args.stagger, barcodeLength, minPhred, asciioffset)
elif args.checkVector == "before":
	barcode_dict, missingBeforeBarcode, badQscore, badLength, badBarcode, tot_reads = parseBarcode_before(inFileNames, args.stagger, barcodeLength, minPhred, asciioffset)


#Writing out barcode and associated phredscore and UMIs to file. 
os.chdir(outFileDirectory)
writeOutFileUMIs(barcode_dict, outFileUMI)             

# If excluded reads are true, bad Barcodes are also written with each issue in a separate file
if args.excludeReads == True:
	writeOutFileBadSeqRecord(missingBeforeBarcode, outFileMissingBeforeBarcode)
	check_file_created(outFileMissingBeforeBarcode)
	writeOutFileBadSeqRecord(badQscore, outFileBadPhred)
	check_file_created(outFileBadPhred)
	writeOutFileBadSeqRecord(badBarcode, outFileBadLength)
	check_file_created(outFileBadLength)
	if args.checkVector == "both":
		writeOutFileBadSeqRecord(missingAfterBarcode, outFileMissingAfterBarcode)
		check_file_created(outFileMissingAfterBarcode)

barcode_counts_dict, UMI_counts = count_read_UMI(barcode_dict)                  
    
#Writing out barcode and associated read/UMI counts to file. 
if args.includeReads:
	writeOutFileBarcodeUMICounts(barcode_counts_dict,outFileUMICounts)
	check_file_created(outFileUMICounts)
	writeOutFileBarcodeCounts(barcode_counts_dict, outFileCounts)
	check_file_created(outFileCounts)
	writeOutFileBarcodeReadCounts(barcode_counts_dict, outFileReadCounts)
	check_file_created(outFileReadCounts)
else:
	writeOutFileBarcodeCounts(barcode_counts_dict, outFileCounts)
	check_file_created(outFileCounts)

print("Writing Summary for {}".format(args.sampleName))
# Writing the summary file for this sample
summary_file = args.sampleName + "_summary.txt"
with open(summary_file, "w") as summary:
	summary.write("Number of reads parsed is {}\n".format(tot_reads))
	summary.write("The number of unique barcode is {}\n".format(len(barcode_dict)))
	summary.write("Number of reads missing sequence (GFP) before barcode is {}\n".format(len(missingBeforeBarcode)))
	if args.checkVector == "both":
		summary.write("Number of reads missing sequence (GFP) after barcode is {}\n".format(len(missingAfterBarcode)))
	summary.write("Number of reads having bad Q Score {}\n".format(len(badQscore)))
	summary.write("Number of reads having bad barcode ('AAAA,''TTTT,' 'GGGG,' 'CCCC' ) {}\n".format(len(badBarcode)))
	summary.write("Total number of reads is:{}\n".format(str(UMI_counts)))

print("Finished parsing Sample {}".format(args.sampleName))




