# This contains the various utility functions used in parseFASTQ_main.py

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import os 
import regex as re #not regular re from python, different package!
from collections import Counter
from itertools import zip_longest
import math


def combine_fastq(inFileNames):
    """
    Combine multiple gzipped FASTQ files into a single gzipped FASTQ file,
    preserving the structure of each FASTQ record.

    Args:
    inFileNames (list): List of input FASTQ file names (gzipped)

    Returns:
    str: Name of the output file
    """
    # Get the absolute path of the first input file
    first_file_path = os.path.abspath(inFileNames[0])
    
    # Get the directory of the first input file
    input_directory = os.path.dirname(first_file_path)
    
    # Get the name of the folder containing the input files
    folder_name = os.path.basename(input_directory)

    # Create the output filename using the folder name
    outFileName = f"{folder_name}_combined.fastq.gz"

    # Create the full path for the output file
    outFilePath = os.path.join(input_directory, outFileName)

    # # Check if the file already exists
    # if os.path.exists(outFilePath):
    #     confirm = input(f"The file {outFileName} already exists. Do you want to overwrite it? (y/n): ")
    #     if confirm.lower() != 'y':
    #         print("Operation cancelled.")
    #         return None  # Return None to indicate cancellation

    # Open all input files and the output file
    input_files = [gzip.open(f, 'rt') for f in inFileNames if os.path.basename(f) != outFileName]
    
    with gzip.open(outFilePath, 'wt') as outfile:
        while True:
            eof = True
            for f in input_files:
                # Read 4 lines (one FASTQ record) from each file
                record = [f.readline() for _ in range(4)]
                if record[0]:  # If we read something (not end of file)
                    eof = False
                    # Write the record to the output file
                    outfile.writelines(record)
            if eof:
                break  # Exit the loop if all files have reached EOF
    
    # Close all input files
    for f in input_files:
        f.close()

    # Print final completion message with folder information
    print(f"All files have been combined into {outFileName}")
    print(f"Combined file is located in folder: {input_directory}")

    # Return the name of the output file
    return outFileName


def find_best_match(sequence, pattern, max_errors, before_or_after):
    """
    This function searches for the best match of the pattern in the sequence, allowing for partial matches
    if a full match is not found. It progressively reduces the pattern length until a match is found or
    the maximum number of errors is reached.
    """
    def find_match(seq, pat, max_err):
        best_match = None
        best_position = None
        best_errors = max_err + 1

        for start in range(len(seq) - len(pat) + 1):
            substring = seq[start:start+len(pat)]
            errors = sum(1 for a, b in zip(substring, pat) if a != b)
            
            if errors <= max_err and errors < best_errors:
                best_match = substring
                best_errors = errors
                best_position = start

        position = [best_position, best_position + len(pat)] if best_position is not None else [None, None]
        return best_match, position, best_errors

    truncation_errors = 0
    pattern_length = len(pattern)

    # Try matching with progressively shorter patterns
    for length in range(pattern_length, pattern_length - max_errors - 1, -1):
        truncated_pattern= pattern[:length]
        # Find the best match for this truncated pattern
        match, position, errors = find_match(sequence, truncated_pattern, max_errors - truncation_errors)
        
        # If a match is found and total errors are within max_errors, return the results
        if match is not None and (errors + truncation_errors) <= max_errors:
            # print(f"Match found for length {length}: {match}")
            return match, position, errors + truncation_errors

        truncation_errors += 1

    # print("No match found even after reducing length")
    return None, [None, None], None


def parseBarcode_both(inFileName, staggerLength, barcodeLength, minQuality_Phred, asciioffset):
	"""Function to parse FASTQ file and create dictionary of key:value pairs. Will only extract reads that contain the vector sequence before and after the barcode, allowing up to 4 or 5 mismatches respectively. 
	Each key is a unique barcode sequence. Each value is a list containing the Phredscore associated with the barcode
	and the UMIs (first 4-6 bases) associated with the read. Function will also output the list of reads that are missing vector sequence before the barcode, missing the vector sequence after the barcode, have a bad Qscore before the barcode or contain a short UMI or primer binding sequence. 
	This is what barcode_dict will look like 
	barcode_dict = {
    "barcode_sequence": [
        "quality_score_of_barcode",
        "sequence_before_stagger",
        "sequence_before_stagger",
        "sequence_before_stagger",
        ... potentially more sequences
    ],
    ... potentially more barcode entries
	}

	"""
	barcode_dict = {}
	# creating dictionaries for issues with identifying vector sequences or low quality barcode sequencing
	missingVectorBefore = []
	missingVectorAfter = []
	badQscore = []
	badLength = []
	badBarcode = []

	vectorBeforeBarcode = "TCGACTAAACGCGCTACTTGAT" #
	vectorAfterBarcode = "ATCCTACTTGTACAGCTCGT"
	tot_reads = 0

	fastQ_file = combine_fastq(inFileName)
	find_after_barcode = staggerLength + len(vectorBeforeBarcode) + barcodeLength
	find_before_barcode = int(math.floor(staggerLength + len(vectorBeforeBarcode) + (barcodeLength / 2)))
	with gzip.open(fastQ_file, 'rt') as fastq:
		print("Started with file:{}".format(fastQ_file))
		for seq_record in FastqGeneralIterator(fastq):
			tot_reads += 1

			# Find best match for vector before barcode
			VBB_match, VBB_position, VBB_errors = find_best_match(seq_record[1][:find_before_barcode], vectorBeforeBarcode, 4, "before")

			# Find best match for vector after barcode
			VBA_match, VBA_position, VBA_errors = find_best_match(seq_record[1][find_after_barcode:], vectorAfterBarcode, 5, "after")

			if VBB_match and VBA_match:
				#  This checks if 5 or more positions in the matched region have a Phred score < minPhred within the GFP primer site.
				if sum([ord(i) - asciioffset < minQuality_Phred for i in seq_record[2][VBB_position[0]:VBB_position[1]]]) >= 5:
					badQscore.append(seq_record) #Skip reads where conditions above are not fulfilled.

				# This checks for homopolymers or unknown nucleotides in the sequence.
				elif(len(re.findall("(AAAA)", seq_record[1])) > 0 or \
				len(re.findall("(TTTT)", seq_record[1])) > 0 or \
				len(re.findall("(GGGG)", seq_record[1])) > 0 or \
				len(re.findall("(CCCC)", seq_record[1])) > 0 or \
				len(re.findall("(NN)", seq_record[1])) > 0):
					badBarcode.append(seq_record)

				# If the sequence passes all checks, this extracts the barcode and updates the barcode_dict.
				else:
					
					barcode = seq_record[1][VBB_position[1]:VBB_position[1]+ barcodeLength] #recording barcode
					if barcode not in barcode_dict:
						barcode_dict[barcode] = [seq_record[2][VBB_position[1]:VBB_position[1]+ barcodeLength]] # record quality score 
						barcode_dict[barcode].append(seq_record[1][0:VBB_position[0]-staggerLength]) # record stagger sequence
					elif barcode in barcode_dict:
						barcode_dict[barcode].append(seq_record[1][0:VBB_position[0]-staggerLength])
			elif VBB_match and not VBA_match:
				missingVectorAfter.append(seq_record)
			elif VBA_match and not VBB_match:
				missingVectorBefore.append(seq_record)

	print("Completed file " + fastQ_file)

	return barcode_dict, missingVectorBefore, missingVectorAfter, badQscore, badLength, badBarcode, tot_reads


def parseBarcode_before(inFileName, staggerLength, barcodeLength, minPhred, asciioffset):
	# Function to parse FASTQ file and create dictionary of key:value pairs. Will only extract reads that contain the vector sequence before the barcode allowing up to 4 mismatches. 
	# Each key is a unique barcode sequence. Each value is a list containing the Phredscore associated with the barcode
	# and the UMIs (first 4-6 bases) associated with the read. Function will also output the list of reads that are missing vector sequence before the barcode, missing the vector sequence after the barcode, have a bad Qscore before the barcode or contain a short UMI or primer binding sequence.
	barcode_dict = {}
	# creating dictionaries for issues with identifying vector sequences or low quality barcode sequencing
	missingVectorBefore = []
	badQscore = []
	badLength = []
	badBarcode = []

	# the GFP sequence as a vector in the primer site before the barcode {read1 from 5'-3'}
	vectorBeforeBarcode = "TCGACTAAACGCGCTACTTGAT" #
	find_before_barcode = int(math.floor(staggerLength + len(vectorBeforeBarcode) + (barcodeLength / 2)))
	# vectorAfterBarcode = re.compile(r'(?e)(?r)(ATCCTACTTGTACAGCTCGT){e<=5}') #vector sequence after barcode as reg expression. Allow up to 5 mismatches and search from end of string first. ***What determines these numbers?
	tot_reads = 0
	
	fastQ_file = combine_fastq(inFileName)
	
	# This starts a loop over all fastQ files associated with each sample.
	with gzip.open(fastQ_file, 'rt') as fastq:

		# This uses BioPython's FastqGeneralIterator to parse each read in the fastQ file.
		for seq_record in FastqGeneralIterator(fastq): #Uses BioPythons FastqGeneralIterator to parse each read.
			tot_reads += 1

			# This searches for the vector before barcode sequence in the current read. Allow up to 4 mismatches. Less than 4 errors - including deletion, insertion and substitution
			# position is the position of the vector before the barcode sequence aka GFP 
			# best_match is the best match string that is closest to GFP both in length and sequence

			best_match, position, error_count= find_best_match(seq_record[1][:find_before_barcode], vectorBeforeBarcode, 4, "before")

			# If the vector before barcode sequence is found: 
			# best_match is None means none is found, and error_count is the count 
			if best_match is not None: 
				
				# This checks if 5 or more positions in the matched region have a Phred score < minPhred.
				if sum([ord(i) - asciioffset < minPhred for i in seq_record[2][position[0]:position[1]]]) >= 5: #Skip reads where >=5 positions in GFP primer site have phredscore < 14.   
					badQscore.append(seq_record) 

				# elif (int(position[1]) - staggerLength) > 30 or (int(position[0]) - staggerLength) < 4: #Skip positions with UMI shorter than 4 bases or UMI+GFP longer than 30 bases. 
				# 	badLength.append(seq_record)

				# This checks for homopolymers or unknown nucleotides in the sequence.
				elif(len(re.findall("(AAAA)", seq_record[1])) > 0 or \
						len(re.findall("(TTTT)", seq_record[1])) > 0 or \
						len(re.findall("(GGGG)", seq_record[1])) > 0 or \
						len(re.findall("(CCCC)", seq_record[1])) > 0 or \
						len(re.findall("(NN)", seq_record[1])) > 0):
						badBarcode.append(seq_record)

				# # This checks for homopolymers or unknown nucleotides in the sequence.
				# elif(len(re.findall("(AAAA)", seq_record[1][position[1]:])) > 0 or \
				# 		len(re.findall("(TTTT)", seq_record[1][position[1]:])) > 0 or \
				# 		len(re.findall("(GGGG)", seq_record[1][position[1]:])) > 0 or \
				# 		len(re.findall("(CCCC)", seq_record[1][position[1]:])) > 0 or \
				# 		len(re.findall("(NN)", seq_record[1][position[1]:])) > 0):
				# 		badBarcode.append(seq_record)
				
				# If the sequence passes all checks, this extracts the barcode and updates the barcode_dict.
				else:
					barcode = seq_record[1][position[1]:position[1] + barcodeLength] #recording barcode
					if barcode not in barcode_dict:
						barcode_dict[barcode] = [seq_record[2][position[1]:position[1] + barcodeLength]] # recording the quality 
						barcode_dict[barcode].append(seq_record[1][0:position[0]-staggerLength]) # recording stagger sequence, demultiplexing
					elif barcode in barcode_dict:
						barcode_dict[barcode].append(seq_record[1][0:position[0]-staggerLength])
	
			# If the vector before barcode sequence is not found, this adds the read to missingVectorBefore.
			else:
				missingVectorBefore.append(seq_record) #Not consider barcode without the GFP tag			
		
		print("Completed file " + fastQ_file)
	return barcode_dict, missingVectorBefore, badQscore, badLength, badBarcode, tot_reads


def writeOutFileUMIs(barcode_dict, outFileName):
	"""Function to write barcode dictionary to gzipped tab delimited file. Each line contains
    a unique barcode, it's phredscore, and all associated UMIs. 
    """
	with gzip.open(outFileName, 'wt') as out_file:
		for barcode in barcode_dict:
			out_file.write(barcode)
			out_file.write("\t" + "\t".join(barcode_dict[barcode]))
			out_file.write("\n")

def writeOutFileBadSeqRecord(badSeqList, outFileName):
    """Function to write a list of bad sequence records to gzipped tab delimited file. 
    """
    with gzip.open(outFileName, 'wt') as out_file:
    	for seqRecord in badSeqList:
            out_file.write("\t".join(seqRecord) + "\n")

def count_read_UMI(barcode_dictionary):
	"""Function to count the number of reads and the number of UMIs associated with each barcode. 
	However, if there is no UMIs, it will simply count total number of reads per barcode. 
    Output is dictionary with the latter as the key and the former as the associated values
    """
	new_dict = {}
	tot_reads = 0
	for i in barcode_dictionary:
		#First component is the number of reads per barcode 
		#Second component is the number of unique UMIs
		new_dict[i] = (sum(Counter(barcode_dictionary[i][1:]).values()), len(Counter(barcode_dictionary[i][1:]).keys()))
		tot_reads += new_dict[i][0]
	return new_dict, tot_reads

def writeOutFileBarcodeCounts(barcode_dict_summary, outFileName):
	"""Function to write reads counts and UMIs associated with each barcode
	to gzipped, tab-delimited file.
	"""
	with gzip.open(outFileName, 'wt') as out_file:
		for barcode in barcode_dict_summary:
			out_file.write(barcode)
			out_file.write("\t" + "\t".join(map(str,barcode_dict_summary[barcode])))
			out_file.write("\n")    

def writeOutFileBarcodeReadCounts(barcode_dict_summary, outFileName):
	"""Function to write reads counts and UMIs associated with each barcode
	to gzipped, tab-delimited file.
	"""
	with gzip.open(outFileName, 'wt') as out_file:
		for barcode in barcode_dict_summary:
			out_file.write(barcode)
			out_file.write("\t" + str(barcode_dict_summary[barcode][0]))
			out_file.write("\n")    

def writeOutFileBarcodeUMICounts(barcode_dict_summary, outFileName):
	"""Function to write UMIs counts associated with each barcode
	to gzipped, tab-delimited file.
	"""
	with gzip.open(outFileName, 'wt') as out_file:
		for barcode in barcode_dict_summary:
			out_file.write(barcode)
			out_file.write("\t" + str(barcode_dict_summary[barcode][1]))
			out_file.write("\n")    

