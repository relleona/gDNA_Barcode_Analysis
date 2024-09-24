# This contains the various utility functions used in parseFASTQ_main.py

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import regex as re #not regular re from python, different package!
from collections import Counter

#String matching function to find the 
def find_best_match(sequence, pattern, max_errors):
	"""
    This function searches for the best match of the pattern in the sequence, where 'best' is defined as the match with 
    the fewest errors (mismatches) up to a maximum allowed number. It uses regex to efficiently find potential starting 
    positions and then performs a character-by-character comparison to count mismatches.

    Parameters:
    sequence (str): The DNA sequence to search within.
    pattern (str): The pattern to search for in the sequence.
    max_errors (int): The maximum number of mismatches allowed for a valid match.

    Returns:
    tuple: A tuple containing four elements:
        1. best_match (str): The best matching substring found, or None if no valid match is found.
        2. position (list): A list containing the [start, end] positions of the best match in the sequence, 
           or [None, None] if no valid match is found.
        3. best_errors (int): The number of mismatches in the best match, or None if no valid match is found.

    This approach ensures that the entire length of the pattern is always considered in the matching process.
	"""
	pattern_length = len(pattern)

	# Find all potential starting positions that match the first character of the pattern.
	start_positions = [match.start() for match in re.finditer(f'(?={pattern[0]})', sequence)]

	#Initialize, if both values stays None, it indicates that there no match that fits the condition of errors <= max_error
	best_match = None
	best_position = None
	#Setting this to max_errors + 1 ensures that any valid match will be considered better than our initial value (error <=max_error)
	best_errors = max_errors + 1


	#This starts a loop that checks each potential starting position.
	for start in start_positions:
		# skips positions that are too close to the end of the sequence to fit the full pattern.
		if start + pattern_length > len(sequence):
			continue
		
		#This extracts a substring from the sequence that's the same length as the pattern.
		substring = sequence[start:start+pattern_length]
		#This counts the number of mismatches between the substring and the pattern.
		errors = sum(1 for a, b in zip(substring, pattern) if a != b)
		
		#It keeps track of the best match encountered (the one with the fewest errors that is still within the maximum allowed errors).
		# if errors <= max_errors (Is this a valid match within our error tolerance?) and errors < best_errors ( Is this match better than the best we've seen so far?)
		if errors <= max_errors and errors < best_errors:
			best_match = substring
			best_errors = errors
			best_position = start

	if best_position is not None:
		#Make a list of start and end of the substring so that it can be used in the later functions 
		position = [best_position, best_position + pattern_length]
	else:
		position = [None, None]


	return best_match, position, best_errors 


def parseBarcode_both(inFilename, staggerLength, barcodeLength, minQuality_Phred, asciioffset):
	"""Function to parse FASTQ file and create dictionary of key:value pairs. Will only extract reads that contain the vector sequence before and after the barcode, allowing up to 4 or 5 mismatches respectively. 
	Each key is a unique barcode sequence. Each value is a list containing the Phredscore associated with the barcode
	and the UMIs (first 4-6 bases) associated with the read. Function will also output the list of reads that are missing vector sequence before the barcode, missing the vector sequence after the barcode, have a bad Qscore before the barcode or contain a short UMI or primer binding sequence. 
	"""
	barcode_dict = {}
	# creating dictionaries for issues with identifying vector sequences or low quality barcode sequencing
	missingVectorBefore = []
	missingVectorAfter = []
	badQscore = []
	badLength = []
	badBarcode = []

	vectorBeforeBarcode = re.compile(r'(?e)(TCGACTAAACGCGCTACTTGAT){e<=4}') #vector primer site as reg expression. Allow up to 4 mismatches. Less than 4 errors - including deletion, insertion and substitution aka GFP primers 
	vectorAfterBarcode = re.compile(r'(?e)(?r)(ATCCTACTTGTACAGCTCGT){e<=5}') #vector sequence after barcode as reg expression. Allow up to 5 mismatches and search from end of string first. ***What determines these numbers?
	tot_reads = 0
	for fastQ_file in inFilename: #Loop over all fastQ files/lanes associated with each sample
		with gzip.open(fastQ_file, 'rt') as fastq:
			print("Started with file:{}".format(fastQ_file))
			for seq_record in FastqGeneralIterator(fastq):
				tot_reads += 1

				# Find best match for vector before barcode
				VBB_match, VBB_position, VBB_errors = find_best_match(seq_record[1], vectorBeforeBarcode, 4)

				# Find best match for vector after barcode
				VBA_match, VBA_position, VBA_errors = find_best_match(seq_record[1], vectorAfterBarcode, 5)

				if VBB_match and VBA_match:
					if sum([ord(i) - asciioffset < minQuality_Phred for i in seq_record[2][VBB_position[0]:VBB_position[1]]]) >= 5:
						badQscore.append(seq_record)
					elif (VBB_position[1] - staggerLength) > 30 or (VBB_position[0] - staggerLength) < 4:
						badLength.append(seq_record)
					elif(len(re.findall("(AAAA)", seq_record[1])) > 0 or \
					len(re.findall("(TTTT)", seq_record[1])) > 0 or \
					len(re.findall("(GGGG)", seq_record[1])) > 0 or \
					len(re.findall("(CCCC)", seq_record[1])) > 0 or \
					len(re.findall("(NN)", seq_record[1])) > 0):
						badBarcode.append(seq_record)
					else:
						barcode = seq_record[1][VBB_position[1]:VBA_position[0]]
						if barcode not in barcode_dict:
							barcode_dict[barcode] = [seq_record[1][VBB_position[1]:VBA_position[0]]]
							barcode_dict[barcode].append(seq_record[1][0:VBB_position[0]-staggerLength])
						elif barcode in barcode_dict:
							barcode_dict[barcode].append(seq_record[1][0:VBB_position[0]-staggerLength])
						elif VBB_match and not VBA_match:
							missingVectorAfter.append(seq_record)
						elif VBA_match and not VBB_match:
							missingVectorBefore.append(seq_record)

	print("Completed file " + fastQ_file)

	return barcode_dict, missingVectorBefore, missingVectorAfter, badQscore, badLength, badBarcode, tot_reads

def parseBarcode_before(inFileNames, staggerLength, barcodeLength, minPhred, asciioffset):
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
	# vectorAfterBarcode = re.compile(r'(?e)(?r)(ATCCTACTTGTACAGCTCGT){e<=5}') #vector sequence after barcode as reg expression. Allow up to 5 mismatches and search from end of string first. ***What determines these numbers?
	tot_reads = 0
	
	
	# This starts a loop over all fastQ files associated with each sample.
	for fastQ_file in inFileNames:
		#Loop over all fastQ files/lanes associated with each sample
		# This opens each gzipped fastQ file in text mode.
		with gzip.open(fastQ_file, 'rt') as fastq:

			# This uses BioPython's FastqGeneralIterator to parse each read in the fastQ file.
			for seq_record in FastqGeneralIterator(fastq): #Uses BioPythons FastqGeneralIterator to parse each read.
				tot_reads += 1

				# This searches for the vector before barcode sequence in the current read. Allow up to 4 mismatches. Less than 4 errors - including deletion, insertion and substitution
				# position is the position of the vector before the barcode sequence aka GFP 
				# best_match is the best match string that is closest to GFP both in length and sequence

				best_match, position, error_count= find_best_match(seq_record[1], vectorBeforeBarcode, 4)
	
				# If the vector before barcode sequence is found: 
				# best_match is None means none is found, and error_count is the count 
				if best_match is not None: 
					
					# This checks if 5 or more positions in the matched region have a Phred score < minPhred.
					if sum([ord(i) - asciioffset < minPhred for i in seq_record[2][position[0]:position[1]]]) >= 5: #Skip reads where >=5 positions in GFP primer site have phredscore < 14.   
						badQscore.append(seq_record) 

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
							barcode_dict[barcode] = [seq_record[2][position[1]:position[1] + barcodeLength]]
							barcode_dict[barcode].append(seq_record[1][0:position[0]-staggerLength])
						elif barcode in barcode_dict:
							barcode_dict[barcode].append(seq_record[1][0:position[0]-staggerLength][0:position[0]-staggerLength])
		
				# If the vector before barcode sequence is not found, this adds the read to missingVectorBefore.
				else:
					missingVectorBefore.append(seq_record) #Not consider barcode without the GFP tag			
		
		print("Completed file " + fastQ_file)
	return barcode_dict, missingVectorBefore, badQscore, badBarcode, tot_reads

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

