import os, subprocess
from argparse import ArgumentParser

#Command line parser
parser = ArgumentParser()
parser.add_argument("pathScript", help = "Specify the path containing the script files.", type = str)
parser.add_argument("pathExperiment", help = "Specify the path to the experiment directory", type = str)
parser.add_argument("pathStaggerFile", help = "Specify the path to the stagger file for the experiment", type = str)
parser.add_argument("checkVector", help = "Option to check vector sequence before or on both sides of the barcode sequence.", default = "both", choices = ["both", "before"])
parser.add_argument("barcodeLength", help = "Desired length of barcode for this experiment")
parser.add_argument("minPhred", help = "Reads with 5 or more bases below this phred score in the primer binding site will not be count")
parser.add_argument("excludeReads", help = "If true, output txt.gz files containing reads excluded from the UMI and count files", default = "False", choices = ["True", "False"])
parser.add_argument("asciioffset", help = "If PhredScore has letters, ascii offset will be 33, otherwise it will be 64. Most recent version of Illumina uses Phred Score of 33. ", default = "33")

args = parser.parse_args()

pathEnvelope = os.path.join(args.pathScript, "Step1_extractBarcode","Envelope.py" )
command = ["python3", pathEnvelope, args.pathExperiment, args.pathScript, "--pathStaggerFile", args.pathStaggerFile, "-r", "-checkVector", args.checkVector,"-barcodeLength", args.barcodeLength, "-Q",args.minPhred ,"-a", args.asciioffset,"-e", args.excludeReads ]
subprocess.call(command)

