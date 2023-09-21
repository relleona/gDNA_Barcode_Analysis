import os, subprocess
from argparse import ArgumentParser

#Command line parser
parser = ArgumentParser()
parser.add_argument("scriptPath", help = "Specify the path containing the script files.", type = str)
parser.add_argument("pathRawFiles", help = "Specify the path containing the raw sample files that need to be reorganised", type = str)
args = parser.parse_args()

#python3 /home/mzo5929/gDNA_Barcode_Extraction/Step0_SampleReorganisation/BaseSpaceSampleOrganisation.py
pathBaseSpaceSampleOrganisation = os.path.join(args.scriptPath, "Step0_SampleReorganisation","BaseSpaceSampleOrganisation.py" )
command = ["python3", pathBaseSpaceSampleOrganisation, args.path]
subprocess.call(command)