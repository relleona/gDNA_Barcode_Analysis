import os, subprocess
from argparse import ArgumentParser

#Command line parser
parser = ArgumentParser()
parser.add_argument("scriptPath", help = "Specify the path containing the script files.", type = str)
parser.add_argument("experimentPath", help = "Specify the path to the experiment directory", type = str)
parser.add_argument("sampleArray", help = "Specify the path to the stagger file for the experiment", type = str)
parser.add_argument("lvHistogramFraction",help = "Specify if the histogram has to be constructed for the full barcode or for a partial part. Options: full, partial", choices=["full", "partial"], default="full")
parser.add_argument("lvHistogramLength", help = "Desired length of barcode for this experiment")
args = parser.parse_args()

pathLvHistogram = os.path.join(args.scriptPath, "Step2_LVHistogram_MultipleSample","LvHistogram.py" )
pathMultipleSamples= os.path.join(args.scriptPath, "Step2_LVHistogram_MultipleSample","multipleSamples.py" )

sampleArray = args.sampleArray.split(',')

# Creating the LV distance histogram for each sample separately
for sample in sampleArray:
    commandLvHistogram = ["python3",pathLvHistogram,args.experimentPath, sample, "-inputFraction",args.lvHistogramFraction, "-inputLength",args.lvHistogramLength]
    subprocess.call(commandLvHistogram)

#To combine multiple samples so that they can be analysed together.

commandMultipleSamples = ["python3", pathMultipleSamples, args.experimentPath] + sampleArray
subprocess.call(commandMultipleSamples)



