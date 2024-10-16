import os, subprocess
from argparse import ArgumentParser
from datetime import datetime

#Command line parser
parser = ArgumentParser()
parser.add_argument("scriptPath", help = "Specify the path containing the script files.", type = str)
parser.add_argument("experimentPath", help = "Specify the path to the experiment directory", type = str)
parser.add_argument("sampleArray", help="Array of samples within the Experiment data.", type=str)
parser.add_argument("lvHistogramFraction",help = "Specify if the histogram has to be constructed for the full barcode or for a partial part. Options: full, partial", choices=["full", "partial"], default="full")
parser.add_argument("lvHistogramLength", help = "Desired length of barcode for this experiment")
args = parser.parse_args()

# pathLvHistogram = os.path.join(args.scriptPath, "Step2_LVHistogram_MultipleSample","LvHistogram.py" )
# pathMultipleSamples= os.path.join(args.scriptPath, "Step2_LVHistogram_MultipleSample","multipleSamples.py" )

# sampleArray = args.sampleArray.split(',')

# # Creating the LV distance histogram for each sample separately
# for sample in sampleArray:
#     commandLvHistogram = ["python3",pathLvHistogram,args.experimentPath, sample, "-inputFraction",args.lvHistogramFraction, "-inputLength",args.lvHistogramLength]
#     subprocess.call(commandLvHistogram)

# #To combine multiple samples so that they can be analysed together.

# commandMultipleSamples = ["python3", pathMultipleSamples, args.experimentPath] + sampleArray +["-inputFraction",args.lvHistogramFraction, "-inputLength",args.lvHistogramLength]
# subprocess.call(commandMultipleSamples)

pathSamples = os.path.join(args.scriptPath, "Step2_LVHistogram_MultipleSample","barcodelength_multiplesamples.py" )
pathLvHistogram= os.path.join(args.scriptPath, "LVHistogram.py" )

# Creating txts of barcode_length and multiple samples
commandMultipleSamples = ["python3", pathSamples, args.experimentPath, args.sampleArray, "-inputFraction",args.lvHistogramFraction, "-inputLength",args.lvHistogramLength]
subprocess.call(commandMultipleSamples)

# Move to the Experiment folder 
os.chdir(args.experimentPath)

# Create txt file in the Experiment folder 
summary_file = "LVmeandistance.txt"

# Generate a timestamp for the output filename
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

# Indicate that it is pre_starcode 
with open(summary_file, "w") as summary:
    summary.write(f"\nSummary of Mean Levenshtein Distances {timestamp}\n")
    summary.write("==========================\n\n")

    summary.write(f"Prestarcode {timestamp}\n")
    summary.write("------------------------------------\n\n")

# Save the path to this summary file 
pathtosummary = os.path.join(args.experimentPath,summary_file)

#Split the data into a list so it can be parsed properly 
sampleArray = args.sampleArray.split(',')
sampleArray.append("Multiple_Samples")

# Creating the LV distance histogram for each sample separately
for sample in sampleArray:

    # Path to the sample that we want to do an LV analysis on and the sample folder it is in
    if args.lvHistogramFraction == "partial" :
        pathtosample = os.path.join(args.experimentPath, "analyzed", sample, "LV_Analysis", "{}_Barcode_{}.txt".format(sample, args.lvHistogramLength))
    else :
        pathtosample = os.path.join(args.experimentPath, "analyzed", sample, "LV_Analysis", "{}_Barcode_full.txt".format(sample))
    pathtofolder = os.path.join(args.experimentPath, "analyzed", sample)

    # Run LV Histogram 
    commandLvHistogram = ["python3", pathLvHistogram, pathtofolder, pathtosample, pathtosummary, sample, "-inputFraction",args.lvHistogramFraction, "-inputLength",args.lvHistogramLength]
    subprocess.call(commandLvHistogram)

