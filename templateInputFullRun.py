import os
import subprocess
#This file contains all the modifiable variables in this analysis pipeline and no other files need to be modified for regular use. Running this file will call all the steps for the analysis. The file structure of the script must be retained to run this pipeline easily. The description for each parameters is in the documentation doc.

parameters = {
    "pathScript": "/home/keerthana/Goyal_Lab/FateMapPipeline_Ubuntu/gDNA_Barcode_Extraction/",
    "pathExperimentFolder": "/home/keerthana/Goyal_Lab/FateMapPipeline_Ubuntu/TimeMachine_Keerthana_code/Subia_Dataset",

    "step1ExtractBarcode": {
        "staggerFile": "/home/keerthana/Goyal_Lab/FateMapPipeline_Ubuntu/TimeMachine_Keerthana_code/Subia_Dataset/Subia_dataset_staggerList.csv",
        "checkVector": "before",
        "barcodeLength": "90",
        "minPhred": "14",
        "excludedReads": "False"
    },

    "step2LvHistogramMultipleSamples": {
        "sampleArray": ["SA2-037_1", "SA2-037_2", "SA2-037_3"],
        "lvHistogramFraction": "partial",
        "lvHistogramLength": "50"
    },

    "pauseBeforeStep3": True,

    "step3StarcodeLvHistogram": {
        "combinedSample": "yes",
        "lengthStarcode": "50",
        "distanceStarcode": "8",
        "threadsStarcode": "4",
        "sampleArrayStarcode": ["Multiple_Samples"]
    }
}
#####################################################################################################################################
#Nothing needs to be modified in the entire pipeline beyond this.
#####################################################################################################################################


pathScript = parameters['pathScript']
pathExperiment = parameters['pathExperimentFolder']

# Step 0
pathStep0 = os.path.join(pathScript,"Step0_SampleReorganisation","Step0.py")
pathRawFiles = os.path.join(pathExperiment, "raw")
commandStep0 = ["python3", pathStep0, pathScript, pathRawFiles]
subprocess.run(commandStep0)

#Step 1
pathStep1 = os.path.join(pathScript,"Step1_extractBarcode","Step1.py")
pathStaggerFile = parameters['step1ExtractBarcode']['staggerFile']
checkVector = parameters['step1ExtractBarcode']["checkVector"]
lengthBarcode = parameters['step1ExtractBarcode']["barcodeLength"]
minPhredScore = parameters['step1ExtractBarcode']["minPhred"]
excludedReads = parameters['step1ExtractBarcode']["excludedReads"]
commandStep1 = ["python3", pathStep1, pathScript, pathExperiment, pathStaggerFile, checkVector,lengthBarcode, minPhredScore, str(excludedReads)]
subprocess.call(commandStep1)

# Step 2
pathStep2 = os.path.join(pathScript,"Step2_LVHistogram_MultipleSample","Step2.py")
sampleArray = parameters["step2LvHistogramMultipleSamples"]["sampleArray"]
sampleArray = ",".join(sampleArray)
lvHistogramFraction = parameters["step2LvHistogramMultipleSamples"]["lvHistogramFraction"]
lvHistogramLength = parameters["step2LvHistogramMultipleSamples"]["lvHistogramLength"]
commandStep2 = ["python3", pathStep2, pathScript, pathExperiment, sampleArray, lvHistogramFraction, lvHistogramLength]
subprocess.call(commandStep2)

#Step 3
if parameters["pauseBeforeStep3"] == False:
    pathStep3 = os.path.join(pathScript,"Step3_Starcode","Step3.py")
    combinedSample = parameters["step3StarcodeLvHistogram"]["combinedSample"]
    lengthStarcode = parameters["step3StarcodeLvHistogram"]["lengthStarcode"]
    distanceStarcode = parameters["step3StarcodeLvHistogram"]["distanceStarcode"]
    threadsStarcode = parameters["step3StarcodeLvHistogram"]["threadsStarcode"]
    sampleArrayStarcode = parameters["step3StarcodeLvHistogram"]["sampleArrayStarcode"]
    sampleArrayStarcode = ",".join(sampleArrayStarcode)
    commandStep3 = ["python3", pathStep3, pathScript,pathExperiment, combinedSample,lengthStarcode, distanceStarcode, threadsStarcode, sampleArrayStarcode]
    subprocess.call(commandStep3)
else:
    print("Step 3 is paused for review of LV plots.")

