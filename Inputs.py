"""
Data Processing Pipeline for FASTQ Reorganization and Job Submission

Description:
This script automates the process of reorganizing FASTQ files and submitting jobs 
for further processing in a high-performance computing environment. It is designed 
to work with the SLURM workload manager on the Quest cluster at Northwestern University.

The script performs the following main tasks:
1. Runs a FASTQ reorganization script (BaseSpaceReorganization.py)
2. Optionally pauses for user verification
3. Generates and submits SLURM job scripts using make_bash.py
4. Monitors the submitted jobs until completion

Requirements:
- Python 3.x
- Access to the Quest cluster at Northwestern University
- SLURM workload manager
- Required Python libraries: datetime, time, os, sys, subprocess

Usage:
1. Modify the 'input_dir' and 'bash_script' dictionaries to set your specific paths and parameters.
2. Set the 'Pause_to_check' variable to True if you want to stop after reorganization, 
   or to a number of minutes if you want to pause before job submission.

Command to run this file 
python3 <path to Step0_2.py> 

Configuration:
- input_dir: Dictionary containing paths to various input files and directories
- bash_script: Dictionary containing SLURM job submission parameters
- Pause_to_check: Control variable for pausing or stopping the script
- max_wait_time1, max_wait_time2: Maximum wait times for job completion

Note:
Modify paths and parameters as necessary for your specific environment.

- created by Aurelia on 2024-08-26 and last modified on  at 
"""

# Import necessary libraries
from datetime import date, timedelta, datetime
import time as time_module  # Import with a different name to avoid conflicts
import os
import sys
import subprocess
from argparse import ArgumentParser


#Variables that needs to be changed 

#Input and output directories needed for the reorganization step
input_dir = {
    #Path to all scripts
    "pathScript": "/projects/b1042/GoyalLab/aleona/gDNA_Barcode_Analysis",
    
    #The path to the downloaded FASTQ files
    "path":"/projects/b1042/GoyalLab/aleona/YG10Xbarcode/Experiments",
    
    #The path to sample sheet, if doesn't exist put "NA"
    "sample_sheet":"NA", #created or taken from ENSEMBLE/GENECODE 
    
    #The path to primer sheet, if doesn't need put "NA"
    "primers_sheet": "NA",
    
    # If stagger file all exists already without the need for Sample Sheet or Primers Sheet put the path, if not place it as "NA"
    "staggerfileall": "/projects/b1042/GoyalLab/aleona/YG10Xbarcode/Exp3/StaggerFileAll.csv", # Can be NA
    
    #The path to folders that will later contain the experiments 
    "exp": "/projects/b1042/GoyalLab/aleona/YG10Xbarcode/Exp3",

    #Depth of subfolders within the Main Experiment folder above to run the job scripts for
    "depth" : "1", #Default is 1 and usually its the same if the format is followed, 

    #If your data is already reorganized with StaggerFile.csv, then it will just initiate the bash script in each of the Experiment folders
    "reorganized" : True
}

parameters = {
    "step1ExtractBarcode": {
        # If they have primers before and after barcode
        "checkVector": "both", # choices : ["both", "before"]

        # Length of the barcode 
        "barcodeLength": "100",
        
        # Minimum phred score 
        "minPhred": "14",
        
        # Specify whether to save unwanted (bad) barcodes to a separate file
        # Set to False to save bad barcodes, True to exclude them   
        "excludedReads": "False", 

        # If PhredScore has letters, ascii offset will be 33, otherwise it will be 64. 
        # Most recent version of Illumina uses Phred Score offset of 33.
        "asciioffset": "33"
    },

    # This step will make a barcode with the length that is specified and conduct an LV analysis on that length. 
    "step2LvHistogramMultipleSamples": {
        # Will you use a part of the barcode for starcode processing or the entire thing 
        "lvHistogramFraction": "partial", # choices : ["partial", "full"]

        # If partial is specified above, what is the length of the barcode?
        # If full is specified above, indicate the same number as barcodeLength above?
        "lvHistogramLength": "90"
    },

    # Do you want to see the LV Histograms and check if you have specified the right length
    "pauseBeforeStep3": False,

    # This step will process starcode according to the length above and the parameters below
    "step3Starcode": {
        # Determine how starcode will process the data
        #   yes: Perform starcode on combined samples from one experiment (Multiple_Samples)
        #   no: Perform starcode on separate samples from the same experiment
        "combinedSample": "yes", # choices : ["yes", "no"]

        # Distance or differences by which starcode will collapse or merge the barcode 
        "distanceStarcode": "8", # Maximum of 8 
    }
}

bash_script ={
    "run_bash": "True",
    "quest" : "True", #Default is TRUE
 ##For creating the bash scripts in Quest
    "account":"b1042",
    "partition": "genomics", 
    "nodes": "1",
    "ntasks":"8",
    "memory":"20GB",
    "time":"12:00:00", 
    "email":"AureliaLeona2028@u.northwestern.edu",
    "mail_type":"END,FAIL" #BEGIN, END, FAIL AND ALL
}



#####################################################################################################################################
#Nothing needs to be modified in the entire pipeline beyond this.
#####################################################################################################################################

###Run Command---------------------------
#'
#' This function executes a shell command and handles its output and errors.
#'
#' @param command A list of strings representing the command to be executed.
#'
#' @return The standard output of the command if successful.
#'
#' @details The function performs the following steps:
#'          1. Prints the command being run
#'          2. Executes the command using subprocess.run()
#'          3. Prints and returns the stdout if successful
#'          4. Prints error details and exits if the command fails
def run_command(command):
    """
    Execute a shell command and handle its output and errors.
    
    Args:
    command (list): A list of strings representing the command to be executed.
    
    Returns:
    str: The standard output of the command if successful.
    
    Raises:
    SystemExit: If the command fails, exits the script with an error message.
    """
    print("Running command:", " ".join(map(str, command)))
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(result.stdout)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}", file=sys.stderr)
        print(f"STDOUT: {e.stdout}", file=sys.stderr)
        print(f"STDERR: {e.stderr}", file=sys.stderr)
        sys.exit(1)

# Printing this to a txt file 
def write_main_folders_dict(output_file, main_folder_dict):
    # print(main_folder_dict)
    with open(output_file, "a") as file:
        for folder_path, samples in main_folder_dict.items():
            # Join samples with comma and space, ensure that there are no duplicates 
            samples_string = ",".join(np.unique(samples))
            # Write the folder path and samples to the file
            file.write(f"{folder_path} = {samples_string}\n")

###Monitor Jobs---------------------------
#'
#' This function monitors the status of submitted jobs and waits for their completion.
#'
#' @param job_ids A list of strings representing the IDs of jobs to monitor.
#' @param max_wait_time An integer representing the maximum time to wait for job completion in seconds.
#'
#' @return None
#'
#' @details The function performs the following steps:
#'          1. Defines an inner function to check if jobs are running
#'          2. Enters a loop to periodically check job status
#'          3. Waits between checks, doubling the interval up to 30 minutes
#'          4. Prints status updates including elapsed time
#'          5. Exits if jobs don't complete within max_wait_time
#'          6. Prints completion message when all jobs are done
#'
def monitor_jobs(job_ids):
    # Define a nested function to check if any of the given jobs are still running
    def are_jobs_running(ids):
        # Run the squeue command to get information about the jobs
        result = subprocess.run(['squeue', '-j', ','.join(ids)], capture_output=True, text=True)
        # Return True if any job ID is found in the squeue output, False otherwise
        return any(id in result.stdout for id in ids)

    # Print the job IDs we're going to monitor
    print(f"Waiting for jobs to complete: {', '.join(job_ids)}")
    
    # Record the start time for calculating elapsed time
    start_time = time_module.time()
    
    # Initial check interval (1 minute)
    check_interval = 60  

    # Continue looping as long as at least one job is still running
    while are_jobs_running(job_ids):
        # Wait for the specified interval before checking again
        time_module.sleep(check_interval)
        
        # Double the check interval, but cap it at 30 minutes (1800 seconds)
        check_interval = min(check_interval * 2, 1800)  
        
        # Calculate and print the elapsed time
        elapsed_time = time_module.time() - start_time
        print(f"Jobs still running... (Elapsed time: {str(timedelta(seconds=int(elapsed_time)))})")

    # If we exit the while loop, it means all jobs have completed
    print("All jobs are done!")



if __name__ == "__main__":
    """
    Main function to orchestrate the data processing pipeline.
    """
    # Assign variables from the configuration dictionaries
    pathScript = input_dir['pathScript']
    path = input_dir['path']
    sample_sheet = input_dir['sample_sheet']
    primers_sheet = input_dir['primers_sheet']
    staggerfileall = input_dir['staggerfileall']

    pathExperiment = input_dir['exp']
    # Make LV_Analysis and matrix folder
    if not os.path.exists(pathExperiment):
        os.makedirs(pathExperiment)
    
    depth = input_dir['depth']
    reorganized = input_dir['reorganized']

    run_bash = bash_script['run_bash']
    quest = bash_script['quest']
    account = bash_script['account']
    partition = bash_script['partition']
    nodes = bash_script['nodes']
    ntasks = bash_script['ntasks']
    memory = bash_script['memory']
    time = bash_script['time']
    email = bash_script['email']
    mail_type = bash_script['mail_type']

    #Step 1
    pathStep1 = os.path.join(pathScript,"Step1_extractBarcode","Envelope.py")
    checkVector = parameters['step1ExtractBarcode']["checkVector"]
    lengthBarcode = parameters['step1ExtractBarcode']["barcodeLength"]
    minPhredScore = parameters['step1ExtractBarcode']["minPhred"]
    excludedReads = parameters['step1ExtractBarcode']["excludedReads"]
    asciioffset = parameters['step1ExtractBarcode']["asciioffset"]

    #Step 2
    pathStep2 = os.path.join(pathScript,"Step2_LVHistogram_MultipleSample","Step2.py")
    lvHistogramFraction = parameters["step2LvHistogramMultipleSamples"]["lvHistogramFraction"]
    lvHistogramLength = parameters["step2LvHistogramMultipleSamples"]["lvHistogramLength"]

    #Step 3
    pathStep3 = os.path.join(pathScript,"Step3_Starcode","Step3.py")
    combinedSample = parameters["step3Starcode"]["combinedSample"]
    distanceStarcode = parameters["step3Starcode"]["distanceStarcode"]


    ## Make input txt file for the new makebash script so there can be an editable bash script 
    # Create a txt file that contains all the necessary inputs for the bash script in the different experiments 
    os.chdir(pathExperiment)

    # Create txt file in the Experiment folder 
    input_file = "Templateinput.txt"

    # Generate a timestamp for the output filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")

    # Open the file in write mode
    with open(input_file, 'w') as summary:
    # Write content to the file
        summary.write(f"\nInputs to make the BASH Scripts\n")
        summary.write("==========================\n\n")

        summary.write(f"Initiate Step 1 {timestamp}\n")
        summary.write("------------------------------------\n")
        summary.write(f"STEP1PATH = {pathStep1}\n")
        summary.write(f"SCRIPT_DIR = {pathScript}\n")
        summary.write(f"CHECKVECTOR = {checkVector}\n")
        summary.write(f"BARCODELENGTH = {lengthBarcode}\n")
        summary.write(f"MINPHRED = {minPhredScore}\n")
        summary.write(f"EXCLUDEREADS = {excludedReads}\n")
        summary.write(f"ASCIIOFFSET = {asciioffset}\n\n")

        summary.write(f"Initiate Step 2 {timestamp}\n")
        summary.write("------------------------------------\n")
        summary.write(f"STEP2PATH = {pathStep2}\n")
        summary.write(f"FRACTION = {lvHistogramFraction}\n")
        summary.write(f"INPUTLENGTH = {lvHistogramLength}\n\n")

        summary.write(f"Initiate Step 3 {timestamp}\n")
        summary.write("------------------------------------\n")
        if parameters["pauseBeforeStep3"] == False:
            summary.write(f"PAUSEBEFORESTEP3 = {False}\n")
        else: 
            summary.write(f"PAUSEBEFORESTEP3 = {True}\n")
        summary.write(f"STEP3PATH = {pathStep3}\n")
        summary.write(f"COMBINESAMPLE = {combinedSample}\n")
        summary.write(f"DISTANCE = {distanceStarcode}\n\n\n")
    
    print(f"Input file '{input_file}' has been created in {pathExperiment}")

    if reorganized == True: 
        main_folders = {}

        for root, folders, files in os.walk(pathExperiment):
            # Get the current depth of the directory
            depth = root[len(pathExperiment):].count(os.sep)

            # We're interested in the second level of subdirectories (depth == 1)
            if depth == 1:
                # The experiment path is the current root
                experiment_path = root
                # The samples are in the 'raw' subfolder
                raw_path = os.path.join(root, 'raw')
                
                if os.path.exists(raw_path):
                    samples = [folder for folder in os.listdir(raw_path) 
                                if os.path.isdir(os.path.join(raw_path, folder))]
                    
                    # Add to our dictionary
                    main_folders[experiment_path] = samples

            # We don't need to go deeper than the first level
            # This condition checks if the current directory is more than one level deep from the starting point.
            if depth > 1:
                # If the condition is true, this line clears the folders list.
                del folders[:]
        
        # Write the list of main_folders path 
        write_main_folders_dict(input_file, main_folders)

        
    else: 
        # If there is a need for sample sheet and primer sheet to make StaggerFileAll
        if staggerfileall == "NA" : 
            # Construct the path to the main reorganization v1 when the data has sample sheet and primers
            Reorg1 = os.path.join(pathScript, "Step0_1reorganization/Reorganizationv1.py")

            # Prepare the command to run the reorganization script
            command_reorg = ["python3", Reorg1, path, "--sample_sheet", sample_sheet, "--primers_sheet",  primers_sheet, "--scriptpath", pathScript]
            # If StaggerFileAll already exists. 
        # elif sample_sheet == "NA" and primers_sheet == "NA": 
        
        else:
            # Construct the path to the main reorganization script when the data already has a staggerfileall
            Reorg2 = os.path.join(pathScript, "Step0_1reorganization/Reorganizationv2.py")

            # Prepare the command to run the reorganization script
            command_reorg = ["python3", Reorg2, path, "--staggerfileall", staggerfileall, "--scriptpath", pathScript, "-exp", pathExperiment]

        # Execute the main reorganization script
        print(f"Executing command: {' '.join(command_reorg)}")
        run_command(command_reorg)


    # Start running the bash script that will initiate the different bash scripts in each experiment folders 

    # Prepare and run the make_bash.py script
    pathmake_bash = os.path.join(pathScript, "Step0_2makebash", "make_bash.py")
    command_makebash = [
        "python3", pathmake_bash, 
        pathExperiment,
        run_bash,
        "-d", depth,
        "-q", quest,
        "-a", account,
        "-p", partition,
        "-n", nodes,
        "-nt", ntasks,
        "-t", time,
        "-m", memory,
        "-e", email,
        "-mt", mail_type
    ]

    # Run the make_bash.py script
    print(f"Executing command: {' '.join(command_makebash)}")
    run_command(command_makebash)


