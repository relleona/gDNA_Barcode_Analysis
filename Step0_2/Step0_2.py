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
from datetime import date, timedelta
import time as time_module  # Import with a different name to avoid conflicts
import os
import sys
import subprocess
from argparse import ArgumentParser


#Variables that needs to be changed 

input_dir = {
    #Input and output directories 
    #Path to Step0_2
    "pathScript": "/projects/b1042/GoyalLab/aleona/gDNA_Barcode_Analysis",
    #The path to the downloaded FASTQ files
    "path":"/projects/b1042/GoyalLab/aleona/YG10Xbarcode/Experiments",
    #The path to sample sheet, if doesn't exist put "NA"
    "sample_sheet":"/projects/b1042/GoyalLab/aleona/gDNA_Barcode_Analysis/Experiment/SAUpdated_GoyalPool3_SampleSheet_Apr2023.xlsx", #created or taken from ENSEMBLE/GENECODE 
    #The path to primer sheet, if doesn't need put "NA"
    "primers_sheet": "/projects/b1042/GoyalLab/aleona/gDNA_Barcode_Analysis/Experiment/20200811_10X_BCSeq_Primers.xlsx",
    # If stagger file all exists already without the need for Sample Sheet or Primers Sheet put the path, if not place it as "NA"
    "staggerfileall": "/projects/b1042/GoyalLab/aleona/YG10Xbarcode/Exp2/StaggerFileAll.csv", # Can be NA
    #The path to Experiments folders
    "exp": "/projects/b1042/GoyalLab/aleona/YG10Xbarcode/Exp2",
    #Depth of subfolders within the Experiment folder above to run the job scripts for
    "depth" : "1" #Default is 1 and usually its the same if the format is followed
}

bash_script ={
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

#Pause to check if the organization is the way you wanted. 
# If you put TRUE, it will stop 
# If you put 10, It will pause for 10 minutes after the job is submitted
Pause_to_check = True 

#Max wait time for job to run -> 10 minutes
max_wait_time = 60*60*10  

#####################################################################################################################################
#Nothing needs to be modified in the entire pipeline beyond this.
#####################################################################################################################################
###Run Command---------------------------
# #'
# # Command line parser
# parser = ArgumentParser(description="Wrapper script for FASTQ file reorganization.")
# parser.add_argument("pathScript", help="Specify the path containing the BaseSpaceReorganization.py, in this case its Step0_2", type=str)
# parser.add_argument("path", help="Specify the path containing the nested sample folders.", type=str)
# parser.add_argument("--sample_sheet", help="Path to the sample sheet Excel file", required=True)
# parser.add_argument("--primers_sheet", help="Path to the primers Excel file", required=True)
# parser.add_argument("-exp", "--experimentfile", help="Specify the output path. If not provided, a 'Experiment' folder will be created in the input path.", type=str, required=False)

# args = parser.parse_args()
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

def extract_job_ids(stdout):
    """
    Extract job IDs from the output of a job submission command.
    
    Args:
    stdout (str): The standard output from a job submission command.
    
    Returns:
    list: A list of strings representing the extracted job IDs.
    """
    return [line.split()[-1] for line in stdout.split('\n') if "Submitted batch job" in line]

def monitor_jobs(job_ids, max_wait_time):
    """
    Monitor the status of submitted jobs and wait for their completion.
    
    Args:
    job_ids (list): A list of strings representing the IDs of jobs to monitor.
    max_wait_time (int): The maximum time to wait for job completion in seconds.
    
    Raises:
    SystemExit: If jobs don't complete within max_wait_time.
    """
    def are_jobs_running(ids):
        result = subprocess.run(['squeue', '-j', ','.join(ids)], capture_output=True, text=True)
        return any(id in result.stdout for id in ids)

    print(f"Waiting for jobs to complete: {', '.join(job_ids)}")
    start_time = time_module.time()
    check_interval = 60  # Start with 1 minute
    while are_jobs_running(job_ids):
        if time_module.time() - start_time > max_wait_time:
            print(f"Jobs didn't complete within {max_wait_time/3600:.2f} hours. Please check their status manually.")
            sys.exit(1)
        time_module.sleep(check_interval)
        check_interval = min(check_interval * 2, 1800)  # Double the interval, max 30 minutes
        elapsed_time = time_module.time() - start_time
        print(f"Jobs still running... (Elapsed time: {str(timedelta(seconds=int(elapsed_time)))})")

    if not are_jobs_running(job_ids):
        print("All jobs are done!")

def pause_for_minutes(minutes):
    """
    Pause the script execution for a specified number of minutes.
    
    Args:
    minutes (int): The number of minutes to pause.
    """
    seconds = minutes * 60
    for remaining in range(seconds, 0, -1):
        sys.stdout.write(f"\rPausing for {remaining} seconds...")
        sys.stdout.flush()
        time_module.sleep(1)
    print("\nPause completed!")

def main():
    """
    Main function to orchestrate the data processing pipeline.
    """
    # Assign variables from the configuration dictionaries
    path_script = input_dir['pathScript']
    path = input_dir['path']
    sample_sheet = input_dir['sample_sheet']
    primers_sheet = input_dir['primers_sheet']
    staggerfileall = input_dir['staggerfileall']
    exp = input_dir['exp']
    depth = input_dir['depth']

    quest = bash_script['quest']
    account = bash_script['account']
    partition = bash_script['partition']
    nodes = bash_script['nodes']
    ntasks = bash_script['ntasks']
    memory = bash_script['memory']
    time = bash_script['time']
    email = bash_script['email']
    mail_type = bash_script['mail_type']

    # # If there is a need for sample sheet and primer sheet to make StaggerFileAll
    # if staggerfileall == "NA" : 
    #     # Construct the path to the main reorganization v1 when the data has sample sheet and primers
    #     Reorg1 = os.path.join(path_script, "Step0_2/Reorganizationv1.py")

    #     # Prepare the command to run the reorganization script
    #     command_reorg = ["python3", Reorg1, path, "--sample_sheet", sample_sheet, "--primers_sheet",  primers_sheet, "--scriptpath", path_script]
    #  # If StaggerFileAll already exists. 
    # # elif sample_sheet == "NA" and primers_sheet == "NA": 
    # else:
    #     # Construct the path to the main reorganization script when the data already has a staggerfileall
    #     Reorg2 = os.path.join(path_script, "Step0_2/Reorganizationv2.py")

    #     # Prepare the command to run the reorganization script
    #     command_reorg = ["python3", Reorg2, path, "--staggerfileall", staggerfileall, "--scriptpath", path_script]

    
    # # Add the experiment folder  if it was provided
    # if exp != "NA":
    #     command_reorg.extend(["-exp", exp])

    # # Execute the main reorganization script
    # print(f"Executing command: {' '.join(command_reorg)}")
    # run_command(command_reorg)

    # # Check if we should pause or stop after reorganization
    # # This pause is so that we can edit the information in the templateInputFullRun.py
    # if Pause_to_check is True: 
    #     print("Actions paused, will not create a job and run them. Rerun the script to create a job and run each job")
    #     return
    # else: 
    #     # Pause for the specified time
    #     pause_for_minutes(Pause_to_check)

    # Prepare and run the make_bash.py script
    pathmake_bash = os.path.join(path_script, "Step0_2", "make_bash.py")
    command_makebash = [
        "python3", pathmake_bash, 
        exp,  
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

    # Run the make_bash.py script and extract job IDs
    stdout = run_command(command_makebash)
    job_ids_step1 = extract_job_ids(stdout)

    # Monitor jobs if any were submitted
    if job_ids_step1:
        print(f"Jobs submitted with IDs: {', '.join(job_ids_step1)}")
        monitor_jobs(job_ids_step1, max_wait_time)
        print("FASTQ reorganization process completed.")
    else:
        print("No jobs were submitted in Step 0_2. Check the script output for errors.")
        print(f"Script output:\n{stdout}")

if __name__ == "__main__":
    main()