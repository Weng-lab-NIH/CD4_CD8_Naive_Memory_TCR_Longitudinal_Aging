"""
FILE NAME:      03_generate_DivE_scripts.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        02_sort_files_for_submission.py; DivE_Calculations_new.R; subsets_parameters_for_dive.csv;
                cd4_cd8_parameters_for_dive.csv; 01_Reformat_for_DivE.py

DESCRIPTION:    This script will generate the submission scripts needed to run DivE as well as the parameters needed
                for DivE Calculations.

INPUT
FILES(S):       1) Reformatted File: a 2 column tab-delineated file where the first column is the TCR (V-CDR3 AA-J)
                and the 2nd column is UMI counts.
                2) Sample Sheet: A CSV file that lists the parameters for every sample needed for DivE.


OUTPUT
FILES(S):       1) None

INPUT
ARGUMENT(S):    1) Input directory: Specify the directory that has the reformatted files.
                2) Sample Sheet: Specify the full path of the CSV file that lists the DivE parameters for every sample
                3) Generate DivE Script: Specify the full path location of the 03_generate_DivE_scripts.py script.
                4) DivE Calculations Script: Specify the full path location of the DivE_Calculations_new.R script



CREATED:        25Aug19

MODIFICATIONS
LOG:



LAST MODIFIED
BY:             Thomas Nguyen

PYTHON
VERSION USED
TO WRITE
SCRIPT:         3.6

VERSION:        1.0

AUTHOR(S):      Thomas Nguyen

STATUS:         Working

TO DO LIST:     None at this moment.

DISCLAIMER(S):

NOTE(S):        for Biowulf (NIH) Implementation only.
"""

import glob, re, os, argparse
import pandas as pd

def make_dir(directory):
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError as exc:
            if exc.errno != os.errno.EEXIST:
                raise
            pass


def obtain_DivE_parameters(parameters_csv,sample_name):
    """
    This function will extract the parameters needed to run DivE from the CSV file.
    """

    param_df = pd.read_csv(parameters_csv,index_col = 1,header = 0)
    rar_iter = param_df.loc[sample_name,'Rarefaction Iteration']
    max_rar_length = param_df.loc[sample_name, 'Maximum Rarefaction Length']
    cell_count = param_df.loc[sample_name, 'Estimated Cell Count']
    norm_cell_count = param_df.loc[sample_name, 'Normalized Cell Count']
    return rar_iter,max_rar_length,cell_count,norm_cell_count


def write_submission_script(current_dir,reformatted_file,rar_iteration,rar_length,cell_count,norm_cell_count):
    """
    This function will write each submission script.
    """
    #Choose which models to use
    models_to_not_use = [7,17,26,35,40,44,48,50]
    models_to_use = [x for x in range(1,59) if x not in models_to_not_use]

    #Find Sample Name

    match = re.search(r"g-P7M1S-(.*)-(.*)_reformatted.txt",reformatted_file)
    barcode,sample_name = match.group(1),match.group(2)
    for model in models_to_use:
        script_name = "Calculate_DivE-{}-{}".format(sample_name,model)
        script_file = open(script_name,'w')
        script = """#!/bin/bash
#SBATCH -J {}-{}-{}
module load R/3.5.0
Rscript DivE_Calculations_new.R {} {} {} {} {} {} {} {} > ./{}/{}-{}-{}.out
module unload R/3.5.0
""".format(barcode,sample_name,model,current_dir,reformatted_file,sample_name,
           rar_iteration, rar_length, cell_count, norm_cell_count,model,
           "Print_Summaries",barcode,sample_name,model)
        script_file.write(script)
        script_file.close()

def write_swarm_file(current_directory):
    """
    This function will write the swarm commands to submit the swarm file.
    """
    submission_files = sorted(glob.glob(current_directory + 'Calculate_DivE*'))
    swarm_filename = 'Submit_DivE.swarm'
    swarm_file = open(swarm_filename,'w')
    for file in submission_files:
        print(file)
        swarm_command = 'bash {}\n'.format(file)
        swarm_file.write(swarm_command)
    swarm_file.close()
    print("Number of Commands in Swarm File: {}".format(len(submission_files)))
    print("Name of Swarm File to submit: {}".format(swarm_filename))



def Main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_directory",
                        help="Specify the sample specific folder that houses the reformatted files, DivE Script, "
                             "and submission-related scripts")
    parser.add_argument("DivE_parameters",help = "Specify the CSV file that contains the DivE parameters")
    args = parser.parse_args()


    plot_dir = os.path.join(args.input_directory,'Plots')
    dive_summary_dir =  os.path.join(args.input_directory,'DivE_Summaries')
    dive_environment_dir = os.path.join(args.input_directory, 'R_objects')
    print_sum_dir = os.path.join(args.input_directory, 'Print_Summaries')

    make_dir(plot_dir)
    make_dir(dive_summary_dir)
    make_dir(dive_environment_dir)
    make_dir(print_sum_dir)


    all_files = sorted(glob.glob(args.input_directory + '*reformatted.txt'))
    for file in all_files:
        match = re.search(r"g-P7M1S-(.*)-(.*)-(alpha|beta)_reformatted.txt",file)
        sample_name = match.group(2)
        rar_iter, max_rar_length, cell_count, norm_cell_count = obtain_DivE_parameters(args.DivE_parameters,sample_name)
        write_submission_script(args.input_directory,file,rar_iter, max_rar_length, cell_count, norm_cell_count)
    write_swarm_file(args.input_directory)

if __name__ == '__main__':
    Main()

