"""
FILE NAME:      02_sort_files_for_submission.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        03_generate_DivE_scripts.py; DivE_Calculations_new.R; subsets_parameters_for_dive.csv;
                cd4_cd8_parameters_for_dive.csv; 01_Reformat_for_DivE.py

DESCRIPTION:    This script will sort all of the donors' files and place them into folders respective of their timepoint,
                cell type, and cell subset. It will also copy over all of the scripts necessary to run DivE.

INPUT
FILES(S):       1) Reformatted File: a 2 column tab-delineated file where the first column is the TCR (V-CDR3 AA-J)
                and the 2nd column is UMI counts.
                2) Sample Sheet: A CSV file that lists the parameters for every sample needed for DivE.

OUTPUT
FILES(S):       1) 50 Submission Scripts to run DivE

INPUT
ARGUMENT(S):    1) Input directory: Specify the directory that has the reformatted files.
                2) Sample Sheet: Specify the full path of the CSV file that lists the DivE parameters for
                the particular samples.



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
import re, sys, os, glob, argparse, shutil, time
from collections import namedtuple, defaultdict
import pandas as pd

file_info = namedtuple("file_info", ["donor", "timepoint",
                                     "cell_type"])


def make_output_directory(directory):
    """
    This function creates the directory if it doesn't already exit. Accounts for concurrent processes.
    """
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError as exc:
            if exc.errno != os.errno.EEXIST:
                raise
        time.sleep(1)


def parse_file_name(filename):
    """
    This function extracts the donor, timepoint, cell type, and cell subset from each file name.
    """
    match = re.search(
        r"g-P7M1S-(.*)-(\d{3,4})([AB])CD(4|8|48)(N|M|NM)-(.*)_reformatted.txt", filename)
    if match:
        donor, timepoint = match.group(2), match.group(3)
        cell_type = 'CD{}{}' .format(match.group(4),match.group(5))
        return donor, timepoint, cell_type
    else:
        print ("Error! File Name does not match conventions, please rename\
         file {}!".format(filename))
        sys.exit(1)


def write_submission_file(output_dir):
    submit_submission = os.path.join(output_dir,'Submit_DivE.sh')
    submission_file = open(submit_submission,'w')
    script = """#!/bin/bash

module load python/3.6
python 03_generate_divE_scripts.py ./ ./sample_sheet.csv
swarm -f Submit_DivE.swarm -g 40 --partition=quick --time=4:00:00 --gres=lscratch:2
module unload python/3.6
"""
    submission_file.write(script)
    submission_file.close()


def Main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_directory",help="Specify the location of the reformatted files.")
    parser.add_argument("sample_sheet",help = "Specify the location of the sample sheet")
    parser.add_argument('generate_dive_file',
                        help = "Specify the location of the Generate_DivE_scripts.py file")
    parser.add_argument('dive_r_script',
                        help = "Specify the location of the DivE_Calculations.R file")



    args = parser.parse_args()

    sorted_files = defaultdict(list)
    input_files = sorted(glob.glob(args.input_directory + '*reformatted.txt'))
    for file in input_files:
        donor, timepoint, cell_type = parse_file_name(file)
        file_breakdown = file_info(donor, timepoint, cell_type)
        sorted_files[file_breakdown].append(file)

    #Read Sample Sheet
    sample_csv = pd.read_csv(args.sample_sheet,index_col = 0,header = 0)
    for specific_file_type, files in sorted_files.items():
        output_directory_name = 'Donor_{}/timepoint-{}-cell-{}'.format(
            specific_file_type.donor, specific_file_type.timepoint,
            specific_file_type.cell_type)
        output_dir_name = os.path.join(
            args.input_directory, output_directory_name)

        make_output_directory(output_dir_name)

        generate_dive_scriptname = os.path.join(output_dir_name, os.path.basename(args.generate_dive_file))
        dive_r_script_scriptname = os.path.join(output_dir_name, os.path.basename(args.dive_r_script))

        #copy the files to the desired directory
        shutil.copy(args.generate_dive_file, generate_dive_scriptname)
        shutil.copy(args.dive_r_script, dive_r_script_scriptname)

        #Create Sample Sheet
        df_index = '{}{}{}'.format(specific_file_type.donor, specific_file_type.timepoint,
            specific_file_type.cell_type)
        rarefaction_iteration = 2000
        rarefaction_length = sample_csv.loc[df_index,'Lowest UMI for DivE']
        cell_count_projection = sample_csv.loc[df_index,'Lymphocyte Count']
        norm_cell_count_projection = sample_csv.loc[df_index, 'Normalized Lymphocyte Count']
        sample_name = df_index
        dive_params_df = pd.DataFrame([
            [sample_name,rarefaction_iteration,
             rarefaction_length,cell_count_projection,
             norm_cell_count_projection]],
        columns = ['Sample Name','Rarefaction Iteration',
                   'Maximum Rarefaction Length','Estimated Cell Count',
                   'Normalized Cell Count'])
        dive_params_name = os.path.join(output_dir_name,'sample_sheet.csv')
        dive_params_df.to_csv(dive_params_name)

        for file in files:
            base_filename = os.path.basename(file)
            new_file_path = os.path.join(output_dir_name, base_filename)
            #Move the reformatted files to the new destination
            shutil.move(file, new_file_path)

        #Write submission script
        write_submission_file(output_dir_name)


if __name__ == '__main__':
    Main()
