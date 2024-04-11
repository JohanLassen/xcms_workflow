import os
import re
import pandas as pd
import argparse
import yaml
import glob


def filter_strings(strings, substring):
    return [s for s in strings if substring in s]


def list_files_recursive(path):
    """Lists all files in a directory and its subdirectories recursively."""

    files = glob.glob('{path}**/*.mzML'.format(path=path), recursive=True)
    return files


if __name__ == "__main__":


    #argparser = argparse.ArgumentParser(description='Process some integers.')
    #argparser.add_argument('--datapath', type=str, help='The path to the parent directory of all mzML files. Subdirectories are allowed.')
    #argparser.add_argument('--outputpath', type=str, help='The path to the output directory. This directory ends up containing all results and temporary files.')
    #argparser.add_argument('--sample_overview', type=str, default=None, help='The path to the sample overview file. Must be a csv file.')
    #args = argparser.parse_args()

    # Import settings
    with open("/faststorage/project/forensics/04johan/katrine_ghb/settings.yaml") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    # Define arguments
    args = argparse.Namespace()
    args.datapath = config["general"]["input_path"]
    args.outputpath = config["general"]["output_path"]
    args.sample_overview = config["general"]["sample_overview_path"]

    print(args.datapath)
    ### Identify files
    #all_files = list_files_recursive(args.datapath)
    mzml_files = glob.glob('{path}**/*.mzML'.format(path=args.datapath), recursive=True)# filter_strings(all_files, "mzML")
    ### Generate directory structure
    tmp_folder         = args.outputpath + "tmp/"
    peaks_folder       = tmp_folder + "peak_picked/"
    integration_folder = tmp_folder + "peak_integrated/"
    results_folder     = args.outputpath + "results/"
    
    ### Create folders
    os.makedirs(tmp_folder, exist_ok=True)
    os.makedirs(peaks_folder, exist_ok=True)
    os.makedirs(integration_folder, exist_ok=True)
    os.makedirs(results_folder, exist_ok=True)

    ### Make a pandas dataframe to orchestrate all steps and their inputs and outputs
    df = pd.DataFrame(mzml_files, columns=["raw_files"])

    # Add sample overview path
    df["sample_overview"] = args.sample_overview 

    ### Generate output file names for each step
    # step 1: peak picking
    outputs = [re.sub(".*{path}".format(path = args.datapath), "", filename) for filename in mzml_files]
    outputs = [re.sub("/", "_", filename) for filename in outputs]
    outputs = [peaks_folder+re.sub("mzML", "rds", filename) for filename in outputs]
    df["peak_picked"] = outputs

    # step 2: file merging
    output_merged = tmp_folder + "merged.rds"
    df["merged"] = output_merged

    # step 3: peak grouping
    output_grouped = tmp_folder + "grouped1.rds"
    df["grouped1"] = output_grouped

    # step 4: peak alignment
    output_aligned = tmp_folder + "aligned.rds"
    df["aligned"] = output_aligned

    # step 5: peak grouping
    output_grouped = tmp_folder + "grouped2.rds"
    df["grouped2"] = output_grouped

    # step 6: peak integration step 1
    output_integrated = integration_folder + "integrated1.rds"
    df["integrated1"] = output_integrated

    # step 7: peak integration step 2
    output_integrated2 = [re.sub(peaks_folder, integration_folder, x) for x in df.peak_picked]
    df["integrated2"] = output_integrated2

    # step 8: peak integration step 3
    output_integration3 = tmp_folder+"peak_integration3_dataset.rds"
    df["integrated3"] = output_integration3

    # step 9: output feature table
    output_feature_table = results_folder+"feature_table.csv"
    df["output_csv"] = output_feature_table

    # step 10: CAMERA
    camera = results_folder+"diffrep.csv"
    df["camera"] = camera

    # Save dataframe to csv
    df.to_csv("run_schedule_{run_id}.csv".format(run_id=config["general"]["run_id"]), index=False)


    