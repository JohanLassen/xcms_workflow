from gwf import Workflow, AnonymousTarget
import os
import re
import pandas as pd
from subprocess import run
import yaml


gwf = Workflow()


### Helper function
def runR(rcall, inputfile, outputfile, extra=None, cores = 1, memory="24g", walltime="01:00:00"):
    inputs = [inputfile]
    outputs = [outputfile]
    options = {
        "nodes": 1,
        'cores': cores,
        'memory': memory,
        "walltime":walltime,
        "account" : "forensics"
    }
    spec = '''
    echo jobinfo $SLURM_JOBID
    {rcall}
    '''.format(rcall = rcall)
    return AnonymousTarget(inputs=inputs, outputs=outputs, spec=spec, options=options)

# Load yaml settings
with open("settings.yaml") as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
scheduler_file = "run_schedule_{run_id}.csv".format(run_id=config["general"]["run_id"])

### Check if run_schedule.csv exists and let user know if not
if not os.path.exists(scheduler_file):
    run(["python", "./scripts/configurate_workflow.py"], check=True)

### Read in run_schedule.csv
targets = pd.read_csv(scheduler_file)


#################### Step 1: Peak Picking ####################
for index in range(len(targets.raw_files)):
    gwf.target_from_template(
        "peak_picking"+str(index),
        runR(
            rcall="Rscript ./scripts/step1_peakcalling.R {inputfile} {outputfile}".format(inputfile=targets.raw_files[index], outputfile=targets.peak_picked[index]),
            inputfile=targets.raw_files[index], 
            outputfile=targets.peak_picked[index])
    )

# #################### Step 2: Merging  ####################
gwf.target_from_template(
    "merging",
    runR(rcall="Rscript ./scripts/step2_merging.R {scheduler_file} {outputfile}".format(scheduler_file=scheduler_file, outputfile=targets.merged[0]),
        inputfile=targets.peak_picked, 
        outputfile=targets.merged[0])
)

# #################### Step 3: Grouping ####################  
gwf.target_from_template(
    "grouping1",
    runR(rcall="Rscript ./scripts/step3_grouping1.R {inputfile} {outputfile}".format(inputfile=targets.merged[0], outputfile=targets.grouped1[0]),
        inputfile=targets.merged[0], 
        outputfile=targets.grouped1[0])
)

# #################### Step 4: Alignment ####################
gwf.target_from_template(
    "alignment",
    runR(rcall="Rscript ./scripts/step4_alignment.R {inputfile} {outputfile}".format(inputfile=targets.grouped1[0], outputfile=targets.aligned[0]),
        inputfile=targets.grouped1[0], 
        outputfile=targets.aligned[0])
)

# #################### Step 5: Grouping2 ####################
gwf.target_from_template(
    "grouping2",
    runR(rcall="Rscript ./scripts/step5_grouping2.R {inputfile} {outputfile}".format(inputfile=targets.aligned[0], outputfile=targets.grouped2[0]),
        inputfile=targets.aligned[0], 
        outputfile=targets.grouped2[0])
)

# #################### Step 6: Integration1 ####################
gwf.target_from_template(
    "integration1",
    runR(rcall="Rscript ./scripts/step6_integration1.R {inputfile} {outputfile}".format(inputfile=targets.grouped2[0], outputfile=targets.integrated1[0]),
        inputfile=targets.grouped2[0], 
        outputfile=targets.integrated1[0])
)

# #################### Step 7: Peak integration 2 ####################
for index in range(len(targets.integrated2)):
    gwf.target_from_template(
        "peak_integration2_"+str(index),
        runR(rcall="Rscript ./scripts/step7_integration2.R {inputfile} {outputfile} {extra}".format(inputfile=targets.integrated1[0], outputfile=targets.integrated2[index], extra=index+1),
            inputfile=targets.integrated1[0], 
            outputfile=targets.integrated2[index], 
            memory = "56g") # R uses indexing from 1.
    )

# #################### Step 8: Integration3 ####################
gwf.target_from_template(
    "peak_integration3",
    runR(rcall="Rscript ./scripts/step8_integration3.R {inputfile} {outputfile} {extra}".format(inputfile=scheduler_file, outputfile=targets.integrated3[0], extra=targets.integrated1[0]),
        inputfile=targets.integrated2, 
        outputfile=targets.integrated3[0],
        extra = targets.integrated1[0]
        )
)

# #################### Step 9: Generate names for features in output ####################
gwf.target_from_template(
    "output_feature_names_pos",
    runR(rcall="Rscript ./scripts/step9_writeoutput.R {inputfile} {outputfile}".format(inputfile=targets.integrated3[0], outputfile=targets.output_csv[0]),
        inputfile=targets.integrated3[0], 
        outputfile=targets.output_csv[0]
        )
)

# #################### Step 10: Generate CAMERA output ####################

gwf.target_from_template(
    "isotope_annotation",
    runR(rcall="Rscript ./scripts/step10_adduct_annotation.R {inputfile} {outputfile}".format(inputfile=targets.integrated3[0], outputfile=targets.camera[0]),
        inputfile=targets.integrated3[0], 
        outputfile=targets.camera[0],
        memory="32g",
        walltime="10:00:00"
        )
)
