from gwf import Workflow, AnonymousTarget
import os
import pandas as pd
from pathlib import Path
from subprocess import run
import yaml


def cxcms(gwf=None, config=None):
    
    if gwf is None:
        gwf = Workflow()

    MODULE_DIR = Path(__file__).parent.resolve()
    SCRIPTS_DIR = MODULE_DIR / "scripts"
    CONTAINER = Path(__file__).parent.parent.parent / "forensics_workflow.sif"
    XCMS_NEW_CONTAINER = MODULE_DIR / "xcms_new.sif"
    ### Helper function
    def runR(rcall, 
             inputfile, 
             outputfile, 
             extra=None, 
             cores = 1, 
             memory="4g", 
             walltime="00:30:00",
             container = None):
        inputs = inputfile if isinstance(inputfile, list) else [inputfile]
        outputs = outputfile if isinstance(outputfile, list) else [outputfile]
        
        cont = container if container is not None else CONTAINER
        # Wrap command with apptainer if container exists
        if cont.exists():
            rcall = f"apptainer exec --no-home {cont} {rcall}"
        else:
            print(f"WARNING: Container not found at {cont}, running without container")

        options = {
            "nodes": 1,
            'cores': cores,
            'memory': memory,
            "walltime":walltime,
            "account" : "forensics"
        }
        spec = '''
        echo jobinfo $SLURM_JOBID
        echo "=== DEBUG INFO ==="
        echo "Current working directory:"
        pwd
        echo "Module directory: {module_dir}"
        echo "Scripts directory: {scripts_dir}"
        echo "=================="
        cd {scripts_parent} && {rcall}
        '''.format(module_dir=MODULE_DIR, scripts_dir=SCRIPTS_DIR, scripts_parent=SCRIPTS_DIR.parent, rcall = rcall)
        return AnonymousTarget(inputs=inputs, outputs=outputs, spec=spec, options=options)

    # Load yaml settings if not provided
    if config is None:
        with open("settings.yaml") as file:
            config = yaml.load(file, Loader=yaml.FullLoader)
    scheduler_file = "{output_path}run_schedule_{run_id}.csv".format(output_path = config["general"]["output_path"], run_id=config["general"]["run_id"])
    
    # Get absolute path to settings.yaml for R scripts
    settings_yaml_path = Path("settings.yaml").resolve()

    ### Check if run_schedule.csv exists and let user know if not
    if not os.path.exists(scheduler_file):
        run(["python", str(SCRIPTS_DIR / "configurate_workflow.py")], check=True)

    ### Save settings for later use and documentation in the output folder
    yaml_save = "{output_path}settings.yaml".format(output_path = config["general"]["output_path"])
    with open(yaml_save, "w") as file:    
        yaml.dump(config, file)

    ### Read in run_schedule.csv
    targets = pd.read_csv(scheduler_file)

    # select first 10 rows of pandas dataframe for testing
    #targets = targets.head(10)


    #################### Step 1: Peak Picking ####################
    for index in range(len(targets.raw_files)):
        if os.path.exists(targets.peak_picked[index]):
            continue

        gwf.target_from_template(
            "peak_picking_{index}_{run_id}".format(index=index, run_id=config["general"]["run_id"]),
            runR(
                rcall="Rscript ./scripts/step1_peakcalling.R '{inputfile}' '{outputfile}' '{settings_path}'".format(inputfile=targets.raw_files[index], outputfile=targets.peak_picked[index], settings_path=settings_yaml_path),
                inputfile=targets.raw_files[index], 
                memory="12g",
                outputfile=targets.peak_picked[index])
        )

    bookkeeping = []

    # Make intermediate txt_completion files for every thousand set of jobs that are completed
    for index in range(0, len(targets.peak_picked), 1000):
        end_index = min(index + 1000, len(targets.peak_picked))

        # Save output file in intmediate folder in tmp folder
        out_folder = config["general"]["output_path"]+"bookkeeping/"
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        output_file = "{output_path}txt_completion_{run_id}_{index}.txt".format(
            output_path=out_folder, 
            run_id=config["general"]["run_id"], 
            index=index // 1000
        )
        
        # Create a target for the completion file
        gwf.target_from_template(
            "txt_completion_{run_id}_{index}".format(run_id=config["general"]["run_id"], index=end_index),
            AnonymousTarget(
                inputs=targets.peak_picked[index:end_index],
                outputs=[output_file],
                spec="echo 'Completed jobs from {start} to {end}' > {output_file}".format(
                    start=index + 1, end=end_index, output_file=output_file
                ),
                options={
                    "nodes": 1,
                    "cores": 1,
                    "memory": "1g",
                    "walltime": "00:05:00",
                    "account": "forensics_TOFscreenings"
                }
            )
        )
        bookkeeping.append(output_file)


    #################### Step 2: Merging ####################
    gwf.target_from_template(
        "merging_{run_id}".format(run_id=config["general"]["run_id"]),
        runR(rcall="Rscript ./scripts/step2_merging.R {scheduler_file} {outputfile}".format(scheduler_file=scheduler_file, outputfile=targets.merged[0], settings_path=settings_yaml_path),
            inputfile=bookkeeping, 
            memory="600g",
            walltime="09:00:00",
            outputfile=targets.merged[0])
    )

    # #################### Step 3, 4, 5: Grouping, Alignment, Grouping2 ####################
    gwf.target_from_template(
        "group_align_group2_{run_id}".format(run_id=config["general"]["run_id"]),
        runR(rcall="Rscript ./scripts/step2_merge_group_align.R {inputfile} {outputfile} '{settings_path}'".format(inputfile=targets.merged[0], outputfile=targets.grouped2[0], settings_path=settings_yaml_path),
            inputfile=targets.merged[0], 
            memory="350g",walltime="09:00:00",
            outputfile=targets.grouped2[0],
            container = XCMS_NEW_CONTAINER)
    )


    #################### Step 6: Integration1 ####################

    integration1_outputs = [targets.integrated2[index].replace(".rds", "_step1.rds") for index in range(len(targets.integrated2))]
    interval_length = 100
    intervals = [[i, i + interval_length] for i in range(0, len(integration1_outputs), interval_length)]
    for index, (start, end) in enumerate(intervals):
        
        x = "Rscript ./scripts/step6_integration1.R {inputfile} '{outputfile}' 'c({start}, {end})' '{settings_path}'".format(
            inputfile=targets.grouped2[0], 
            outputfile=" ".join(integration1_outputs[start:end]), 
            start=start+1, 
            end=end,
            settings_path=settings_yaml_path)
        
        # Create the target for each interval
        gwf.target_from_template(
            "integration1_{run_id}_{index}".format(run_id=config["general"]["run_id"], index=index),
            runR(rcall=x,
                inputfile=targets.grouped2[0], 
                memory="200g",
                cores=2,
                walltime="03:00:00",
                outputfile=integration1_outputs[start:end])  # This will create a list of output files for each target
        )

    step2_outputs = list()

    #################### Step 7: Peak integration 2 ####################
    for index in range(len(integration1_outputs)):
        inputs = integration1_outputs[index]
        step2_outputs.append(targets.integrated2[index])
        
        x = "Rscript ./scripts/step7_integration2.R '{inputfile}' '{outputfile}'".format(
            inputfile=inputs, 
            outputfile=targets.integrated2[index])
        
        gwf.target_from_template(
            "peak_integration2_{index}_{run_id}".format(index=index, run_id=config["general"]["run_id"]),
            runR(rcall=x,
                inputfile=inputs, 
                outputfile=targets.integrated2[index],
                cores=1, 
                walltime="00:20:00",
                memory = "10g")
        )

    integration2_bookkeeping = []
    # Make intermediate txt_completion files for every thousand set of integration2 jobs
    for index in range(0, len(step2_outputs), 1000):
        end_index = min(index + 1000, len(step2_outputs))

        # Save output file in intermediate folder
        out_folder = config["general"]["output_path"]+"bookkeeping/"
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        output_file = "{output_path}integration2_completion_{run_id}_{index}.txt".format(
            output_path=out_folder,
            run_id=config["general"]["run_id"],
            index=index // 1000
        )

        # Create a target for the completion file
        gwf.target_from_template(
            "integration2_completion_{run_id}_{index}".format(run_id=config["general"]["run_id"], index=index // 1000),
            AnonymousTarget(
                inputs=step2_outputs[index:end_index],
                outputs=[output_file],
                spec="echo 'Completed integration2 jobs from {start} to {end}' > {output_file}".format(
                    start=index + 1, end=end_index, output_file=output_file
                ),
                options={
                    "nodes": 1,
                    "cores": 1,
                    "memory": "1g",
                    "walltime": "00:05:00",
                    "account": "forensics_TOFscreenings"
                }
            )
        )
        integration2_bookkeeping.append(output_file)


    #################### Step 8: Integration3 ####################
    gwf.target_from_template(
        "peak_integration3_{run_id}".format(run_id=config["general"]["run_id"]),
        runR(rcall="Rscript ./scripts/step8_integration3.R {inputfile} {outputfile} {extra} {settings_path}".format(inputfile=scheduler_file, outputfile=targets.output_csv[0], extra=targets.grouped2[0], settings_path=settings_yaml_path),
            inputfile=integration2_bookkeeping,
            outputfile=[targets.output_csv[0], targets.object[0], 
                        targets.feature_list[0], targets.assay[0], 
                        targets.feature_definitions[0], targets.col_data[0]],
            walltime="20:00:00",
            memory="400g",
            container = Path("basenev")
            ))

    return gwf

