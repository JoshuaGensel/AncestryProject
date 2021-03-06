import subprocess
import os
from tree_processing import tree_processing
from tree_construction import tree_construction
from tree_analysis import tree_analysis
from datetime import datetime
from random import *
import getopt,sys

n_runs = None

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))

def main():
    
    start_time = datetime.now()
    Td_values = [100,300,1000]
    N_sample_values = [10, 20]
    seeds = [randint(10**8,10**9-1), randint(0**8,10**9-1),randint(0**8,10**9-1)]
    file_name_starts = []
    for i in seeds:
        file_name_starts.append(f"ID_{i}")

    slim_drift = f"-d txt_directory='{os.path.join(ROOT_DIR, 'data', 'drift/')}'".replace('\\','/')
    slim_burn_in = f"-d burnin_file='{os.path.join(ROOT_DIR, 'data', 'burn_in', 'burnin.trees')}'".replace('\\','/')
    slim_ts_raw = f"-d trees_directory='{os.path.join(ROOT_DIR, 'data', 'ts_raw/')}'".replace('\\','/')


    print("\t Running SLiM_model ...")
    for s,Td in zip(seeds,Td_values):
        command = f"slim -s {s} -d Td={Td} {slim_drift} {slim_burn_in} {slim_ts_raw} SLiM_model.slim"
        subprocess.run(command.split(),stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    model_time = datetime.now()
    print(f"\t \t model runs finished after: {model_time-start_time}")

    print("\t Running tree_processing.py ...")
    ts_raw = os.fsencode(os.path.join(ROOT_DIR, 'data', 'ts_raw'))
    for file in os.listdir(ts_raw):
        filename = os.fsdecode(file)
        for start in file_name_starts:
            if filename.startswith(start):
                for n in N_sample_values:
                    tree_processing(filename,n)
            else:
                continue
    processing_time = datetime.now()
    print(f"\t \t tree processing finished after: {processing_time-model_time}")

    print("\t Running tree_construction.py ...")
    fasta = os.fsencode(os.path.join(ROOT_DIR, 'data', 'fasta'))
    for file in os.listdir(fasta):
        filename = os.fsdecode(file)
        for start in file_name_starts:
            if filename.startswith(start): 
                tree_construction(filename)
            else:
                continue
    construction_time = datetime.now()
    print(f"\t \t tree construction finished after: {construction_time-processing_time}")

    print("\t Running tree_analysis.py ...")
    parsed_trees = os.fsencode(os.path.join(ROOT_DIR, 'data', 'tree_genealogy'))
    for file in os.listdir(parsed_trees):
        filename = os.fsdecode(file)
        for start in file_name_starts:
            if filename.startswith(start): 
                tree_analysis(filename)
            else:
                continue
    analysis_time = datetime.now()
    print(f"\t \t tree analysis finished after: {analysis_time-construction_time}")
        
    end_time = datetime.now()
    print(f'\t Finished runs in: {end_time - start_time}')

full_cmd_arguments = sys.argv
argument_list = full_cmd_arguments[1:]
short_options = "hn:"
long_options = ["help", "n_runs="]
try:
    arguments, values = getopt.getopt(argument_list, short_options, long_options)
except getopt.error as err:
    print (str(err))
    sys.exit(2)
for current_argument, current_value in arguments:
    if current_argument in ("-h", "--help"):
        print("-n/--n_runs defines how many full parameter runs you want to execute (will run on one core)")
    elif current_argument in ("-p", "--path"):
        path = str(current_value)
    elif current_argument in ("-n", "--n_runs"):
        n_runs = int(current_value)



if n_runs != None:
    print(f"Starting {n_runs} runs.")
    for _ in range(n_runs):
        main()
    print(f"Finished all {n_runs} runs :)")