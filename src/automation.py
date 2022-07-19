import subprocess
import os
from tree_processing import tree_processing
from tree_construction import tree_construction
from tree_analysis import tree_analysis

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))

n_runs = 1
Td_values = [100,300,1000]
N_sample_values = [10,20,30]

slim_drift = f"-d txt_directory={os.path.join(ROOT_DIR, 'data', 'drift')}"
slim_burn_in = f"-d burnin_file={os.path.join(ROOT_DIR, 'data', 'burn_in')}"
slim_ts_raw = f"-d trees_directory={os.path.join(ROOT_DIR, 'data', 'ts_raw')}"


# for run in range(n_runs):
#     for Td in Td_values:
#         command = f"slim -d Td={Td} {slim_drift} {slim_burn_in} {slim_ts_raw} SLiM_model.slim"
#         subprocess.run(command.split())

ts_raw = os.fsencode(os.path.join(ROOT_DIR, 'data', 'ts_raw'))

for file in os.listdir(ts_raw):
     filename = os.fsdecode(file)
     if filename.endswith(".trees"):
        for n in N_sample_values:
            tree_processing(filename,n)
     else:
         continue

fasta = os.fsencode(os.path.join(ROOT_DIR, 'data', 'fasta'))
 
for file in os.listdir(fasta):
     filename = os.fsdecode(file)
     if filename.endswith(".fasta"): 
        tree_construction(filename)
     else:
        continue
     
parsed_trees = os.fsencode(os.path.join(ROOT_DIR, 'data', 'tree_genealogy'))
 
for file in os.listdir(parsed_trees):
     filename = os.fsdecode(file)
     if filename.endswith(".tree"): 
         tree_analysis(filename)
     else:
         continue