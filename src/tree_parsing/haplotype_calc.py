from logging import FATAL
#from numpy.lib.financial import rate
from numpy.lib.function_base import append
import pyslim
import json
import tqdm
import tskit
import msprime
import tsinfer
import numpy as np
import pandas as pd
from pandas import DataFrame
import getopt, sys
import re

path = "D:\Daten\programming_projects\AncestryProject\output\other\Tdiff_100_SB_50_T_100200_run_3473623206.trees"
nSamples = 10

#Commandline options

full_cmd_arguments = sys.argv
argument_list = full_cmd_arguments[1:]
short_options = "hp:n:"
long_options = ["help", "path=", "nSamples="]
try:
    arguments, values = getopt.getopt(argument_list, short_options, long_options)
except getopt.error as err:
    print (str(err))
    sys.exit(2)
for current_argument, current_value in arguments:
    if current_argument in ("-h", "--help"):
        print("-p/--path specifies the path for the .trees-file\n-n/--nSamples defines the sampled number of individuals per subpop")
    elif current_argument in ("-p", "--path"):
        path = str(current_value)
    elif current_argument in ("-n", "--nSamples"):
        nSamples = int(current_value)


def haplotypeCalc(input):

    #Importing and simplifying tree sequence

    treeseq_raw = pyslim.load(input)
    all_nodes = treeseq_raw.samples()
    nodes_F = []
    nodes_M = []

    for u in all_nodes:
        if treeseq_raw.mutation_at(u,0) != -1 and treeseq_raw.individual(treeseq_raw.node(u).individual).metadata["sex"] == pyslim.INDIVIDUAL_TYPE_FEMALE:
            nodes_F.append(u)
        if treeseq_raw.mutation_at(u,916567) != -1 and treeseq_raw.individual(treeseq_raw.node(u).individual).metadata["sex"] == pyslim.INDIVIDUAL_TYPE_MALE:
            nodes_M.append(u)
        else:
            pass
        
    
    treeseq_mtDNA = treeseq_raw.simplify(samples=nodes_F)
    treeseq_YChrom = treeseq_raw.simplify(samples=nodes_M)

    if(treeseq_mtDNA.num_trees != 1): raise ValueError("more than one tree!")
    if(treeseq_YChrom.num_trees != 1): raise ValueError("more than one tree!")

    if(treeseq_mtDNA.first().num_roots > 1): raise ValueError("more than one root!")
    if(treeseq_YChrom.first().num_roots > 1): raise ValueError("more than one root!")

    nodes_F_p1 = []
    nodes_F_p2 = []
    nodes_F_p3 = []
    for i in treeseq_mtDNA.individuals_alive_at(0):
        if treeseq_mtDNA.individual(i).metadata["subpopulation"] == 1:
            nodes_F_p1.extend(treeseq_mtDNA.individual(i).nodes)
        elif treeseq_mtDNA.individual(i).metadata["subpopulation"] == 2:
            nodes_F_p2.extend(treeseq_mtDNA.individual(i).nodes)
        elif treeseq_mtDNA.individual(i).metadata["subpopulation"] == 3:
            nodes_F_p3.extend(treeseq_mtDNA.individual(i).nodes)
        else:
            raise ValueError("Individual without population!")
    keep_nodes_F = np.concatenate((np.random.choice(nodes_F_p1, nSamples, replace=False), np.random.choice(nodes_F_p2, nSamples, replace=False), np.random.choice(nodes_F_p3, nSamples, replace=False)))
    sts_mtDNA = treeseq_mtDNA.simplify(keep_nodes_F)

    nodes_M_p1 = []
    nodes_M_p2 = []
    nodes_M_p3 = []
    for i in treeseq_YChrom.individuals_alive_at(0):
        if treeseq_YChrom.individual(i).metadata["subpopulation"] == 1:
            nodes_M_p1.extend(treeseq_YChrom.individual(i).nodes)
        elif treeseq_YChrom.individual(i).metadata["subpopulation"] == 2:
            nodes_M_p2.extend(treeseq_YChrom.individual(i).nodes)
        elif treeseq_YChrom.individual(i).metadata["subpopulation"] == 3:
            nodes_M_p3.extend(treeseq_YChrom.individual(i).nodes)
        else:
            raise ValueError("Individual without population!")
    keep_nodes_M = np.concatenate((np.random.choice(nodes_M_p1, nSamples, replace=False), np.random.choice(nodes_M_p2, nSamples, replace=False), np.random.choice(nodes_M_p3, nSamples, replace=False)))
    sts_YChrom = treeseq_YChrom.simplify(keep_nodes_M)

    mts_mtDNA = msprime.mutate(sts_mtDNA, rate=6.85*10**-7, random_seed=1, keep=False)
    mts_YChrom = msprime.mutate(sts_YChrom, rate=3.01*10**-8, random_seed=1, keep=False)

    mts_mtDNA = mts_mtDNA.delete_sites([i.id for i in mts_mtDNA.sites() if i.position < 900000])  
    mts_YChrom = mts_YChrom.delete_sites([i.id for i in mts_YChrom.sites() if i.position >= 900000])
    
    def convert_to_json_metadata(ts):
    # make a new ts with json metadata
        new_tables = ts.dump_tables()
        # iterate through (nearly) all the tables
        for table_name, table in new_tables.name_map.items():
            # provenance table doesn't have metadata
            if table_name not in ["provenances"]:
                # packset_metadata doesn't validate, so dump json in here and switch schema after
                table.packset_metadata([json.dumps(row.metadata).encode() for row in table])
                table.metadata_schema = tskit.MetadataSchema({'codec': 'json'})
        # May also need to convert top level metadata?
        return new_tables.tree_sequence()

    mtDNA_samplets = convert_to_json_metadata(mts_mtDNA)
    
    mtDNA_sample_data = tsinfer.SampleData.from_tree_sequence(mtDNA_samplets, path = "mtDNAtest.samples")
    inferred_tree_mtDNA = tsinfer.infer(sample_data= mtDNA_sample_data)
 
    print("SLiM-Tree mtDNA:")
    print(mts_mtDNA.draw_text())
    print("Inferred tsinfer-Tree mtDNA:")
    print(inferred_tree_mtDNA.draw_text())
    
    # outputFile = inputFile.replace(".trees", "_mtDNA.txt")
    # output = open(outputFile, "w")
    # for h,v in zip(mts_mtDNA.haplotypes(), mts_mtDNA.samples()):
    #     pop=sts_mtDNA.individual(v).metadata["subpopulation"]
    #     output.write(f"Pop{pop}_Ind{v}\t"+h+"\n")
    # output.close()

    # outputFile = inputFile.replace(".trees", "_YChrom.txt")
    # output = open(outputFile, "w")
    # for h,v in zip(mts_YChrom.haplotypes(), mts_YChrom.samples()):
    #     pop=sts_YChrom.individual(v).metadata["subpopulation"]
    #     output.write(f"Pop{pop}_Ind{v}\t"+h+"\n")
    # output.close()
    


haplotypeCalc(path)
