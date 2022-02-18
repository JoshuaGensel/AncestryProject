from logging import FATAL
from numpy.lib.function_base import append
import pyslim
import json
import tskit
import msprime
import tsinfer
import numpy as np
import getopt, sys
import re
import os

path = "D:/Daten/programming_projects/AncestryProject/output/ts_raw/Tdiff_100_SB_50_T_200_run_641774310.trees"
nSamples = 10
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..','..'))

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


def haplotypeCalc(path):

    head, tail = os.path.split(path)
    inputFileName = tail
    #Importing and simplifying tree sequence

    ts_raw = pyslim.load(path)
    all_nodes = ts_raw.samples()
    nodes_F = []
    nodes_M = []

    for u in all_nodes:
        if ts_raw.mutation_at(u,0) != -1 and ts_raw.individual(ts_raw.node(u).individual).metadata["sex"] == pyslim.INDIVIDUAL_TYPE_FEMALE:
            nodes_F.append(u)
        if ts_raw.mutation_at(u,916567) != -1 and ts_raw.individual(ts_raw.node(u).individual).metadata["sex"] == pyslim.INDIVIDUAL_TYPE_MALE:
            nodes_M.append(u)
        else:
            pass
        
    
    ts_M_raw = ts_raw.simplify(samples=nodes_F)
    ts_Y_raw = ts_raw.simplify(samples=nodes_M)

    if(ts_M_raw.num_trees != 1): raise ValueError("more than one tree!")
    if(ts_Y_raw.num_trees != 1): raise ValueError("more than one tree!")

    if(ts_M_raw.first().num_roots > 1): raise ValueError("more than one root!")
    if(ts_Y_raw.first().num_roots > 1): raise ValueError("more than one root!")

    nodes_F_p1 = []
    nodes_F_p2 = []
    nodes_F_p3 = []
    for i in ts_M_raw.individuals_alive_at(0):
        if ts_M_raw.individual(i).metadata["subpopulation"] == 1:
            nodes_F_p1.extend(ts_M_raw.individual(i).nodes)
        elif ts_M_raw.individual(i).metadata["subpopulation"] == 2:
            nodes_F_p2.extend(ts_M_raw.individual(i).nodes)
        elif ts_M_raw.individual(i).metadata["subpopulation"] == 3:
            nodes_F_p3.extend(ts_M_raw.individual(i).nodes)
        else:
            raise ValueError("Individual without population!")
    keep_nodes_F = np.concatenate((np.random.choice(nodes_F_p1, nSamples, replace=False), 
                                   np.random.choice(nodes_F_p2, nSamples, replace=False), 
                                   np.random.choice(nodes_F_p3, nSamples, replace=False)))
    sts_M = ts_M_raw.simplify(keep_nodes_F)

    nodes_M_p1 = []
    nodes_M_p2 = []
    nodes_M_p3 = []
    for i in ts_Y_raw.individuals_alive_at(0):
        if ts_Y_raw.individual(i).metadata["subpopulation"] == 1:
            nodes_M_p1.extend(ts_Y_raw.individual(i).nodes)
        elif ts_Y_raw.individual(i).metadata["subpopulation"] == 2:
            nodes_M_p2.extend(ts_Y_raw.individual(i).nodes)
        elif ts_Y_raw.individual(i).metadata["subpopulation"] == 3:
            nodes_M_p3.extend(ts_Y_raw.individual(i).nodes)
        else:
            raise ValueError("Individual without population!")
    keep_nodes_M = np.concatenate((np.random.choice(nodes_M_p1, nSamples, replace=False), 
                                   np.random.choice(nodes_M_p2, nSamples, replace=False), 
                                   np.random.choice(nodes_M_p3, nSamples, replace=False)))
    sts_Y = ts_Y_raw.simplify(keep_nodes_M)

    mts_M = msprime.mutate(sts_M, rate=6.85*10**-7, random_seed=1, keep=False)
    mts_Y = msprime.mutate(sts_Y, rate=3.01*10**-8, random_seed=1, keep=False)

    mts_M = mts_M.delete_sites([i.id for i in mts_M.sites() if i.position < 900000])  
    mts_Y = mts_Y.delete_sites([i.id for i in mts_Y.sites() if i.position >= 900000])
    
    mts_M.dump(os.path.join(ROOT_DIR, 'output', 'ts_parsed', inputFileName.replace('.trees', '_mtDNA.trees')))
    mts_Y.dump(os.path.join(ROOT_DIR, 'output', 'ts_parsed', inputFileName.replace('.trees', '_YChrom.trees')))
    
    def convert_to_json_metadata(ts):
        # make a new ts with json metadata
        new_tables = ts.dump_tables()
        # iterate through (nearly) all the tables
        for table_name, table in new_tables.name_map.items():
            # provenance table doesn't have metadata
            if table_name not in ["provenances"]:
                metadata = []
                for row in table:
                    try:
                        row_metadata = row.metadata or {}
                        metadata.append(json.dumps(row_metadata).encode())
                    except TypeError:
                        raise TypeError(f"Can't convert {row.metadata} to JSON")
                # packset_metadata doesn't validate, so dump json in here and switch schema after
                table.packset_metadata(metadata)
                table.metadata_schema = tskit.MetadataSchema({'codec': 'json'})
        # May also need to convert top level metadata?
        return new_tables.tree_sequence()

    sample_ts_M = convert_to_json_metadata(mts_M)
    sample_ts_Y = convert_to_json_metadata(mts_Y)
    
    with tsinfer.SampleData(
        path=os.path.join(ROOT_DIR, 'output', 'tsinfer_samples', inputFileName.replace('.trees', '_mtDNA.samples')), 
        sequence_length=sample_ts_M.sequence_length, 
        num_flush_threads=2
    ) as sample_data_M:
        for var in sample_ts_M.variants():
            sample_data_M.add_site(var.site.position, var.genotypes, var.alleles)
    its_M = tsinfer.infer(sample_data= sample_data_M)
    print(its_M.draw_text())
    
    with tsinfer.SampleData(
        path=os.path.join(ROOT_DIR, 'output', 'tsinfer_samples', inputFileName.replace('.trees', '_YChrom.samples')), 
        sequence_length=sample_ts_Y.sequence_length, 
        num_flush_threads=2
    ) as sample_data_Y:
        for var in sample_ts_Y.variants():
            sample_data_Y.add_site(var.site.position, var.genotypes, var.alleles)
    its_Y = tsinfer.infer(sample_data= sample_data_Y)
    print(its_Y.draw_text())
    
    its_M.dump(os.path.join(ROOT_DIR, 'output', 'ts_inferred', inputFileName.replace('.trees', '_mtDNA.trees')))
    its_Y.dump(os.path.join(ROOT_DIR, 'output', 'ts_inferred', inputFileName.replace('.trees', '_YChrom.trees')))
    
    sample_data_Y.close()
    sample_data_M.close()
    os.remove(os.path.join(ROOT_DIR, 'output', 'tsinfer_samples', inputFileName.replace('.trees', '_mtDNA.samples')))
    os.remove(os.path.join(ROOT_DIR, 'output', 'tsinfer_samples', inputFileName.replace('.trees', '_YChrom.samples')))
    
 
    output = open(os.path.join(ROOT_DIR, 'output', 'haplotypes', inputFileName.replace('.trees', '_mtDNA.txt')), "w")
    for h,v in zip(mts_M.haplotypes(), mts_M.samples()):
        pop=sts_M.individual(v).metadata["subpopulation"]
        output.write(f"Pop{pop}_Ind{v}\t"+h+"\n")
    output.close()

    output = open(os.path.join(ROOT_DIR, 'output', 'haplotypes', inputFileName.replace('.trees', '_YChrom.txt')), "w")
    for h,v in zip(mts_Y.haplotypes(), mts_Y.samples()):
        pop=sts_Y.individual(v).metadata["subpopulation"]
        output.write(f"Pop{pop}_Ind{v}\t"+h+"\n")
    output.close()
    
    
haplotypeCalc(path)