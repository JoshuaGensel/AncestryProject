from logging import FATAL
#from numpy.lib.financial import rate
from numpy.lib.function_base import append
import pyslim
from tqdm import tqdm
import tskit
import msprime
import tsinfer
import numpy as np
import pandas as pd
from pandas import DataFrame
import getopt, sys
import re

path = "D:\Daten\programming_projects\AncestryProject\output\other\Tdiff_100_SB_50_T_200_run_2252115279.trees"
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

    treeseq_raw = tskit.load(input)
    all_nodes = treeseq_raw.samples()
    nodes_F = []
    nodes_M = []

    for u in all_nodes:
        if treeseq_raw.mutation_at(u,0) != -1 and treeseq_raw.individual(treeseq_raw.node(u).individual).metadata["sex"] == pyslim.INDIVIDUAL_TYPE_FEMALE and treeseq_raw.node(u).time == 0:
            nodes_F.append(u)
        if treeseq_raw.mutation_at(u,916568) != -1 and treeseq_raw.individual(treeseq_raw.node(u).individual).metadata["sex"] == pyslim.INDIVIDUAL_TYPE_MALE and treeseq_raw.node(u).time == 0:
            nodes_M.append(u)
        else:
            pass

    treeseq_mtDNA = treeseq_raw.simplify(samples=nodes_F, keep_input_roots=True)
    treeseq_YChrom = treeseq_raw.simplify(samples=nodes_M, keep_input_roots=True)

    if(treeseq_mtDNA.num_trees != 1): raise ValueError("more than one tree!")
    if(treeseq_YChrom.num_trees != 1): raise ValueError("more than one tree!")

    rts_mtDNA = treeseq_mtDNA.recapitate(recombination_rate=0, Ne=5000, random_seed=1)
    rts_YChrom = treeseq_YChrom.recapitate(recombination_rate=0, Ne=5000, random_seed=1)

    if(rts_mtDNA.first().num_roots > 1): raise ValueError("more than one root!")
    if(rts_YChrom.first().num_roots > 1): raise ValueError("more than one root!")

    nodes_F_p1 = []
    nodes_F_p2 = []
    nodes_F_p3 = []
    for i in rts_mtDNA.individuals_alive_at(0):
        if rts_mtDNA.individual(i).metadata["subpopulation"] == 1:
            nodes_F_p1.extend(rts_mtDNA.individual(i).nodes)
        elif rts_mtDNA.individual(i).metadata["subpopulation"] == 2:
            nodes_F_p2.extend(rts_mtDNA.individual(i).nodes)
        elif rts_mtDNA.individual(i).metadata["subpopulation"] == 3:
            nodes_F_p3.extend(rts_mtDNA.individual(i).nodes)
        else:
            raise ValueError("Individual without population!")
    keep_nodes_F = np.concatenate((np.random.choice(nodes_F_p1, nSamples, replace=False), np.random.choice(nodes_F_p2, nSamples, replace=False), np.random.choice(nodes_F_p3, nSamples, replace=False)))
    sts_mtDNA = rts_mtDNA.simplify(keep_nodes_F)

    nodes_M_p1 = []
    nodes_M_p2 = []
    nodes_M_p3 = []
    for i in rts_YChrom.individuals_alive_at(0):
        if rts_YChrom.individual(i).metadata["subpopulation"] == 1:
            nodes_M_p1.extend(rts_YChrom.individual(i).nodes)
        elif rts_YChrom.individual(i).metadata["subpopulation"] == 2:
            nodes_M_p2.extend(rts_YChrom.individual(i).nodes)
        elif rts_YChrom.individual(i).metadata["subpopulation"] == 3:
            nodes_M_p3.extend(rts_YChrom.individual(i).nodes)
        else:
            raise ValueError("Individual without population!")
    keep_nodes_M = np.concatenate((np.random.choice(nodes_M_p1, nSamples, replace=False), np.random.choice(nodes_M_p2, nSamples, replace=False), np.random.choice(nodes_M_p3, nSamples, replace=False)))
    sts_YChrom = rts_YChrom.simplify(keep_nodes_M)

    mts_mtDNA = msprime.mutate(sts_mtDNA, rate=6.85*10**-7, random_seed=1, keep=False)
    mts_YChrom = msprime.mutate(sts_YChrom, rate=3.01*10**-8, random_seed=1, keep=False)

    mts_mtDNA = mts_mtDNA.delete_sites([i.id for i in mts_mtDNA.sites() if i.position < 900000])  
    mts_YChrom = mts_YChrom.delete_sites([i.id for i in mts_YChrom.sites() if i.position >= 900000])

    progress = tqdm.tqdm(total = mts_mtDNA.num_sites)
    with tsinfer.SampleData.from_tree_sequence(
        mts_mtDNA, path = "mtDNAtest.samples", use_times = False
    ) as mtDNA_sample_data:
        for var in mts_mtDNA.variants():
            mtDNA_sample_data.add_site(var.site.position, var.genotypes, var.alleles)
            progress.update()
        progress.close()
    
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
