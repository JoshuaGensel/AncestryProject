from logging import FATAL
from msilib.schema import Error
from unicodedata import name
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
import textwrap

#default values for important variables

path = "D:/Daten/programming_projects/AncestryProject/data/ts_raw/Td_100_SB_50_T_5_run_1214366183.trees"
nSamples = 10
ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))

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

    #saving the name of the .trees input file
    head, tail = os.path.split(path)
    inputFileName = tail
    
    #Importing and simplifying tree sequence to nodes that resemble mtDNA of females or Y-chromosome of males
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

    #Checking that there is now only 1 tree with 1 root for mtDNA and Y-chromosome
    if(ts_M_raw.num_trees != 1): raise ValueError("more than one tree!")
    if(ts_Y_raw.num_trees != 1): raise ValueError("more than one tree!")

    if(ts_M_raw.first().num_roots > 1): raise ValueError("more than one root!")
    if(ts_Y_raw.first().num_roots > 1): raise ValueError("more than one root!")

    #randomly sampling n = nSamples number of female individuals per subpopulation
    #than simplifying down to only consider these samples
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

    #now same sampling procedure for males
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

    #output tree as newick tree file
    #for mtDNA tree
    output = open(os.path.join(ROOT_DIR, 'data', 'tree_genealogy', inputFileName.replace('.trees', '_M.tree')), "w")
    output.write(sts_M.first().as_newick())
    output.close()
    #for Y-chromosome tree
    output = open(os.path.join(ROOT_DIR, 'data', 'tree_genealogy', inputFileName.replace('.trees', '_Y.tree')), "w")
    output.write(sts_Y.first().as_newick())
    output.close()
    
    #overlaying mutations on the tree sequences with msprime
    mts_M = msprime.mutate(sts_M, rate=6.85*10**-7, random_seed=1, keep=False)
    mts_Y = msprime.mutate(sts_Y, rate=3.01*10**-8, random_seed=1, keep=False)

    #deleting mutations in the wrong region (i.e. mutations in the mtDNA of a node that is supposed to be a Y-chromosome)
    mts_M = mts_M.delete_sites([i.id for i in mts_M.sites() if i.position < 900000])  
    mts_Y = mts_Y.delete_sites([i.id for i in mts_Y.sites() if i.position >= 900000])  
    
    #output haplotypes of the sampled individuals as FASTA
    #first for mtDNA
    output = open(os.path.join(ROOT_DIR, 'data', 'fasta', inputFileName.replace('.trees', '_M.fasta')), "w")
    for h,v in zip(mts_M.haplotypes(), mts_M.samples()):
        seq = h.replace('0','A').replace('1','C')
        seq = textwrap.fill(seq, 80)
        pop=sts_M.individual(v).metadata["subpopulation"]
        name = f">Pop{pop}_Ind{v}\n"
        if pop == 3:
            if sts_M.mutation_at(v,10) != -1:
                name = f">Pop{pop}_Ind{v}_SourceP1\n"
            elif sts_M.mutation_at(v,5) != -1:
                name = f">Pop{pop}_Ind{v}_SourceP2\n"
            else:
                raise Error("P3 ind without Sourcepopulation!")
        output.write(name+seq+"\n")
    output.close()

    #here for Y-chromosome
    output = open(os.path.join(ROOT_DIR, 'data', 'fasta', inputFileName.replace('.trees', '_Y.fasta')), "w")
    for h,v in zip(mts_Y.haplotypes(), mts_Y.samples()):
        seq = h.replace('0','A').replace('1','C')
        seq = textwrap.fill(seq, 80)
        pop = sts_Y.individual(v).metadata["subpopulation"]
        name = f">Pop{pop}_Ind{v}\n"
        if pop == 3:
            if sts_Y.mutation_at(v,916558) != -1:
                name = f">Pop{pop}_Ind{v}_SourceP1\n"
            elif sts_Y.mutation_at(v,916563) != -1:
                name = f">Pop{pop}_Ind{v}_SourceP2\n"
            else:
                raise Error("P3 ind without Sourcepopulation!")
        output.write(name+seq+"\n")
    output.close()
    
    


    
haplotypeCalc(path)