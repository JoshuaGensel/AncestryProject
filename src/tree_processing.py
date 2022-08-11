import pyslim
import tskit
import msprime
import numpy as np
import os
import textwrap

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))

def tree_processing(inputFileName,nSamples):

    #saving the name of the .trees input file
    path = os.path.join(ROOT_DIR, 'data', 'ts_raw',inputFileName)
    
    #Importing and simplifying tree sequence to nodes that resemble mtDNA of females or Y-chromosome of males
    ts_raw = pyslim.load(path)
    
    #picking nodes that resemble mtDNA of female individuals or Y-Chromosome of a male indivduals in the last generation
    nodes_M_p1 = []
    nodes_M_p2 = []
    nodes_M_p3 = []
    nodes_Y_p1 = []
    nodes_Y_p2 = []
    nodes_Y_p3 = []

    for u in ts_raw.samples(time=0):
        if ts_raw.mutation_at(u,0) != -1 and ts_raw.individual(ts_raw.node(u).individual).metadata["sex"] == pyslim.INDIVIDUAL_TYPE_FEMALE:
            if ts_raw.individual(ts_raw.node(u).individual).metadata["subpopulation"] == 1:
                nodes_M_p1.append(u)
            elif ts_raw.individual(ts_raw.node(u).individual).metadata["subpopulation"] == 2:
                nodes_M_p2.append(u)
            elif ts_raw.individual(ts_raw.node(u).individual).metadata["subpopulation"] == 3:
                nodes_M_p3.append(u)
            else:
                raise ValueError("Female Individual without population!")
        elif ts_raw.mutation_at(u,916567) != -1 and ts_raw.individual(ts_raw.node(u).individual).metadata["sex"] == pyslim.INDIVIDUAL_TYPE_MALE:
            if ts_raw.individual(ts_raw.node(u).individual).metadata["subpopulation"] == 1:
                nodes_Y_p1.append(u)
            elif ts_raw.individual(ts_raw.node(u).individual).metadata["subpopulation"] == 2:
                nodes_Y_p2.append(u)
            elif ts_raw.individual(ts_raw.node(u).individual).metadata["subpopulation"] == 3:
                nodes_Y_p3.append(u)
            else:
                raise ValueError("Male Individual without population!")
        else:
            pass

    #fucntion for creating newick and fasta output from a simplified tree sequence for mtDNA and MSY, n is just for naming the output files

    def createTreeOutput(sts_M,sts_Y,n):
        #output tree as newick tree file
        #for mtDNA tree
        node_labels_M = {}
        for s in sts_M.samples():
            pop = sts_M.individual(s).metadata["subpopulation"]
            name = f"P{pop}_n{s}_SP{pop}"
            if pop == 3:
                if sts_M.mutation_at(s,10) != -1:
                    name = f"P{pop}_n{s}_SP1"
                elif sts_M.mutation_at(s,5) != -1:
                    name = f"P{pop}_n{s}_SP2"
                else:
                    raise ValueError("P3 ind without Sourcepopulation!")
            node_labels_M[s] = name
        output = open(os.path.join(ROOT_DIR, 'data', 'tree_genealogy', inputFileName.replace('.trees', f'_NS_{n}_M.tree')), "w")
        output.write(sts_M.first().as_newick(node_labels=node_labels_M))
        output.close()
        #for Y-chromosome tree
        node_labels_Y = {}
        for s in sts_Y.samples():
            pop = sts_Y.individual(s).metadata["subpopulation"]
            name = f"P{pop}_n{s}_SP{pop}"
            if pop == 3:
                if sts_Y.mutation_at(s,916558) != -1:
                    name = f"P{pop}_n{s}_SP1"
                elif sts_Y.mutation_at(s,916563) != -1:
                    name = f"P{pop}_n{s}_SP2"
                else:
                    raise ValueError("P3 ind without Sourcepopulation!")
            node_labels_Y[s] = name
        output = open(os.path.join(ROOT_DIR, 'data', 'tree_genealogy', inputFileName.replace('.trees', f'_NS_{n}_Y.tree')), "w")
        output.write(sts_Y.first().as_newick(node_labels=node_labels_Y))
        output.close()

        #overlaying mutations on the tree sequences with msprime
        mts_M = msprime.mutate(sts_M, rate=6.85*10**-7, random_seed=1, keep=False)
        mts_Y = msprime.mutate(sts_Y, rate=3.01*10**-8, random_seed=1, keep=False)

        #deleting mutations in the wrong region (i.e. mutations in the mtDNA of a node that is supposed to be a Y-chromosome)
        mts_M = mts_M.delete_sites([i.id for i in mts_M.sites() if i.position < 900000])  
        mts_Y = mts_Y.delete_sites([i.id for i in mts_Y.sites() if i.position >= 900000])  
        
        #output haplotypes of the sampled individuals as FASTA
        #first for mtDNA
        output = open(os.path.join(ROOT_DIR, 'data', 'fasta', inputFileName.replace('.trees', f'_NS_{n}_M.fasta')), "w")
        for h,v in zip(mts_M.haplotypes(), mts_M.samples()):
            seq = h.replace('0','A').replace('1','C')
            seq = textwrap.fill(seq, 80)
            pop=sts_M.individual(v).metadata["subpopulation"]
            name = f">P{pop}_n{v}_SP{pop}\n"
            if pop == 3:
                if sts_M.mutation_at(v,10) != -1:
                    name = f">P{pop}_n{v}_SP1\n"
                elif sts_M.mutation_at(v,5) != -1:
                    name = f">P{pop}_n{v}_SP2\n"
                else:
                    raise ValueError("P3 ind without Sourcepopulation!")
            output.write(name+seq+"\n")
        output.close()

        #here for Y-chromosome
        output = open(os.path.join(ROOT_DIR, 'data', 'fasta', inputFileName.replace('.trees', f'_NS_{n}_Y.fasta')), "w")
        for h,v in zip(mts_Y.haplotypes(), mts_Y.samples()):
            seq = h.replace('0','A').replace('1','C')
            seq = textwrap.fill(seq, 80)
            pop = sts_Y.individual(v).metadata["subpopulation"]
            name = f">P{pop}_n{v}_SP{pop}\n"
            if pop == 3:
                if sts_Y.mutation_at(v,916558) != -1:
                    name = f">P{pop}_n{v}_SP1\n"
                elif sts_Y.mutation_at(v,916563) != -1:
                    name = f">P{pop}_n{v}_SP2\n"
                else:
                    raise ValueError("P3 ind without Sourcepopulation!")
            output.write(name+seq+"\n")
        output.close()
        
    #simplifying down for maximum samples provided
    
    maxSample_value = max(nSamples)
    
    keep_nodes_M_p1 = np.random.choice(nodes_M_p1, maxSample_value, replace=False)
    keep_nodes_M_p2 = np.random.choice(nodes_M_p2, maxSample_value, replace=False)
    keep_nodes_M_p3 = np.random.choice(nodes_M_p3, maxSample_value, replace=False)
    keep_nodes_M = np.concatenate((keep_nodes_M_p1, keep_nodes_M_p2, keep_nodes_M_p3))
    sts_M1 = ts_raw.simplify(keep_nodes_M)

    keep_nodes_Y_p1 = np.random.choice(nodes_Y_p1, maxSample_value, replace=False)
    keep_nodes_Y_p2 = np.random.choice(nodes_Y_p2, maxSample_value, replace=False)
    keep_nodes_Y_p3 = np.random.choice(nodes_Y_p3, maxSample_value, replace=False)
    keep_nodes_Y = np.concatenate((keep_nodes_Y_p1, keep_nodes_Y_p2, keep_nodes_Y_p3))
    sts_Y1 = ts_raw.simplify(keep_nodes_Y)

    #Checking that there is now only 1 tree with 1 root for mtDNA and Y-chromosome for max sample trees
    if(sts_M1.num_trees != 1): raise ValueError("more than one tree!")
    if(sts_Y1.num_trees != 1): raise ValueError("more than one tree!")

    if(sts_M1.first().num_roots > 1): raise ValueError("more than one root!")
    if(sts_Y1.first().num_roots > 1): raise ValueError("more than one root!")
    
    #creating output for max sample tree
    createTreeOutput(sts_M=sts_M1, sts_Y=sts_Y1, n=maxSample_value)

    #simplifying from max sample ts for every samplesize instead of raw to save time, then create output
    for i in nSamples:
        if i == max(nSamples):
            pass
        keep_nodes_M_i = np.concatenate((np.random.choice(sts_M1.samples(population=1), i, replace=False), 
                                    np.random.choice(sts_M1.samples(population=2), i, replace=False), 
                                    np.random.choice(sts_M1.samples(population=3), i, replace=False)))
        sts_M_i = sts_M1.simplify(keep_nodes_M_i)

        keep_nodes_Y_i = np.concatenate((np.random.choice(sts_Y1.samples(population=1), i, replace=False), 
                                        np.random.choice(sts_Y1.samples(population=2), i, replace=False), 
                                        np.random.choice(sts_Y1.samples(population=3), i, replace=False)))
        sts_Y_i = sts_Y1.simplify(keep_nodes_Y_i)

        #Checking that there is now only 1 tree with 1 root for mtDNA and Y-chromosome
        if(sts_M_i.num_trees != 1): raise ValueError("more than one tree!")
        if(sts_Y_i.num_trees != 1): raise ValueError("more than one tree!")

        if(sts_M_i.first().num_roots > 1): raise ValueError("more than one root!")
        if(sts_Y_i.first().num_roots > 1): raise ValueError("more than one root!")

        createTreeOutput(sts_M=sts_M_i, sts_Y=sts_Y_i, n=i)
    
# directory = os.fsencode(os.path.join(ROOT_DIR, 'data', 'ts_raw'))

# N_sample_values = [10,20,50]

# for file in os.listdir(directory):
#      filename = os.fsdecode(file)
#      if filename.endswith(".trees"):
#         tree_processing(filename,N_sample_values)
#      else:
#         continue