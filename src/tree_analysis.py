from msilib.schema import Error
from Bio import Phylo
import pandas as pd
import os
import ete3
import re

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))

filename = 'Td_50_SB_50_T_5_run_3482603321_M.tree'

parameter_vector = re.findall(r'\d+', filename)
TD = int(parameter_vector[0])
SB = int(parameter_vector[1])
T = int(parameter_vector[2])
ID = int(parameter_vector[3])

test_phylo_tree = Phylo.read(os.path.join(ROOT_DIR, 'data', 'tree_genealogy', filename), "newick")
test_ete_tree = ete3.Tree(os.path.join(ROOT_DIR, 'data', 'tree_genealogy', filename))

non_informative_nodes = 0
for n in test_ete_tree.traverse():
    leafpops = set()
    for l in n.get_leaf_names():
        leafpops.add(l[:2])
    if("P1" in leafpops and "P2" in leafpops):
        non_informative_nodes += 1

non_informative_nodes = (non_informative_nodes-1)/(len([n for n in test_ete_tree.traverse()])-len([n for n in test_ete_tree.get_leaves()]))

def get_SP(node):
    sisternodes = set()
    for n in node.up.get_leaves():
        sisternodes.add(n.name[:2])
    sisternodes.remove("P3")
    if("P1" in sisternodes and "P2" in sisternodes):
        return "unknown"
    elif("P1" in sisternodes):
        return "SP1"
    elif("P2" in sisternodes):
        return "SP2"
    else:
        return get_SP(node.up)

P3_nodes = []
for i in test_ete_tree.iter_leaves():
    if i.name.startswith("P3"):
        P3_nodes.append(i)
        
N_SAMPLES = len(P3_nodes)

p1_correct = 0
p1_false = 0
p2_correct = 0
p2_false = 0
unknown_SP1 = 0
unknown_SP2 = 0


for i in P3_nodes:
    #print(i.name, ": ", get_SP(i))
    if get_SP(i) == "SP1":
        if i.name[-3:] == get_SP(i):
            p1_correct += 1
        else:
            p1_false += 1
    elif get_SP(i) == "SP2":
        if i.name[-3:] == get_SP(i):
            p2_correct += 1
        else:
            p2_false += 1
    elif get_SP(i) == "unknown":
        if i.name[-3:] == "SP1":
            unknown_SP1 += 1
        else:
            unknown_SP2 += 1
    else:
        raise Error("Nothing inferred!")
    

p1_correct = p1_correct/N_SAMPLES
p1_false = p1_false/N_SAMPLES
unknown_SP1 = unknown_SP1/N_SAMPLES
p2_correct = p2_correct/N_SAMPLES
p2_false = p2_false/N_SAMPLES
unknown_SP2 = unknown_SP2/N_SAMPLES

p1_proportion = p1_correct + p1_false
p2_proportion = p2_correct + p2_false
unknown_proportion = unknown_SP1 + unknown_SP2
false_proportion = p1_false + p2_false

#print(false_proportion, unknown_proportion)
#Phylo.draw(test_phylo_tree)