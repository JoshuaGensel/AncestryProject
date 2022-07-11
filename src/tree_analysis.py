import csv
from fileinput import filename
from msilib.schema import Error
import os
import ete3
import re


ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))

def tree_analysis(filename):
    parameter_vector = re.findall(r'\d+', filename)
    TD = int(parameter_vector[0])
    SB = int(parameter_vector[1])
    T = int(parameter_vector[2])
    ID = int(parameter_vector[3])
    SOURCE = filename[-6]

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

    p1_proportion = (p1_correct + p1_false)/(N_SAMPLES-unknown_SP1-unknown_SP2)
    p2_proportion = (p2_correct + p2_false)/(N_SAMPLES-unknown_SP1-unknown_SP2)
    unknown_proportion = (unknown_SP1 + unknown_SP2)/N_SAMPLES
    false_proportion = (p1_false + p2_false)/N_SAMPLES

    with open(os.path.join(ROOT_DIR, 'data', 'tree_analysis_data', 'treedata.csv'),"a") as outputfile:
        values_writer = csv.writer(outputfile, delimiter=',')
        values_writer.writerow([ID,SOURCE,TD,T,SB,N_SAMPLES,non_informative_nodes,p1_proportion,p2_proportion,unknown_proportion,false_proportion])
        outputfile.close()

directory = os.fsencode(os.path.join(ROOT_DIR, 'data', 'tree_genealogy'))
    
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith(".tree"): 
         tree_analysis(filename)
     else:
         continue