import csv
from msilib.schema import Error
import os
import ete3
import re


ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))

def tree_analysis(filename):
    
    #getting run parameters from filename
    parameter_vector = re.findall(r'\d+', filename)
    ID = int(parameter_vector[0])
    TD = int(parameter_vector[1])
    TA = int(parameter_vector[2])
    NS = int(parameter_vector[3])
    SOURCE = filename[-6]

    #loading genealogy
    genealogy = ete3.Tree(os.path.join(ROOT_DIR, 'data', 'tree_genealogy', filename))

    #function for calculating noninformative nodes
    def calc_NO_INFO(inputTree):
        non_informative_nodes = 0
        for n in inputTree.traverse():
            leafpops = set()
            for l in n.get_leaf_names():
                leafpops.add(l[:2])
            if("P1" in leafpops and "P2" in leafpops):
                non_informative_nodes += 1
        return(non_informative_nodes-1)/(len([n for n in inputTree.traverse()])-len([n for n in inputTree.iter_leaves()]))

    #Setting NO_INFO metrics
    G_NO_INFO = calc_NO_INFO(genealogy)

    #function for infering the sourcepopulation of a nodebased on the tree clustering
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

    #function for calculating P1-proportion, porportion of unknown nodes (which cluster with both sourcepops) 
    # and proportion of nodes which sourcepopulation get inferred wrongly by looking at the clustering in the tree 
    
    def calc_inference_metrics(inputTree):
    
        P3_nodes = []
        for i in inputTree.iter_leaves():
            if i.name.startswith("P3"):
                P3_nodes.append(i)

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

        return([(p1_correct + p1_false)/(NS-unknown_SP1-unknown_SP2),(unknown_SP1 + unknown_SP2)/NS,(p1_false + p2_false)/NS])

    G_metrics = calc_inference_metrics(genealogy)
    G_P1 = G_metrics[0]
    G_UNKNOWN = G_metrics[1]
    G_FALSE = G_metrics[2]

    def findDriftfile():
        directory = os.fsencode(os.path.join(ROOT_DIR, 'data', 'drift'))
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename == f"ID_{ID}_TD_{TD}.txt": 
                return os.path.join(ROOT_DIR, 'data', 'drift',filename)
            else:
                continue
        return "no driftfile found!"
    
    driftfile = open(findDriftfile(),'r')
    lines = driftfile.readlines()
    for l in lines:
        values = l.split(sep="	")
        if values[2] == str(TD+100000):
            initialFrequencies = values
        elif values[2] == str(TA+TD+100000):
            trueFrequencies = values
        else:
            continue

    if(SOURCE == "M"):
        INIT_P1 = float(initialFrequencies[5])
        TRUE_P1 = float(trueFrequencies[5])
        with open(os.path.join(ROOT_DIR, 'data', 'tree_analysis_data', 'treedata_M.csv'),"a") as outputfile:
            values_writer = csv.writer(outputfile, delimiter=',')
            values_writer.writerow([ID,TD,TA,NS,INIT_P1,TRUE_P1,G_NO_INFO,G_P1,G_FALSE,G_UNKNOWN])
            outputfile.close()
    elif(SOURCE == "Y"):
        INIT_P1 = float(initialFrequencies[3])
        TRUE_P1 = float(trueFrequencies[3])
        with open(os.path.join(ROOT_DIR, 'data', 'tree_analysis_data', 'treedata_Y.csv'),"a") as outputfile:
            values_writer = csv.writer(outputfile, delimiter=',')
            values_writer.writerow([ID,TD,TA,NS,INIT_P1,TRUE_P1,G_NO_INFO,G_P1,G_FALSE,G_UNKNOWN])
            outputfile.close()
    else:
        raise Error("No genetic source in filename!")

directory = os.fsencode(os.path.join(ROOT_DIR, 'data', 'tree_genealogy'))
 
for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith(".tree"): 
         tree_analysis(filename)
     else:
         continue