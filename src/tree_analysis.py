from fileinput import filename
from tkinter import constants
from tkinter.tix import Tree
from unicodedata import name
from Bio import Phylo
from matplotlib.pyplot import draw, get
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

P3_nodes = []
for i in test_ete_tree.iter_leaves():
    if i.name.startswith("P3"):
        P3_nodes.append(i)

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

for i in P3_nodes:
    print(i.name, ": ", get_SP(i))

Phylo.draw(test_phylo_tree)