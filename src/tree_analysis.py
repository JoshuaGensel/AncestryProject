from Bio import Phylo
from matplotlib.pyplot import draw
import pandas as pd
import os

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))

testtree = Phylo.read(os.path.join(ROOT_DIR, 'data', 'tree_genealogy', 'Td_100_SB_50_T_5_run_1214366183_M.tree'), "newick")
Phylo.draw(testtree)