import os
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import *

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))

def tree_construction(filename):

    aln = AlignIO.read(os.path.join(ROOT_DIR, 'data', 'fasta', filename), 'fasta')

    #constructing distance-matrix/nj tree
    constructor = DistanceTreeConstructor(DistanceCalculator("identity"), 'nj')
    tree_nj = constructor.build_tree(aln)
    tree_nj.root_at_midpoint()
    Phylo.write(tree_nj, os.path.join(ROOT_DIR, 'data', 'tree_nj', filename).replace('.fasta','.tree'),'newick')


    
# directory = os.fsencode(os.path.join(ROOT_DIR, 'data', 'fasta'))
 
# for file in os.listdir(directory):
#      filename = os.fsdecode(file)
#      if filename.endswith(".fasta"): 
#          tree_construction(filename)
#      else:
#          continue
