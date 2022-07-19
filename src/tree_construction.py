import os
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import *

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))

def tree_construction(filename):

    aln = AlignIO.read(os.path.join(ROOT_DIR, 'data', 'fasta', filename), 'fasta')

    #constructing distance-matrix/UPGMA tree
    constructor_upgma = DistanceTreeConstructor(DistanceCalculator("identity"), 'upgma')
    tree_upgma = constructor_upgma.build_tree(aln)
    Phylo.write(tree_upgma, os.path.join(ROOT_DIR, 'data', 'tree_upgma', filename).replace('.fasta','.tree'),'newick')

    #constructing maximum-parsimony tree based on neighbor joining tree
    constructor_nj = DistanceTreeConstructor(DistanceCalculator("identity"), 'nj')
    tree_nj = constructor_nj.build_tree(aln)
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor = ParsimonyTreeConstructor(searcher, tree_nj)
    tree_mp = constructor.build_tree(aln)
    Phylo.write(tree_mp, os.path.join(ROOT_DIR, 'data', 'tree_mp', filename).replace('.fasta','.tree'),'newick')
    
# directory = os.fsencode(os.path.join(ROOT_DIR, 'data', 'fasta'))
 
# for file in os.listdir(directory):
#      filename = os.fsdecode(file)
#      if filename.endswith(".fasta"): 
#          tree_construction(filename)
#      else:
#          continue
