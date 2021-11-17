import os
import dendropy
import argparse
import glob

genes_per_model_condition = glob.glob('avian_dataset/*.f200')

for gene_tree_path in genes_per_model_condition:
    print(gene_tree_path)
    gene_trees = dendropy.TreeList.get(path=gene_tree_path, schema='newick', suppress_edge_lengths=True)
    gene_trees.write_to_path(dest='data/' + gene_tree_path + '.stripped.tre', schema='newick', suppress_edge_lengths=True,
                                        suppress_internal_node_labels=True)
