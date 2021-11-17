import os
import dendropy
import argparse
import glob




gene_tree_path = "avian_dataset/avian-0_5X-1000-500-all.f200"

gene_trees = dendropy.TreeList.get(path=gene_tree_path, schema='newick', suppress_edge_lengths=True)
gene_trees.write_to_path(dest='data/avian-0_5X-1000-500-all.f200.stripped.tre', schema='newick', suppress_edge_lengths=True,
                                    suppress_internal_node_labels=True)
