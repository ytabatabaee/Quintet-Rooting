"""
Extract taxa subsets from a species tree and a corresponding set of gene trees
All rights reserved.
License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""

import dendropy
import os
import argparse
import random
import itertools


def main(args):
    n = args.testnum
    m = args.taxanum
    species_tree_path = args.speciestree
    gene_trees_path = args.genespath
    dataset_path = args.datapath
    taxa_list = args.taxalist

    tns = dendropy.TaxonNamespace()
    species_tree = dendropy.Tree.get(path=species_tree_path, schema='newick', taxon_namespace=tns)
    species_taxa = [t.label for t in tns]
    gene_trees = dendropy.TreeList.get(path=gene_trees_path, schema='newick', taxon_namespace=tns)

    if taxa_list: # take input list of taxa
        samples = list([t for t in taxa_list if t in species_taxa])
        n = 1
    else:
        samples = random.sample(list(itertools.combinations(species_taxa, m)), n)

    if not os.path.exists(dataset_path):
        os.makedirs(dataset_path)

    for i in range(1, n + 1):
        smpl_taxa = samples[i-1]
        test_path = dataset_path + '/' + str(i)
        if not os.path.exists(test_path):
            os.makedirs(test_path)
        with open(test_path + '/taxa_labels.txt', 'w') as fp:
            for t in smpl_taxa:
                fp.write(t + "\n")
        subtree = species_tree.extract_tree_with_taxa_labels(labels=smpl_taxa, suppress_unifurcations=True)
        subtree.write_to_path(dest=test_path + '/model-species.tre', schema='newick')
        all_gene_trees = dendropy.TreeList()
        for g in gene_trees:
            g_taxa = set([x.taxon.label for x in g.leaf_nodes()])
            smpl_g_taxa = [t for t in smpl_taxa if t in g_taxa]
            try:
                gene_subtree = g.extract_tree_with_taxa_labels(labels=smpl_g_taxa, suppress_unifurcations=True)
                if len(gene_subtree.leaf_nodes()) < 2:
                    continue
                all_gene_trees.append(gene_subtree)
            except:
                continue
        all_gene_trees.write_to_path(dest=test_path + '/gene-trees.tre', schema='newick')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--speciestree", type=str, help="Path to file containing species tree in newick",
                        required=True, default=None)
    parser.add_argument("-g", "--genespath", type=str, help="Path to file containing gene trees in newick",
                        required=True, default=None)
    parser.add_argument("-m", "--taxanum", type=int, help="number of taxa to select",
                        required=False, default=5)
    parser.add_argument("-n", "--testnum", type=int, help="number of sample tests to create",
                        required=False, default=1)
    parser.add_argument("-d", "--datapath", type=str, help="Output path (where dataset should be created)",
                        required=False, default='extracted_data')
    parser.add_argument("-tl", "--taxalist", type = str, help="list of input taxa",
                        required=False, nargs='*', default=[])
    main(parser.parse_args())
