import os
import dendropy
import argparse

def map_taxonnamespace(s, map):
    # map a taxon namespace
    a = '' + s
    for i in range(len(map)):
        a = a.replace(map[i], str(i+1))
    return a

def extract_quintet_tree(indices):
    true_species_tree_path = "avian_dataset/avian-model-species.tre"
    gene_tree_path = 'data/avian-0_5X-1000-500-all.f200.stripped.tre'
    str_indices = ''.join([str(i) for i in indices])

    if os.path.exists('data/species_tree_mapped' + str_indices + '.tre'):
        return

    tns = dendropy.TaxonNamespace()
    species_tree = dendropy.Tree.get(path=true_species_tree_path, schema='newick',
                                     taxon_namespace=tns, rooting="default-rooted")
    subtree_taxa = [tns[i].label for i in indices]
    species_subtree = species_tree.extract_tree_with_taxa_labels(labels=subtree_taxa, suppress_unifurcations=True)
    s = map_taxonnamespace(species_subtree.as_string(schema='newick'), subtree_taxa)
    species_subtree_mapped = dendropy.Tree.get(data=s, schema='newick')
    species_subtree_mapped.write_to_path(dest='data/species_tree_mapped' + str_indices + '.tre', schema='newick', suppress_edge_lengths=True,
                                    suppress_internal_node_labels=True)

    # striping genes
    gene_trees = dendropy.TreeList.get(path=gene_tree_path, schema='newick', taxon_namespace=tns)
    induced_trees = dendropy.TreeList()

    for g in gene_trees:
        subtree = g.extract_tree_with_taxa_labels(labels=subtree_taxa, suppress_unifurcations=True)
        s = map_taxonnamespace(subtree.as_string(schema='newick'), subtree_taxa)
        induced_trees.read(data=s, schema='newick')
    induced_trees.write_to_path(dest='data/avian_genes_mapped' + str_indices + '.tre', schema='newick')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indices", type=int, nargs="+", help="indices of taxa",
                        required=True, default=None)
    return parser.parse_args()


if __name__ == "__main__":
    extract_quintet_tree(parse_args().indices)
