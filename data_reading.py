import os
import dendropy

def map_taxonnamespace(s, map):
    # map a taxon namespace to
    a = '' + s
    for i in range(len(map)):
        a = a.replace(subtree_taxa[i], str(i+1))
    return a

true_species_tree_path = "avian_dataset/avian-model-species.tre"
gene_tree_path = "avian_dataset/avian-0_5X-1000-500-all.f200"

tns = dendropy.TaxonNamespace()
species_tree = dendropy.Tree.get(path=true_species_tree_path, schema='newick',
                                 taxon_namespace=tns, rooting="default-rooted")
subtree_taxa = [t.label for t in tns[16:21]]
#print(tns)
#print(subtree_taxa)

species_subtree = species_tree.extract_tree_with_taxa_labels(labels=subtree_taxa, suppress_unifurcations=True)

#print(species_tree.as_string(schema='newick'))
#print(species_subtree.as_string(schema='newick'))
species_subtree.write_to_path(dest='avian_species.tre', schema='newick', suppress_edge_lengths=True,
                                suppress_internal_node_labels=True)
s = map_taxonnamespace(species_subtree.as_string(schema='newick'), subtree_taxa)
species_subtree_mapped = dendropy.Tree.get(data=s, schema='newick')
species_subtree_mapped.write_to_path(dest='species_tree_mapped.tre', schema='newick', suppress_edge_lengths=True,
                                suppress_internal_node_labels=True)

# striping genes

gene_trees = dendropy.TreeList.get(path=gene_tree_path, schema='newick', taxon_namespace=tns)
induced_trees = dendropy.TreeList()

for g in gene_trees:
    subtree = g.extract_tree_with_taxa_labels(labels=subtree_taxa, suppress_unifurcations=True)
    induced_trees.read(data=subtree.as_string(schema='newick'), schema='newick')
induced_trees.write_to_path(dest='avian_genes.tre', schema='newick', suppress_edge_lengths=True,
                                suppress_internal_node_labels=True)


quintets = dendropy.TreeList.get(path='topologies/quintets.tre', schema='newick')
genes = dendropy.TreeList.get(path='avian_genes.tre', schema='newick', rooting='default-unrooted')

induced_trees_mapped = dendropy.TreeList()

for g in genes:
    s = map_taxonnamespace(g.as_string(schema='newick'), subtree_taxa)
    induced_trees_mapped.read(data=s, schema='newick')

induced_trees_mapped.write_to_path(dest='avian_genes_mapped.tre', schema='newick', suppress_edge_lengths=True,
                                suppress_internal_node_labels=True)

s = genes[0].as_string(schema='newick')
print(s)
s = map_taxonnamespace(s, subtree_taxa)
print(s)
