import os
import glob
import subprocess
import dendropy

gene_tree_list = glob.glob('data/avian_genes_mapped*.tre')

count = 0

for gene_tree_path in gene_tree_list:
    indices_string = ''.join(c for c in gene_tree_path if c.isdigit())
    species_tree_path = 'data/species_tree_mapped' + indices_string + '.tre'
    output_path = 'data/species_tree_astral' + indices_string + '.tre'
    print(species_tree_path)
    #print(output_path)
    cmd = 'java -jar Astral/astral.5.7.8.jar -i' + gene_tree_path + ' -o ' + output_path + ' 2>out.log'
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    tns = dendropy.TaxonNamespace()
    species_tree_toplogy = dendropy.Tree.get(path=species_tree_path, schema='newick',
                                        rooting="force-unrooted", taxon_namespace=tns)
    astral_tree_toplogy = dendropy.Tree.get(path=output_path, schema='newick',
                                        rooting="force-unrooted", taxon_namespace=tns, suppress_edge_lengths=True)
    if dendropy.calculate.treecompare.symmetric_difference(species_tree_toplogy, astral_tree_toplogy) == 0:
        count += 1

print(count/len(gene_tree_list))
