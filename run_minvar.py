import os
import glob
import subprocess
import dendropy

species_tree_list = glob.glob('data/species_tree_mapped_with_lengths*.tre')


for species_tree in species_tree_list:
    indices_string = ''.join(c for c in species_tree if c.isdigit())
    true_species_tree_path = 'data/species_tree_mapped' + indices_string + '.tre'
    output_path = 'mp_species_tree' + indices_string + '.tre'
    cmd = 'python ../MinVar-Rooting/FastRoot.py -m MP -i ' + species_tree + ' -o ' + output_path
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

print(count/len(gene_tree_list))'''
