import os
import glob
import subprocess
import dendropy

species_tree_list = glob.glob('data/species_tree_mapped_with_lengths*.tre')
print(len(species_tree_list))

count = 0
avg_distance = 0

for species_tree_with_lengths in species_tree_list:
    indices_string = ''.join(c for c in species_tree_with_lengths if c.isdigit())
    true_species_tree_path = 'data/species_tree_mapped' + indices_string + '.tre'
    output_path = 'data/mp_species_tree' + indices_string + '.tre'
    cmd = 'python ../MinVar-Rooting/FastRoot.py -m MP -i ' + species_tree_with_lengths + ' -o ' + output_path
    #print(true_species_tree_path)

    #print(output_path)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    tns = dendropy.TaxonNamespace()
    true_species_tree = dendropy.Tree.get(path=true_species_tree_path, schema='newick',
                                        rooting="force-rooted", taxon_namespace=tns)
    mp_species_tree = dendropy.Tree.get(path=output_path, schema='newick',
                                        rooting="force-rooted", taxon_namespace=tns, suppress_edge_lengths=True)
    d = dendropy.calculate.treecompare.symmetric_difference(true_species_tree, mp_species_tree)
    avg_distance += d
    count += int(d == 0)

print(count/len(species_tree_list))
print(avg_distance/len(species_tree_list))
