import os
import glob
import subprocess
import dendropy

model_condition = 'avian-1X-1000-500-all'
dataset_path = '../data/avian_dataset/extracted_quintets/'
model_list = glob.glob(dataset_path + model_condition + '/gene_trees_mapped*.tre')
method = 'MV' # only this and the model condition need to be changed afterwards
output_dir = '../data/avian_dataset/' + method + '_trees/'

count = 0
avg_clade_distance = 0

for item in model_list[:10]:
    indices_string = ''.join(c for c in item.split('/')[-1] if c.isdigit())
    species_tree_with_lengths = dataset_path + 'species_tree_mapped_with_lengths' + indices_string + '.tre'
    true_species_tree_path = dataset_path + 'species_tree_mapped' + indices_string + '.tre'

    if not os.path.exists(output_dir + model_condition):
        os.makedirs(output_dir + model_condition)

    output_path = output_dir + model_condition + '/estimated_species_tree' + indices_string + '.tre'
    cmd = 'python ../../MinVar-Rooting/FastRoot.py -m ' + method +' -i ' + species_tree_with_lengths + ' -o ' + output_path

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    tns = dendropy.TaxonNamespace()
    true_species_tree = dendropy.Tree.get(path=true_species_tree_path, schema='newick',
                                        rooting="force-rooted", taxon_namespace=tns)
    mp_species_tree = dendropy.Tree.get(path=output_path, schema='newick',
                                        rooting="force-rooted", taxon_namespace=tns, suppress_edge_lengths=True)
    d = dendropy.calculate.treecompare.symmetric_difference(true_species_tree, mp_species_tree)
    avg_clade_distance += d
    count += int(d == 0)

data_size = len(model_list[:10])

print("Test count")
print(data_size)
print("Percentage of tests correctly rooted:")
print(count/data_size*100)
print("Average RF distance (not normalized, i.e. fp+fn)")
print(avg_clade_distance/data_size)
