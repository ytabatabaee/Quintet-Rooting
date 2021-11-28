import os
import glob
import subprocess
import dendropy

model_condition = 'avian-2X-1000-500-all' # only change this
dataset_path = 'data/avian_dataset_final/extracted_quintets/'
model_list = glob.glob(dataset_path + model_condition + '/gene_trees_mapped*.tre')
output_dir = 'data/avian_dataset_final/astral_trees/'

correct_topology_count = 0
avg_rf_dist = 0

for item in model_list:
    indices_string = ''.join(c for c in item.split('/')[-1] if c.isdigit())
    true_species_tree_path = dataset_path + 'species_tree_mapped' + indices_string + '.tre'
    gene_tree_path = dataset_path + model_condition + '/gene_trees_mapped' + indices_string + '.tre'
    output_path = output_dir + model_condition + '/astral_species_tree' + indices_string + '.tre'

    if not os.path.exists(output_dir + model_condition):
        os.makedirs(output_dir + model_condition)

    cmd = 'java -jar Astral/astral.5.7.8.jar -i' + gene_tree_path + ' -o ' + output_path + ' 2>out.log'
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    tns = dendropy.TaxonNamespace()
    species_tree_toplogy = dendropy.Tree.get(path=true_species_tree_path, schema='newick',
                                        rooting="force-unrooted", taxon_namespace=tns)
    astral_tree_toplogy = dendropy.Tree.get(path=output_path, schema='newick',
                                        rooting="force-unrooted", taxon_namespace=tns, suppress_edge_lengths=True)
    d = dendropy.calculate.treecompare.symmetric_difference(species_tree_toplogy, astral_tree_toplogy)
    avg_rf_dist += d
    correct_topology_count += 1 * (d == 0)


data_size = len(model_list)

print("Percentage of tests where the infered tree had the correct topology:")
print(correct_topology_count/data_size*100)
print("Average RF distance (not normalized, i.e. fp+fn)")
print(avg_rf_dist/data_size)
