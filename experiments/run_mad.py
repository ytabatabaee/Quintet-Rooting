import os
import glob
import subprocess
import dendropy

def rf_distance(t1, t2):
    t1.encode_bipartitions()
    t2.encode_bipartitions()
    return dendropy.calculate.treecompare.symmetric_difference(t1, t2)

model_condition = 'avian-0_5X-1000-500-all' # only change this
dataset_path = '../data/avian_dataset/extracted_quintets/'
model_list = glob.glob(dataset_path + model_condition + '/gene_trees_mapped*.tre')
output_dir = '../data/avian_dataset/mad_trees/'

count = 0
avg_clade_distance = 0

for item in model_list:
    indices_string = ''.join(c for c in item.split('/')[-1] if c.isdigit())
    species_tree_with_lengths = dataset_path + 'species_tree_mapped_with_lengths' + indices_string + '.tre'
    true_species_tree_path = dataset_path + 'species_tree_mapped' + indices_string + '.tre'

    if not os.path.exists(output_dir + model_condition):
        os.makedirs(output_dir + model_condition)

    output_path = output_dir + model_condition + '/estimated_species_tree' + indices_string + '.tre'
    st_with_length = dendropy.Tree.get(path=species_tree_with_lengths, schema='newick')
    with open('temp.tre', 'w') as f:
        f.write(str(st_with_length)[:str(st_with_length).rindex(')')+1]+';')
    cmd = 'python ../../mad/mad.py  temp.tre -n'

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    print(out.decode("utf-8"))
    tns = dendropy.TaxonNamespace()
    true_species_tree = dendropy.Tree.get(path=true_species_tree_path, schema='newick',
                                        rooting="force-rooted", taxon_namespace=tns)
    mp_species_tree = dendropy.Tree.get(path='temp.tre.rooted', schema='newick',
                                        rooting="force-rooted", taxon_namespace=tns, suppress_edge_lengths=True)
    d = rf_distance(true_species_tree, mp_species_tree)
    avg_clade_distance += d
    count += int(d == 0)

data_size = len(model_list)


print("Test count")
print(data_size)
print("Percentage of tests correctly rooted:")
print(count/data_size*100)
print("Average RF distance (not normalized, i.e. fp+fn)")
print(avg_clade_distance/data_size)
