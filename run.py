import numpy as np
import itertools
import data_reading
import rooting
import os
import time
import subprocess
import random
import glob

def generate_data(indices_list):
    for indices in indices_list:
        string_indices = ' '.join([str(i) for i in indices])
        cmd = 'python3 data_reading.py -i ' + string_indices
        print(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        out, err = p.communicate()

def generate_indices(taxa_count, n):
    indices_list = set()
    while len(indices_list) < n:
        indices_list.add(tuple(sorted(random.sample(range(taxa_count), 5))))
    return list(indices_list)


taxa_count = 48
n = 100

topk_count = 0
correct_topology_count = 0
correct_tree_count = 0
avg_rf_dist = 0

'''indices_list = generate_indices(taxa_count, n)
with open('test_indices_with_lengths.txt', 'a') as fp:
    fp.write(str(indices_list))
generate_data(indices_list)'''

gene_tree_list = glob.glob('data/avian_genes_mapped*.tre')
print(len(gene_tree_list))

for gene_tree_path in gene_tree_list[:10]:
    start_time = time.time()
    indices_string = ''.join(c for c in gene_tree_path if c.isdigit())
    species_tree_path = 'data/species_tree_mapped' + indices_string + '.tre'
    cmd = 'python3 rooting.py -i ' + species_tree_path + ' -o ' + gene_tree_path
    print(cmd)

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    print(out.decode("utf-8"))
    lis = list(out.decode("utf-8").split("\n"))
    length = len(lis)

    topk_count += int(lis[length-5])
    correct_tree_count += int(lis[length-4])
    correct_topology_count += int(lis[length-3])
    avg_rf_dist += int(lis[length-2])

    print(time.time() - start_time)

data_size = len(gene_tree_list[:10])

print("Percentage of tests where the true species tree is among the top 3 (of 105) rooted candidates:")
print(topk_count/data_size*100)
print("Percentage of tests where the infered tree had the correct topology:")
print(correct_topology_count/data_size*100)
print("Percentage of tests where the inferred tree was the true rooted species tree :")
print(correct_tree_count/data_size*100)
print("Average RF distance (rooted, not normalized, i.e. fp+fn)")
print(avg_rf_dist/data_size)
