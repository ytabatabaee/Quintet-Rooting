import numpy as np
import itertools
import data_reading
import rooting
import os
import time
import subprocess
import random
import glob
import re

def generate_data(indices_list, true_species_tree_path, dataset_path):
    for indices in indices_list:
        string_indices = ' '.join([str(i) for i in indices])
        cmd = 'python3 data_reading.py -i ' + string_indices + ' -t ' + true_species_tree_path + ' -d ' + dataset_path
        print(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        print(out.decode("utf-8"))


def generate_indices(taxa_count, n):
    indices_list = set()
    while len(indices_list) < n:
        indices_list.add(tuple(sorted(random.sample(range(taxa_count), 5))))
    return list(indices_list)


def read_indicies_from_file():
    with open('test_indices_with_lengths.txt', 'r') as fp:
        indices_list_str = fp.read()
    indices_list_str = indices_list_str.split(')')
    indices_list = []
    for item in indices_list_str:
        indices_list.append([int(s) for s in item.replace('(', ' ').replace(',', ' ').split() if s.isdigit()])
    return indices_list[:-1]


taxa_count = 48
n = 10

model_condition = 'avian-0_5X-1000-500-all' # only this need to be changed afterward, perhaps take as input?
true_species_tree_path = 'avian_dataset/avian-model-species.tre'
#gene_tree_path = 'data/avian_dataset/' + model_condition + '-all.f200.stripped.tre'
dataset_path = 'data/avian_dataset/extracted_quintets/'

start_time = time.time()

#indices_list = generate_indices(taxa_count, n)
#indices_list = read_indicies_from_file()
#generate_data(indices_list, true_species_tree_path, dataset_path)

#print(time.time() - start_time)

model_list = glob.glob(dataset_path + model_condition + '/gene_trees_mapped*.tre')
#species_tree_list = glob.glob(dataset_path + 'species_tree_mapped_with_lengths*.tre')
#print(len(model_list))

topk_count = 0
correct_topology_count = 0
correct_tree_count = 0
avg_rf_dist = 0
avg_rf_dist_unrooted = 0

types_to_num = ['c', 'b', 'p']
class_confusion_matrix = np.zeros((3, 3))


for item in model_list[:10]:
    start_time = time.time()
    indices_string = ''.join(c for c in item.split('/')[-1] if c.isdigit())
    #species_tree = dataset_path + 'species_tree_mapped_with_lengths' + indices_string + '.tre'
    species_tree_path = dataset_path + 'species_tree_mapped' + indices_string + '.tre'
    gene_tree_path = dataset_path + model_condition + '/gene_trees_mapped' + indices_string + '.tre'
    cmd = 'python3 rooting.py -i ' + species_tree_path + ' -o ' + gene_tree_path
    print(cmd)

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    print(out.decode("utf-8"))
    lis = list(out.decode("utf-8").split("\n"))
    length = len(lis)

    true_type = lis[length-8]
    avg_rf_dist_unrooted += int(lis[length-7])
    type = lis[length-6]
    topk_count += int(lis[length-5])
    correct_tree_count += int(lis[length-4])
    correct_topology_count += int(lis[length-3])
    avg_rf_dist += int(lis[length-2])

    class_confusion_matrix[types_to_num.index(type)][types_to_num.index(true_type)] += 1

    print(time.time() - start_time)

data_size = len(model_list[:10])

print("class confusion matrix of types (rows are predicted, columns are true) (c, b, p)")
print(class_confusion_matrix)
print("Percentage of tests where the true species tree is among the top 5 (of 105) rooted candidates:")
print(topk_count/data_size*100)
print("Percentage of tests where the infered tree had the correct topology:")
print(correct_topology_count/data_size*100)
print("Percentage of tests where the inferred tree was the true rooted species tree :")
print(correct_tree_count/data_size*100)
print("Average Clade distance (rooted, not normalized, i.e. fp+fn)")
print(avg_rf_dist/data_size)
print("Average RF distance (not normalized, i.e. fp+fn)")
print(avg_rf_dist_unrooted/data_size)
