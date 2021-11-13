import numpy as np
import itertools
import data_reading
import rooting
import os
import time
import subprocess

taxa_count = 48
indices_list = list(itertools.combinations(np.arange(0, taxa_count), 5))[:100]
topk_count = 0
correct_topology_count = 0
correct_tree_count = 0
avg_rf_dist = 0

for indices in indices_list:
    start_time = time.time()
    string_indices = ''
    for i in indices:
        string_indices += str(i) + ' '

    cmd = 'python3 data_reading.py -i ' + string_indices
    print(cmd)

    #cmd = 'pwd'
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    print(out)
    print(err)

    cmd = 'python3 rooting.py -i random.txt -o random.txt'
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

    print(err)

    #os.system("python3 data_reading.py -i 1 2 3 4 5")
    print(time.time() - start_time)

print("percentage of the times where the true species tree is among the top k candidates:")
print(topk_count/len(indices_list))
print("percentage of the times where the infered tree had the correct topology:")
print(correct_topology_count/len(indices_list))
print("percentage of the times where the inferred tree was the true rooted species tree :")
print(correct_tree_count/len(indices_list))
print("avg rf distance (rooted)")
print(avg_rf_dist/len(indices_list))
