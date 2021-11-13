import numpy as np
import itertools
import data_reading
import rooting
import os
import time
import subprocess

taxa_count = 48
indices_list = list(itertools.combinations(np.arange(0, taxa_count), 5))[:10]
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
    print(err)


    #os.system("python3 data_reading.py -i 1 2 3 4 5")
    print(time.time() - start_time)
