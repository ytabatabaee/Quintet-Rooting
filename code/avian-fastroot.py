import os
import dendropy

base_path = "F:/CS581/Project/Data/Astral-III/avian/Simulated/Test_Cases"
#astral_path = "F:/Astral.5.7.8/Astral/astral.5.7.8.jar"

#test_cases = os.listdir(base_path)
#total_cases = len(test_cases)
test_case_start = 1
test_case_end = 2

conditions = ['0_5X-1000-500', '1X-1000-250', '1X-1000-500', '1X-1000-1000', '1X-1000-1500', '2X-1000-500']

"""
for case in range(test_case_start, test_case_end + 1):
	case_path = base_path + "/" + str(case)
	print(str(case))
	for condition in conditions:
		input = case_path + '/avian-' + condition + '-all.f200'
		output = case_path + '/' + condition + '.astral'
		command = 'java -jar ' + astral_path + ' -i ' + input + ' -o ' + output
		os.system(command)
		print(condition)
print(" ******* ")
"""

fastroot_path = "F:/MinVar-Rooting-master/MinVar-Rooting-master/FastRoot.py"



methods = ['MP', 'MV']

for condition in conditions:
	correct = [0.0, 0.0, 0.0]
	dist = [0.0, 0.0, 0.0]
	log_file = open(base_path + '/' + condition + '.log', 'w')
	for case in range(test_case_start, test_case_end + 1):
		print('case ', str(case))
		case_path = base_path + "/" + str(case)
		input = case_path + '/' + condition + '.astral'
		for i in range(len(methods)):
			print('method ', methods[i])
			output = case_path + '/' + condition + '.' + methods[i]
			os.system('python ' + fastroot_path + ' -i ' + input + ' -m ' + methods[i] + ' -o ' + output)
			taxa = dendropy.TaxonNamespace()
			in_tree = dendropy.Tree.get(path=input, rooting='force-rooted', schema='newick', taxon_namespace=taxa)
			out_tree = dendropy.Tree.get(path=output, rooting='force-rooted', schema='newick', taxon_namespace=taxa)
			in_tree.encode_bipartitions()
			out_tree.encode_bipartitions()
			method_dist = dendropy.calculate.treecompare.symmetric_difference(in_tree, out_tree)
			dist[i] += method_dist
			if(method_dist == 0):
				correct[i] += 1.0
				log_file.write('case ' + str(case) + '\t' + methods[i] + '\trf\t' + str(method_dist) + '\tT')
			else:
				log_file.write('case ' + str(case) + '\t' + methods[i] + '\trf\t' + str(method_dist) + '\tF')
		log_file.write('\n')
		log_file.flush()

	log_file.write('\n')
	total_cases = test_case_end - test_case_start + 1
	for i in range(len(methods)):
		correct_pct = (correct[i] * 100.0) / total_cases
		dist_avg = dist[i] / total_cases
		log_file.write(methods[i] + '\tcorrect_percentage\t' + str(correct_pct) + '\taverage_dist_distance\t' + str(dist_avg) + '\n')
	log_file.flush()
	log_file.close()