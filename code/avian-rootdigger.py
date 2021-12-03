import dendropy
import os
import argparse

def extract_model_subtree(test_case_no, test_case_folder, model_tree_path, output_path):
	print("subtree extraction\n")
	model_tree = dendropy.Tree.get(path=model_tree_path, rooting='force-rooted', schema='newick')
	
	subtree_taxa = []
	infile = open(test_case_folder + '/' + str(test_case_no) + '/' + str(test_case_no) + '.txt', "r")
	for line in infile:
		subtree_taxa.append(line.rstrip())
	subtree = model_tree.extract_tree_with_taxa_labels(labels=subtree_taxa, suppress_unifurcations=True)
	outfile = open(output_path, "w")
	outfile.write(subtree.as_string('newick'))
	outfile.flush()
	outfile.close()
	command = "sed -e \"s/\[&R\] //g\" -i " + output_path
	os.system(command)
	print("extracted\n")

def run_rootdigger(test_case_folder, test_case_count, model_tree_path, rootdigger_path):
	
	rootdigger_log = open(test_case_folder + "/avian-rootdigger.log", "w")
	
	correct_count = 0.0
	total_sym_diff = 0.0
	
	for test_case_no in range(1, test_case_count + 1):
		print("Test case ", test_case_no)

		model_subtree_path = test_case_folder + '/' + str(test_case_no) + '/' + str(test_case_no) + '.model.tre'
		alignment_path = test_case_folder + '/' + str(test_case_no) + '/' + str(test_case_no) + '.fasta'
		rooted_path = test_case_folder + '/' + str(test_case_no) + '/' + str(test_case_no) + '.root.tre'
		log_path = test_case_folder + '/' + str(test_case_no) + '/' + str(test_case_no) + '.log'

		extract_model_subtree(test_case_no, test_case_folder, model_tree_path, model_subtree_path)

		#command = rootdigger_path + ' --msa ' + alignment_path + ' --tree ' + model_subtree_path + ' > ' + log_path
		command = rootdigger_path + ' --msa ' + alignment_path + ' --tree ' + model_subtree_path
		os.system(command)
		
		command = "mv " + model_subtree_path + ".rooted.tree " + rooted_path

		#command = "grep -v -E \"^\[\" " + log_path + " | grep -v -E \"^Inf\" > " + rooted_path
		os.system(command)

		tns = dendropy.TaxonNamespace()
		ref_tree = dendropy.Tree.get(path=model_subtree_path, rooting='force-rooted', schema='newick', taxon_namespace=tns)
		est_tree = dendropy.Tree.get(path=rooted_path, rooting='force-rooted', schema='newick', taxon_namespace=tns)
		ref_tree.encode_bipartitions()
		est_tree.encode_bipartitions()
		sym_diff = dendropy.calculate.treecompare.symmetric_difference(ref_tree, est_tree)

		rootdigger_log.write("******* Test case - " + str(test_case_no) + " *******\n")
		rootdigger_log.write("Reference: " + ref_tree.as_string('newick') + "\n")
		rootdigger_log.write("Estimated: " + est_tree.as_string('newick') + "\n")
		rootdigger_log.write("Symmetric difference: " + str(sym_diff) + "\n\n")
		rootdigger_log.flush()
		if(sym_diff == 0):
			correct_count += 1.0
		else:
			total_sym_diff += sym_diff
		print("Done\n")

	rootdigger_log.close()		

	rootdigger_stat = open(test_case_folder + "/avian-rootdigger-stat.txt", "w")

	correct_pct = (correct_count / test_case_count) * 100.0
	avg_sym_diff = (total_sym_diff / test_case_count)

	rootdigger_stat.write("Correct rooting percentage = " + str(correct_pct) + "\n")
	rootdigger_stat.write("Average difference = " + str(avg_sym_diff) + "\n")
	rootdigger_stat.flush()
	rootdigger_stat.close()

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-f", "--test_case_folder", type=str, help="Folder with all test cases (subfolders per test case)", required=True)
    parser.add_argument("-c", "--test_case_count", type=int, help="No. of test cases", required=True)
    parser.add_argument("-m", "--model_tree_path", type=str, help="Model tree path", required=True)
    parser.add_argument("-r", "--rootdigger_path", type=str, help="Rootdigger binary path", required=True)
    
    args = parser.parse_args()
    
    test_case_folder = args.test_case_folder
    test_case_count = args.test_case_count
    model_tree_path = args.model_tree_path
    rootdigger_path = args.rootdigger_path
    
    run_rootdigger(test_case_folder, test_case_count, model_tree_path, rootdigger_path)

if __name__ == "__main__":
    main()
