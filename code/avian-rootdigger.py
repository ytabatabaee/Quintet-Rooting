import dendropy
import os
import argparse

def extract_model_subtree(test_case_no, test_case_folder, model_tree_path, output_path):
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
	print("subtree extracted\n")

def run_rootdigger(test_case_folder, test_case_count, model_tree_path, rootdigger_path, output_folder, log_filename, stat_filename):
	
	rootdigger_log = open(output_folder + "/" + log_filename, "w")
	
	correct_count = 0.0
	total_sym_diff = 0.0
	
	for test_case_no in range(1, test_case_count + 1):
		print("Test case ", test_case_no)
        
		case_output_folder = output_folder + '/' + str(test_case_no)
		os.mkdir(case_output_folder)

		model_subtree_path = test_case_folder + '/' + str(test_case_no) + '/' + str(test_case_no) + '.model.tre'
                
		extract_model_subtree(test_case_no, test_case_folder, model_tree_path, model_subtree_path)
        
		alignment_path = test_case_folder + '/' + str(test_case_no) + '/' + str(test_case_no) + '.fasta'
		log_path = case_output_folder + '/' + str(test_case_no) + '.log'

		command = rootdigger_path + ' --msa ' + alignment_path + ' --tree ' + model_subtree_path + ' --exhaustive > ' + log_path
		os.system(command)
        
		default_rooted_path = model_subtree_path + '.rooted.tree'
		default_lwr_path = model_subtree_path + '.lwr.tree'
		default_ckp_path = model_subtree_path + '.ckp'
        
		output_rooted_path = case_output_folder + '/' + str(test_case_no) + '.rooted.tree'
		output_lwr_path = case_output_folder + '/' + str(test_case_no) + '.lwr.tree'
		output_ckp_path = case_output_folder + '/' + str(test_case_no) + '.ckp'
        
		command = "mv " + default_rooted_path + " " + output_rooted_path
		os.system(command)
		command = "mv " + default_lwr_path + " " + output_lwr_path
		os.system(command)
		command = "mv " + default_ckp_path + " " + output_ckp_path
		os.system(command)

		tns = dendropy.TaxonNamespace()
		ref_tree = dendropy.Tree.get(path=model_subtree_path, rooting='force-rooted', schema='newick', taxon_namespace=tns)
		est_tree = dendropy.Tree.get(path=output_rooted_path, rooting='force-rooted', schema='newick', taxon_namespace=tns)
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

	rootdigger_stat = open(output_folder + "/" + stat_filename, "w")

	correct_pct = (correct_count / test_case_count) * 100.0
	avg_sym_diff = (total_sym_diff / test_case_count)

	rootdigger_stat.write("Correct rooting percentage = " + str(correct_pct) + "\n")
	rootdigger_stat.write("Average symmetric difference = " + str(avg_sym_diff) + "\n")
	rootdigger_stat.flush()
	rootdigger_stat.close()

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-f", "--test_case_folder", type=str, help="Folder with all test cases (subfolders per test case)", required=True)
    parser.add_argument("-c", "--test_case_count", type=int, help="No. of test cases", required=True)
    parser.add_argument("-m", "--model_tree_path", type=str, help="Path to model tree containing 48 species", required=True)
    parser.add_argument("-r", "--rootdigger_path", type=str, help="Rootdigger binary path", required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="Subfolders per test case will be created", required=True)
    parser.add_argument("-l", "--log_filename", type=str, help="Filename for per test case log", required=False, default="rootdigger.log")
    parser.add_argument("-s", "--stat_filename", type=str, help="Filename for aggregated statistics", required=False, default="rootdigger.stat")

    args = parser.parse_args()
    
    test_case_folder = args.test_case_folder
    test_case_count = args.test_case_count
    model_tree_path = args.model_tree_path
    rootdigger_path = args.rootdigger_path
    output_folder = args.output_folder
    log_filename = args.log_filename
    stat_filename = args.stat_filename
    
    run_rootdigger(test_case_folder, test_case_count, model_tree_path, rootdigger_path, output_folder, log_filename, stat_filename)

if __name__ == "__main__":
    main()
