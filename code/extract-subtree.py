import dendropy
import os
import sys
import argparse

def extract_subtrees(input_path, output_path, subtree_taxa, input_filename, output_filename):
    model_conditions = os.listdir(input_path)
    for condition in model_conditions:
        #print(condition)
        condition_output_path = os.path.join(output_path, condition)
        if not os.path.exists(condition_output_path):
            os.mkdir(condition_output_path)
        replicates = os.listdir(os.path.join(input_path, condition))
        for replicate in replicates:
            replicate_output_path = os.path.join(output_path, condition, replicate)
            if not os.path.exists(replicate_output_path):
                os.mkdir(replicate_output_path)
            input_file_path = os.path.join(input_path, condition, replicate, input_filename)
            output_file_path = os.path.join(output_path, condition, replicate, output_filename)
            output_file = open(output_file_path, "a")
            input_trees = dendropy.TreeList.get(path=input_file_path, schema='newick')
            print(len(input_trees))
            for tree in input_trees:
                subtree = tree.extract_tree_with_taxa_labels(labels=subtree_taxa, suppress_unifurcations=True)
                output_file.write(subtree.as_string('newick'))
            output_file.close()


def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=str, help="Input folder path (subfolders per model condition)", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output foler path (subfolders per model condition will be created)", required=True)
    parser.add_argument("-if", "--input_filename", type=str, help="Name of file containing input trees in each replicate", required=True)
    parser.add_argument("-of", "--output_filename", type=str, help="Name of file containing output subtrees in each replicate", required=False, default="subtrees.tre")
    parser.add_argument("-t", "--taxa", type=str, nargs='+', help="Taxa in extracted trees", required=True)
    """
    parser.add_argument("-rc", "--replicate_count", type=int, help="No. of replicates per model condition", required=False)
    parser.add_argument("-cc", "--condition_count", type=int, help="No. of model conditions", required=False, default=None)
    """
    
    args = parser.parse_args()
    input_path = args.input
    output_path = args.output
    input_filename = args.input_filename
    output_filename = args.output_filename
    subtree_taxa = args.taxa
    
    extract_subtrees(input_path, output_path, subtree_taxa, input_filename, output_filename)    
if __name__ == "__main__":
    main()
