import dendropy
import os

taken_taxa = ["1", "2", "3", "4", "5"]

base_data_path = "F:/CS581/Project/Data"
base_out_path = base_data_path + "/5-Taxa"

def extract_D1_true_species_trees():
    D1_true_path = base_data_path + "/D1-true/simphy-uyen"
    D1_out_path = base_out_path + "/D1-true"
    os.mkdir(D1_out_path)
    replicates = 100
    sub_dirs = os.listdir(D1_true_path)
    for sub_dir in sub_dirs:
        if sub_dir.startswith("outgroup"):
            model_out_path = D1_out_path + "/" + sub_dir
            os.mkdir(model_out_path)
            for i in range (1, replicates + 1):
                replicate = 'f{r:03}'
                species_tree = D1_true_path + "/" + sub_dir + "/" + replicate + "/s_tree.trees"
                species_quintet = tree.extract_tree_with_taxa_labels(labels=taken_taxa, suppress_unifurcations=True)
                out_file = open(model_out_path + "/" + replicate + ".species.tre", "a")
                out_file.write(species_quintet.as_string('newick'))
                out_file.close


def extract_D2_true_species_trees():
    taxa_count = ["sp-2000", "sp-5000"]
    D2_true_path = base_data_path + "/D2-true/simphy-uyen"
    D2_out_path = base_out_path + "/D2-true"
    os.mkdir(D2_out_path)
    replicates = 20
    for taxa_count_ in taxa_count:
        tc_out_path = D2_out_path + "/" + taxa_count_
        os.mkdir(tc_out_path)
        sub_dirs = os.listdir(D1_true_path)
        for sub_dir in sub_dirs:
            if sub_dir.startswith("outgroup"):
                model_out_path = tc_out_path + "/" + sub_dir
                os.mkdir(model_out_path)
                for i in range (1, replicates + 1):
                    replicate = 'f{r:02}'
                    species_tree = D1_true_path + "/" + sub_dir + "/" + replicate + "/s_tree.trees"
                    species_quintet = tree.extract_tree_with_taxa_labels(labels=taken_taxa, suppress_unifurcations=True)
                    out_file = open(model_out_path + "/" + replicate + ".species.tre", "a")
                    out_file.write(species_quintet.as_string('newick'))
                    out_file.close

                
            
def extract_D1_estimated_gene_trees():
    D1_est_path = base_data_path + "/D1-estimated"
    D1_out_path = base_out_path + "/D1-Estimated"
    os.mkdir(D1_out_path)
    total_taxa_count = 100
    for model_condition in os.listdir(D1_est_path):
        print(model_condition)
        model_out_path = D1_out_path + "/" + model_condition
        os.mkdir(model_out_path)
        for r in range (1, total_taxa_count+1):
            replicate = f'{r:03}'
            gene_tree_path = D1_est_path + "/" + model_condition + "/" + replicate + "/estimatedgenetrees/estimatedgenetre.gtr"
            gene_trees = dendropy.TreeList.get(path=gene_tree_path, schema='newick')
            #print(len(gene_trees))
            out_file = open(model_out_path + "/" + replicate + ".quintet.tre", "a")
            for tree in gene_trees:
                #print(tree.as_string('newick') + "\n")
                quintet = tree.extract_tree_with_taxa_labels(labels=taken_taxa, suppress_unifurcations=True)
                #print(quintet.as_string('newick') + "\n\n")
                out_file.write(quintet.as_string('newick'))          
            out_file.close()

def extract_D2_estimated_gene_trees():
    D2_est_path = base_data_path + "/D2-estimated/sp-5000"
    D2_out_path = base_out_path + "/D2-Estimated/sp-5000"
    os.mkdir(D2_out_path)

    total_taxa_count = 20

    for model_condition in os.listdir(D2_est_path):
        print(model_condition)
        model_out_path = D2_out_path + "/" + model_condition
        os.mkdir(model_out_path)
        for r in range (1, total_taxa_count+1):
            replicate = f'{r:02}'
            gene_tree_path = D2_est_path + "/" + model_condition + "/" + replicate + "/estimatedgenetrees/estimatedgenetre.gtr"
            gene_trees = dendropy.TreeList.get(path=gene_tree_path, schema='newick')
            #print(len(gene_trees))
            out_file = open(model_out_path + "/" + replicate + ".quintet.tre", "a")
            for tree in gene_trees:
                #print(tree.as_string('newick') + "\n")
                quintet = tree.extract_tree_with_taxa_labels(labels=taken_taxa, suppress_unifurcations=True)
                #print(quintet.as_string('newick') + "\n\n")
                out_file.write(quintet.as_string('newick'))          
            out_file.close()



extract_D1_true_species_trees()
extract_D2_true_species_trees()
