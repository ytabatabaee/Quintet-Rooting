import sys
import os
import argparse
import dendropy
import numpy as np
import time
import itertools


def main(args):
    species_tree_path = args.speciestree
    gene_tree_path = args.genetrees
    output_path = args.outputtree

    # reading gene tree and species tree topology files
    tns = dendropy.TaxonNamespace()
    species_tree_toplogy = dendropy.Tree.get(path=species_tree_path, schema='newick',
                                         taxon_namespace=tns, rooting="force-unrooted", suppress_edge_lengths=True)
    gene_trees = dendropy.TreeList.get(path=gene_tree_path, schema='newick', taxon_namespace=tns, rooting="force-unrooted")
    unrooted_quintets_base = dendropy.TreeList.get(path='topologies/quintets.tre', schema='newick')
    rooted_quintets_base = dendropy.TreeList()
    rooted_quintets_base.read(path='topologies/caterpillar.tre', schema='newick', rooting="default-rooted")
    rooted_quintets_base.read(path='topologies/pseudo_caterpillar.tre', schema='newick', rooting="default-rooted")
    rooted_quintets_base.read(path='topologies/balanced.tre', schema='newick', rooting="default-rooted")
    rooted_quintet_indices = np.load('rooted_quintet_indices.npy')

    rooted_candidates = get_all_rooted_trees(species_tree_toplogy)
    r_score = np.zeros(len(rooted_candidates))
    print(len(rooted_candidates))

    taxon_set = [t.label for t in tns]
    all_quintet_taxa = list(itertools.combinations(taxon_set, 5))
    print(len(all_quintet_taxa))
    for i in range(len(rooted_candidates)):
        r = rooted_candidates[i]
        print(r.as_string(schema='newick'))
        for q_taxa in all_quintet_taxa:
            print(q_taxa)
            start_time = time.time()
            cost = extract_quintets(r, gene_trees, q_taxa, unrooted_quintets_base, rooted_quintets_base, tns, rooted_quintet_indices)
            print(time.time() - start_time)
            print(cost)
            r_score[i] += cost
        print(r_score[i])
    print(len(all_quintet_taxa))
    min_idx = np.argmin(r_score)
    for i in range(len(rooted_candidates)):
        print(rooted_candidates[i].as_string(schema='newick'))
        print(r_score[i])
    print("The best tree is:")
    print(min_idx)
    print(rooted_candidates[min_idx])
    return



def extract_quintets(species_tree, gene_trees, q_taxa, unrooted_quintets_base, rooted_quintets_base, tns, rooted_quintet_indices):
    subtree = species_tree.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
    print(subtree.as_string(schema='newick'))
    quintets = [dendropy.Tree.get(data=map_taxonnamespace(str(q), q_taxa)+';', schema='newick', rooting='force-unrooted', taxon_namespace=tns) for q in unrooted_quintets_base]
    rooted_quintets = [dendropy.Tree.get(data=map_taxonnamespace(str(q), q_taxa)+';', schema='newick', rooting='force-rooted', taxon_namespace=tns) for q in rooted_quintets_base]
    # compute ui values
    u_count = np.zeros(len(quintets))
    for g in gene_trees:
        g_subtree = g.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
        for i in range(len(quintets)):
            if dendropy.calculate.treecompare.symmetric_difference(quintets[i], g_subtree) == 0:
                u_count[i] += 1
                break
    u_distribution = u_count / len(gene_trees)
    # get the index of the rooted tree
    idx = -1
    for i in range(len(rooted_quintets)):
        if dendropy.calculate.treecompare.symmetric_difference(rooted_quintets[i], subtree) == 0:
            print(rooted_quintets[i].as_string(schema='newick'))
            idx = i
            break
    print("rooted tree index, " + str(idx+1))
    tree_type = None
    if idx < 60:
        tree_type = 'c'
    elif idx >= 60 and idx < 75:
        tree_type = 'p'
    elif idx >= 75 and idx < 105:
        tree_type = 'b'
    indices = rooted_quintet_indices[idx]
    return cost(u_distribution, indices, tree_type)


def map_taxonnamespace(s, map):
    a = '' + s
    for i in range(len(map)):
        a = a.replace(str(i+1), map[i])
    return a


def get_all_rooted_trees(unrooted_tree):
    # generate all the 2n-3 possible rooted trees
    rooted_candidates = []
    tree = dendropy.Tree(unrooted_tree)
    for edge in tree.preorder_edge_iter():
        try:
            tree.reroot_at_edge(edge, update_bipartitions=True)
            rooted_candidates.append(dendropy.Tree(tree))
        except:
            continue
    # removing duplicates (by default, it returns 2n-2 trees)
    for i1 in range(len(rooted_candidates)):
        for i2 in range(i1+1, len(rooted_candidates)):
            if dendropy.calculate.treecompare.symmetric_difference(rooted_candidates[i1], rooted_candidates[i2]) == 0:
                rooted_candidates.pop(i2)
                break
    return rooted_candidates


def invariant_metric(a, b):
    return np.abs(a - b)# / (a + b)


def inequality_metric(a, b): # it should be a < b
    return (a-b) * (a > b)#a - b# penalty when a > b


def cost(u, indices, type):
    invariant_score = 0
    inequality_score = 0
    equivalence_classes = []
    inequality_classes = []
    #print(indices)
    #for i in indices:
    #    print(u[i])
    if type == 'c':
        equivalence_classes = [[0], [1], [2], [3, 12], [4, 11], [5, 8], [6, 7, 9, 10, 13, 14]]
        inequality_classes = [[0, 1], [0, 3], [1, 4], [3, 4], [4, 6], [2, 1], [2, 5], [5, 4]] # [a, b] -> a > b, a and b are cluster indices
    elif type == 'b':
        equivalence_classes = [[0], [1, 2], [3, 12], [4, 5, 8, 11], [6, 7, 9, 10, 13, 14]]
        inequality_classes = [[0, 1], [0, 2], [1, 3], [2, 3], [3, 4]]
    elif type == 'p':
        equivalence_classes = [[0], [1, 2], [3, 12], [7, 10], [4, 5, 6, 8, 9, 11, 13, 14]]
        inequality_classes = [[0, 1], [0, 2], [0, 3], [1, 4], [2, 4], [3, 4]]
    for c in equivalence_classes:
        #print([indices[i]+1 for i in c]) #** important
        inclass_distance = 0
        for i in range(len(c)):
            for j in range(len(c)):
                inclass_distance += invariant_metric(u[indices[c[i]]], u[indices[c[j]]])
                #print(str(u[indices[c[i]]]) + ' - ' + str(u[indices[c[j]]]))
        invariant_score += inclass_distance/(len(c))#(len(c)-1)
    #invariant_score = invariant_score# / len(equivalence_classes)

    for ineq in inequality_classes: #distance between clusters
        outclass_distance = 0
        #print("check " + str(ineq[0]) + " > " + str(ineq[1]))
        for i in equivalence_classes[ineq[0]]:
            for j in equivalence_classes[ineq[1]]:
                #print("check " + str(u[indices[i]]) + " > " + str(u[indices[j]]) + " with score " +  str(inequality_metric(u[indices[j]], u[indices[i]])))
                outclass_distance += inequality_metric(u[indices[j]], u[indices[i]])
        inequality_score += outclass_distance/(len(equivalence_classes[ineq[0]]))#*len(equivalence_classes[ineq[1]]))
    #inequality_score = 0# / len(inequality_classes)

    return invariant_score + inequality_score


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--speciestree", type=str, help="Input unrooted species tree",
                        required=True, default=None)

    parser.add_argument("-g", "--genetrees", type=str, help="Input gene trees",
                        required=True, default=None)

    parser.add_argument("-o", "--outputtree", type=str, help="Output rooted species tree",
                        required=True, default=None)

    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
