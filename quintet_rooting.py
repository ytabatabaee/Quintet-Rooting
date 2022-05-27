import sys
import os
import argparse
import dendropy
import numpy as np
import time
import itertools
import random


def main(args):
    species_tree_path = args.speciestree
    gene_tree_path = args.genetrees
    output_path = args.outputtree
    sampling_method = args.samplingmethod
    random.seed(args.seed)

    # reading gene tree and species tree topology files
    tns = dendropy.TaxonNamespace()
    tns_base = dendropy.TaxonNamespace()
    unrooted_species_tree_topology = dendropy.Tree.get(path=species_tree_path, schema='newick',
                                         taxon_namespace=tns, rooting="force-unrooted", suppress_edge_lengths=True)
    gene_trees = dendropy.TreeList.get(path=gene_tree_path, schema='newick', taxon_namespace=tns, rooting="force-unrooted")
    unrooted_quintets_base = dendropy.TreeList.get(path='topologies/quintets.tre', taxon_namespace=tns_base, schema='newick')
    rooted_quintets_base = dendropy.TreeList(taxon_namespace=tns_base)
    rooted_quintets_base.read(path='topologies/caterpillar.tre', schema='newick', rooting="default-rooted")
    rooted_quintets_base.read(path='topologies/pseudo_caterpillar.tre', schema='newick', rooting="default-rooted")
    rooted_quintets_base.read(path='topologies/balanced.tre', schema='newick', rooting="default-rooted")
    rooted_quintet_indices = np.load('rooted_quintet_indices.npy')

    # search space of rooted trees
    rooted_candidates = get_all_rooted_trees(unrooted_species_tree_topology)
    r_score = np.zeros(len(rooted_candidates))

    # set of sampled quintets
    taxon_set = [t.label for t in tns]
    sample_quintet_taxa = []
    if len(taxon_set) == 5 or sampling_method == 'd':
        sample_quintet_taxa = list(itertools.combinations(taxon_set, 5))
    elif sampling_method == 'tc':
        sample_quintet_taxa = triplet_cover_sample(taxon_set)
    elif sampling_method == 'le':
        sample_quintet_taxa = linear_quintet_encoding_sample(unrooted_species_tree_topology, taxon_set)
    elif sampling_method == 'rl':
        sample_quintet_taxa = random_linear_sampling(taxon_set)


    start = time.time()

    # preprocessing
    quintet_scores = np.zeros((len(sample_quintet_taxa), 7))
    quintet_unrooted_indices = np.zeros(len(sample_quintet_taxa), dtype=int)
    quintets_r_all = []

    for j in range(len(sample_quintet_taxa)):
        q_taxa = sample_quintet_taxa[j]
        quintets_u = [dendropy.Tree.get(data=map_taxonnamespace(str(q), q_taxa)+';', schema='newick', rooting='force-unrooted', taxon_namespace=tns) for q in unrooted_quintets_base]
        quintets_r = [dendropy.Tree.get(data=map_taxonnamespace(str(q), q_taxa)+';', schema='newick', rooting='force-rooted', taxon_namespace=tns) for q in rooted_quintets_base] # you don't need all of this
        subtree_u = unrooted_species_tree_topology.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
        quintet_tree_dist = gene_tree_distribution(gene_trees, q_taxa, quintets_u)
        quintet_unrooted_indices[j] = get_quintet_unrooted_index(subtree_u, q_taxa, quintets_u)
        quintet_scores[j] = compute_cost_rooted_quintets(quintet_tree_dist, quintet_unrooted_indices[j], rooted_quintet_indices)
        quintets_r_all.append(quintets_r)

    #print('preprocessing time', time.time() - start)

    start = time.time()

    # computing scores
    for i in range(len(rooted_candidates)):
        r = rooted_candidates[i]
        for j in range(len(sample_quintet_taxa)):
            q_taxa = sample_quintet_taxa[j]
            subtree_r = r.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
            r_idx = get_quintet_rooted_index(subtree_r, q_taxa, quintets_r_all[j], quintet_unrooted_indices[j])
            r_score[i] += quintet_scores[j][r_idx]

    min_idx = np.argmin(r_score)
    rooted_candidates[min_idx].write_to_path(dest=output_path, schema='newick')
    #print(r_score)

    #print('scoring time', time.time() - start)

    if args.confidencescore:
        confidence_scores = (np.max(r_score) - r_score)/np.sum(np.max(r_score) - r_score)
        tree_ranking_indices = np.argsort(r_score)
        with open(output_path+".rank.cfn", 'w') as fp:
            for i in tree_ranking_indices:
                fp.write(str(rooted_candidates[i]) + ';\n')
                fp.write(str(confidence_scores[i]) + '\n')


def random_linear_sampling(taxon_set):
    n = len(taxon_set)
    return random.sample(list(itertools.combinations(taxon_set, 5)), 2*n-3)


def triplet_cover_sample(taxon_set):
    sample_quintet_taxa = []
    all_triplet_taxa = list(itertools.combinations(taxon_set, 3))
    for trip in all_triplet_taxa:
        trip_prime = random.sample([x for x in taxon_set if x not in trip], 2)
        sample_quintet_taxa.append(tuple(list(trip) + list(trip_prime)))
    return sample_quintet_taxa


def linear_quintet_encoding_sample(unrooted_tree, taxon_set):
    sample_quintet_taxa = []
    tree = dendropy.Tree(unrooted_tree)
    for edge in tree.preorder_edge_iter():
        seed_node = tree.seed_node
        try:
            if edge.is_leaf():
                quintet = []
                tri_node = edge.tail_node
                tree.reroot_at_node(tri_node, update_bipartitions=True)
                tri_partition_taxa = []
                for child in tree.seed_node.child_node_iter():
                    partition = []
                    for c in child.leaf_nodes():
                        partition.append(c.taxon.label)
                    tri_partition_taxa.append(partition)
                if len(tri_partition_taxa) > 2:
                    for partition in tri_partition_taxa:
                        quintet.extend(random.sample(partition, 1))
                    quintet.extend(random.sample([x for x in taxon_set if x not in quintet], 5 - len(quintet)))
                    sample_quintet_taxa.append(tuple(quintet))
                tree.reroot_at_node(seed_node, update_bipartitions=True)

            elif edge.is_internal():
                quintet = []
                adj_edges  = edge.get_adjacent_edges()
                tree.reroot_at_edge(edge, update_bipartitions=True)
                four_partition_taxa = []
                for e in adj_edges:
                    partition = []
                    for c in e.head_node.leaf_nodes():
                        partition.append(c.taxon.label)
                    four_partition_taxa.append(partition)
                for i in range(len(four_partition_taxa)):
                    for j in range(i+1, len(four_partition_taxa)):
                        if set(four_partition_taxa[i]).issubset(set(four_partition_taxa[j])):
                            four_partition_taxa[j] = list(set(four_partition_taxa[j]) - set(four_partition_taxa[i]))
                four_partition_taxa = [p for p in four_partition_taxa if p != []]
                if len(four_partition_taxa) > 2:
                    for partition in four_partition_taxa:
                        quintet.extend(random.sample(partition, 1))
                    quintet.extend(random.sample([x for x in taxon_set if x not in quintet], 5 - len(quintet)))
                    sample_quintet_taxa.append(tuple(quintet))
                tree.reroot_at_node(seed_node, update_bipartitions=True)
        except Exception as ex:
            continue
    return sample_quintet_taxa

# ADR Table 5
u2r_mapping = np.array([[0, 1, 58, 59, 66, 75, 104],
                       [2, 3, 52, 53, 63, 78, 103],
                       [4, 5, 46, 47, 60, 87, 102],
                       [6, 7, 56, 57, 69, 76, 101],
                       [8, 9, 40, 41, 64, 81, 100],
                       [10, 11, 34, 35, 61, 90, 99],
                       [12, 13, 50, 51, 70, 79, 95],
                       [14, 15, 38, 39, 67, 82, 94],
                       [16, 17, 28, 29, 62, 93, 96],
                       [18, 19, 44, 45, 71, 86, 88],
                       [20, 21, 32, 33, 68, 85, 91],
                       [22, 23, 26, 27, 65, 84, 97],
                       [24, 25, 54, 55, 72, 77, 98],
                       [30, 31, 48, 49, 73, 80, 92],
                       [36, 37, 42, 43, 74, 83, 89]])


# score the 7 rooted variants of a quintet
def compute_cost_rooted_quintets(u_distribution, u_idx, rooted_quintet_indices):
    rooted_tree_indices = u2r_mapping[u_idx]
    costs = np.zeros(7)
    for i in range(7):
        idx = rooted_tree_indices[i]
        unlabeled_topology = idx_2_unlabeled_topology(idx)
        indices = rooted_quintet_indices[idx]
        costs[i] = cost(u_distribution, indices, unlabeled_topology)
    return costs


def idx_2_unlabeled_topology(idx):
    if idx < 60:
        return 'c'
    elif idx >= 60 and idx < 75:
        return 'p'
    elif idx >= 75 and idx < 105:
        return 'b'
    return None


def gene_tree_distribution(gene_trees, q_taxa, quintets_u):
    u_count = np.zeros(len(quintets_u))
    for g in gene_trees:
        g_subtree = g.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
        for i in range(len(quintets_u)):
            if dendropy.calculate.treecompare.symmetric_difference(quintets_u[i], g_subtree) == 0:
                u_count[i] += 1
                break
    u_distribution = u_count / len(gene_trees)
    return u_distribution


def get_quintet_unrooted_index(subtree_u, q_taxa, quintets_u):
    idx_u = -1
    for i in range(len(quintets_u)):
        if dendropy.calculate.treecompare.symmetric_difference(quintets_u[i], subtree_u) == 0:
            idx_u = i
            break
    return idx_u


def get_quintet_rooted_index(subtree_r, q_taxa, quintets_r, u_idx):
    idx_r = -1
    for i in range(7):
        idx = u2r_mapping[u_idx][i]
        if dendropy.calculate.treecompare.symmetric_difference(quintets_r[idx], subtree_r) == 0:
            idx_r = i
            break
    return idx_r


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
    # removing duplicates
    for i1 in range(len(rooted_candidates)):
        for i2 in range(i1+1, len(rooted_candidates)):
            if dendropy.calculate.treecompare.symmetric_difference(rooted_candidates[i1], rooted_candidates[i2]) == 0:
                rooted_candidates.pop(i2)
                break
    return rooted_candidates


def invariant_metric(a, b):
    return np.abs(a - b)


def inequality_metric(a, b):
    return (a - b) * (a > b)


def cost(u, indices, type):
    invariant_score = 0
    inequality_score = 0
    equivalence_classes = []
    inequality_classes = []
    if type == 'c':
        equivalence_classes = [[0], [1], [2], [3, 12], [4, 11], [5, 8], [6, 7, 9, 10, 13, 14]]
        inequality_classes = [[0, 1], [0, 3], [1, 4], [3, 4], [4, 6], [2, 1], [2, 5], [5, 4]] # [a, b] -> a > b, a and b are cluster indices
    elif type == 'b':
        equivalence_classes = [[0], [1, 2], [3, 12], [4, 5, 8, 11], [6, 7, 9, 10, 13, 14]]
        inequality_classes = [[0, 1], [0, 2], [1, 3], [2, 3], [3, 4]]
    elif type == 'p':
        equivalence_classes = [[0], [1, 2], [3, 12], [7, 10], [4, 5, 6, 8, 9, 11, 13, 14]]
        inequality_classes = [[0, 1], [0, 2], [0, 3], [1, 4], [2, 4], [3, 4]]
    for c in equivalence_classes: #distance inside clusters
        inclass_distance = 0
        for i in range(len(c)):
            for j in range(len(c)):
                inclass_distance += invariant_metric(u[indices[c[i]]], u[indices[c[j]]])
        invariant_score += inclass_distance/(len(c))

    for ineq in inequality_classes: #distance between clusters
        outclass_distance = 0
        for i in equivalence_classes[ineq[0]]:
            for j in equivalence_classes[ineq[1]]:
                outclass_distance += inequality_metric(u[indices[j]], u[indices[i]])
        inequality_score += outclass_distance/(len(equivalence_classes[ineq[0]]))
    return invariant_score + inequality_score


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--speciestree", type=str, help="input unrooted species tree in newick format",
                        required=True, default=None)

    parser.add_argument("-g", "--genetrees", type=str, help="input gene trees in newick format",
                        required=True, default=None)

    parser.add_argument("-o", "--outputtree", type=str, help="output file containing a rooted species tree",
                        required=True, default=None)

    parser.add_argument("-sm", "--samplingmethod", type=str, help="quintet sampling method (TC for triplet cover, LE for linear encoding)",
                        required=False, default='d')

    parser.add_argument("-cfs", "--confidencescore",  action='store_true', help="output confidence scores for each possible rooted tree")

    parser.add_argument("-c", "--cost", type=str, help="cost function (MC for minimal constraints)",
                        required=False, default='d')

    parser.add_argument("-rs", "--seed", type=int, help="random seed", required=False, default=1234)

    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
