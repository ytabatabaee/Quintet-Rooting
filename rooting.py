import sys
import os
import argparse
import dendropy
import numpy as np

MAX_NUM_INVARIANT = 20
CANDIDATE_SIZE = 7
MAX_VAL = 100


def main(args):
    input_path = args.input
    output_path = args.output
    mode = args.mode

    tns = dendropy.TaxonNamespace()
    quintets = dendropy.TreeList.get(path='topologies/quintets.tre', schema='newick', taxon_namespace=tns)

    true_species_tree = dendropy.Tree.get(path=input_path, schema='newick', taxon_namespace=tns)
    species_tree_toplogy = dendropy.Tree.get(data=true_species_tree.as_string(schema='newick'), schema='newick',
                                            rooting="force-unrooted", taxon_namespace=tns)

    gene_trees = dendropy.TreeList.get(path=output_path, schema='newick', taxon_namespace=tns)
    u_count = np.zeros(len(quintets))

    for g in gene_trees:
        for i in range(len(quintets)):
            d = dendropy.calculate.treecompare.symmetric_difference(quintets[i], g)
            if d == 0:
                u_count[i] += 1
                break

    u_distribution = u_count / len(gene_trees)
    print("estimated gene tree distribution")
    print(u_distribution)
    #print(np.sum(u_distribution))

    score_v = np.zeros(105)

    caterpillars = dendropy.TreeList.get(path='topologies/caterpillar.tre', schema='newick', taxon_namespace=tns)
    for i in range(len(caterpillars)):
        unrooted_cater = dendropy.Tree.get(data=caterpillars[i].as_string(schema='newick'), schema='newick',
                                        rooting="force-unrooted", taxon_namespace=tns)
        #if dendropy.calculate.treecompare.symmetric_difference(unrooted_cater, species_tree_toplogy) != 0:
        #    score_v[i] = MAX_VAL
        #    continue
        score_v[i] = score(caterpillars[i], caterpillars[0], u_distribution, tns, quintets, "c")
        #print(i, caterpillars[i], score_v[i],
        #      dendropy.calculate.treecompare.symmetric_difference(caterpillars[i], true_species_tree), "c")

    pseudo_caterpillars = dendropy.TreeList.get(path='topologies/pseudo_caterpillar.tre', schema='newick', taxon_namespace=tns)
    for i in range(len(pseudo_caterpillars)):
        unrooted_pseudo = dendropy.Tree.get(data=pseudo_caterpillars[i].as_string(schema='newick'), schema='newick',
                                        rooting="force-unrooted", taxon_namespace=tns)
        #if dendropy.calculate.treecompare.symmetric_difference(unrooted_pseudo, species_tree_toplogy) != 0:
        #    score_v[len(caterpillars) + i] = MAX_VAL
        #    continue
        score_v[len(caterpillars) + i] = score(pseudo_caterpillars[i], pseudo_caterpillars[6], u_distribution, tns, quintets, "p")
        #print(i+len(caterpillars), pseudo_caterpillars[i], score_v[len(caterpillars) + i],
        #      dendropy.calculate.treecompare.symmetric_difference(pseudo_caterpillars[i], true_species_tree), "p")

    balanced = dendropy.TreeList.get(path='topologies/balanced.tre', schema='newick', taxon_namespace=tns)
    for i in range(len(balanced)):
        unrooted_balanced = dendropy.Tree.get(data=balanced[i].as_string(schema='newick'), schema='newick',
                                        rooting="force-unrooted", taxon_namespace=tns)
        #if dendropy.calculate.treecompare.symmetric_difference(unrooted_balanced, species_tree_toplogy) != 0:
        #    score_v[i + len(caterpillars) + len(pseudo_caterpillars)] = MAX_VAL
        #    continue
        score_v[i + len(caterpillars) + len(pseudo_caterpillars)] = score(balanced[i], balanced[0], u_distribution, tns, quintets, "b")
        #print(i + len(caterpillars) + len(pseudo_caterpillars), balanced[i], score_v[i + len(caterpillars) + len(pseudo_caterpillars)],
        #      dendropy.calculate.treecompare.symmetric_difference(balanced[i], true_species_tree), "b")

    min_indices = sorted(range(len(score_v)), key = lambda sub: score_v[sub])[:CANDIDATE_SIZE]

    rooted_candidates = []
    rooted_candidate_types = []
    r_base = []

    for idx in min_indices:
        if idx < 60:
            rooted_candidates.append(caterpillars[idx])
            rooted_candidate_types.append('c')
            r_base.append(caterpillars[0])
        elif idx >= 60 and idx < 75:
            rooted_candidates.append(pseudo_caterpillars[idx-60])
            rooted_candidate_types.append('p')
            r_base.append(pseudo_caterpillars[6])
        elif idx >= 75 and idx < 105:
            rooted_candidates.append(balanced[idx-75])
            rooted_candidate_types.append('b')
            r_base.append(balanced[0])

    print("best rooted tree", rooted_candidates[0])
    print("true species tree:", true_species_tree)

    unrooted_tree = dendropy.Tree.get(data=rooted_candidates[0].as_string(schema='newick'), schema='newick',
                                    rooting="force-unrooted", taxon_namespace=tns)

    correct_topology_flag = False
    correct_tree_flag = False
    if dendropy.calculate.treecompare.symmetric_difference(rooted_candidates[0], true_species_tree) == 0:
        correct_tree_flag = True
    if dendropy.calculate.treecompare.symmetric_difference(unrooted_tree, species_tree_toplogy) == 0:
        correct_topology_flag = True

    top_five_flag = False
    for i in range(len(rooted_candidates)):
        print(rooted_candidates[i], score(rooted_candidates[i], r_base[i], u_distribution, tns, quintets, rooted_candidate_types[i]),
              dendropy.calculate.treecompare.symmetric_difference(rooted_candidates[i], true_species_tree), rooted_candidate_types[i])
        if dendropy.calculate.treecompare.symmetric_difference(rooted_candidates[i], true_species_tree) == 0:
              top_five_flag = True


    map = taxon_set_map(r_base[0], rooted_candidates[0], tns)
    mapped_indices = quintets_map(quintets, tns, map)
    #for i in mapped_indices:
    #    print(u_distribution[i])

    print(rooted_candidate_types[0])
    print(int(top_five_flag))
    print(int(correct_tree_flag))
    print(int(correct_topology_flag))
    print(dendropy.calculate.treecompare.symmetric_difference(rooted_candidates[0], true_species_tree))
    return


def quintets_map(quintets, tns, map):
    q_mapped = quintets.clone() # this mean q_mapped is different from quintets!
    mapped_indices = np.zeros(len(quintets), dtype=int)
    for i in range(len(quintets)):
        q_mapped_str = str(quintets[i])
        for j in range(len(map)):
            idx = str(quintets[i]).index(str(tns[j]).replace('\'', ''))
            q_mapped_str = q_mapped_str[:idx] + str(map[j]).replace('\'', '') + q_mapped_str[idx+1:]

        q_mapped[i] = dendropy.Tree.get(data=q_mapped_str+';', schema="newick")
        for k in range(len(quintets)):
            d = dendropy.calculate.treecompare.symmetric_difference(quintets[k], q_mapped[i])
            if d == 0:
                mapped_indices[i] = k
                break
    return mapped_indices

def taxon_set_map(t1, t2, tns):
    map = ['0']*5
    for i in range(len(tns)):
        idx  = str(t1).index(str(tns[i]).replace('\'', ''))
        map[i] = str(t2)[idx]
    return map

def invariant_metric(a, b):
    return np.abs(a - b)# / (a + b)
    #return np.power(a-b, 2)
    #return np.log(max(a, b) / min(a, b))

def inequality_metric(a, b): # it should be a < b
    return (a-b) * (a > b)#a - b# penalty when a > b


# returns the set of invariants of this rooted tree using table 5 of allman's paper
# this doesn't have to be a function, but can be a lookup table as well
'''def invariants(u, indices, type):
    invariants = np.zeros(MAX_NUM_INVARIANT)
    if type == 'c':
        invariants[0] = invariant_metric(u[indices[13]], u[indices[14]])
        invariants[1] = invariant_metric(u[indices[10]], u[indices[14]])
        invariants[2] = invariant_metric(u[indices[9]], u[indices[14]])
        invariants[3] = invariant_metric(u[indices[7]], u[indices[14]])
        invariants[4] = invariant_metric(u[indices[6]], u[indices[14]])
        invariants[5] = invariant_metric(u[indices[5]], u[indices[8]])
        invariants[6] = invariant_metric(u[indices[4]], u[indices[11]])
        invariants[7] = invariant_metric(u[indices[3]], u[indices[12]])
        invariants[8] = invariant_metric(u[indices[1]] + u[indices[8]], u[indices[2]] + u[indices[11]])
        invariants[9] = inequality_metric(u[indices[1]], u[indices[0]])
        invariants[10] = inequality_metric(u[indices[3]], u[indices[0]])
        invariants[11] = inequality_metric(u[indices[4]], u[indices[3]])
        invariants[12] = inequality_metric(u[indices[4]], u[indices[1]])
        invariants[13] = inequality_metric(u[indices[6]], u[indices[4]])
        invariants[14] = inequality_metric(u[indices[1]], u[indices[2]])
        invariants[15] = inequality_metric(u[indices[5]], u[indices[2]])
        invariants[16] = inequality_metric(u[indices[4]], u[indices[5]])
    elif type == 'b':
        invariants[0] = invariant_metric(u[indices[13]], u[indices[14]])
        invariants[1] = invariant_metric(u[indices[10]], u[indices[14]])
        invariants[2] = invariant_metric(u[indices[9]], u[indices[14]])
        invariants[3] = invariant_metric(u[indices[8]], u[indices[11]])
        invariants[4] = invariant_metric(u[indices[7]], u[indices[14]])
        invariants[5] = invariant_metric(u[indices[6]], u[indices[14]])
        invariants[6] = invariant_metric(u[indices[5]], u[indices[11]])
        invariants[7] = invariant_metric(u[indices[4]], u[indices[11]])
        invariants[8] = invariant_metric(u[indices[3]], u[indices[12]])
        invariants[9] = invariant_metric(u[indices[1]], u[indices[2]])
        invariants[10] = inequality_metric(u[indices[1]], u[indices[0]])
        invariants[11] = inequality_metric(u[indices[2]], u[indices[0]])
        #invariants[12] = inequality_metric(u[indices[3]], u[indices[0]])
        invariants[13] = inequality_metric(u[indices[4]], u[indices[1]])
        invariants[14] = inequality_metric(u[indices[4]], u[indices[3]])
        invariants[15] = inequality_metric(u[indices[6]], u[indices[4]])
    elif type == 'p':
        invariants[0] = invariant_metric(u[indices[13]], u[indices[14]])
        invariants[1] = invariant_metric(u[indices[11]], u[indices[14]])
        invariants[2] = invariant_metric(u[indices[9]], u[indices[14]])
        invariants[3] = invariant_metric(u[indices[8]], u[indices[14]])
        invariants[4] = invariant_metric(u[indices[7]], u[indices[10]])
        invariants[5] = invariant_metric(u[indices[6]], u[indices[14]])
        invariants[6] = invariant_metric(u[indices[5]], u[indices[14]])
        invariants[7] = invariant_metric(u[indices[4]], u[indices[14]])
        invariants[8] = invariant_metric(u[indices[3]], u[indices[12]])
        invariants[9] = invariant_metric(u[indices[1]], u[indices[2]])
        invariants[10] = inequality_metric(u[indices[1]], u[indices[0]])
        invariants[11] = inequality_metric(u[indices[3]], u[indices[0]])
        invariants[12] = inequality_metric(u[indices[7]], u[indices[0]])
        invariants[13] = inequality_metric(u[indices[4]], u[indices[1]])
        invariants[14] = inequality_metric(u[indices[4]], u[indices[3]])
        invariants[15] = inequality_metric(u[indices[4]], u[indices[7]])
    else:
        return None
    return invariants'''


def invariants(u, indices, type):
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
        #print([u[indices[i]] for i in c]) #** important
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


def score(r, r_base, u_distribution, tns, quintets, type):
    map = taxon_set_map(r_base, r, tns)
    mapped_indices = quintets_map(quintets, tns, map)
    #return np.sum(invariants(u_distribution, mapped_indices, type))
    return invariants(u_distribution, mapped_indices, type)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, help="Input file name (5-taxon gene trees)",
                        required=True, default=None)

    parser.add_argument("-o", "--output", type=str, help="Output file name (5-taxon species tree)",
                        required=True, default=None)

    parser.add_argument("-m", "--mode", type=str, choices=["n", "c"], help="Scoring algorithm mode, 'n' for naive scoring, 'c' for clustering",
                        required=False, default="n")

    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
