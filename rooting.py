import sys
import os
import argparse
import dendropy
import numpy as np

MAX_NUM_INVARIANT = 10
CANDIDATE_SIZE = 5


def main(args):
    input_path = args.input
    output_path = args.output
    mode = args.mode

    #gene_trees = dendropy.TreeList.get(path=input_path, schema='newick', rooting="default-unrooted",
    #                                    suppress_edge_lengths=False)

    #gene_trees.write_to_path(dest='test.tre', schema='newick', suppress_edge_lengths=True,
    #                            suppress_internal_node_labels=True)

    tns = dendropy.TaxonNamespace()
    quintets = dendropy.TreeList.get(path='topologies/quintets.tre', schema='newick', taxon_namespace=tns)

    true_species_tree = dendropy.Tree.get(path='species_tree_mapped.tre', schema='newick', taxon_namespace=tns)
    gene_trees = dendropy.TreeList.get(path='avian_genes_mapped.tre', schema='newick', taxon_namespace=tns)
    u_count = np.zeros(len(quintets))

    for g in gene_trees:
        for i in range(len(quintets)):
            d = dendropy.calculate.treecompare.symmetric_difference(quintets[i], g)
            if d == 0:
                u_count[i] += 1
                break

    u_distribution = u_count / len(gene_trees)
    print("estimated gene tree distribution")
    # print(u_distribution)
    #print(np.sum(u_distribution))

    score_v = np.zeros(105)

    caterpillars = dendropy.TreeList.get(path='topologies/caterpillar.tre', schema='newick', taxon_namespace=tns)
    for i in range(len(caterpillars)):
        score_v[i] = score(caterpillars[i], caterpillars[0], u_distribution, tns, quintets, "c")

    pseudo_caterpillars = dendropy.TreeList.get(path='topologies/pseudo_caterpillar.tre', schema='newick', taxon_namespace=tns)
    for i in range(len(pseudo_caterpillars)):
        score_v[len(caterpillars) + i] = score(pseudo_caterpillars[i], pseudo_caterpillars[0], u_distribution, tns, quintets, "p")

    balanced = dendropy.TreeList.get(path='topologies/balanced.tre', schema='newick', taxon_namespace=tns)
    for i in range(len(balanced)):
        score_v[i + len(caterpillars) + len(pseudo_caterpillars)] = score(balanced[i], balanced[0], u_distribution, tns, quintets, "b")

    idx = np.argpartition(score_v, CANDIDATE_SIZE)
    #print(idx)
    min_indices = idx[:CANDIDATE_SIZE]
    # min_indices = np.where(score_v == score_v.min())[0]
    #print(min_indices)

    rooted_candidates = []
    rooted_candidate_types = []
    r_base = []

    for idx in min_indices:
        if idx < 60:
            rooted_candidates.append(caterpillars[idx])
            rooted_candidate_types.append('c')
            r_base.append(caterpillars[0])
        elif idx > 60 and idx < 75:
            rooted_candidates.append(pseudo_caterpillars[idx-60])
            rooted_candidate_types.append('p')
            r_base.append(pseudo_caterpillars[0])
        elif idx > 75 and idx < 105:
            rooted_candidates.append(balanced[idx-75])
            rooted_candidate_types.append('b')
            r_base.append(balanced[0])

    for i in range(len(rooted_candidates)):
        print(rooted_candidates[i], score(rooted_candidates[i], r_base[i], u_distribution, tns, quintets, rooted_candidate_types[i]))
        print("d:", dendropy.calculate.treecompare.symmetric_difference(rooted_candidates[i], true_species_tree))


    print("real species tree:", true_species_tree)
    #print("min score", np.min(score_v))
    #print(score_v)
    # print(sorted(score_v))

    print("infered rooted tree:")
    inffered_rooted_tree = []
    for i in range(len(rooted_candidates)):
        if check_inequalities(rooted_candidates[i], r_base[i], u_distribution, tns, quintets, rooted_candidate_types[i]):
            inffered_rooted_tree.append(rooted_candidates[i])
            print(rooted_candidates[i])
            print(dendropy.calculate.treecompare.symmetric_difference(rooted_candidates[i], true_species_tree))
            print(score(rooted_candidates[i], r_base[i], u_distribution, tns, quintets, rooted_candidate_types[i]))

    print(len(inffered_rooted_tree))

    return


def quintets_map(quintets, tns, map):
    q_mapped = quintets.clone()
    mapped_indices = np.zeros(len(quintets), dtype=int)

    for i in range(len(quintets)):
        q_mapped_str = str(quintets[i])
        for j in range(len(map)):
            idx = str(q_mapped[i]).index(str(tns[j]).replace('\'', ''))
            q_mapped_str = q_mapped_str[:idx] + str(map[j]).replace('\'', '') + q_mapped_str[idx+1:]

        q_mapped[i] = dendropy.Tree.get(data=q_mapped_str+';', schema="newick")
        for k in range(len(quintets)):
            d = dendropy.calculate.treecompare.symmetric_difference(quintets[k], q_mapped[i])
            if d == 0:
                mapped_indices[i] = k
                break

    return mapped_indices

def check_inequalities(r, r_base, u_distribution, tns, quintets, type):
    map = taxon_set_map(r_base, r, tns)
    indices = quintets_map(quintets, tns, map)
    eq_holds = True
    if type == 'c':
        eq_holds = (u_distribution[indices[0]] > u_distribution[indices[1]]) & (u_distribution[indices[2]] > u_distribution[indices[1]])
    elif type == 'b':
        eq_holds = (u_distribution[indices[0]] > u_distribution[indices[1]])
    elif type == 'p':
        eq_holds = (u_distribution[indices[0]] > u_distribution[indices[1]])
    return eq_holds


def taxon_set_map(t1, t2, tns):
    map = ['0']*5
    for i in range(len(tns)):
        idx  = str(t1).index(str(tns[i]).replace('\'', ''))
        map[i] = str(t2)[idx]
    return map

def root(u_distribution):
    return

# returns the set of invariants of this rooted tree using table 5 of allman's paper
# this doesn't have to be a function, but can be a lookup table as well
def invariants(u, indices, type):
    invariants = np.zeros(MAX_NUM_INVARIANT)
    if type == 'c':
        invariants[0] = u[indices[13]] - u[indices[14]]
        invariants[1] = u[indices[10]] - u[indices[14]]
        invariants[2] = u[indices[9]] - u[indices[14]]
        invariants[3] = u[indices[7]] - u[indices[14]]
        invariants[4] = u[indices[6]] - u[indices[14]]
        invariants[5] = u[indices[5]] - u[indices[8]]
        invariants[6] = u[indices[4]] - u[indices[11]]
        invariants[7] = u[indices[3]] - u[indices[12]]
        invariants[8] = u[indices[1]] - u[indices[2]] + u[indices[8]] - u[indices[11]]
    elif type == 'b':
        invariants[0] = u[indices[13]] - u[indices[14]]
        invariants[1] = u[indices[10]] - u[indices[14]]
        invariants[2] = u[indices[9]] - u[indices[14]]
        invariants[3] = u[indices[8]] - u[indices[11]]
        invariants[4] = u[indices[7]] - u[indices[14]]
        invariants[5] = u[indices[6]] - u[indices[14]]
        invariants[6] = u[indices[5]] - u[indices[11]]
        invariants[7] = u[indices[4]] - u[indices[11]]
        invariants[8] = u[indices[3]] - u[indices[12]]
        invariants[9] = u[indices[1]] - u[indices[2]]
    elif type == 'p':
        invariants[0] = u[indices[13]] - u[indices[14]]
        invariants[1] = u[indices[11]] - u[indices[14]]
        invariants[2] = u[indices[9]] - u[indices[14]]
        invariants[3] = u[indices[8]] - u[indices[14]]
        invariants[4] = u[indices[7]] - u[indices[10]]
        invariants[5] = u[indices[6]] - u[indices[14]]
        invariants[6] = u[indices[5]] - u[indices[14]]
        invariants[7] = u[indices[4]] - u[indices[14]]
        invariants[8] = u[indices[3]] - u[indices[12]]
        invariants[9] = u[indices[1]] - u[indices[2]]
    else:
        return None
    return invariants


'''def invariants(u, indices, type):
    invariants = np.zeros(0)
    equivalence_classes = []
    if type == 'c':
        equivalence_classes = [[0], [1], [2], [3, 12], [4, 11], [5, 8], [6, 7, 9, 10, 13, 14]]
    elif type == 'b':
        equivalence_classes = [[0], [1, 2], [3, 12], [4, 5, 8, 11], [6, 7, 9, 10, 13, 14]]
    elif type == 'p':
        equivalence_classes = [[0], [1, 2], [3, 12], [7, 10], [4, 5, 6, 8, 9, 11, 13, 14]]
    for c in equivalence_classes:
        for i in range(len(c)):
            for j in range(i+1, len(c)):
                invariants = np.append(invariants, u[indices[c[i]]] - u[indices[c[j]]])
    return invariants'''


def score(r, r_base, u_distribution, tns, quintets, type):
    map = taxon_set_map(r_base, r, tns)
    mapped_indices = quintets_map(quintets, tns, map)
    return np.sum(np.absolute(invariants(u_distribution, mapped_indices, type)))


def tree_shape(u_distribution):
    # sort the u values
    # clustering
    # find the clusters in the u-distriution and the sizes of these clusters
    # if smallest class size == 8, then the tree is pseudo_caterpillar
    # if the class of the second smallest prob was 4 then balanced otherwise caterpillar
    return

def branch_lengths():
    # todo: this function doesn't seem to be working for complicated equations
    from sympy import solve, Poly, Eq, Function, exp, Symbol
    from sympy.abc import x, y, z, a, b
    #x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    u = gen_caterpillar()
    print(u*100)
    print(np.sum(u))
    print(solve([1 - 2/3*x - 2/3*y + 1/3*x*y + 1/18*x*y**3 + 1/90*x*y**3*z**6 - u[1],
                1/3*y - 1/6*x*y - 1/9*x*y**3 + 1/90*x*y**3*z**6 - u[2],
                1/3*y - 1/6*x*y - 1/18*x*y**3 - 2/45*x*y**3*z**6 - u[3],
                1/3*x - 1/3*x*y + 1/18*x*y**3 + 1/90*x*y**3*z**6 - u[4],
                1/6*x*y - 1/9*x*y**3 + 1/90*x*y**3*z**6 - u[5],
                1/6*x*y - 1/18*x*y**3 - 2/45*x*y**3*z**6 - u[6],
                1/18*x*y**3 + 1/90*x*y**3*z**6 - u[7]]))
    #print(solve([x**2 - 1,
    #            y-x + 3]))
    #solve([x**2 - 3, y - 1], set=True)
    #print(x, y)
    return

# can use these functions to know how far ui values are from
# where they should be

def gen_caterpillar(x=0.1, y=0.2, z=0.3):
    u = np.zeros(16)
    u[1] = 1 - 2/3*x - 2/3*y + 1/3*x*y + 1/18*x*y**3 + 1/90*x*y**3*z**6
    u[2] = 1/3*y - 1/6*x*y - 1/9*x*y**3 + 1/90*x*y**3*z**6
    u[3] = 1/3*y - 1/6*x*y - 1/18*x*y**3 - 2/45*x*y**3*z**6
    u[4] = u[13] = 1/3*x - 1/3*x*y + 1/18*x*y**3 + 1/90*x*y**3*z**6
    u[5] = u[12] = 1/6*x*y - 1/9*x*y**3 + 1/90*x*y**3*z**6
    u[6] = u[9] = 1/6*x*y - 1/18*x*y**3 - 2/45*x*y**3*z**6
    u[7] = u[8] = u[10] = u[11] = u[14] = u[15] = 1/18*x*y**3 + 1/90*x*y**3*z**6
    return u

def gen_balanced(x=0.1, y=0.2, z=0.3):
    u = np.zeros(16)
    u[1] = 1 - 2/3*x - 2/3*y*z + 1/3*x*y*z + 1/15*x*y**3*z
    u[2] = u[3] = 1/3*y*z - 1/6*x*y*z - 1/10*x*y**3*z
    u[4] = u[13] = 1/3*x - 1/3*x*y*z + 1/15*x*y**3*z
    u[5] = u[6] = u[9] = u[12] = 1/6*x*y*z - 1/10*x*y**3*z
    u[7] = u[8] = u[10] = u[11] = u[14] = u[15] = 1/15*x*y**3*z
    return u

def gen_pseudo_caterpillar(x=0.1, y=0.2, z=0.3):
    u = np.zeros(16)
    u[1] = 1 - 2/3*x - 2/3*y + 4/9*x*y - 2/45*x*y*z**6
    u[2] = u[3] = 1/3*y - 5/18*x*y + 1/90*x*y*z**6
    u[4] = u[13] = 1/3*x - 5/18*x*y + 1/90*x*y*z**6
    u[5] = u[6] = u[7] = u[9] = u[10] = u[12] = u[14] = u[15] = 1/18*x*y + 1/90*x*y*z**6
    u[8] = u[11] = 1/9*x*y - 2/45*x*y*z**6
    return u


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
