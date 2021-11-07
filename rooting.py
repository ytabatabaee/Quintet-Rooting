import sys
import os
import argparse
import dendropy
import numpy as np


def main(args):
    input_path = args.input
    output_path = args.output
    mode = args.mode

    gene_trees = dendropy.TreeList.get(path=input_path, schema='newick', rooting="default-unrooted",
                                        suppress_edge_lengths=False)

    gene_trees.write_to_path(dest='test.tre', schema='newick', suppress_edge_lengths=True,
                                suppress_internal_node_labels=True) # todo this should be done in the data generation file?
                                #or maybe not? since other methods need the rooted version


    #gene_trees[0].encode_bipartitions()
    #gene_trees[1].encode_bipartitions()
    # print(gene_trees[1].as_string(schema='newick'))
    #print(gene_trees[1].as_string(schema='newick'))
    #  = dendropy.calculate.treecompare.find_missing_bipartitions(gene_trees[0], gene_trees[1])
    # d = dendropy.calculate.treecompare.symmetric_difference(gene_trees[0], gene_trees[1])
    #for c1 in c:
    #    print(c1)

    tns = dendropy.TaxonNamespace()
    quintets = dendropy.TreeList.get(path='quintets.tre', schema='newick', taxon_namespace=tns)
    gene_trees = dendropy.TreeList.get(path='test.tre', schema='newick', taxon_namespace=tns)
    u_count = np.zeros(len(quintets))

    for g in gene_trees:
        #print(g)
        for i in range(len(quintets)):
            d = dendropy.calculate.treecompare.symmetric_difference(quintets[i], g)
            if d == 0:
                #print(i)
                #print('found quintet', quintets[i])
                u_count[i] += 1
                break

    u_distribution = u_count / len(gene_trees)
    print(u_distribution)
    print(np.sum(u_distribution))
    # print(d)
    # print(gene_trees[0].as_string(schema='newick'))

    #print(gene_trees.read_from_path(input_path, schema="newick"))
    # print(trees)
    #tree.deroot()
    #leaf_dict = tree.label_to_node(selection='leaves')

    return
    # read data
    # compute gene tree distribution
    # score the distribution
    # return the rooted tree


def root(u_distribution):
    return

def invariants(r):
    # todo returns the set of invariants of this rooted tree using table 5 of allman's paper
    # this doesn't have to be a function, but can be a lookup table as well
    return

def score(r, u_distribution):
    return

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
