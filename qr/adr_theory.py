import dendropy
import numpy as np
from qr.utils import *
import os


def gen_unrooted_gene_dist(x_c, y_c, z_c, tree_shape):
    """
    Given branch lengths of a 5-taxon rooted model tree in coalescent units, generate the true unrooted
    gene tree probability distribution corresponding to it according to Appendix B of ADR paper
    :param float x_c, y_c, z_c: internal branch lengths of the model tree in coalescent units
    :param str tree_shape: topological shape of tree the model tree
    :rtype: np.ndarray
    """
    x = np.exp(-x_c)
    y = np.exp(-y_c)
    z = np.exp(-z_c)
    u = np.zeros(16)

    if tree_shape == 'c':  # ((((a, b):x, c):y, d):z, e)
        u[1] = 1 - 2 / 3 * x - 2 / 3 * y + 1 / 3 * x * y + 1 / 18 * x * (y ** 3) + 1 / 90 * x * (y ** 3) * (z ** 6)
        u[2] = 1 / 3 * y - 1 / 6 * x * y - 1 / 9 * x * (y ** 3) + 1 / 90 * x * (y ** 3) * (z ** 6)
        u[3] = 1 / 3 * y - 1 / 6 * x * y - 1 / 18 * x * (y ** 3) - 2 / 45 * x * (y ** 3) * (z ** 6)
        u[4] = u[13] = 1 / 3 * x - 1 / 3 * x * y + 1 / 18 * x * (y ** 3) + 1 / 90 * x * (y ** 3) * (z ** 6)
        u[5] = u[12] = 1 / 6 * x * y - 1 / 9 * x * (y ** 3) + 1 / 90 * x * (y ** 3) * (z ** 6)
        u[6] = u[9] = 1 / 6 * x * y - 1 / 18 * x * (y ** 3) - 2 / 45 * x * (y ** 3) * (z ** 6)
        u[7] = u[8] = u[10] = u[11] = u[14] = u[15] = 1 / 18 * x * (y ** 3) + 1 / 90 * x * (y ** 3) * (z ** 6)
    elif tree_shape == 'b':  # (((a, b):x, c):y, (d, e):z)
        u[1] = 1 - 2 / 3 * x - 2 / 3 * y * z + 1 / 3 * x * y * z + 1 / 15 * x * (y ** 3) * z
        u[2] = u[3] = 1 / 3 * y * z - 1 / 6 * x * y * z - 1 / 10 * x * (y ** 3) * z
        u[4] = u[13] = 1 / 3 * x - 1 / 3 * x * y * z + 1 / 15 * x * (y ** 3) * z
        u[5] = u[6] = u[9] = u[12] = 1 / 6 * x * y * z - 1 / 10 * x * (y ** 3) * z
        u[7] = u[8] = u[10] = u[11] = u[14] = u[15] = 1 / 15 * x * (y ** 3) * z
    elif tree_shape == 'p':  # (((a, b):x, (d, e):y):z, c)
        u[1] = 1 - 2 / 3 * x - 2 / 3 * y + 4 / 9 * x * y - 2 / 45 * x * y * (z ** 6)
        u[2] = u[3] = 1 / 3 * y - 5 / 18 * x * y + 1 / 90 * x * y * (z ** 6)
        u[4] = u[13] = 1 / 3 * x - 5 / 18 * x * y + 1 / 90 * x * y * (z ** 6)
        u[5] = u[6] = u[7] = u[9] = u[10] = u[12] = u[14] = u[15] = 1 / 18 * x * y + 1 / 90 * x * y * (z ** 6)
        u[8] = u[11] = 1 / 9 * x * y - 2 / 45 * x * y * (z ** 6)

    return u


"""
A mapping between the indices of 15 unrooted binary trees and the 105 rooted binary
trees according to Table 5 in Allman, Degnan, and Rhodes (J Math Biol, 2011)
"""
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


def draw_hasse_diagram(indices, tree_shape, file_name):
    """
    Given the partial order indices and topological shape of a model tree R, write its hasse diagram to a file
    :param list indices: partial order on tree R in the form of a list of indices
    :param str tree_shape: topological shape of R
    :param str file_name: name of the output file
    :rtype: list
    """
    import graphviz
    h = []
    equiv_classes, inequalities = get_partial_order(tree_shape)
    for c1, c2 in inequalities:
        h.append([['u' + str(indices[e] + 1) for e in equiv_classes[c1]],
                  ['u' + str(indices[e] + 1) for e in equiv_classes[c2]]])
    g = graphviz.Digraph('G', filename=file_name)
    for e1, e2 in h:
        g.edge(",".join(e1), ",".join(e2))
    g.view()
    return h


def get_partial_order(tree_shape):
    """
    Given a rooted model tree topological shape, return the basic form of the tree's partial order
    :param str tree_shape: topological shape of the model tree
    :rtype: list, list
    """
    equiv_classes = []
    inequalities = []
    if tree_shape == 'c':
        equiv_classes = [[0], [1], [2], [3, 12], [4, 11], [5, 8], [6, 7, 9, 10, 13, 14]]
        inequalities = [[0, 1], [0, 3], [1, 4], [3, 4], [4, 6], [2, 1], [2, 5], [5, 4]]
    elif tree_shape == 'b':
        equiv_classes = [[0], [1, 2], [3, 12], [4, 5, 8, 11], [6, 7, 9, 10, 13, 14]]
        inequalities = [[0, 1], [0, 2], [1, 3], [2, 3], [3, 4]]
    elif tree_shape == 'p':
        equiv_classes = [[0], [1, 2], [3, 12], [7, 10], [4, 5, 6, 8, 9, 11, 13, 14]]
        inequalities = [[0, 1], [0, 2], [0, 3], [1, 4], [2, 4], [3, 4]]
    return equiv_classes, inequalities


def gen_rooted_quintet_indices():
    tns = dendropy.TaxonNamespace()
    script_path = os.path.realpath(__file__).rsplit("/", 1)[0]
    quintets = dendropy.TreeList.get(path=script_path+'/topologies/quintets.tre', schema='newick', taxon_namespace=tns)
    caterpillars = dendropy.TreeList.get(path=script_path+'/topologies/caterpillar.tre', schema='newick', taxon_namespace=tns)
    pseudo_caterpillars = dendropy.TreeList.get(path=script_path+'/topologies/pseudo_caterpillar.tre', schema='newick',
                                                taxon_namespace=tns)
    balanced = dendropy.TreeList.get(path=script_path+'/topologies/balanced.tre', schema='newick', taxon_namespace=tns)

    rooted_quintet_indices = np.zeros((105, 15), dtype=int)
    for i in range(len(caterpillars)):
        rooted_quintet_indices[i] = get_indices(caterpillars[i], caterpillars[0], tns, quintets, "c")
    for i in range(len(pseudo_caterpillars)):
        rooted_quintet_indices[i + len(caterpillars)] = get_indices(pseudo_caterpillars[i], pseudo_caterpillars[6], tns, quintets, "p")
    for i in range(len(balanced)):
        rooted_quintet_indices[i + len(caterpillars) + len(pseudo_caterpillars)] = get_indices(balanced[i], balanced[0], tns, quintets, "b")
    np.save(script_path+'/rooted_quintet_indices.npy', rooted_quintet_indices)

    return rooted_quintet_indices


def quintets_map(quintets, tns, taxon_map):
    q_mapped = quintets.clone()
    mapped_indices = np.zeros(len(quintets), dtype=int)
    for i in range(len(quintets)):
        q_mapped_str = str(quintets[i])
        for j in range(len(taxon_map)):
            idx = str(quintets[i]).index(str(tns[j]).replace('\'', ''))
            q_mapped_str = q_mapped_str[:idx] + str(taxon_map[j]).replace('\'', '') + q_mapped_str[idx + 1:]

        q_mapped[i] = dendropy.Tree.get(data=q_mapped_str + ';', schema="newick")
        for k in range(len(quintets)):
            d = dendropy.calculate.treecompare.symmetric_difference(quintets[k], q_mapped[i])
            if d == 0:
                mapped_indices[i] = k
                break
    return mapped_indices


def print_partial_order(indices, tree_shape):
    equiv_classes, inequalities = get_partial_order(tree_shape)
    for c in equiv_classes:
        print([indices[i] + 1 for i in c])
    return


def get_indices(r, r_base, tns, quintets, tree_shape):
    taxon_map = taxon_set_map(r_base, r, tns)
    mapped_indices = quintets_map(quintets, tns, taxon_map)
    print_partial_order(mapped_indices, tree_shape)
    return mapped_indices
