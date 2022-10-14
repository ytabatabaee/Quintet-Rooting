import numpy as np
from qr.adr_theory import *


def cost(u, indices, tree_shape, cost_func, k, q_size, shape_coef, abratio):
    """
    Given the probability distribution of unrooted quintet trees u,
    the partial order of a tree R, and the type of the fitness function,
    returns the cost Cost(R, u)
    :param shape_coef: coefficient for shape penalty term
    :param k: number of genes
    :param n: number of taxa in the species tree
    :param np.ndarray u: unrooted quintet tree probability distribution
    :param list indices: partial order on tree R in the form of a list of indices
    :param str tree_shape: topological shape of tree R
    :param str cost_func: type of the fitness function
    :rtype: float
    """
    invariant_score = 0
    inequality_score = 0
    equiv_classes, inequalities = get_partial_order(tree_shape)
    est_shape = topological_shape(u, k, q_size)
    # similarity inside equiv classes
    for c in equiv_classes:
        intraclass_sim = 0
        for i in range(len(c)):
            for j in range(len(c)):
                intraclass_sim += invariant_metric(u[indices[c[i]]], u[indices[c[j]]])
        if cost_func == 'star':
            invariant_score += intraclass_sim
        else:
            invariant_score += intraclass_sim / (len(c))

    # distance between equiv classes
    for ineq in inequalities:
        interclass_distance = 0
        for i in equiv_classes[ineq[0]]:
            for j in equiv_classes[ineq[1]]:
                interclass_distance += inequality_metric(u[indices[j]], u[indices[i]])
        if cost_func == 'star':
            inequality_score += interclass_distance
        else:
            inequality_score += interclass_distance / (len(equiv_classes[ineq[0]]))

    return invariant_score * abratio + inequality_score + shape_coef * int(est_shape != tree_shape)


def topological_shape(u, k, q_size):
    u_sorted = np.sort(u)
    threshold = A(k, q_size) # this could correspond to a lower bound on f(r) of the species tree
    if u_sorted[6] - u_sorted[5] < threshold:
        return 'p'
    elif (u_sorted[6] - u_sorted[5] >= threshold) and (u_sorted[8] - u_sorted[7] < threshold):
        return 'b'
    else:
        return 'c'


def A(k, q_size):
    return 2 * np.sqrt(np.log(30 * k * q_size) / (2 * k))


def invariant_metric(a, b):
    return np.abs(a - b)


def inequality_metric(a, b):
    return (a - b) * (a > b)
