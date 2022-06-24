import numpy as np
from qr.adr_theory import *


def cost(u, indices, tree_shape, cost_func):
    """
    Given the probability distribution of unrooted quintet trees u,
    the partial order of a tree R, and the type of the fitness function,
    returns the cost Cost(R, u)
    :param np.ndarray u: unrooted quintet tree probability distribution
    :param list indices: partial order on tree R in the form of a list of indices
    :param str tree_shape: topological shape of tree R
    :param str cost_func: type of the fitness function
    :rtype: float
    """
    invariant_score = 0
    inequality_score = 0
    equiv_classes, inequalities = get_partial_order(tree_shape)
    # similarity inside equiv classes
    if cost_func == 'inq':
        invariant_score = 0
    else:
        for c in equiv_classes:
            intraclass_sim = 0
            for i in range(len(c)):
                for j in range(len(c)):
                    intraclass_sim += invariant_metric(u[indices[c[i]]], u[indices[c[j]]])
            invariant_score += intraclass_sim / (len(c))
    # distance between equiv classes
    for ineq in inequalities:
        interclass_distance = 0
        for i in equiv_classes[ineq[0]]:
            for j in equiv_classes[ineq[1]]:
                interclass_distance += inequality_metric(u[indices[j]], u[indices[i]])
        if cost_func == 'inq':
            inequality_score += interclass_distance
        else:
            inequality_score += interclass_distance / (len(equiv_classes[ineq[0]]))

    return invariant_score + inequality_score


def invariant_metric(a, b):
    """
    Invariant penalty term
    :param float a, b: probabilities of two trees
    :rtype: float
    """
    return np.abs(a - b)


def inequality_metric(a, b):
    """
    Inequality penalty term
    :param float a, b: probabilities of two trees
    :rtype: float
    """
    return (a - b) * (a > b)
