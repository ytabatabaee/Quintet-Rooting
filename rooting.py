import sys
import os
import argparse
import treeswift
import numpy as np


def main(args):
    return
  # compute gene tree distribution
  # score the distribution
  # return the rooted tree


def root(u_distribution):
    return

def score(r, u_distribution):
    return

def tree_shape(u_distribution):
    # sort the u values
    # kernel density estimation
    return

def branch_lengths():
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
    return


if __name__ == "__main__":
    #main(parse_args())
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
