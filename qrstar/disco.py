"""
Integrated DISCO tree decomposition support.

This module adapts the core tree-level DISCO v1.4.1 implementation from:
https://github.com/JSdoubleL/DISCO

Original DISCO software copyright (c) James Willson and contributors.
The upstream DISCO project is distributed under the MIT License. QR-STAR is
distributed under GPLv3; this adapted module is included under QR-STAR's GPLv3
distribution with attribution to the upstream DISCO authors.
"""

import warnings

import treeswift


def read_tree_newick(newick, family_id=None):
    stripped = newick.strip()
    location = "" if family_id is None else " for gene family %s" % family_id
    if not stripped or stripped == ";":
        raise ValueError("Empty gene-family tree%s" % location)
    if not stripped.endswith(";"):
        raise ValueError("Malformed Newick input%s: missing terminating semicolon" % location)
    if stripped.count("(") != stripped.count(")"):
        raise ValueError("Malformed Newick input%s: unbalanced parentheses" % location)
    try:
        tree = treeswift.read_tree_newick(stripped)
    except Exception as exc:
        raise ValueError("Malformed Newick input%s: %s" % (location, exc)) from exc
    if isinstance(tree, list):
        raise ValueError("Empty gene-family tree%s" % location)
    return tree


def leaf_labels(tree):
    return [leaf.get_label() for leaf in tree.traverse_leaves()]


def species_set(tree, gene_to_species):
    return {gene_to_species(label) for label in leaf_labels(tree)}


def unroot(tree):
    """
    Unroots a treeswift tree. Adapted from DISCO's version of treeswift deroot.
    This preserves two-leaf subtrees instead of contracting (A,B) to A.
    """
    if tree.root is None:
        return tree
    if tree.root.num_children() == 2:
        left, right = tree.root.child_nodes()
        if not right.is_leaf():
            right.contract()
        elif not left.is_leaf():
            left.contract()
    tree.is_rooted = False
    return tree


def reroot_on_edge(tree, node):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if not node.is_root():
            if (
                not hasattr(node, "edge_length")
                or node.edge_length is None
                or node.edge_length == 0
            ):
                node.edge_length = 1
            tree.reroot(node, length=node.edge_length / 2)


def get_min_root(tree, gene_to_species=lambda x: x):
    """
    Calculates the DISCO/ASTRAL-Pro minimum duplication-loss-score root.
    Returns the treeswift node corresponding to the best edge, its score, and ties.
    """

    def score(total_set, set1, set2):
        if len(set1.intersection(set2)) != 0:
            if total_set == set1 or total_set == set2:
                if set1 == set2:
                    return 1
                return 2
            return 3
        return 0

    if tree.root.num_children() == 0:
        tree.root.s = {gene_to_species(tree.root.get_label())}
        return tree.root, 0, []

    if tree.root.num_children() != 2:
        reroot_on_edge(tree, tree.root.child_nodes()[0])

    tree.resolve_polytomies()

    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.down = {gene_to_species(node.get_label())}
            node.d_score = 0
        else:
            if node.num_children() != 2:
                raise ValueError("DISCO rooting requires binary nodes after resolving polytomies")
            left, right = node.child_nodes()
            node.down = left.down.union(right.down)
            node.d_score = left.d_score + right.d_score + score(node.down, left.down, right.down)

    min_score, best_root, ties = float("inf"), None, []
    for node in tree.traverse_preorder():
        node.skip = node.is_root()
        if node.is_root():
            root = node

    left, right = root.child_nodes()
    left.up = right.down
    left.u_score = right.d_score
    right.up = left.down
    right.u_score = left.d_score
    left.skip = True
    right.skip = True
    min_score = left.u_score + left.d_score + score(left.up.union(left.down), left.up, left.down)

    if not left.is_leaf():
        best_root = left
    elif not right.is_leaf():
        best_root = right
    else:
        best_root = root
    ties = [best_root]

    for node in tree.traverse_preorder(leaves=False):
        if node.skip:
            continue
        parent = node.get_parent()
        other = parent.child_nodes()[0] if parent.child_nodes()[0] != node else parent.child_nodes()[1]
        node.up = parent.up.union(other.down)
        node.u_score = parent.u_score + other.d_score + score(node.up, parent.up, other.down)
        total_score = node.u_score + node.d_score + score(node.up.union(node.down), node.up, node.down)
        if total_score == min_score:
            ties.append(node)
        if total_score < min_score:
            min_score = total_score
            best_root = node
            ties = [node]
    return best_root, min_score, ties


def tag(tree, gene_to_species=lambda x: x):
    """
    Tags internal nodes as speciation (S) or duplication (D) according to rooting.
    """
    tree.suppress_unifurcations()
    tree.resolve_polytomies()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.s = {gene_to_species(node.get_label())}
            node.n_dup = 0
        else:
            left, right = node.child_nodes()
            node.s = left.s.union(right.s)
            node.n_dup = left.n_dup + right.n_dup
            if len(left.s.intersection(right.s)) == 0:
                node.tag = "S"
            else:
                node.tag = "D"
                node.n_dup += 1
    tree.n_dup = tree.root.n_dup


def decompose(tree):
    """
    Decomposes a rooted, tagged DISCO tree at duplication nodes.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = []
        for node in tree.traverse_postorder(leaves=False):
            if node.tag == "D":
                left, right = node.child_nodes()
                delete = left if len(left.s) < len(right.s) else right
                out.append(tree.extract_subtree(delete))
                out[-1].suppress_unifurcations()
                node.remove_child(delete)
        tree.suppress_unifurcations()
        out.append(tree)
        return out


def relabel(tree, gene_to_species=lambda x: x):
    for leaf in tree.traverse_postorder(internal=False):
        leaf.set_label(gene_to_species(leaf.get_label()))
    return tree


def decompose_gene_family(newick, gene_to_species, family_id=None):
    """
    Runs core DISCO on one gene-family tree and returns decomposed trees.
    """
    tree = read_tree_newick(newick, family_id=family_id)
    root, score, ties = get_min_root(tree, gene_to_species)
    reroot_on_edge(tree, root)
    tag(tree, gene_to_species)
    out = decompose(tree)
    for subtree in out:
        unroot(subtree)
        relabel(subtree, gene_to_species)
        subtree.suppress_unifurcations()
    return out
