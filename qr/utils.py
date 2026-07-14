import numpy as np
import dendropy
import re


def plot_unrooted_gene_dist(u_distribution, ax, title):
    """
    Plot a given probability distribution of unrooted quintet trees
    :param np.ndarray u_distribution: unrooted quintet tree probability distribution
    :param plt.axes ax : plot axis object
    :param str title: plot title
    :rtype: plt.axes
    """
    from matplotlib.pyplot import plt
    labels = ['$T_1$', '$T_2$', '$T_3$', '$T_4$', '$T_5$', '$T_6$', '$T_7$', '$T_8$', '$T_9$', '$T_{10}$', '$T_{11}$',
              '$T_{12}$', '$T_{13}$', '$T_{14}$', '$T_{15}$']
    ax.stem(range(1, 16), u_distribution, markerfmt=' ', linefmt='black', basefmt=" ")
    ax.set_ylabel('Probability')
    ax.set_title(title)
    ax.set_xticks(np.arange(1, 16, 1))
    ax.set_xticklabels(labels)
    return ax


def taxon_set_map(t1, t2, tns):
    """
    Returns a mapping of taxa from tree t2 to tree t1 (base)
    :param dendropy.Tree t1, t2: two trees on with the same taxon namespace
    :param dendropy.TaxonNamespace tns: taxon namespace of trees
    :rtype: list
    """
    taxon_map = ['0'] * len(tns)
    for i in range(len(tns)):
        idx = str(t1).index(str(tns[i]).replace('\'', ''))
        taxon_map[i] = str(t2)[idx]
    return taxon_map


# the code of this function is taken from https://gist.github.com/bgusach/a967e0587d6e01e889fd1d776c5f3729
# thanks to @bgusach
def multireplace(string, replacements, ignore_case=False):
    """
    Given a string and a replacement map, it returns the replaced string.
    :param str string: string to execute replacements on
    :param dict replacements: replacement dictionary {value to find: value to replace}
    :param bool ignore_case: whether the match should be case insensitive
    :rtype: str
    """
    if not replacements:
        # Edge case that'd produce a funny regex and cause a KeyError
        return string

    # If case insensitive, we need to normalize the old string so that later a replacement
    # can be found. For instance with {"HEY": "lol"} we should match and find a replacement for "hey",
    # "HEY", "hEy", etc.
    if ignore_case:
        def normalize_old(s):
            return s.lower()

        re_mode = re.IGNORECASE
    else:
        def normalize_old(s):
            return s

        re_mode = 0

    replacements = {normalize_old(key): val for key, val in replacements.items()}

    # Place longer ones first to keep shorter substrings from matching where the longer ones should take place
    # For instance given the replacements {'ab': 'AB', 'abc': 'ABC'} against the string 'hey abc', it should produce
    # 'hey ABC' and not 'hey ABc'
    rep_sorted = sorted(replacements, key=len, reverse=True)
    rep_escaped = map(re.escape, rep_sorted)

    # Create a big OR regex that matches any of the substrings to replace
    pattern = re.compile("|".join(rep_escaped), re_mode)

    # For each match, look up the new string in the replacements, being the key the normalized old string
    return pattern.sub(lambda match: replacements[normalize_old(match.group(0))], string)


def map_taxon_namespace(string, taxa_labels):
    """
    Given a quintet tree with taxa labels 1-5 as a string and a set of taxa,
    returns a mapped version of strings with taxa in set taxa_labels
    :param str string: string of a quintet tree with taxa 1-5
    :param tuple taxa_labels: taxa labels to be mapped to
    :rtype: str
    """
    taxa_map_dict = dict()
    for i in range(len(taxa_labels)):
        taxa_map_dict[str(i + 1)] = taxa_labels[i]
    return multireplace(string, taxa_map_dict)


def idx_2_unlabeled_topology(idx):
    """
    Given an index of a rooted binary tree (1-105), returns its topological shape
    :param int idx: index of rooted binary tree
    :rtype: str
    """
    if idx < 60:
        return 'c'
    elif 60 <= idx < 75:
        return 'p'
    elif 75 <= idx < 105:
        return 'b'
    return None


def get_quintet_unrooted_index(subtree_u, quintets_u):
    """
    Returns the index of an unrooted quintet tree on a set of taxa
    :param dendropy.Tree subtree_u: an unrooted 5-taxon tree
    :param list quintets_u: list of 15 unrooted quintet trees
    :rtype: int
    """
    idx_u = -1
    for i in range(len(quintets_u)):
        if dendropy.calculate.treecompare.symmetric_difference(quintets_u[i], subtree_u) == 0:
            idx_u = i
            break
    return idx_u


def get_quintet_rooted_index(subtree_r, quintets_r, u_idx):
    """
    Returns the index of a rooted quintet tree
    :param dendropy.Tree subtree_r: an unrooted 5-taxon tree on q_taxa
    :param list quintets_r: list of 105 rooted quintet trees on q_taxa
    :param int u_idx: index of unrooted quintet tree u
    :rtype: int
    """
    from qr.adr_theory import u2r_mapping
    idx_r = -1
    for i in range(7):
        idx = u2r_mapping[u_idx][i]
        if dendropy.calculate.treecompare.symmetric_difference(quintets_r[idx], subtree_r) == 0:
            idx_r = i
            break
    return idx_r


def _quintet_label_map(q_taxa):
    return {q_taxa[i]: str(i + 1) for i in range(5)}


def _leaf_labels(node, label_map=None):
    if label_map is None:
        return frozenset(leaf.taxon.label for leaf in node.leaf_iter())
    return frozenset(label_map[leaf.taxon.label] for leaf in node.leaf_iter())


def unrooted_quintet_signature(tree, label_map=None):
    """
    Returns a canonical split signature for an unrooted 5-taxon tree.
    """
    all_taxa = _leaf_labels(tree.seed_node, label_map)
    splits = set()
    for edge in tree.postorder_edge_iter():
        if edge.head_node is None:
            continue
        side = _leaf_labels(edge.head_node, label_map)
        other = all_taxa - side
        if len(side) in (0, len(all_taxa)) or min(len(side), len(other)) <= 1:
            continue
        canon = side if (len(side), sorted(side)) <= (len(other), sorted(other)) else other
        splits.add(tuple(sorted(canon)))
    return tuple(sorted(splits))


def rooted_quintet_signature(tree, label_map=None):
    """
    Returns a canonical clade signature for a rooted 5-taxon tree.
    """
    n_taxa = len(list(tree.leaf_node_iter()))
    clades = []
    for node in tree.postorder_node_iter():
        if node.is_leaf() or node is tree.seed_node:
            continue
        clade = _leaf_labels(node, label_map)
        if 1 < len(clade) < n_taxa:
            clades.append(tuple(sorted(clade)))
    return tuple(sorted(clades, key=lambda c: (len(c), c)))


def build_unrooted_quintet_lookup(quintets_u):
    return {unrooted_quintet_signature(q): i for i, q in enumerate(quintets_u)}


def build_rooted_quintet_lookup(quintets_r):
    return {rooted_quintet_signature(q): i for i, q in enumerate(quintets_r)}


def precompute_unrooted_split_sets(tree):
    all_taxa = _leaf_labels(tree.seed_node)
    split_sets = []
    for edge in tree.postorder_edge_iter():
        if edge.head_node is None:
            continue
        side = _leaf_labels(edge.head_node)
        if len(side) not in (0, len(all_taxa)):
            split_sets.append(side)
    return split_sets


def precompute_rooted_clade_sets(tree):
    n_taxa = len(list(tree.leaf_node_iter()))
    clade_sets = []
    for node in tree.postorder_node_iter():
        if node.is_leaf() or node is tree.seed_node:
            continue
        clade = _leaf_labels(node)
        if 1 < len(clade) < n_taxa:
            clade_sets.append(clade)
    return clade_sets


def _canonical_split(side, all_taxa):
    other = all_taxa - side
    return side if (len(side), sorted(side)) <= (len(other), sorted(other)) else other


def _root_split(tree, all_taxa):
    children = list(tree.seed_node.child_node_iter())
    if len(children) < 2:
        return None
    return _canonical_split(_leaf_labels(children[0]), all_taxa)


def precompute_rooting_candidate_data(unrooted_tree):
    """
    Returns root-edge split identifiers and rooted clade sets in the same order
    as get_all_rooted_trees(), without materializing every rooted tree.
    """
    tree = dendropy.Tree(unrooted_tree)
    all_taxa = _leaf_labels(tree.seed_node)
    root_splits = []
    rooted_clades = []
    raw_indices = []
    for raw_idx, edge in enumerate(tree.preorder_edge_iter()):
        try:
            tree.reroot_at_edge(edge, update_bipartitions=False)
            root_split = _root_split(tree, all_taxa)
            if root_split is None:
                continue
            root_splits.append(root_split)
            rooted_clades.append(precompute_rooted_clade_sets(tree))
            raw_indices.append(raw_idx)
        except:
            continue

    if root_splits:
        root_splits.pop(0)
        rooted_clades.pop(0)
        raw_indices.pop(0)

    return root_splits, rooted_clades, raw_indices


def materialize_rooted_candidate(unrooted_tree, raw_index):
    tree = dendropy.Tree(unrooted_tree)
    for idx, edge in enumerate(tree.preorder_edge_iter()):
        try:
            tree.reroot_at_edge(edge, update_bipartitions=True)
        except:
            continue
        if idx == raw_index:
            return dendropy.Tree(tree)
    raise ValueError("Root candidate not found in unrooted tree")


def unrooted_quintet_signature_from_splits(split_sets, q_taxa):
    q_set = set(q_taxa)
    label_map = _quintet_label_map(q_taxa)
    splits = set()
    for side in split_sets:
        q_side = side & q_set
        q_other = q_set - q_side
        if len(q_side) in (0, 5) or min(len(q_side), len(q_other)) <= 1:
            continue
        side_labels = tuple(sorted(label_map[t] for t in q_side))
        other_labels = tuple(sorted(label_map[t] for t in q_other))
        splits.add(side_labels if (len(side_labels), side_labels) <= (len(other_labels), other_labels)
                   else other_labels)
    return tuple(sorted(splits))


def rooted_quintet_signature_from_clades(clade_sets, q_taxa):
    q_set = set(q_taxa)
    label_map = _quintet_label_map(q_taxa)
    clades = set()
    for clade in clade_sets:
        q_clade = clade & q_set
        if 1 < len(q_clade) < 5:
            clades.add(tuple(sorted(label_map[t] for t in q_clade)))
    return tuple(sorted(clades, key=lambda c: (len(c), c)))


def get_quintet_unrooted_index_from_splits(split_sets, q_taxa, unrooted_lookup):
    return unrooted_lookup[unrooted_quintet_signature_from_splits(split_sets, q_taxa)]


def get_quintet_rooted_index_from_clades(clade_sets, q_taxa, u_idx, rooted_lookup):
    from qr.adr_theory import u2r_mapping
    rooted_idx = rooted_lookup[rooted_quintet_signature_from_clades(clade_sets, q_taxa)]
    for i in range(7):
        if u2r_mapping[u_idx][i] == rooted_idx:
            return i
    return -1


def gene_tree_distribution(gene_trees, q_taxa, quintets_u, normalized):
    """
    Given a set of gene trees, labels of 5 taxa 'q_taxa' and the set of unrooted
    quintet trees, estimates the quintet distribution on the induced gene subtrees
    on this 5 taxa
    :param list gene_trees: a set of unrooted gene trees
    :param tuple q_taxa: labels of 5 taxa
    :param list quintets_u: list of 15 unrooted quintet trees on q_taxa
    :param normalized: normalization by the number of gene trees having a quintet rather than all gene trees
    :rtype: np.ndarray
    """
    u_count = np.zeros(len(quintets_u))
    for g in gene_trees:
        g_taxa = set([x.taxon.label for x in g.leaf_nodes()])
        if not set(q_taxa).issubset(g_taxa):
            continue
        g_subtree = g.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
        for i in range(len(quintets_u)):
            if dendropy.calculate.treecompare.symmetric_difference(quintets_u[i], g_subtree) == 0:
                u_count[i] += 1
                break
    if normalized and sum(u_count) != 0:
        u_distribution = u_count / sum(u_count)
    else:
        u_distribution = u_count / len(gene_trees)
    return u_distribution
