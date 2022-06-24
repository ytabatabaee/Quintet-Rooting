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


def gene_tree_distribution(gene_trees, q_taxa, quintets_u):
    """
    Given a set of gene trees, labels of 5 taxa 'q_taxa' and the set of unrooted
    quintet trees, estimates the quintet distribution on the induced gene subtrees
    on this 5 taxa
    :param list gene_trees: a set of unrooted gene trees
    :param tuple q_taxa: labels of 5 taxa
    :param list quintets_u: list of 15 unrooted quintet trees on q_taxa
    :rtype: np.ndarray
    """
    u_count = np.zeros(len(quintets_u))
    for g in gene_trees:
        g_subtree = g.extract_tree_with_taxa_labels(labels=q_taxa, suppress_unifurcations=True)
        for i in range(len(quintets_u)):
            if dendropy.calculate.treecompare.symmetric_difference(quintets_u[i], g_subtree) == 0:
                u_count[i] += 1
                break
    u_distribution = u_count / len(gene_trees)
    return u_distribution
