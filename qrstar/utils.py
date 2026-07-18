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


def collect_newick_leaf_labels(path):
    """
    Returns the union of leaf labels in a Newick file without constructing trees.
    """
    labels = set()
    expect_leaf = True
    with open(path) as fp:
        data = fp.read()

    i = 0
    while i < len(data):
        ch = data[i]
        if ch == '[':
            depth = 1
            i += 1
            while i < len(data) and depth:
                if data[i] == '[':
                    depth += 1
                elif data[i] == ']':
                    depth -= 1
                i += 1
            continue
        if ch in ' \t\r\n':
            i += 1
            continue
        if ch in '(,':
            expect_leaf = True
            i += 1
            continue
        if ch == ')':
            expect_leaf = False
            i += 1
            continue
        if ch == ';':
            expect_leaf = True
            i += 1
            continue
        if not expect_leaf:
            i += 1
            continue

        if ch in '\'"':
            quote = ch
            label = []
            i += 1
            while i < len(data):
                ch = data[i]
                i += 1
                if ch == quote:
                    if i < len(data) and data[i] == quote:
                        label.append(quote)
                        i += 1
                        continue
                    break
                label.append(ch)
            if label:
                labels.add(''.join(label))
        elif ch not in ':,();':
            start = i
            while i < len(data) and data[i] not in ':,();[] \t\r\n':
                i += 1
            if i > start:
                labels.add(data[start:i])
        else:
            i += 1
            continue
        expect_leaf = False
    return labels


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
    from qrstar.adr_theory import u2r_mapping
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


def _precompute_node_masks(tree, taxon_bit_map):
    node_masks = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node_masks[node] = taxon_bit_map[node.taxon.label]
        else:
            mask = 0
            for child in node.child_node_iter():
                mask |= node_masks[child]
            node_masks[node] = mask
    return node_masks


def precompute_unrooted_split_masks(tree, taxon_bit_map, all_taxa_mask):
    node_masks = _precompute_node_masks(tree, taxon_bit_map)
    split_masks = []
    for edge in tree.postorder_edge_iter():
        if edge.head_node is None:
            continue
        side = node_masks[edge.head_node]
        if side not in (0, all_taxa_mask):
            split_masks.append(side)
    return split_masks


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


def precompute_rooting_candidate_splits(unrooted_tree):
    """
    Returns root-edge split identifiers in the same order as get_all_rooted_trees().
    """
    all_taxa = _leaf_labels(unrooted_tree.seed_node)
    root_splits = []
    raw_indices = []
    for raw_idx, edge in enumerate(unrooted_tree.preorder_edge_iter()):
        if edge.head_node is None:
            continue
        side = _leaf_labels(edge.head_node)
        if len(side) in (0, len(all_taxa)):
            continue
        root_splits.append(_canonical_split(side, all_taxa))
        raw_indices.append(raw_idx)

    if root_splits:
        root_splits.pop(0)
        raw_indices.pop(0)

    return root_splits, raw_indices


def precompute_rooting_candidate_masks(unrooted_tree, taxon_bit_map, all_taxa_mask):
    """
    Returns root-edge split masks in the same order as get_all_rooted_trees().
    """
    node_masks = _precompute_node_masks(unrooted_tree, taxon_bit_map)
    root_masks = []
    raw_indices = []
    for raw_idx, edge in enumerate(unrooted_tree.preorder_edge_iter()):
        if edge.head_node is None:
            continue
        side = node_masks[edge.head_node]
        if side in (0, all_taxa_mask):
            continue
        root_masks.append(side)
        raw_indices.append(raw_idx)

    if root_masks:
        root_masks.pop(0)
        raw_indices.pop(0)

    return root_masks, raw_indices


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


_QUINTET_MASK_SIZE = [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
                      1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5]


def _quintet_mask_signature(tree):
    label_map = {str(i + 1): 1 << i for i in range(5)}
    clades = set()
    for node in tree.postorder_node_iter():
        if node.is_leaf() or node is tree.seed_node:
            continue
        mask = 0
        for leaf in node.leaf_iter():
            mask |= label_map[leaf.taxon.label]
        if 1 < _QUINTET_MASK_SIZE[mask] < 5:
            clades.add(mask)
    return tuple(sorted(clades, key=lambda m: (_QUINTET_MASK_SIZE[m], m)))


def _quintet_mask_signature_bits(tree):
    bits = 0
    for mask in _quintet_mask_signature(tree):
        bits |= 1 << mask
    return bits


def build_rooted_quintet_mask_lookup(quintets_r):
    return {_quintet_mask_signature(q): i for i, q in enumerate(quintets_r)}


def build_rooted_quintet_mask_bits_lookup(quintets_r):
    return {_quintet_mask_signature_bits(q): i for i, q in enumerate(quintets_r)}


def build_rooted_quintet_local_index(rooted_quintet_mask_lookup):
    from qrstar.adr_theory import u2r_mapping
    local_index = np.full((len(u2r_mapping), 105), -1, dtype=np.int8)
    for u_idx in range(len(u2r_mapping)):
        for i in range(7):
            local_index[u_idx][u2r_mapping[u_idx][i]] = i
    return local_index


def build_taxon_bit_map(taxon_labels):
    return {label: 1 << i for i, label in enumerate(taxon_labels)}


def taxa_mask(taxa, taxon_bit_map):
    mask = 0
    for taxon in taxa:
        mask |= taxon_bit_map[taxon]
    return mask


def split_set_masks(split_sets, taxon_bit_map):
    return [taxa_mask(split_set, taxon_bit_map) for split_set in split_sets]


def root_split_indices(root_splits, split_sets, all_taxa):
    split_indices = {}
    for i, split_set in enumerate(split_sets):
        split_indices.setdefault(_canonical_split(split_set, all_taxa), i)
    return [split_indices[root_split] for root_split in root_splits]


def root_split_indices_from_masks(split_masks, root_masks, all_taxa_mask):
    split_indices = {}
    for i, split_mask in enumerate(split_masks):
        split_indices.setdefault(split_mask, i)
        split_indices.setdefault(all_taxa_mask ^ split_mask, i)
    return [split_indices[root_mask] for root_mask in root_masks]


def root_split_index_positions(root_split_idxs, split_count):
    positions = [[] for _ in range(split_count)]
    for root_idx, split_idx in enumerate(root_split_idxs):
        positions[split_idx].append(root_idx)
    return [np.asarray(pos, dtype=int) for pos in positions]


def build_root_in_split_matrix(split_masks, root_masks, all_taxa_mask):
    root_in_split = np.zeros((len(split_masks), len(root_masks)), dtype=bool)
    for split_idx, split_mask in enumerate(split_masks):
        root_in_split[split_idx] = [bool(split_mask & root_mask) and
                                    bool(split_mask & (all_taxa_mask ^ root_mask))
                                    for root_mask in root_masks]
    return root_in_split


def quintet_split_mask_info(q_taxa, split_masks, taxon_bit_map):
    full_q_mask = taxa_mask(q_taxa, taxon_bit_map)
    local_bits = [taxon_bit_map[taxon] for taxon in q_taxa]
    info = []
    for split_idx, split_mask in enumerate(split_masks):
        q_intersection = split_mask & full_q_mask
        if q_intersection == 0 or q_intersection == full_q_mask:
            continue
        local_mask = 0
        for i in range(5):
            if q_intersection & local_bits[i]:
                local_mask |= 1 << i
        info.append((split_idx, local_mask))
    return info


def rooted_quintet_indices_for_all_roots(q_split_info, root_in_split, root_split_positions, rooted_quintet_mask_bits_lookup,
                                         rooted_quintet_local_index, u_idx):
    signatures = np.zeros(root_in_split.shape[1], dtype=np.uint64)
    for split_idx, q_mask in q_split_info:
        complement = 31 ^ q_mask
        false_bit = np.uint64(1 << q_mask) if 1 < _QUINTET_MASK_SIZE[q_mask] < 5 else np.uint64(0)
        true_bit = np.uint64(1 << complement) if 1 < _QUINTET_MASK_SIZE[complement] < 5 else np.uint64(0)
        signatures |= np.where(root_in_split[split_idx], true_bit, false_bit)
        root_edge_bit = false_bit | true_bit
        if root_edge_bit:
            signatures[root_split_positions[split_idx]] |= root_edge_bit

    unique_signatures, inverse = np.unique(signatures, return_inverse=True)
    unique_indices = np.fromiter((rooted_quintet_local_index[u_idx][rooted_quintet_mask_bits_lookup[int(sig)]]
                                  for sig in unique_signatures), dtype=np.int8, count=len(unique_signatures))
    return unique_indices[inverse]


def rooted_quintet_index_from_split_masks(q_split_info, split_masks, all_taxa_mask, root_mask, root_split_idx,
                                          rooted_quintet_mask_lookup, rooted_quintet_local_index, u_idx):
    signature = set()
    root_complement = all_taxa_mask ^ root_mask
    for split_idx, q_mask in q_split_info:
        if split_idx == root_split_idx:
            masks = (q_mask, 31 ^ q_mask)
        else:
            split_mask = split_masks[split_idx]
            if (split_mask & root_mask) and (split_mask & root_complement):
                masks = (31 ^ q_mask,)
            else:
                masks = (q_mask,)
        for mask in masks:
            if 1 < _QUINTET_MASK_SIZE[mask] < 5:
                signature.add(mask)

    rooted_idx = rooted_quintet_mask_lookup[tuple(sorted(signature, key=lambda m: (_QUINTET_MASK_SIZE[m], m)))]
    return rooted_quintet_local_index[u_idx][rooted_idx]


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
    from qrstar.adr_theory import u2r_mapping
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
