import argparse
import time
import dendropy
import numpy as np
import sys
from table_five import TreeSet

from qr.adr_theory import *
from qr.fitness_cost import *
from qr.quintet_sampling import *
from qr.utils import *
from qr.version import __version__


def main(args):
    st_time = time.time()
    script_path = os.path.realpath(__file__).rsplit("/", 1)[0]

    # input args
    species_tree_path = args.speciestree
    gene_tree_path = args.genetrees
    output_path = args.outputtree
    sampling_method = args.samplingmethod.lower()
    random.seed(args.seed)
    cost_func = 'd' if args.legacyqr else args.cost.lower()
    shape_coef = args.coef
    mult_le = args.multiplicity
    abratio = args.abratio

    header = """*********************************
*     QR-STAR & QR v""" + __version__ + """    *
*********************************"""
    sys.stdout.write(header + '\n')


    tns = dendropy.TaxonNamespace()
    unrooted_species = dendropy.Tree.get(path=species_tree_path, schema='newick',
                                         taxon_namespace=tns, rooting="force-unrooted", suppress_edge_lengths=True)
    if len(tns) < 5:
        raise Exception("Species tree " + species_tree_path + " has less than 5 taxa!\n")
    set_recursion_limit_for_taxa(len(tns))
    gene_tree_taxa = collect_newick_leaf_labels(gene_tree_path)
    missing_gene_taxa = set(t.label for t in tns) - gene_tree_taxa
    if missing_gene_taxa:
        sys.stdout.write("Warning: %d species-tree taxa are absent from all gene trees; "
                         "quintets containing them will receive zero counts.\n" % len(missing_gene_taxa))
    gene_trees = TreeSet(gene_tree_path)

    # reading fixed quintet topology files
    tns_base = dendropy.TaxonNamespace()
    unrooted_quintets_base = dendropy.TreeList.get(path=script_path + '/qr/topologies/quintets.tre',
                                                   taxon_namespace=tns_base, schema='newick')
    rooted_quintets_base = dendropy.TreeList(taxon_namespace=tns_base)
    rooted_quintets_base.read(path=script_path + '/qr/topologies/caterpillar.tre', schema='newick',
                              rooting="default-rooted")
    rooted_quintets_base.read(path=script_path + '/qr/topologies/pseudo_caterpillar.tre', schema='newick',
                              rooting="default-rooted")
    rooted_quintets_base.read(path=script_path + '/qr/topologies/balanced.tre', schema='newick',
                              rooting="default-rooted")
    rooted_quintet_indices = np.load(script_path + '/qr/rooted_quintet_indices.npy')
    unrooted_quintet_lookup = build_unrooted_quintet_lookup(unrooted_quintets_base)
    rooted_quintet_lookup = build_rooted_quintet_lookup(rooted_quintets_base)
    rooted_quintet_mask_lookup = build_rooted_quintet_mask_lookup(rooted_quintets_base)
    rooted_quintet_mask_bits_lookup = build_rooted_quintet_mask_bits_lookup(rooted_quintets_base)
    rooted_quintet_local_index = build_rooted_quintet_local_index(rooted_quintet_mask_lookup)
    taxon_set = [t.label for t in tns]

    sys.stdout.write('Loading time: %.2f sec\n' % (time.time() - st_time))
    ss_time = time.time()

    # search space of rooted trees
    rooted_candidate_splits, rooted_candidate_raw_indices = precompute_rooting_candidate_splits(unrooted_species)
    unrooted_species_splits = precompute_unrooted_split_sets(unrooted_species)
    taxon_bit_map = build_taxon_bit_map(taxon_set)
    all_taxa_mask = taxa_mask(taxon_set, taxon_bit_map)
    split_masks = split_set_masks(unrooted_species_splits, taxon_bit_map)
    root_masks = [taxa_mask(split_set, taxon_bit_map) for split_set in rooted_candidate_splits]
    root_split_idxs = np.asarray(root_split_indices(rooted_candidate_splits, unrooted_species_splits,
                                                    frozenset(taxon_set)))
    root_in_split = build_root_in_split_matrix(split_masks, root_masks, all_taxa_mask)
    r_score = np.zeros(len(rooted_candidate_splits))

    sys.stdout.write('Creating search space time: %.2f sec\n' % (time.time() - ss_time))
    sm_time = time.time()

    # set of sampled quintets
    sample_quintet_taxa = []
    if len(taxon_set) == 5 or sampling_method == 'exh':
        sample_quintet_taxa = list(itertools.combinations(taxon_set, 5))
    elif sampling_method == 'tc':
        sample_quintet_taxa = triplet_cover_sample(taxon_set)
    elif sampling_method == 'le':
        sample_quintet_taxa = linear_quintet_encoding_sample(unrooted_species, taxon_set, mult_le)
    elif sampling_method == 'rl':
        sample_quintet_taxa = random_linear_sample(taxon_set)

    sys.stdout.write('Quintet sampling time: %.2f sec\n' % (time.time() - sm_time))
    proc_time = time.time()

    sys.stdout.write("Number of taxa (n): %d\n" % len(tns))
    sys.stdout.write("Number of gene trees (k): %d\n" % len(gene_trees))
    sys.stdout.write("Size of search space (|R|): %d\n" % len(rooted_candidate_splits))
    sys.stdout.write("Size of sampled quintets set (|Q*|): %d\n" % len(sample_quintet_taxa))

    # preprocessing
    quintet_scores = np.zeros((len(sample_quintet_taxa), 7))
    quintet_unrooted_indices = np.zeros(len(sample_quintet_taxa), dtype=int)
    quintet_split_info = []

    for j in range(len(sample_quintet_taxa)):
        q_taxa = sample_quintet_taxa[j]
        quintet_split_info.append(quintet_split_mask_info(q_taxa, split_masks, taxon_bit_map))
        if set(q_taxa).issubset(gene_tree_taxa):
            quintet_counts = np.asarray(gene_trees.tally_single_quintet(q_taxa))
        else:
            quintet_counts = np.zeros(15)
        quintet_normalizer = sum(quintet_counts) if args.normalized else len(gene_trees)
        quintet_tree_dist = quintet_counts
        if quintet_normalizer != 0:
            quintet_tree_dist = quintet_tree_dist / quintet_normalizer
        quintet_unrooted_indices[j] = get_quintet_unrooted_index_from_splits(unrooted_species_splits, q_taxa,
                                                                              unrooted_quintet_lookup)
        quintet_scores[j] = compute_cost_rooted_quintets(quintet_tree_dist, quintet_unrooted_indices[j],
                                                         rooted_quintet_indices, cost_func, len(gene_trees),
                                                         len(sample_quintet_taxa), shape_coef, abratio)

    sys.stdout.write('Preprocessing time: %.2f sec\n' % (time.time() - proc_time))
    sc_time = time.time()

    # computing scores
    for j in range(len(sample_quintet_taxa)):
        r_indices = rooted_quintet_indices_for_all_roots(quintet_split_info[j], root_in_split, root_split_idxs,
                                                         rooted_quintet_mask_bits_lookup,
                                                         rooted_quintet_local_index, quintet_unrooted_indices[j])
        r_score += quintet_scores[j][r_indices]

    min_idx = np.argmin(r_score)
    best_rooted_candidate = materialize_rooted_candidate(unrooted_species, rooted_candidate_raw_indices[min_idx])
    with open(output_path, 'w') as fp:
        fp.write(str(best_rooted_candidate) + ';\n')

    sys.stdout.write('Scoring time: %.2f sec\n' % (time.time() - sc_time))
    sys.stdout.write('Best rooting: \n%s \n' % str(best_rooted_candidate))

    # computing confidence scores
    if args.confidencescore:
        sys.stdout.write('Scores of all rooted trees:\n %s \n' % str(r_score))
        confidence_scores = (np.max(r_score) - r_score) / np.sum(np.max(r_score) - r_score)
        tree_ranking_indices = np.argsort(r_score)
        with open(output_path + ".rank.cfn", 'w') as fp:
            for i in tree_ranking_indices:
                fp.write(str(materialize_rooted_candidate(unrooted_species, rooted_candidate_raw_indices[i])) + ';\n')
                fp.write(str(confidence_scores[i]) + '\n')

    sys.stdout.write('Total execution time: %.2f sec\n' % (time.time() - st_time))


def compute_cost_rooted_quintets(u_distribution, u_idx, rooted_quintet_indices, cost_func, k, q_size, shape_coef, abratio):
    """
    Scores the 7 possible rootings of an unrooted quintet
    :param np.ndarray u_distribution: unrooted quintet tree probability distribution
    :param int u_idx: index of unrooted binary tree
    :param np.ndarray rooted_quintet_indices: indices of partial orders for all rooted quintet trees
    :param str cost_func: type of the fitness function
    :rtype: np.ndarray
    """
    rooted_tree_indices = u2r_mapping[u_idx]
    costs = np.zeros(7)
    for i in range(7):
        idx = rooted_tree_indices[i]
        unlabeled_topology = idx_2_unlabeled_topology(idx)
        indices = rooted_quintet_indices[idx]
        costs[i] = cost(u_distribution, indices, unlabeled_topology, cost_func, k, q_size, shape_coef, abratio)
    return costs


def get_all_rooted_trees(unrooted_tree):
    """
    Generates all the possible rooted trees with a given unrooted topology
    :param dendropy.Tree unrooted_tree: an unrooted tree topology
    :rtype: list
    """
    rooted_candidates = []
    tree = dendropy.Tree(unrooted_tree)
    for edge in tree.preorder_edge_iter():
        try:
            tree.reroot_at_edge(edge, update_bipartitions=True)
            rooted_candidates.append(dendropy.Tree(tree))
        except:
            continue
    # removing duplicates
    rooted_candidates[0].resolve_polytomies(update_bipartitions=True)
    for i in range(1, len(rooted_candidates)):
        if dendropy.calculate.treecompare.symmetric_difference(rooted_candidates[0], rooted_candidates[i]) == 0:
            rooted_candidates.pop(0)
            break
    return rooted_candidates


def set_recursion_limit_for_taxa(n_taxa):
    """
    DendroPy uses recursive copying/traversal internally and can exceed
    Python's default recursion limit on large or highly unbalanced trees.
    """
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 10 * n_taxa + 1000))


def parse_args():
    parser = argparse.ArgumentParser(description=str('== QR-STAR & QR v' + __version__ + ' =='))

    parser.add_argument("-t", "--speciestree", type=str,
                        help="input unrooted species tree in newick format",
                        required=True, default=None)

    parser.add_argument("-g", "--genetrees", type=str,
                        help="input gene trees in newick format",
                        required=True, default=None)

    parser.add_argument("-o", "--outputtree", type=str,
                        help="output file containing a rooted species tree",
                        required=True, default=None)

    parser.add_argument("-sm", "--samplingmethod", type=str,
                        help="quintet sampling method (LE for linear encoding (default), EXH for exhaustive",
                        required=False, default='LE')

    parser.add_argument("-c", "--cost", type=str,
                        help="cost function (STAR for QR-STAR, D for legacy QR)",
                        required=False, default='STAR')

    parser.add_argument("--legacyqr", action='store_true',
                        help="run the original QR method (equivalent to -c D)")

    parser.add_argument("-cfs", "--confidencescore", action='store_true',
                        help="output confidence scores for each possible rooted tree as well as a ranking")

    parser.add_argument("-mult", "--multiplicity", type=int,
                        help="multiplicity (number of quintets mapped to each edge) in QR-LE",
                        required=False, default=1)

    parser.add_argument("-norm", "--normalized", action='store_true',
                        help="normalization for unresolved gene trees or missing taxa",
                        required=False, default=False)

    parser.add_argument("-coef", "--coef", type=float,
                        help="coefficient for shape penalty term in QR-STAR", required=False, default=0)

    parser.add_argument("-abratio", "--abratio", type=float,
                        help="Ratio between invariant and inequality penalties used in QR-STAR", required=False, default=1)

    parser.add_argument("-rs", "--seed", type=int,
                        help="random seed", required=False, default=1234)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main(parse_args())
