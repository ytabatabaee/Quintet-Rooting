import random
import dendropy
import itertools


def linear_quintet_encoding_sample(unrooted_tree, taxon_set):
    """
    Given an unrooted species tree T and its taxa labels, returns a list of quintets
    of taxa corresponding to a linear encoding of T
    :param Tree unrooted_tree: unrooted species tree topology T
    :param list taxon_set: labels of taxa of L(T)
    :rtype: list
    """
    sample_quintet_taxa = []
    tree = dendropy.Tree(unrooted_tree)
    for edge in tree.preorder_edge_iter():
        seed_node = tree.seed_node
        try:
            if edge.is_leaf():
                quintet = []
                tri_node = edge.tail_node
                tree.reroot_at_node(tri_node, update_bipartitions=True)
                tri_partition_taxa = []
                for child in tree.seed_node.child_node_iter():
                    partition = []
                    for c in child.leaf_nodes():
                        partition.append(c.taxon.label)
                    tri_partition_taxa.append(partition)
                if len(tri_partition_taxa) > 2:
                    for partition in tri_partition_taxa:
                        quintet.extend(random.sample(partition, 1))
                    quintet.extend(random.sample([x for x in taxon_set if x not in quintet], 5 - len(quintet)))
                    sample_quintet_taxa.append(tuple(quintet))
                tree.reroot_at_node(seed_node, update_bipartitions=True)

            elif edge.is_internal():
                quintet = []
                adj_edges = edge.get_adjacent_edges()
                tree.reroot_at_edge(edge, update_bipartitions=True)
                four_partition_taxa = []
                for e in adj_edges:
                    partition = []
                    for c in e.head_node.leaf_nodes():
                        partition.append(c.taxon.label)
                    four_partition_taxa.append(partition)
                for i in range(len(four_partition_taxa)):
                    for j in range(i + 1, len(four_partition_taxa)):
                        if set(four_partition_taxa[i]).issubset(set(four_partition_taxa[j])):
                            four_partition_taxa[j] = list(set(four_partition_taxa[j]) - set(four_partition_taxa[i]))
                four_partition_taxa = [p for p in four_partition_taxa if p != []]
                if len(four_partition_taxa) > 2:
                    for partition in four_partition_taxa:
                        quintet.extend(random.sample(partition, 1))
                    quintet.extend(random.sample([x for x in taxon_set if x not in quintet], 5 - len(quintet)))
                    sample_quintet_taxa.append(tuple(quintet))
                tree.reroot_at_node(seed_node, update_bipartitions=True)
        except Exception as ex:
            continue
    return sample_quintet_taxa


def random_linear_sample(taxon_set):
    """
    Samples 2n-3 random quintets of taxa from given labels
    :param list taxon_set: labels of taxa
    :rtype: list
    """
    n = len(taxon_set)
    return random.sample(list(itertools.combinations(taxon_set, 5)), 2 * n - 3)


def triplet_cover_sample(taxon_set):
    """
    Samples C(n, 3) quintets covering all triplets of taxa in the leafset
    :param list taxon_set: labels of taxa
    :rtype: list
    """
    sample_quintet_taxa = []
    all_triplet_taxa = list(itertools.combinations(taxon_set, 3))
    for trip in all_triplet_taxa:
        trip_prime = random.sample([x for x in taxon_set if x not in trip], 2)
        sample_quintet_taxa.append(tuple(list(trip) + list(trip_prime)))
    return sample_quintet_taxa
