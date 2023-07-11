# Quintet Rooting

**Quintet Rooting (QR)** is a polynomial-time method for rooting species trees from multi-locus datasets. QR is designed based on the theoretical work by [Allman, Degnan, and Rhodes (J Math Biol, 2011)](https://link.springer.com/article/10.1007/s00285-010-0355-7) that prove the identifiability of rooted 5-taxon trees from unrooted gene trees under the multi-species coalescent (MSC) model. QR is especially useful for multi-locus datasets with gene tree discordance due to incomplete lineage sorting (ILS). QR scores different rootings of a given unrooted species tree according to the distribution of unrooted quintets (i.e. 5-leaf trees) induced by a given set of gene trees, and returns the best rooting as well as a ranking over all rooted trees in the search space with a confidence score assigned to each.

**QR-STAR** is a variant of QR that has an additional step for determining the topological shape of each quintet and a different cost function, and is statistically consistent under the multi-species coalescent (MSC) model. See usage instructions for commands to run each variant.

## Dependencies
QR is implemented in Python 3. It was developed and tested in Python version 3.7.0 and has the following dependencies:
- [Python 3.x](https://www.python.org)
- [Dendropy 4.x](https://dendropy.org/index.html)
- [Numpy](https://numpy.org)
- [table-five](https://github.com/RuneBlaze/fifteen)

If you have Python 3 and pip, you can use `pip install -r requirements.txt` to install all dependencies.

## Usage Instructions
We recommend that you clone the repository and run `quintet_rooting.py` in the base directory.

### Rooting an unrooted species tree
**Input:** A file containing an unrooted species tree (with at least 5 taxa) and a file containing a set of unrooted single-copy gene trees, both in newick format (with or without branch lengths).

**Output:** A file containing the rooted species tree in newick format, and when run with `-cfs`, an additional file containing a ranking over all rooted trees in the search space sorted according to their confidence scores.
```
$ python3 quintet_rooting.py -t <species-topology.tre> -g <input-genes.tre> -o <output-tree.tre> -c STAR
```
**Arguments**
- **Required**
```
 -t,  --speciestree        input unrooted species tree in newick format
 -g,  --genetrees          input single-copy gene trees in newick format
 -o,  --output             output file containing a rooted species tree
```
- **Optional**
```
 -h,  --help               show this help message and exit
 -sm, --samplingmode       TC for triplet cover, LE for linear encoding, EXH for exhaustive
 -c,  --cost               cost function (STAR for QR-STAR)
 -cfs, --confidencescore   output confidence scores for each possible rooted tree
 -mult, --multiplicity     multiplicity (number of quintets mapped to each edge) in QR-LE
 -norm, --normalized       using normalization for unresolved gene trees or missing taxa
 -coef, --coef             shape coefficient in QR-STAR
 -abratio, --abratio       ratio of invariants to inequalities in QR-STAR
 -rs,  --seed              random seed
```
**Example**
The `example` directory contains one example set, with a 10-taxon avian species tree and 1000 gene trees (without branch lengths). The commands below show examples of different modes of running QR on this dataset.

QR-STAR in default mode (*recommended*):
```
$ python3 quintet_rooting.py -t ./example/avian-species-10.tre -g ./example/avian-genes-10.tre -o ./example/avian-rooted-10.tre -c STAR
```
QR in exhaustive mode:
```
$ python3 quintet_rooting.py -t ./example/avian-species-10.tre -g ./example/avian-genes-10.tre -o ./example/avian-rooted-10.tre -sm EXH
```

### Quintet Sampling Method
QR and QR-STAR can run with different sampling methods. The default version (called linear encoding) runs in O(nk), where n is the number of taxa and k is the number of gene trees. The exhaustive version runs in O(n<sup>5</sup>k) and scores all quintets. When the number of taxa is small enough (< 30), the exhaustive version could be run and may provide better accuracy compared to the LE sampling.

## Additional Files
The basic topology of all rooted and unrooted binary 5-leaf trees are provided in the `./qr/topologies` directory (taxa are simply shown with numbers 1-5). The `./qr/rooted_quintet_indices.npy` file contains the set of equivalence classes for the distribution of unrooted gene trees for each of the 105 5-taxon rooted species tree. The `./qr/adr_theory.py` file provides useful functions related to the ADR theory, such as functions for visualizing the partial order of each 5-taxon rooted tree with a hasse diagram. These scripts have an additional dependency on [Graphviz](https://pypi.org/project/graphviz/) and [Matplotlib 3.x](https://matplotlib.org).

## Publications
Y. Tabatabaee, K. Sarkar, and T. Warnow (2022). Quintet Rooting: rooting species trees under the multi-species coalescent model, Bioinformatics, Volume 38, Issue Supplement_1, Pages i109–i117, https://doi.org/10.1093/bioinformatics/btac224

Y. Tabatabaee, S. Roch and T. Warnow (2023). Statistically consistent rooting of species trees under the multispecies coalescent model. International Conference on Research in Computational Molecular Biology, Pages 41-57, https://doi.org/10.1101/2022.10.26.513897

### Data Availability
Datasets used in these papers are available in the following repositories: [QR datasets](https://github.com/ytabatabaee/QR-paper) and [QR-STAR datasets](https://github.com/ytabatabaee/QR-STAR-paper)

## Acknowledgements
The algorithm was originally designed by Tandy Warnow and Yasamin Tabatabaee. The code is contributed by Yasamin Tabatabaee, Baqiao Liu and Kowshika Sarker.
