# Quintet Rooting

Quintet Rooting (QR) is a method for rooting species trees from multi-locus datasets that uses phylogenetic invariants and inequalities established by [Allman, Degnan, and Rhodes (J Math Biol, 2011)](https://link.springer.com/article/10.1007/s00285-010-0355-7) under the multi-species coalescent (MSC) model to identify the position of the root in an unrooted species tree. Quintet Rooting is polynomial-time and is especially useful for multi-locus datasets with gene tree discordance due to incomplete lineage sorting (ILS). QR scores different rootings of a given unrooted species tree according to the distribution of unrooted quintets (i.e. 5-leaf trees) induced by a given set of gene trees, and returns the best rooting as well as a ranking over all rooted trees in the search space with a confidence score assigned to each.

## Dependencies
Quintet Rooting is implemented in Python 3. It was developed and tested in Python version 3.6.8 and has the following dependencies:
- [Python 3.x](https://www.python.org)
- [Dendropy 4.x](https://dendropy.org/index.html)
- [Numpy](https://numpy.org)

If you have Python 3 and pip, you can use `pip install -r requirements.txt` to install all dependencies.

## Usage Instructions
Quintet Rooting must be run in a directory containing files in the `./qr` directory. We recommend that you clone the repository and run `quintet_rooting.py` in the base directory.

### Rooting an unrooted species tree
**Input:** A file containing an unrooted species tree (with at least 5 taxa) and a file containing a set of unrooted gene trees, both in newick format (may or may not include branch lengths).

**Output:** A file containing the rooted species tree in newick format, and when run with `-cfs`, an additional file containing a ranking over all rooted trees in the search space sorted according to their confidence scores.
```
$ python3 quintet_rooting.py -t <species-topology.tre> -g <input-genes.tre> -o <output-tree.tre>
```
**Arguments**
- **Required**
```
 -t,  --speciestree        input unrooted species tree in newick format
 -g,  --genetrees          input gene trees in newick format
 -o,  --output             output file containing a rooted species tree
```
- **Optional**
```
 -h,  --help               show this help message and exit
 -sm, --samplingmode       TC for triplet cover, LE for linear encoding
 -c,  --cost               cost function (INQ for inequalities only)
 -cfs,--confidencescore    output confidence scores for each possible rooted tree
 -rs,  --seed              random seed
```
**Example**
The `example` directory contains one example set, containing a 10-taxon avian species tree with 1000 genes (without branch lengths). The commands below show different modes of running Quintet Rooting on this data:
```
$ python3 quintet_rooting.py -t ./example/avian-species-10.tre -g ./example/avian-genes-10.tre -o ./example/avian-rooted-10.tre -sm LE
```
```
$ python3 quintet_rooting.py -t ./example/avian-species-10.tre -g ./example/avian-genes-10.tre -o ./example/avian-rooted-10.tre -sm TC -c INQ
```
### Quintet Sampling Method
Quintet Rooting can run with three different sampling methods. The default version runs in O(n<sup>5</sup>k), where n is the number of taxa and k is the number of gene trees, and exhaustively scores all quintets. The `TC` version scores a subset of quintets that guarantee to cover all triplets and runs in O(n<sup>3</sup>k). The `LE` version uses a linear mapping of edges in the unrooted topology to quintets and runs in O(nk). The default version is generally more accurate, although we recommend using the `LE` version for datasets with more than 20 species, as the others will take longer to run. All three sampling methods have provided reasonably close accuracy under simulations.

## Additional Files
The basic topology of all rooted and unrooted binary 5-leaf trees are provided in the `./qr/topologies` directory (taxa are simply shown with numbers 1-5). The `./qr/rooted_quintet_indices.npy` file contains the set of equivalence classes for the distribution of unrooted gene trees for each of the 105 5-taxon rooted species tree. The `./qr/adr_theory.py` file provides useful functions related to the ADR theory, such as functions for visualizing the partial order of each 5-taxon rooted tree with a hasse diagram. These scripts have an additional dependency on [Graphviz](https://pypi.org/project/graphviz/) and [Matplotlib 3.x](https://matplotlib.org).

## Publication
Yasamin Tabatabaee, Kowshika Sarker, Tandy Warnow, Quintet Rooting: rooting species trees under the multi-species coalescent model, Bioinformatics, Volume 38, Issue Supplement_1, July 2022, Pages i109–i117, https://doi.org/10.1093/bioinformatics/btac224

## Bug Reports
For questions or to report bugs and errors, please contact Yasamin Tabatabaee [(syt3@illinois.edu)](mailto:syt3@illinois.edu).
