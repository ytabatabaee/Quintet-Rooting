# Quintet Rooting

Quintet Rooting (QR) is a method for rooting species trees from multi-locus datasets that uses phylogenetic invariants and inequalities established by Allman, Degnan, and Rhodes (J Math Biol, 2011) under the multi-species coalescent (MSC) model to identify the position of the root in an unrooted tree. Quintet Rooting is polynomial-time and is especially useful for multi-locus datasets with gene tree discordance due to incomplete lineage sorting (ILS). Quintet Rooting scores different rootings of a given unrooted species tree according to the distribution of unrooted quintets (i.e. 5-leaf trees) induced by a given set of gene trees, and returns the best rooting as well as a ranking over all rooted trees in the search space with a confidence score assigned to each.

## Contents
- [Dependencies](#dependencies)
- [Usage Instructions](#usage-instructions)
  * [Rooting an unrooted species tree](#rooting-an-unrooted-species-tree)
  * [Quintet Sampling Method](#quintet-sampling-method)
- [Additional Files](#additional-files)

## Dependencies
Quintet Rooting is implemented in Python 3. It was developed and tested in Python version 3.6.8 and has the following dependencies:
- [Python 3.x](https://www.python.org)
- [Dendropy 4.x](https://dendropy.org/index.html)
- [Numpy](https://numpy.org)

If you have Python 3 and pip, you can use `pip install -r requirements.txt` to install the other dependencies.

## Usage Instructions
Quintet Rooting must be run in a directory containing the `rooted_quintet_indices.npy` file and files in the `./topologies` directory. We recommend that you clone the repository and run `quintet_rooting.py` in the base directory.

### Rooting an unrooted species tree
**Input:** A file containing a species tree topology and a file containing a set of gene trees, both in newick format

**Output:** A file containing the rooted species tree in newick format, and when run with `-cfs`, a file containing a ranking over all rooted trees in the search space sorted according to their confidence scores.
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
 -c,  --cost               MC for minimal constraints
 -cfs,--confidencescore    output confidence scores for each possible rooted tree
 -rs,  --seed              random seed
```
**Example**
The `example` directory contains two example sets, one containing a 40-taxon avian species tree with 1000 genes and the other containing a 10-taxon mammalian species tree with 200 genes. The commands below show different modes of running Quintet Rooting on these data:
```
$ python3 quintet_rooting.py -t ./example/avian-species-40.tre -g ./example/avian-genes-40.tre -o ./example/avian-rooted-40.tre -sm LE
```
```
$ python3 quintet_rooting.py -t ./example/mammalian-species-10.tre -g ./example/mammalian-genes-10.tre -o ./example/mammalian-rooted-10.tre -sm TC -c MC
```
### Quintet Sampling Method
Quintet Rooting can run with three different sampling methods. The default version runs in O(n<sup>5</sup>k), where n is the number of taxa and k is the number of gene trees, and exhaustively scores all quintets. The `TC` version scores a subset of quintets that guarantee to cover all triplets and runs in O(n<sup>3</sup>k). The `LE` version uses a linear mapping of edges in the unrooted topology to quintets and and runs in O(nk). The default version is generally more accurate, although we recommond using the `LE` version for datasets with more than 20 species, as the others will take longer to run. All three sampling methods are statistically consistent under the multi-species coalescent model and have provided reasonably close accuracy under simulations.

## Additional Files
The basic topology of all rooted and unrooted 5-leaf trees are provided in the `./topologies` direcotory (taxa are simply shown with numbers 1-5). The `rooted_quintet_indices.npy` file is in standard binary format in Numpy (and therefore can be used with `numpy.load()` function) and contains the set of equivalence classes for the distribution of unrooted gene trees under each 5-taxon rooted species tree.
