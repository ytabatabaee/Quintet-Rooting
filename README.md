# QR-STAR

**QR-STAR** is a statistically consistent method for rooting species trees given unrooted gene trees under the multispecies coalescent (MSC) model. It is designed for datasets with gene tree discordance due to incomplete lineage sorting (ILS). QR-STAR scores candidate rootings of an unrooted species tree using the distribution of unrooted quintet gene trees  and returns the highest-scoring rooted species tree.

This repository contains the reference implementation of QR-STAR introduced in [Tabatabaee et al., RECOMB & J. Comp. Biol. (2023)](https://doi.org/10.1101/2022.10.26.513897), as well as the original **Quintet Rooting (QR)** algorithm introduced in [Tabatabaee et al., Bioinformatics (2022)](https://doi.org/10.1093/bioinformatics/btac224). QR-STAR is the recommended method for all new analyses.

## Dependencies
QR-STAR is implemented in Python 3 and has the following dependencies:
- [Python 3.x](https://www.python.org)
- [Dendropy 4.x](https://dendropy.org/index.html)
- [Numpy](https://numpy.org)
- [table-five](https://github.com/RuneBlaze/fifteen)

If you have Python 3 and pip, you can use `pip install -r requirements.txt` to install all dependencies.

## Usage Instructions

**Input:** A file containing an unrooted species tree (with at least 5 taxa) and a file containing a set of unrooted single-copy gene trees, both in newick format (with or without branch lengths, may contain missing taxa).

**Output:** A file containing the rooted species tree in newick format, and when run with `-cfs`, an additional file containing a ranking over all rooted trees in the search space sorted according to their confidence scores.
```
$ python3 quintet_rooting.py -t <species-topology.tre> -g <input-genes.tre> -o <output-tree.tre>
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
 -c,  --cost               cost function (STAR for QR-STAR default, D for legacy QR)
      --legacyqr           run the original QR method (equivalent to -c D)
 -cfs, --confidencescore   output confidence scores for each possible rooted tree
 -mult, --multiplicity     multiplicity (number of quintets mapped to each edge) in QR-LE
 -norm, --normalized       using normalization for unresolved gene trees or missing taxa
 -coef, --coef             shape coefficient in QR-STAR
 -abratio, --abratio       ratio of invariants to inequalities in QR-STAR
 -rs,  --seed              random seed
```
## Example
The `example` directory contains a 10-taxon avian species tree with 1000 gene trees. The commands below show examples of different modes of running QR-STAR and QR on this dataset.

QR-STAR in default mode (*recommended*):
```
$ python3 quintet_rooting.py -t ./example/avian-species-10.tre -g ./example/avian-genes-10.tre -o ./example/avian-rooted-10.tre
```
QR-STAR with exhaustive sampling:
```
$ python3 quintet_rooting.py -t ./example/avian-species-10.tre -g ./example/avian-genes-10.tre -o ./example/avian-rooted-10.tre -sm EXH
```
Original QR:
```
$ python3 quintet_rooting.py -t ./example/avian-species-10.tre -g ./example/avian-genes-10.tre -o ./example/avian-rooted-10.tre --legacyqr
```

## Publications
Y. Tabatabaee, K. Sarkar, and T. Warnow (2022). Quintet Rooting: rooting species trees under the multi-species coalescent model, Bioinformatics, Volume 38, Issue Supplement_1, Pages i109–i117, https://doi.org/10.1093/bioinformatics/btac224

Y. Tabatabaee, S. Roch and T. Warnow (2023). Statistically consistent rooting of species trees under the multispecies coalescent model. International Conference on Research in Computational Molecular Biology, Pages 41-57, https://doi.org/10.1101/2022.10.26.513897

Y. Tabatabaee, S. Roch and T. Warnow (2023). QR-STAR: A polynomial-time statistically consistent method for rooting species trees under the coalescent. Journal of Computational Biology 30.11 (2023): 1146-1181, https://doi.org/10.1089/cmb.2023.0185

### Data Availability
Datasets used in these papers are available in the following repositories: [QR datasets](https://github.com/ytabatabaee/QR-paper) and [QR-STAR datasets](https://github.com/ytabatabaee/QR-STAR-paper)

## Acknowledgements
The algorithm was originally designed by Tandy Warnow and Yasamin Tabatabaee. The code is contributed by Yasamin Tabatabaee, Baqiao Liu and Kowshika Sarker.
