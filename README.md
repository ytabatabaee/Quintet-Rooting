# QR-STAR: rooting species trees under the multi-species coalescent
[![PyPI version](https://img.shields.io/pypi/v/qrstar)](https://pypi.org/project/qrstar/)
[![Python versions](https://img.shields.io/pypi/pyversions/qrstar)](https://pypi.org/project/qrstar/)
[![License](https://img.shields.io/github/license/ytabatabaee/QR-STAR)](https://github.com/ytabatabaee/QR-STAR/blob/main/LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1089%2Fcmb.2023.0185-blue)](https://doi.org/10.1089/cmb.2023.0185)
[![Downloads](https://img.shields.io/pepy/dt/qrstar?label=downloads)](https://pepy.tech/project/qrstar)

**QR-STAR** is a statistically consistent method for rooting species trees given unrooted gene trees under the multispecies coalescent (MSC) model. It can be applied to datasets with gene tree discordance due to incomplete lineage sorting and gene duplication and loss (with [DISCO](https://github.com/jsdoublel/DISCO) integration). QR-STAR scores candidate rootings of an unrooted species tree using the distribution of unrooted quintet gene subtrees and returns the highest-scoring rooted species tree.

> **QR-STAR is scalable to large phylogenomic datasets:** In our benchmark, it rooted species trees with **10,000 species given 1,000 gene trees in on average 21 minutes**.

This repository provides the reference implementation of QR-STAR, introduced in [Tabatabaee et al., RECOMB and *Journal of Computational Biology* (2023)](https://doi.org/10.1089/cmb.2023.0185), as well as the original **Quintet Rooting (QR)** algorithm [Tabatabaee et al., *Bioinformatics* (2022)](https://doi.org/10.1093/bioinformatics/btac224) and the **DISCO+QR** pipeline for multi-copy input [Willson et al., *Bioinformatics Advances* (2023)](https://doi.org/10.1093/bioadv/vbad015). *QR-STAR is the recommended method for all new analyses.*

## Installation
QR-STAR is implemented in Python 3 and can be installed from PyPI:
```
$ python3 -m pip install qrstar
```

To install the development version from this repository:
```
$ git clone https://github.com/ytabatabaee/QR-STAR.git
$ cd QR-STAR
$ python3 -m pip install .
```
To verify successful installation and view command-line options:
```
$ qrstar -h
```

## Usage

**Input:** A file containing a resolved unrooted species tree with at least 5 taxa and a file containing a set of unrooted single- or muli-copy gene trees (may contain missing taxa or polytomies), both in newick format.

**Output:** The rooted species tree in newick format. If `-o/--outputtree` is provided, the tree is written to that file; otherwise, it is printed to standard output. When run with `-cfs` and `-o`, an additional file contains a ranking over all rooted trees in the search space sorted according to their confidence scores.
```
$ qrstar -t <species-topology.tre> -g <input-genes.tre> [-o <output-tree.tre>]
```
**Arguments**
- **Required**
```
 -t,  --speciestree        input unrooted species tree in newick format
 -g,  --genetrees          input single-copy gene trees in newick format
```
- **Optional**
```
 -h,  --help               show this help message and exit
 -o,  --outputtree         output file containing a rooted species tree; stdout if omitted
 -sm, --samplingmethod     TC for triplet cover, LE for linear encoding, EXH for exhaustive
 -c,  --cost               cost function (STAR for QR-STAR default, D for legacy QR)
      --legacyqr           run the original QR method (equivalent to -c D)
 -cfs, --confidencescore   output confidence scores for each possible rooted tree
 -mult, --multiplicity     multiplicity (number of quintets mapped to each edge) in QR-LE
 -norm, --normalized       using normalization for unresolved gene trees or missing taxa
      --multicopy          enable integrated DISCO decomposition of multi-copy gene-family trees
      --delimiter          delimiter used to map gene-copy labels to species labels
      --nth-delimiter      use the species label before the nth delimiter
      --gene-species-map   two-column gene-copy to species mapping file
      --save-disco-trees   write retained DISCO decomposed single-copy trees to this file
 -coef, --coef             shape coefficient in QR-STAR
 -abratio, --abratio       ratio of invariants to inequalities in QR-STAR
 -rs,  --seed              random seed
```

### Multi-copy input

QR-STAR can analyze multi-copy gene-family trees by running the integrated [DISCO](https://github.com/jsdoublel/DISCO) decomposition before the QR-STAR rooting analysis. DISCO roots and decomposes each multi-copy gene family tree, and QR-STAR runs on the resulting single-copy trees with at least 5 taxa. 

When an explicit gene-to-species mapping file is available, use the command:

```
$ qrstar -t <species_tree.tre> -g <gene_families.tre> --multicopy --gene-species-map <gene_to_species.tsv> -o <rooted_tree.tre>
```

The mapping file is whitespace- or tab-delimited with two columns:

```
gene_copy_1    species_A
gene_copy_2    species_A
gene_copy_3    species_B
```

For delimiter-based gene-copy to species mapping, use the command:

```
$ qrstar -t <species_tree.tre> -g <gene_families.tre> --multicopy --delimiter "|"  [-o rooted_tree.tre]
```

For labels such as `Homo_sapiens|ENSG001`, `--delimiter "|"` maps the gene copy to species `Homo_sapiens`. Use `--nth-delimiter` to preserve DISCO's nth-delimiter interpretation when species names contain the delimiter.

If both `--gene-species-map` and `--delimiter` are provided, the explicit mapping file takes precedence. Use `--save-disco-trees decomposed.tre` to retain the decomposed single-copy gene trees produced by DISCO.

## Example
The `example` directory contains three example sets with 10, 100, and 1000 taxon species trees, each with 1000 gene trees. The commands below show examples of running QR-STAR, QR, and DISCO+QR-STAR on these datasets.

QR-STAR in default mode (*recommended*):
```
$ qrstar -t ./example/avian-species-10.tre -g ./example/avian-genes-10.tre -o ./example/avian-rooted-10.tre -cfs
$ qrstar -t ./example/s_tree.trees -g ./example/truegenetrees -o ./example/qrstar_truegenetrees.tre > ./example/qrstar_truegenetrees.log
```
QR-STAR with exhaustive sampling (*$O(n^5)$*):
```
$ qrstar -t ./example/avian-species-10.tre -g ./example/avian-genes-10.tre -o ./example/avian-rooted-10.tre -sm EXH
```
Original QR:
```
$ qrstar -t ./example/avian-species-10.tre -g ./example/avian-genes-10.tre -o ./example/avian-rooted-10.tre --legacyqr
```

DISCO+QR-STAR (*for multi-copy input*):
```
$ qrstar -t ./example/discoqr-species-100.tre -g ./example/discoqr-gene-families-100.tre --multicopy --delimiter "|" -o ./example/discoqr-rooted-100.tre
```

## Publications

Please cite the paper corresponding to the method used in your analysis:

* If you use **QR-STAR**, including the recommended default command, please cite:

  > Y. Tabatabaee, S. Roch, and T. Warnow (2023).
  > “QR-STAR: A polynomial-time statistically consistent method for rooting species trees under the coalescent.”
  > *Journal of Computational Biology*, 30(11): 1146–1181.
  > https://doi.org/10.1089/cmb.2023.0185

* If you use the original **Quintet Rooting (QR)** algorithm, please cite:

  > Y. Tabatabaee, K. Sarkar, and T. Warnow (2022).
  > “Quintet Rooting: Rooting species trees under the multi-species coalescent model.”
  > *Bioinformatics*, 38(Supplement 1): i109–i117.
  > https://doi.org/10.1093/bioinformatics/btac224

* If you use **DISCO + QR-STAR** with `--multicopy`, please cite both the QR-STAR paper above and the following paper:

  > J. Willson, Y. Tabatabaee, B. Liu, and T. Warnow (2023).
  > “DISCO+QR: rooting species trees in the presence of GDL and ILS.”
  > *Bioinformatics Advances*, 3(1): vbad015.
  > https://doi.org/10.1093/bioadv/vbad015

An earlier version of the QR-STAR work appeared at RECOMB 2023:

> Y. Tabatabaee, S. Roch, and T. Warnow (2023).
> “Statistically consistent rooting of species trees under the multispecies coalescent model.”
> *International Conference on Research in Computational Molecular Biology*, pages 41–57.
> Preprint: https://doi.org/10.1101/2022.10.26.513897


### Data Availability
Datasets used in these papers are available in the following repositories: [QR datasets](https://github.com/ytabatabaee/QR-paper),  [QR-STAR datasets](https://github.com/ytabatabaee/QR-STAR-paper), and [DISCO+QR datasets](https://databank.illinois.edu/datasets/IDB-5748609).

## Acknowledgements
The algorithm was originally designed by Tandy Warnow and Yasamin Tabatabaee. The code is contributed by Yasamin Tabatabaee, Baqiao Liu and Kowshika Sarker. Multi-copy preprocessing adapts core DISCO decomposition functionality from the [DISCO](https://github.com/jsdoublel/DISCO) software written by James Willson.
