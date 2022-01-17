# Quintet Rooting

Quintet Rooting is a method for rooting species trees from multi-locus datasets, which is based on a proof of identifiability of the rooted species tree under the multi-species coalescent (MSC) model established by Allman, Degnan, and Rhodes (J Math Biol, 2011).

## Dependencies
- Python 3.x
- [Dendropy 4.x](https://dendropy.org/index.html)
- [Numpy](https://numpy.org)

To install the above dependencies, you can use `pip install -r requirements.txt`.

## Usage
### Rooting an unrooted species tree
```
$ python3 root.py -t <species-topology.tre> -g <input-genes.tre> -o <output.tre>
```
**Arguments**
```
 -h, --help             show this help message and exit
 -t, --speciestree      input unrooted species tree in newick format
 -g, --genetrees        input gene trees in newick format
 -o, --output           output rooted species tree
```
