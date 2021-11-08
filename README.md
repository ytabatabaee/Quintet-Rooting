# Rooting Species Trees using Phylogenetic Invariants

## Dependencies
- Python 3.x
- [Dendropy 4.x](https://dendropy.org/index.html)
- [Numpy](https://numpy.org)

To install the above dependencies, you can use `pip install -r requirements.txt`. 

## Usage
```
$ python3 rooting.py -i input_file -o output_file -m n
```
**Arguments**
```
 -h, --help       show this help message and exit
 -i, --input      input gene trees in newick format
 -o, --output     inferred rooted species tree in newick format
 -m, --mode       'n' stands for naive (scoring function approach) and 'c' stands for clustering
```

