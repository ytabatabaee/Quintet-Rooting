import statistics
from dataclasses import dataclass, field
from pathlib import Path

import dendropy


MIN_DISCO_SPECIES = 5


class MultiCopyInputError(ValueError):
    pass


@dataclass
class DiscoTreeRecord:
    family_id: str
    subtree_index: int
    newick: str
    species_count: int


@dataclass
class MultiCopySummary:
    input_gene_families: int = 0
    single_copy_input_families: int = 0
    multi_copy_input_families: int = 0
    disco_trees_generated: int = 0
    disco_trees_retained: int = 0
    disco_trees_discarded: int = 0
    species_counts: list = field(default_factory=list)
    records: list = field(default_factory=list)
    mapping_message: str = ""

    @property
    def median_species(self):
        return statistics.median(self.species_counts) if self.species_counts else None

    @property
    def min_species(self):
        return min(self.species_counts) if self.species_counts else None

    @property
    def max_species(self):
        return max(self.species_counts) if self.species_counts else None


def load_species_labels(species_tree_path):
    tns = dendropy.TaxonNamespace()
    try:
        dendropy.Tree.get(
            path=species_tree_path,
            schema="newick",
            taxon_namespace=tns,
            rooting="force-unrooted",
            suppress_edge_lengths=True,
        )
    except Exception as exc:
        raise MultiCopyInputError("Malformed species-tree Newick input: %s" % exc) from exc
    return {taxon.label for taxon in tns}


def parse_gene_species_map(path):
    mapping = {}
    try:
        with open(path) as fp:
            for line_no, line in enumerate(fp, 1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                fields = stripped.split()
                if len(fields) != 2:
                    raise MultiCopyInputError(
                        "Mapping file %s line %d must contain exactly two columns" % (path, line_no)
                    )
                gene, species = fields
                if gene in mapping:
                    raise MultiCopyInputError(
                        "Duplicate mapping for gene copy %s in %s line %d" % (gene, path, line_no)
                    )
                mapping[gene] = species
    except OSError as exc:
        raise MultiCopyInputError("Could not read gene-species mapping file %s: %s" % (path, exc)) from exc
    if not mapping:
        raise MultiCopyInputError("Gene-species mapping file %s is empty" % path)
    return mapping


def delimiter_mapper(delimiter, nth_delimiter):
    if nth_delimiter < 1:
        raise MultiCopyInputError("--nth-delimiter must be at least 1")

    def mapper(label):
        parts = label.split(delimiter)
        if len(parts) <= nth_delimiter - 1:
            raise MultiCopyInputError(
                "Gene label %s does not contain delimiter %r enough times for --nth-delimiter %d"
                % (label, delimiter, nth_delimiter)
            )
        species = delimiter.join(parts[:nth_delimiter])
        if species == "":
            raise MultiCopyInputError("Gene label %s maps to an empty species name" % label)
        return species

    return mapper


def build_gene_to_species(args):
    if args.gene_species_map:
        mapping = parse_gene_species_map(args.gene_species_map)

        def mapper(label):
            if label not in mapping:
                raise MultiCopyInputError("Gene label %s is missing from mapping file %s" % (label, args.gene_species_map))
            return mapping[label]

        if args.delimiter is not None:
            message = (
                "Gene-copy mapping: using explicit mapping file %s; ignoring --delimiter %r"
                % (args.gene_species_map, args.delimiter)
            )
        else:
            message = "Gene-copy mapping: using explicit mapping file %s" % args.gene_species_map
        return mapper, message

    if args.delimiter is not None:
        return delimiter_mapper(args.delimiter, args.nth_delimiter), (
            "Gene-copy mapping: using delimiter %r with --nth-delimiter %d" % (args.delimiter, args.nth_delimiter)
        )

    return lambda label: label, "Gene-copy mapping: using gene labels as species labels"


def _iter_gene_family_lines(path):
    seen = 0
    try:
        with open(path) as fp:
            for line_no, line in enumerate(fp, 1):
                stripped = line.strip()
                if not stripped:
                    continue
                seen += 1
                yield line_no, stripped
    except OSError as exc:
        raise MultiCopyInputError("Could not read gene-family tree file %s: %s" % (path, exc)) from exc
    if seen == 0:
        raise MultiCopyInputError("Gene-family tree file %s contains no trees" % path)


def _validate_tree_labels(labels, mapped_species, species_tree_labels, mapping_is_identity, family_id):
    if not labels:
        raise MultiCopyInputError("Empty gene-family tree %s" % family_id)
    missing_species = sorted(set(mapped_species) - species_tree_labels)
    if missing_species:
        raise MultiCopyInputError(
            "Gene family %s maps to species absent from the species tree: %s"
            % (family_id, ", ".join(missing_species))
        )
    if mapping_is_identity and len(set(mapped_species)) != len(mapped_species):
        raise MultiCopyInputError(
            "Gene family %s contains duplicate species labels but no --delimiter or --gene-species-map was supplied"
            % family_id
        )


def process_gene_families(
    gene_family_path,
    output_path,
    gene_to_species,
    species_tree_labels,
    mapping_message="",
    minimum_species=MIN_DISCO_SPECIES,
    mapping_is_identity=False,
):
    from qrstar import disco

    summary = MultiCopySummary(mapping_message=mapping_message)
    output_path = Path(output_path)
    with open(output_path, "w") as out:
        for input_index, newick in _iter_gene_family_lines(gene_family_path):
            family_id = "family_%d" % input_index
            try:
                parsed_tree = disco.read_tree_newick(newick, family_id=family_id)
                labels = disco.leaf_labels(parsed_tree)
                mapped_species = [gene_to_species(label) for label in labels]
                _validate_tree_labels(labels, mapped_species, species_tree_labels, mapping_is_identity, family_id)
                input_species_count = len(set(mapped_species))
                if len(labels) == input_species_count:
                    summary.single_copy_input_families += 1
                else:
                    summary.multi_copy_input_families += 1
                decomposed = disco.decompose_gene_family(newick, gene_to_species, family_id=family_id)
            except MultiCopyInputError:
                raise
            except Exception as exc:
                raise MultiCopyInputError("Could not process gene family %s: %s" % (family_id, exc)) from exc

            summary.input_gene_families += 1
            retained_for_family = 0
            for subtree_index, subtree in enumerate(decomposed, 1):
                subtree_species = {leaf.get_label() for leaf in subtree.traverse_leaves()}
                species_count = len(subtree_species)
                summary.disco_trees_generated += 1
                if species_count < minimum_species:
                    summary.disco_trees_discarded += 1
                    continue
                retained_for_family += 1
                summary.disco_trees_retained += 1
                summary.species_counts.append(species_count)
                subtree_newick = subtree.newick()
                summary.records.append(
                    DiscoTreeRecord(
                        family_id=family_id,
                        subtree_index=retained_for_family,
                        newick=subtree_newick,
                        species_count=species_count,
                    )
                )
                out.write(subtree_newick + "\n")
            if retained_for_family == 0:
                raise MultiCopyInputError(
                    "Gene family %s produced no DISCO trees with at least %d species"
                    % (family_id, minimum_species)
                )

    if summary.disco_trees_retained == 0:
        raise MultiCopyInputError("DISCO produced no usable trees with at least %d species" % minimum_species)
    return summary


def write_summary(summary, stream):
    stream.write("Input mode: multi-copy (DISCO)\n")
    stream.write("Input gene families: %d\n" % summary.input_gene_families)
    stream.write("Single-copy input families: %d\n" % summary.single_copy_input_families)
    stream.write("Multi-copy input families: %d\n" % summary.multi_copy_input_families)
    stream.write("DISCO trees generated: %d\n" % summary.disco_trees_generated)
    stream.write("DISCO trees retained (>=5 species): %d\n" % summary.disco_trees_retained)
    stream.write("DISCO trees discarded (<5 species): %d\n" % summary.disco_trees_discarded)
    stream.write("Gene trees supplied to QR-STAR: %d\n" % summary.disco_trees_retained)
    stream.write("Quintet normalization: enabled\n")
    if summary.median_species is not None:
        stream.write("Median species per decomposed tree: %s\n" % summary.median_species)
        stream.write("Minimum species per decomposed tree: %d\n" % summary.min_species)
        stream.write("Maximum species per decomposed tree: %d\n" % summary.max_species)
