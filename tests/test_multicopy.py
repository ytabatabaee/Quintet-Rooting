import argparse
import os
import subprocess
import sys
from pathlib import Path

import dendropy
import pytest

from qrstar import disco
from qrstar import cli
from qrstar.multicopy import (
    MultiCopyInputError,
    build_gene_to_species,
    parse_gene_species_map,
    process_gene_families,
)


SPECIES_TREE = "((A,B),(C,(D,(E,F))));\n"
MULTICOPY_TREE = "(((A|1,B|1),(C|1,(D|1,(E|1,F|1)))),((A|2,B|2),(C|2,(D|2,(E|2,F|2)))));\n"
SINGLECOPY_TREE = "((A,B),(C,(D,(E,F))));\n"


def species_labels():
    return {"A", "B", "C", "D", "E", "F"}


def split_signature(newick):
    tns = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(
        data=newick,
        schema="newick",
        taxon_namespace=tns,
        rooting="force-unrooted",
        suppress_edge_lengths=True,
    )
    all_taxa = frozenset(t.label for t in tns)
    splits = set()
    for edge in tree.postorder_edge_iter():
        if edge.head_node is None or edge.head_node.is_leaf():
            continue
        side = frozenset(leaf.taxon.label for leaf in edge.head_node.leaf_iter())
        if 0 < len(side) < len(all_taxa):
            splits.add(frozenset(side if len(side) <= len(all_taxa - side) else all_taxa - side))
    return splits


def write(path, text):
    path.write_text(text)
    return path


def multicopy_args(**kwargs):
    values = {
        "gene_species_map": None,
        "delimiter": None,
        "nth_delimiter": 1,
    }
    values.update(kwargs)
    return argparse.Namespace(**values)


def test_delimiter_mapping_and_nth_delimiter():
    mapper, message = build_gene_to_species(multicopy_args(delimiter="|", nth_delimiter=2))
    assert mapper("clade|A|copy1") == "clade|A"
    assert "delimiter" in message


def test_mapping_file_precedence_and_conflicts(tmp_path):
    mapping = write(tmp_path / "map.tsv", "A|1 X\n")
    mapper, message = build_gene_to_species(
        multicopy_args(gene_species_map=str(mapping), delimiter="|", nth_delimiter=1)
    )
    assert mapper("A|1") == "X"
    assert "ignoring --delimiter" in message

    duplicate = write(tmp_path / "duplicate.tsv", "A|1 A\nA|1 B\n")
    with pytest.raises(MultiCopyInputError, match="Duplicate mapping"):
        parse_gene_species_map(str(duplicate))


def test_missing_mapping_entry_is_error(tmp_path):
    gene_file = write(tmp_path / "genes.tre", MULTICOPY_TREE)
    out_file = tmp_path / "decomp.tre"
    mapping = {"A|1": "A"}

    def mapper(label):
        if label not in mapping:
            raise MultiCopyInputError("missing %s" % label)
        return mapping[label]

    with pytest.raises(MultiCopyInputError, match="missing"):
        process_gene_families(str(gene_file), out_file, mapper, species_labels())


def test_integrated_disco_decomposition_matches_expected_topology():
    mapper = lambda label: label.split("|")[0]
    decomposed = disco.decompose_gene_family(MULTICOPY_TREE, mapper, family_id="family_1")
    assert len(decomposed) == 2
    expected = split_signature("((A,B),C,(D,(E,F)));")
    assert [sorted(leaf.get_label() for leaf in t.traverse_leaves()) for t in decomposed] == [
        ["A", "B", "C", "D", "E", "F"],
        ["A", "B", "C", "D", "E", "F"],
    ]
    assert [split_signature(t.newick()) for t in decomposed] == [expected, expected]


def test_process_gene_families_summary_and_save_file(tmp_path):
    gene_file = write(tmp_path / "genes.tre", MULTICOPY_TREE)
    out_file = tmp_path / "decomp.tre"
    mapper = lambda label: label.split("|")[0]
    summary = process_gene_families(str(gene_file), out_file, mapper, species_labels())
    assert summary.input_gene_families == 1
    assert summary.multi_copy_input_families == 1
    assert summary.single_copy_input_families == 0
    assert summary.disco_trees_generated == 2
    assert summary.disco_trees_retained == 2
    assert summary.disco_trees_discarded == 0
    assert [record.family_id for record in summary.records] == ["family_1", "family_1"]
    assert len(out_file.read_text().strip().splitlines()) == 2


def test_single_copy_families_work_in_multicopy_mode(tmp_path):
    gene_file = write(tmp_path / "single.tre", SINGLECOPY_TREE)
    out_file = tmp_path / "decomp.tre"
    summary = process_gene_families(str(gene_file), out_file, lambda label: label, species_labels())
    assert summary.single_copy_input_families == 1
    assert summary.multi_copy_input_families == 0
    assert summary.disco_trees_retained == 1


def test_small_family_and_no_usable_tree_errors(tmp_path):
    gene_file = write(tmp_path / "small.tre", "((A|1,B|1),(A|2,B|2));\n")
    out_file = tmp_path / "decomp.tre"
    mapper = lambda label: label.split("|")[0]
    with pytest.raises(MultiCopyInputError, match="produced no DISCO trees"):
        process_gene_families(str(gene_file), out_file, mapper, species_labels())


def test_species_absent_from_species_tree_is_error(tmp_path):
    gene_file = write(tmp_path / "bad_species.tre", "((A|1,B|1),(Z|1,(D|1,(E|1,F|1))));\n")
    out_file = tmp_path / "decomp.tre"
    mapper = lambda label: label.split("|")[0]
    with pytest.raises(MultiCopyInputError, match="absent from the species tree"):
        process_gene_families(str(gene_file), out_file, mapper, species_labels())


def test_identity_mapping_rejects_duplicate_copy_labels(tmp_path):
    gene_file = write(tmp_path / "duplicates.tre", "((A,B),(A,(D,(E,F))));\n")
    out_file = tmp_path / "decomp.tre"
    with pytest.raises(MultiCopyInputError, match="duplicate species labels"):
        process_gene_families(
            str(gene_file),
            out_file,
            lambda label: label,
            species_labels(),
            mapping_is_identity=True,
        )


def test_malformed_newick_is_error(tmp_path):
    gene_file = write(tmp_path / "bad.tre", "((A|1,B|1),(C|1,D|1)\n")
    out_file = tmp_path / "decomp.tre"
    mapper = lambda label: label.split("|")[0]
    with pytest.raises(MultiCopyInputError, match="Malformed Newick|Could not process"):
        process_gene_families(str(gene_file), out_file, mapper, species_labels())


def test_cli_multicopy_end_to_end_and_temp_cleanup(tmp_path, monkeypatch):
    species = write(tmp_path / "species.tre", SPECIES_TREE)
    genes = write(tmp_path / "genes.tre", MULTICOPY_TREE)
    rooted = tmp_path / "rooted.tre"
    temp_disco = tmp_path / "temp_disco.tre"

    def fake_mkstemp(prefix, suffix):
        fd = os.open(temp_disco, os.O_CREAT | os.O_RDWR)
        return fd, str(temp_disco)

    monkeypatch.setattr(cli.tempfile, "mkstemp", fake_mkstemp)
    args = argparse.Namespace(
        speciestree=str(species),
        genetrees=str(genes),
        outputtree=str(rooted),
        samplingmethod="EXH",
        seed=1234,
        legacyqr=False,
        cost="STAR",
        coef=0,
        multiplicity=1,
        abratio=1,
        normalized=False,
        confidencescore=False,
        multicopy=True,
        delimiter="|",
        nth_delimiter=1,
        gene_species_map=None,
        save_disco_trees=None,
    )
    cli.main(args)
    assert rooted.exists()
    assert args.normalized is True
    assert not temp_disco.exists()


def test_cli_save_disco_trees_and_two_step_equivalence(tmp_path):
    species = write(tmp_path / "species.tre", SPECIES_TREE)
    genes = write(tmp_path / "genes.tre", MULTICOPY_TREE)
    direct_rooted = tmp_path / "direct.tre"
    two_step_rooted = tmp_path / "two_step.tre"
    decomp = tmp_path / "decomp.tre"
    repo = Path(__file__).resolve().parents[1]

    subprocess.run(
        [
            sys.executable,
            "-m",
            "qrstar.cli",
            "-t",
            str(species),
            "-g",
            str(genes),
            "--multicopy",
            "--delimiter",
            "|",
            "--save-disco-trees",
            str(decomp),
            "-o",
            str(direct_rooted),
            "-sm",
            "EXH",
        ],
        cwd=repo,
        check=True,
    )
    assert decomp.exists()
    subprocess.run(
        [
            sys.executable,
            "-m",
            "qrstar.cli",
            "-t",
            str(species),
            "-g",
            str(decomp),
            "-norm",
            "-o",
            str(two_step_rooted),
            "-sm",
            "EXH",
        ],
        cwd=repo,
        check=True,
    )
    assert direct_rooted.read_text() == two_step_rooted.read_text()


def test_cli_help_documents_multicopy_options():
    result = subprocess.run(
        [sys.executable, "-m", "qrstar.cli", "-h"],
        cwd=Path(__file__).resolve().parents[1],
        text=True,
        capture_output=True,
        check=True,
    )
    assert "--multicopy" in result.stdout
    assert "--delimiter" in result.stdout
    assert "--gene-species-map" in result.stdout
    assert "--save-disco-trees" in result.stdout


def test_single_copy_regression_against_prechange_cli(tmp_path):
    repo = Path(__file__).resolve().parents[1]
    old_cli = tmp_path / "old_cli.py"
    old_cli.write_text(subprocess.check_output(["git", "show", "HEAD:qrstar/cli.py"], cwd=repo, text=True))
    output_old = tmp_path / "old.tre"
    output_new = tmp_path / "new.tre"
    cmd_args = [
        "-t",
        str(repo / "example/avian-species-10.tre"),
        "-g",
        str(repo / "example/avian-genes-10.tre"),
        "-sm",
        "EXH",
        "-o",
    ]
    subprocess.run([sys.executable, str(old_cli), *cmd_args, str(output_old), "-cfs"], cwd=repo, check=True)
    subprocess.run(
        [sys.executable, "-m", "qrstar.cli", *cmd_args, str(output_new), "-cfs"],
        cwd=repo,
        check=True,
    )
    assert output_new.read_text() == output_old.read_text()
    assert (tmp_path / "new.tre.rank.cfn").read_text() == (tmp_path / "old.tre.rank.cfn").read_text()
