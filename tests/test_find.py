import os
import shutil
import tempfile

import subprocess
from time import sleep

import pytest

from deduce_uces.Logger import NullLogger
from deduce_uces.cli import make_argparser, parse_args
from deduce_uces.io.bed import read_bed, BedItem
from deduce_uces.io.fasta import read_fasta_sequences
from deduce_uces.main import run_find
import sys


def run_deduce(inputs, raw_args, tmp_path):
    logger = NullLogger()

    parser = make_argparser()

    for input in inputs:
        input_path = os.path.join(os.path.split(__file__)[0], "data", input)
        shutil.copy(input_path, tmp_path)

        # bedtools getfasta is very finicky about invisible '\r' characters...
        if os.name == "nt":
            subprocess.run(
                ["dos2unix", os.path.join(tmp_path, os.path.basename(input_path))],
                check=True,
            )

    current_wd = os.getcwd()
    os.chdir(tmp_path)

    try:
        args = parse_args(
            parser, ["find"] + raw_args + ["--output", "deduce_uces", tmp_path]
        )

        run_find(args, logger)
    finally:
        os.chdir(current_wd)


def run_deduce_fasta(inputs, raw_args):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_deduce(inputs, raw_args, tmp_path)

        uces = [
            uce.seq
            for uce in read_fasta_sequences(
                os.path.join(tmp_path, "deduce_uces.deduce.fa")
            )
        ]

    return uces


def run_deduce_bed(inputs, raw_args):
    with tempfile.TemporaryDirectory() as tmp_path:
        run_deduce(inputs, raw_args + ["--output-format", "bed"], tmp_path)
        beds = {}

        for input in inputs:
            genome_name = os.path.splitext(os.path.basename(input))[0]
            output_path = os.path.join(tmp_path, f"deduce_uces-{genome_name}.bed")

            if os.path.exists(output_path):
                with open(output_path, "r") as f:
                    beds[genome_name] = [b for b in read_bed(f)]
            else:
                beds[genome_name] = []

    positions = {g: [b.start for b in bs] for g, bs in beds.items()}
    scores = {g: [b.tags["score"] for b in bs] for g, bs in beds.items()}

    return beds, positions, scores


def test_small_scale_strict():
    raw_args = [
        "--core-kmer-mapping-limit",
        "1",
        "--uce-min-length",
        "40",
        "--core-kmer-threshold",
        "1",
        "--core-kmer-size",
        "20",
        "--mapper",
        "bowtie",
        "--threads",
        "1",
    ]

    assert run_deduce_fasta(["small/a.fa", "small/b.fa", "small/c.fa"], raw_args) in [
        ["GGAGCCTTATGGCATAGTCGTCCGCGGAGCACTCTGGTAA"],
        ["TTACCAGAGTGCTCCGCGGACGACTATGCCATAAGGCTCC"],
    ]
    assert run_deduce_fasta(
        ["small/a_rc.fa", "small/b.fa", "small/c.fa"], raw_args
    ) in [
        ["GGAGCCTTATGGCATAGTCGTCCGCGGAGCACTCTGGTAA"],
        ["TTACCAGAGTGCTCCGCGGACGACTATGCCATAAGGCTCC"],
    ]


def test_small_scale_strict_bed():
    inputs = ["small/a.fa", "small/b.fa", "small/c.fa"]

    raw_args = [
        "--uce-max-occurrences",
        "1",
        "--uce-min-length",
        "40",
        "--core-kmer-threshold",
        "1",
        "--core-kmer-size",
        "20",
        "--mapper",
        "bowtie",
    ]

    beds, _, _ = run_deduce_bed(inputs, raw_args)

    assert beds == {
        "a": [
            BedItem(
                chromosome="a",
                start=20,
                end=60,
                tags={"name": "0", "score": "100.0"},
            )
        ],
        "b": [
            BedItem(
                chromosome="b",
                start=20,
                end=60,
                tags={"name": "0", "score": "100.0"},
            )
        ],
        "c": [
            BedItem(
                chromosome="c",
                start=18,
                end=58,
                tags={"name": "0", "score": "100.0"},
            )
        ],
    }


def test_small_scale_duplicate_core_relaxed():
    inputs = [
        "small/b.fa",
        "small/c.fa",
        "small/a_duplicate_core.fa",
    ]

    raw_args = [
        "--uce-max-occurrences",
        "1",
        "--uce-min-length",
        "40",
        # RELAXED CORE KMER THRESHOLD
        "--core-kmer-threshold",
        "0.5",
        "--core-kmer-size",
        "20",
        "--mapper",
        "bowtie",
    ]

    uces = run_deduce_fasta(inputs, raw_args)

    assert uces == ["GGAGCCTTATGGCATAGTCGTCCGCGGAGCACTCTGGTAA"]


# Single SNP in A's UCE
def test_small_scale_snp_strict():
    inputs = [
        "small/b.fa",
        "small/c.fa",
        "small/a_snp.fa",
    ]

    raw_args = [
        "--uce-max-occurrences",
        "1",
        "--uce-min-length",
        "40",
        "--core-kmer-threshold",
        "1",
        "--core-kmer-size",
        "20",
        "--mapper",
        "bowtie",
    ]

    uces = run_deduce_fasta(inputs, raw_args)

    assert uces == []


def test_small_scale_snp_relaxed():
    inputs = [
        "small/b.fa",
        "small/c.fa",
        "small/a_snp.fa",
    ]

    raw_args = [
        "--uce-max-occurrences",
        "1",
        "--uce-min-length",
        "40",
        "--core-kmer-threshold",
        "1",
        "--core-kmer-size",
        "15",
        "--uce-core-min-homology",
        "0.9",
        "--reference",
        "b",
        "--mapper",
        "bowtie",
    ]

    uces = run_deduce_fasta(inputs, raw_args)

    assert len(uces) == 1
    assert uces[0] in [
        "GGAGCCTTATGGCATAGTCGTCCGCGGAGCACTCTGGTAAGGG",
        "CCCTTACCAGAGTGCTCCGCGGACGACTATGCCATAAGGCTCC",
    ]


@pytest.mark.skip(reason="non-deterministic, needs fixing")
def test_small_scale_snp_relaxed_bed():
    inputs = [
        "small/b.fa",
        "small/c.fa",
        "small/a_snp.fa",
    ]

    raw_args = [
        "--uce-max-occurrences",
        "1",
        "--uce-min-length",
        "40",
        "--core-kmer-threshold",
        "1",
        "--core-kmer-size",
        "15",
        "--uce-core-min-homology",
        "0.9",
        "--reference",
        "b",
        "--mapper",
        "bowtie",
    ]
    beds, _, _ = run_deduce_bed(inputs, raw_args)

    assert beds == {
        "a_snp": [
            BedItem(
                chromosome="a",
                start=20,
                end=60,
                tags={"name": "0", "score": "97.5"},
            )
        ],
        "b": [
            BedItem(
                chromosome="b",
                start=20,
                end=63,
                tags={"name": "0", "score": "100.0"},
            )
        ],
        "c": [
            BedItem(
                chromosome="c",
                start=18,
                end=61,
                tags={"name": "0", "score": "93.02325581395348"},
            )
        ],
    }


@pytest.mark.skip(reason="non-deterministic, needs fixing")
def test_small_scale_with_extension():
    inputs = [
        "small/a.fa",
        "small/b.fa",
        "small/c.fa",
        "small/d.fa",
        "small/e.fa",
    ]

    raw_args = [
        "--uce-max-occurrences",
        "3",
        "--uce-min-length",
        "40",
        "--core-kmer-threshold",
        "0.5",
        "--core-kmer-size",
        "20",
        "--uce-core-min-homology",
        "0.9",
        "--reference",
        "a",
        "--mapper",
        "bowtie",
    ]

    beds, positions, scores = run_deduce_bed(inputs, raw_args)

    assert len(beds) == 5

    assert positions == {
        "a": [20],
        "b": [20],
        "c": [18],
        "d": [20],
        "e": [18],
    }

    assert scores == {
        "a": ["100.0"],
        "b": ["100.0"],
        "c": ["100.0"],
        "d": ["95.0"],
        "e": ["92.5"],
    }

    # Now with junk
    beds, positions, scores = run_deduce_bed(inputs + ["small/junk.fa"], raw_args)

    assert len(beds) == 6

    assert positions == {
        "a": [20],
        "b": [20],
        "c": [18],
        "d": [18],
        "e": [20],
        "junk": [],
    }

    assert scores == {
        "a": ["100.0"],
        "b": ["100.0"],
        "c": ["100.0"],
        "d": ["95.0"],
        "e": ["92.5"],
        "junk": [],
    }


def test_lower_homology_100():
    inputs = [
        "homology/ref.fa",
        "homology/1.fa",
    ]

    raw_args = [
        "--uce-max-occurrences",
        "1000",
        "--uce-min-length",
        "1000",
        "--core-kmer-threshold",
        "1",
        "--core-kmer-size",
        "15",
        "--uce-core-min-homology",
        "1",
        "--reference",
        "ref",
        "--mapper",
        "bowtie",
    ]

    uces = run_deduce_fasta(inputs, raw_args)

    assert len(uces) == 1
    assert uces[0] in ["G" * 1000, "C" * 1000]


def test_lower_homology_0999():
    inputs = [
        "homology/ref.fa",
        "homology/0999.fa",
    ]

    raw_args = [
        "--uce-max-occurrences",
        "1000",
        "--uce-min-length",
        "1000",
        "--core-kmer-threshold",
        "1",
        "--core-kmer-size",
        "15",
        "--uce-core-min-homology",
        "0.999",
        "--reference",
        "ref",
        "--mapper",
        "bowtie",
    ]

    uces = run_deduce_fasta(inputs, raw_args)

    assert len(uces) == 1
    assert uces[0] in ["G" * 1000, "C" * 1000]


# 100 incorrect bases in a single run
def test_lower_homology_0900_run():
    inputs = [
        "homology/ref.fa",
        "homology/0900.fa",
    ]

    raw_args = [
        "--uce-max-occurrences",
        "1000",
        "--uce-min-length",
        "1000",
        "--core-kmer-threshold",
        "1",
        "--core-kmer-size",
        "15",
        "--uce-core-min-homology",
        "0.900",
        "--reference",
        "ref",
        "--mapper",
        "bowtie",
    ]

    uces = run_deduce_fasta(inputs, raw_args)
    assert len(uces) == 1
    assert uces[0] in ["G" * 1000, "C" * 1000]


# 100 incorrect bases, spread throughout the UCE
def test_lower_homology_0900_dispersed():
    inputs = [
        "homology/ref.fa",
        "homology/0900_dispersed.fa",
    ]

    # Sanity check: should find no UCEs with homology > 90%
    raw_args = [
        "--uce-max-occurrences",
        "1000",
        "--uce-min-length",
        "1000",
        "--core-kmer-threshold",
        "1",
        "--core-kmer-size",
        "15",
        "--uce-core-min-homology",
        "0.901",
        "--reference",
        "ref",
        "--mapper",
        "bowtie",
    ]

    uces = run_deduce_fasta(inputs, raw_args)

    assert uces == []

    # Should find one UCE with homology == 90%
    raw_args = [
        "--uce-max-occurrences",
        "1000",
        "--uce-min-length",
        "1000",
        "--core-kmer-threshold",
        "1",
        "--core-kmer-size",
        "5",
        "--uce-core-min-homology",
        "0.900",
        "--reference",
        "ref",
        "--mapper",
        "bowtie",
    ]

    uces = run_deduce_fasta(inputs, raw_args)
    assert len(uces) == 1
    assert uces[0] in ["G" * 1000, "C" * 1000]
