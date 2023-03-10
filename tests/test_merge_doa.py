import csv
import tempfile
from deduce_uces.algs.indexed_seq import IndexedSeq
from deduce_uces.algs.name_mapping import NameMapping, create_name_mapping
from deduce_uces.algs.uce_params import UCEParams, get_uce_max_gap
from deduce_uces.commands.gen_synteny_map import sort_and_group
from deduce_uces.stages.candidates import read_core_kmers_in_regions
import json
import os
import random
from functools import wraps

import pytest

from deduce_uces.algs.merge_doa import (
    CoreKmer,
    GenomeInfo,
    core_kmers_to_uces,
    get_candidate_regions,
)
from deduce_uces.Logger import NullLogger
from deduce_uces.run import ProgramContext
from deduce_uces.stages.trie import create_trie, save_trie

NAME_MAPPING = NameMapping(
    contig_to_uchar={
        "A": {"chr1": 0},
        "B": {"chr1": 0},
        "C": {"chr1": 0},
    },
    sequence_to_uchar={"A": 0, "B": 1, "C": 2},
    uchar_to_contig={
        0: {0: "chr1"},
        1: {0: "chr1"},
        2: {0: "chr1"},
    },
    uchar_to_sequence={0: "A", 1: "B", 2: "C"},
)


def data(dir, name):
    return os.path.abspath(os.path.join("tests", "data", dir, name))


# Based on https://tbenthompson.com/post/how_i_test/
def golden(test_name: str):
    def decorator(sut):
        save = os.environ.get("SAVE_GOLDEN") == "*"

        @wraps(sut)
        def wrapper(*args, **kwargs):
            actual = sut(*args, **kwargs)
            filename = os.path.join("tests", "golden", test_name + ".golden")
            if save:
                with open(filename, "w") as f:
                    f.write(json.dumps(actual))

            with open(filename, "r") as f:
                expected = json.loads(f.read())

            if actual != expected:
                print("=== GOT ===")
                print(actual)

                print("=== WANT ===")
                print(expected)
            assert actual == expected

        return wrapper

    return decorator


@pytest.fixture
def randomseed():
    random.seed(0)


def get_core_kmer_dataset(filename, random_offset=True):
    input_path = os.path.join(os.path.split(__file__)[0], "data", "merge", filename)

    genome_to_int = {"A": 0, "B": 1, "C": 2}
    chr_to_int = {"chr1": 0}

    # To avoid confusing indices and positions in the genome, throw in some randomness here
    # Now any code that mixes up the two will be in all sorts of bother
    random_offset = random.randint(0, 10000) if random_offset else 0
    with open(input_path, "r") as f:
        reader = csv.reader(f)
        return [
            CoreKmer(
                line[0],
                genome_to_int[line[1]],
                chr_to_int[line[2]],
                int(line[3]) + random_offset,
            )
            for line in reader
        ]


# A: GCCACATGGCTTTCTTGTTCTGGTCGGATCCATCGTT
#     =========   ========================
# B: ACCACATGGCAAACTTGTTCTGGTCGGATCCATCGTT
#     =========   ============ ==========
# C: GCCACATGGCGTTCTTGTTCTGGTCCGATCCATCGTA

CORE_KMER_TEST_CASES = [
    (25, 1, 3, 0, []),
    (20, 1, 2, 0, ["CTTGTTCTGGTCGGATCCATCGTT"]),
    (20, 0.95, 3, 0, ["TCTTGTTCTGGTCGGATCCATCGT"]),
    (
        12,
        1,
        3,
        0,
        ["CTTGTTCTGGTC"],
    ),
    (
        9,
        1,
        3,
        0,
        ["CTTGTTCTGGTC", "CCACATGGC", "GATCCATCGT"],
    ),
    (
        9,
        0.90,
        3,
        0,
        ["CTTGTTCTGGTCGGATCCATCGT", "CCACATGGC"],
    ),
]


@pytest.mark.parametrize(
    "uce_min_length,uce_min_homology,uce_min_sequences,region_max_gap,expected",
    CORE_KMER_TEST_CASES,
)
def test_core_kmers_to_uces(
    uce_min_length,
    uce_min_homology,
    uce_min_sequences,
    region_max_gap,
    expected,
):
    core_kmers = get_core_kmer_dataset("core_kmers.csv", random_offset=False)

    with tempfile.TemporaryDirectory() as tmp_path:
        core_kmer_tries = {}

        for genome, cks in sort_and_group(core_kmers, lambda x: x.genome):
            trie_filename = os.path.join(tmp_path, f"{genome}.trie")

            t = create_trie(cks)
            save_trie(t, trie_filename)

            core_kmer_tries[genome] = trie_filename

        context = ProgramContext(
            threads=1, working_dir=os.getcwd(), logger=NullLogger()
        )

        params = UCEParams(
            uce_min_length=uce_min_length,
            uce_min_homology=uce_min_homology,
            core_kmer_len=3,
            uce_min_support=uce_min_sequences,
            uce_max_gap=get_uce_max_gap(
                uce_min_homology, uce_min_length, region_max_gap
            ),
            core_kmer_mismatches=0,
        )

        print(params)

        assert (
            set(
                core_kmers_to_uces(
                    core_kmer_tries,
                    GenomeInfo("A", "a.bam", data("simple", "a.fa")),
                    [
                        GenomeInfo("B", "b.bam", data("simple", "b.fa")),
                        GenomeInfo("C", "c.bam", data("simple", "c.fa")),
                    ],
                    NAME_MAPPING,
                    params,
                    context=context,
                )
            )
            == set(expected)
        )


@pytest.mark.slow
@golden("merge_bej_chrX")
def test_core_kmers_to_uces_bej_chrX():
    context = ProgramContext(
        threads=1, working_dir=data("bej04_chrX", ""), logger=NullLogger()
    )

    genomes = [
        GenomeInfo(
            "hg16", data("bej04_chrX", "hg16.bam"), data("bej04_chrX", "hg16.fa")
        ),
        GenomeInfo("mm3", data("bej04_chrX", "mm3.bam"), data("bej04_chrX", "mm3.fa")),
        GenomeInfo("rn3", data("bej04_chrX", "rn3.bam"), data("bej04_chrX", "rn3.fa")),
    ]

    params = UCEParams(
        uce_min_length=200,
        uce_min_homology=1,
        core_kmer_len=50,
        uce_min_support=3,
        uce_max_gap=get_uce_max_gap(1, 200, 0),
        core_kmer_mismatches=0,
    )

    name_mapping = create_name_mapping(
        {g.name: IndexedSeq(g.fa).get_contigs() for g in genomes}
    )

    core_kmers = read_core_kmers_in_regions(genomes, params, name_mapping, context)

    return sorted(
        list(
            set(
                core_kmers_to_uces(
                    core_kmers,
                    GenomeInfo("hg16", "hg16.bam", data("bej04_chrX", "hg16.fa")),
                    [
                        GenomeInfo("mm3", "mm3.bam", data("bej04_chrX", "mm3.fa")),
                        GenomeInfo("rn3", "rn3.bam", data("bej04_chrX", "rn3.fa")),
                    ],
                    name_mapping,
                    params,
                    context=context,
                )
            )
        )
    )


CANDIDATE_REGIONS_TEST_CASES = [
    ("core_kmer_regions.csv", 3, 0, [(25, 32)]),
    ("core_kmer_regions.csv", 8, 0, [(25, 32)]),
    ("core_kmer_regions.csv", 10, 0, [(25, 32)]),
    ("core_kmer_regions.csv", 12, 0, []),
    ("core_kmer_regions_gap.csv", 2, 0, [(25, 28), (50, 54)]),
    ("core_kmer_regions_gap.csv", 3, 0, [(25, 28), (50, 54)]),
    ("core_kmer_regions_gap.csv", 4, 1, [(25, 28), (50, 54)]),
    ("core_kmer_regions_adjacent.csv", 4, 0, [(25, 27), (50, 54), (58, 62)]),
    ("core_kmer_regions_adjacent.csv", 5, 0, [(25, 27), (50, 54), (58, 62)]),
    ("core_kmer_regions_adjacent.csv", 6, 0, [(50, 54), (58, 62)]),
]


@pytest.mark.parametrize(
    "filename,uce_min_length,region_max_gap,expected",
    CANDIDATE_REGIONS_TEST_CASES,
)
def test_get_candidate_regions(filename, uce_min_length, region_max_gap, expected):
    core_kmers = get_core_kmer_dataset(filename, random_offset=False)

    params = UCEParams(
        uce_min_length=uce_min_length,
        uce_min_homology=1,
        core_kmer_len=3,
        uce_min_support=1,
        uce_max_gap=get_uce_max_gap(1, uce_min_length, region_max_gap),
        core_kmer_mismatches=0,
    )

    regions = [
        (region[0].position, region[-1].position)
        for region in get_candidate_regions(core_kmers, params)
    ]

    assert regions == expected
