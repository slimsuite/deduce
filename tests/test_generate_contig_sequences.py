import csv
import os
import pytest
from deduce_uces.algs.indexed_seq import IndexedSeq
import tempfile

from hypothesis import given
from hypothesis.strategies import composite, integers, lists, text

from deduce_uces.algs.merge_doa import CoreKmer, generate_contig_sequences_indexed
from deduce_uces.algs.uce_params import UCEParams, get_uce_max_gap


@composite
def contig_data(draw, fasta=text(min_size=3, alphabet=["A", "C", "G", "T"])):
    fa = draw(fasta)
    core_kmer_len = draw(integers(min_value=3, max_value=len(fa)))

    core_kmer_positions = draw(
        lists(
            integers(min_value=0, max_value=len(fa) - core_kmer_len), unique=True
        ).map(sorted)
    )

    core_kmers = sorted(
        [
            CoreKmer(
                id=0,
                genome=0,
                contig="chr",
                position=pos,
            )
            for pos in core_kmer_positions
        ]
    )

    min_len = draw(integers(min_value=core_kmer_len, max_value=len(fa)))

    return fa, core_kmers, min_len, core_kmer_len


# Ensure that the sequences generated are always present in the FASTA
@given(contig_data())
def test_generate_contig_sequences_seqs_in_fasta(data):
    fa, core_kmers, min_len, core_kmer_len = data

    params = UCEParams(
        uce_min_length=min_len,
        uce_min_homology=1,
        core_kmer_len=core_kmer_len,
        uce_min_support=1,
        uce_max_gap=get_uce_max_gap(1, min_len, 0),
        core_kmer_mismatches=0,
    )

    with tempfile.NamedTemporaryFile() as fasta:
        fasta.writelines([b">chr\n", fa.encode()])
        fasta.flush()

        index = IndexedSeq(fasta.name, build_index=True)

        assert all(
            s in fa
            for s in generate_contig_sequences_indexed(
                (ck for ck in core_kmers),
                min_len=min_len,
                uce_params=params,
                index=index,
            )
        )


def get_contig_sequence_dataset(filename):
    input_path = os.path.join(os.path.split(__file__)[0], "data", "merge", filename)
    with open(input_path, "r") as f:
        reader = csv.reader(f)
        return [CoreKmer(int(line[0]), 0, "chr1", int(line[0])) for line in reader]


GENERATE_CONTIG_SEQUENCE_TEST_CASES = [
    (11, ["ATCCCCAAATA"]),
    (10, ["ATCCCCAAAT", "TCCCCAAATA"]),
    (
        6,
        [
            "AGGGTT",
            "GGGTTC",
            "GGTTCA",
            "GTTCAG",
            "CCGTAG",
            "CGTAGT",
            "GTAGTA",
            "TAGTAC",
            "ATCCCC",
            "TCCCCA",
            "CCCCAA",
            "CCCAAA",
            "CCAAAT",
            "CAAATA",
        ],
    ),
]

TEST_SEQ = "AGGGTTCAGCGAAAAAGGCCGTAGTACCCCATCCCCAAATA"


@pytest.mark.parametrize(
    "min_len,expected",
    GENERATE_CONTIG_SEQUENCE_TEST_CASES,
)
def test_generate_contig_sequences(min_len, expected):
    contig_sequences = get_contig_sequence_dataset("contig_sequences.csv")

    with tempfile.NamedTemporaryFile() as fasta:
        fasta.writelines([b">chr1\n", TEST_SEQ.encode()])
        fasta.flush()

        index = IndexedSeq(fasta.name, build_index=True)

        params = UCEParams(
            uce_min_length=min_len,
            uce_min_homology=1,
            core_kmer_len=5,
            uce_min_support=1,
            uce_max_gap=get_uce_max_gap(1, min_len, 0),
            core_kmer_mismatches=0,
        )

        assert (
            list(
                generate_contig_sequences_indexed(
                    (x for x in contig_sequences),
                    min_len=min_len,
                    uce_params=params,
                    index=index,
                )
            )
            == expected
        )
