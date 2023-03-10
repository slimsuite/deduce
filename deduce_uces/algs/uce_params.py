import math
from typing import NamedTuple


UCEParams = NamedTuple(
    "UCEParams",
    [
        ("uce_min_length", int),
        ("uce_min_homology", float),
        ("uce_min_support", int),
        ("uce_max_gap", int),
        ("core_kmer_len", int),
        ("core_kmer_mismatches", int),
    ],
)


def get_uce_max_gap(
    uce_min_homology: float, uce_min_length: int, region_extra_gap: int
) -> int:
    return (
        uce_min_length
        - math.ceil(uce_min_homology * float(uce_min_length))
        + region_extra_gap
    )


def calculate_uce_params(number_of_genomes: int, args) -> UCEParams:
    uce_min_support = math.ceil(args.core_kmer_threshold * number_of_genomes)

    uce_max_gap = get_uce_max_gap(
        args.uce_core_min_homology, args.uce_min_length, args.core_kmer_max_gap
    )

    return UCEParams(
        uce_min_length=args.uce_min_length,
        uce_min_homology=args.uce_core_min_homology,
        core_kmer_len=args.core_kmer_size,
        uce_min_support=uce_min_support,
        uce_max_gap=uce_max_gap,
        core_kmer_mismatches=args.core_kmer_mismatches,
    )
