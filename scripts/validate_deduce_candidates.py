import argparse
import csv

from itertools import groupby

from Bio.Seq import Seq

from deduce_uces.io.fasta import read_fasta_sequences
from deduce_uces.io.sam import parse_sam_lazy, is_unmapped


# https://stackoverflow.com/a/3844832
def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate a dedUCE candidates file against a known set of UCEs"
    )

    parser.add_argument(
        "candidates",
        metavar="CANDIDATES",
        type=str,
    )

    parser.add_argument(
        "expected_uces",
        metavar="EXPECTED_UCES",
        type=str,
    )

    return parser.parse_args()


def main():
    args = parse_args()
    expected_uces = list(read_fasta_sequences(args.expected_uces))

    n_matched = 0
    n_unmatched = 0

    candidates = [candidate for candidate in parse_sam_lazy(args.candidates)]

    uces = set(u.seq for u in expected_uces)
    uces_matched = set()

    for candidate in candidates:
        did_match = False

        for uce in expected_uces:
            if (
                candidate.sequence in uce.seq
                or Seq(candidate.sequence).reverse_complement() in uce.seq
            ):
                did_match = True
                n_matched += 1
                uces_matched.add(uce.seq)
                break

        if not did_match:
            n_unmatched += 1

    print(f"Candidates matched to expected UCEs: {n_matched}/{n_matched + n_unmatched}")
    print(f"Candidates not matched: {n_unmatched}/{n_matched + n_unmatched}")
    print(f"UCEs not matched: {len(uces - uces_matched)}/{len(uces)}")

    for u in uces - uces_matched:
        print(u)


if __name__ == "__main__":
    main()
