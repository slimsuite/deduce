import argparse
import csv

from itertools import groupby

from deduce_uces.io.fasta import read_fasta_sequences
from deduce_uces.io.sam import parse_sam_lazy


# https://stackoverflow.com/a/3844832
def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate a dedUCE UCE file against a known set of UCEs"
    )

    parser.add_argument(
        "actual",
        metavar="ACTUAL",
        type=str,
    )

    parser.add_argument(
        "expected",
        metavar="EXPECTED",
        type=str,
    )

    parser.add_argument("--save-matches", default=False, action="store_true")

    return parser.parse_args()


def seqs_match(actual, expected):
    if len(actual) > len(expected):
        return expected in actual or expected.reverse_complement() in actual

    return actual in expected or actual.reverse_complement() in expected


def main():
    args = parse_args()
    actual = list(read_fasta_sequences(args.actual))
    expected = list(read_fasta_sequences(args.expected))

    mapping = []

    n_actual_unmatched = 0
    actual_but_not_expected = []
    for a in actual:
        did_match = False

        for e in expected:
            if seqs_match(a.seq.upper(), e.seq.upper()):
                did_match = True
                mapping.append([a.name, e.name])
                break

        if not did_match:
            mapping.append([a.name, "-"])
            n_actual_unmatched += 1
            actual_but_not_expected.append(a.seq)

    if args.save_matches:
        with open("matched_candidates.csv", "w") as f:
            w = csv.writer(f)

            w.writerow(["candidate", "uce"])
            w.writerows(mapping)

    expected_names_that_were_mapped = set(x[1] for x in mapping) - {"-"}
    expected_names_that_were_unmapped = (
        set(u.name for u in expected) - expected_names_that_were_mapped
    )

    print(
        f"Expected and present: {len(expected_names_that_were_mapped)}/{len(expected)}"
    )
    print(
        f"Expected, but not present: {len(expected_names_that_were_unmapped)}/{len(expected)}"
    )
    for x in expected_names_that_were_unmapped:
        print(x)
    print(f"Present, but not expected: {n_actual_unmatched }/{len(actual)}")

    # for i, x in enumerate(actual_but_not_expected):
    #     print(f">actual_uce.{i}")
    #     print(x)


if __name__ == "__main__":
    main()
