import os
from typing import List, NamedTuple

from Bio import SeqIO

from deduce_uces.ext.bedtools import intersect
from deduce_uces.io.bed import read_bed
from deduce_uces.Logger import Logger
from deduce_uces.run import ProgramContext

FastaSeq = NamedTuple(
    "FastaSeq", [("id", str), ("len", int), ("forward", str), ("reverse", str)]
)


def seqs_match(actual: FastaSeq, expected: FastaSeq):
    if actual.len > expected.len:
        smaller = expected
        bigger = actual
    else:
        smaller = actual
        bigger = expected

    return (
        smaller.forward in bigger.forward
        or smaller.reverse in bigger.reverse
        or smaller.forward in bigger.reverse
        or smaller.reverse in bigger.forward
    )


def read_fasta(filename: str) -> List[FastaSeq]:
    return [
        FastaSeq(
            record.id,
            len(record.seq),
            record.seq.upper(),
            record.seq.reverse_complement().upper(),
        )
        for record in SeqIO.parse(filename, "fasta")
    ]


def minimise(seqs: List[FastaSeq]) -> List[FastaSeq]:
    minimal = []

    print("\tsorting...")
    sorted_seqs = list(sorted(seqs, key=lambda s: s.len))

    print("\tminimising...")
    for i in range(len(sorted_seqs)):
        for j in range(i + 1, len(sorted_seqs)):
            if (
                sorted_seqs[i].forward in sorted_seqs[j].forward
                or sorted_seqs[i].reverse in sorted_seqs[j].reverse
            ):
                break
        else:
            minimal.append(sorted_seqs[i])

    print(f"\tminimised: was {len(seqs)} UCEs, now {len(minimal)} UCEs")
    return minimal


def run_validate(args, logger: Logger):
    context = ProgramContext(threads=args.threads, working_dir=os.curdir, logger=logger)

    if args.format == "fasta":
        expected_uces = read_fasta(args.published_uces)
        actual_uces = read_fasta(args.deduce_uces)

        if args.minimise:
            print("Minimising expected UCEs...")
            expected_uces = minimise(expected_uces)

            print("Minimising actual UCEs...")
            actual_uces = minimise(actual_uces)

        print(f"Total UCEs expected: {len(expected_uces)}")
        print(f"Total UCEs found: {len(actual_uces)}")

        expected = set(u.id for u in expected_uces)
        expected_and_not_found = set(u.id for u in expected_uces)
        found_and_not_expected = set(u.id for u in actual_uces)

        i = 0
        total_actual = len(actual_uces)
        report_step = total_actual // 10

        for a in actual_uces:
            if i % report_step == 0:
                print(f"\tChecking {i}/{total_actual}...")

            for e in expected_uces:
                if seqs_match(a, e):
                    expected_and_not_found.discard(e.id)
                    found_and_not_expected.discard(a.id)

            i += 1

        print(f"Expected and found: {len(expected) - len(expected_and_not_found)}")
        print(f"Expected and not found: {len(expected_and_not_found)}")
        print(f"Not expected and found: {len(found_and_not_expected)}")

        if args.show_missing:
            print("\nMissing UCEs\n===========")
            expected_seqs = {u.id: u.forward for u in expected_uces}
            for uce in expected_and_not_found:
                print(f">{uce}")
                print(f"{expected_seqs[uce]}")

            print("\nUnexpected UCEs\n===========")
            for uce in found_and_not_expected:
                print(uce)

    else:
        intersection_file = intersect(
            args.published_uces, args.deduce_uces, [], context
        )

        with open(intersection_file, "r") as f:
            for feature in read_bed(f, ["name"]):
                print(feature)
