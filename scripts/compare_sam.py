import argparse

from deduce_uces.io.sam import parse_sam_lazy


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare two SAM files to find discrepancies"
    )

    parser.add_argument("a", metavar="A", type=str)
    parser.add_argument("b", metavar="B", type=str)

    return parser.parse_args()


def main():
    args = parse_args()

    a_mapped_alignments = [
        alignment
        for alignment in parse_sam_lazy(args.a)
        if "NM" in alignment.tags and alignment.tags["NM"] == 0
    ]
    b_mapped_alignments = [
        alignment
        for alignment in parse_sam_lazy(args.b)
        if "NM" in alignment.tags and alignment.tags["NM"] == 0
    ]

    a_mapped_names = set(alignment.query_name for alignment in a_mapped_alignments)
    b_mapped_names = set(alignment.query_name for alignment in b_mapped_alignments)

    shared_mapped_names = a_mapped_names.intersection(b_mapped_names)

    print(len(shared_mapped_names))


if __name__ == "__main__":
    main()
