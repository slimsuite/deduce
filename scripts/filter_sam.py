import argparse

from deduce_uces.io.sam import parse_sam_lazy, is_unmapped


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter a SAM to only include mapped sequences"
    )

    parser.add_argument("sam", metavar="SAM", type=str)

    return parser.parse_args()


def main():
    args = parse_args()

    mapped = set()
    partially_mapped = set()
    unmapped = set()

    for alignment in parse_sam_lazy(args.sam):
        if is_unmapped(alignment):
            unmapped.add(alignment.query_name)
        elif alignment.tags["NM"] > 0:
            partially_mapped.add(alignment.query_name)
        else:
            mapped.add(alignment.query_name)

    print(f"Mapped: {len(mapped)}")
    print(f"Partially mapped: {len(partially_mapped)}")
    print(f"Unmapped: {len(unmapped)}")
    print(list(unmapped)[:5])


if __name__ == "__main__":
    main()
