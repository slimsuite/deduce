import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        description="Fix a BED file which is using 1-indexing"
    )

    parser.add_argument(
        "bed",
        metavar="BED",
        type=str,
    )

    args = parser.parse_args()

    with open(args.bed) as f:
        reader = csv.reader(f, delimiter="\t")
        lines = [[l[0], int(l[1]) - 1, int(l[2]) - 1, l[3], l[4]] for l in reader]

    with open(args.bed + ".fixed", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(lines)


if __name__ == "__main__":
    main()
