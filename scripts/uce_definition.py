import argparse
import itertools
from typing import List, Mapping

import seaborn as sns

# Apply the default theme
sns.set_theme()

from deduce_uces.io.bed import read_bed, BedItem

AllUCEs = Mapping[str, List[BedItem]]
ParameterSpace = List[Mapping[str, int]]

ALL_PARAMS = {"min_length", "min_appearances", "min_identity"}
PARAM_MINS = {"min_length": 50, "min_appearances": 1, "min_identity": 70}


def generate_parameter_space(all_uces: AllUCEs) -> ParameterSpace:
    max_length = max([a.end - a.start for u in all_uces.values() for a in u])

    min_lengths = range(50, max_length, 10)
    min_appearances = range(1, 2)
    min_identity = range(70, 100)

    return [
        {"min_length": l, "min_appearances": a, "min_identity": i}
        for l, a, i in itertools.product(min_lengths, min_appearances, min_identity)
    ]


def count_uces_by_parameter_set(all_uces: AllUCEs, p):
    matching = 0
    for appearances in all_uces.values():
        # Apply length constraint
        # Apply identity constraint
        # Apply number of appearances
        if (
            len(
                [
                    a
                    for a in appearances
                    if a.end - a.start >= p["min_length"]
                    and a.tags["score"] >= p["min_identity"]
                ]
            )
            > p["min_appearances"]
        ):
            matching += 1

    return {**p, "uces": matching}


def plot_single_parameter(all_uces: AllUCEs, pspace: ParameterSpace, param: str):
    constants = ALL_PARAMS - {param}

    filtered = []
    for p in pspace:
        if all(p[c] == PARAM_MINS[c] for c in constants):
            filtered.append(p)

    data = [count_uces_by_parameter_set(all_uces, p) for p in filtered]

    sns.relplot(
        data=data,
        x="min_length",
        y="uces",
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Take in a big ol buncha BED files and make some graphs"
    )

    parser.add_argument(
        "alignment",
        metavar="ALIGNMENT",
        nargs="+",
        type=str,
    )

    return parser.parse_args()


def main():
    args = parse_args()

    all_uces = {}
    for a in args.alignment:
        with open(a) as f:
            for item in read_bed(f):
                if item.tags["name"] not in all_uces:
                    all_uces[item.tags["name"]] = []

                all_uces[item.tags["name"]].append(item)

    parameter_space = generate_parameter_space(all_uces)

    plot_single_parameter(all_uces, parameter_space, "min_length")


if __name__ == "__main__":
    main()
