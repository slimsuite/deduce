import csv
from typing import NamedTuple, Union, Mapping, Collection, Iterator, List

from deduce_uces.run import run_subprocess, ProgramContext

BedItem = NamedTuple(
    "BedItem",
    [
        ("chromosome", str),
        ("start", int),
        ("end", int),
        ("tags", Mapping[str, Union[int, str]]),
    ],
)


def write_bed(f, items: Collection[BedItem]):
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(
        [
            [
                item.chromosome,
                item.start,
                item.end,
                item.tags["name"] if "name" in item.tags else "",
                item.tags["score"] if "score" in item.tags else "",
            ]
            + [v for k, v in item.tags.items() if k not in ["name", "score"]]
            for item in items
        ]
    )


def read_bed(f, tags: List[str] = None) -> Iterator[BedItem]:
    if tags is None:
        tags = ["name", "score"]

    reader = csv.reader(f, delimiter="\t")

    for row in reader:
        yield BedItem(
            chromosome=row[0],
            start=int(row[1]),
            end=int(row[2]),
            tags={k: v for k, v in zip(tags, row[3:])},
        )
