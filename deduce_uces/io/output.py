import os
from typing import Iterable, List, Dict, NamedTuple, Optional

from deduce_uces.io.bed import BedItem, write_bed
from deduce_uces.io.fasta import write_fasta
from deduce_uces.run import ProgramContext

UCEInstance = NamedTuple(
    "UCEInstance",
    [
        ("genome", str),
        ("contig", str),
        ("start", int),
        ("end", int),
        ("homology", Optional[float]),
    ],
)
UCE = NamedTuple(
    "UCE",
    [
        ("id", int),
        ("consensus_sequence", str),
        ("instances", List[UCEInstance]),
    ],
)


def output_uces_fasta(output_prefix, uces: Iterable[str], context: ProgramContext):
    fasta_filename = os.path.join(context.working_dir, output_prefix + ".deduce.fa")

    write_fasta(fasta_filename, [(i, uce) for i, uce in enumerate(uces)])


def output_uces_bed(output_prefix, uces: Iterable[UCE], context: ProgramContext):
    genome_items: Dict[str, List[BedItem]] = {}

    for uce in uces:
        for instance in uce.instances:
            genome_name, chromosome, start, end, homology = instance
            if genome_name not in genome_items:
                genome_items[genome_name] = []

            genome_items[genome_name].append(
                BedItem(
                    chromosome=chromosome,
                    start=start,
                    end=end,
                    tags={"name": uce.id, "score": str(100 * homology)},
                )
            )

    for genome, items in genome_items.items():
        with open(
            os.path.join(context.working_dir, f"{output_prefix}-{genome}.bed"), "w"
        ) as f:
            write_bed(f, items)
