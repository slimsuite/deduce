import csv
import os
import sys
from typing import Counter

import numpy as np
import pysam
import scipy.stats

from deduce_uces.Genome import Genome
from deduce_uces.io.fasta import read_fasta_sequences
from deduce_uces.io.sam import count_mismatches
from deduce_uces.Logger import Logger
from deduce_uces.mappers.mapper import Mapper, MappingOutputType, MappingParameters
from deduce_uces.mappers.minimap_mapper import MinimapMapper
from deduce_uces.run import ProgramContext

STAT_COLS = ["mean", "median", "sd", "min", "max", "iqr"]


def get_stats(arr):
    if len(arr) == 0:
        return ["-"] * 6
    return [
        np.mean(arr),
        np.median(arr),
        np.std(arr),
        np.amin(arr),
        np.amax(arr),
        scipy.stats.iqr(arr),
    ]


def output_stats(stats):
    for s, col in zip(stats, STAT_COLS):
        print(f"{col}\t{s}")


def output_human(results):
    for res in results:
        print("\n\nSUMMARY\n=======")
        print(f"Aln\t{res[0]}")
        print(f"Hit\t{res[1][0]}")
        print(f"(exact)\t{res[1][1]}")
        print(f"Miss\t{res[1][2]}")

        print("\n\nMISMATCHES\n==========")
        output_stats(res[2])

        print("\n\nALIGNMENT SCORES\n================")
        output_stats(res[3])


def output_csv_header():
    csv_writer = csv.writer(sys.stdout)
    csv_writer.writerow(
        ["alignment", "hit", "exact_hit", "miss"]
        + ["mm_" + c for c in STAT_COLS]
        + ["score_" + c for c in STAT_COLS]
    )


def output_csv_row(row):
    csv_writer = csv.writer(sys.stdout)
    csv_writer.writerow([row[0]] + list(row[1]) + row[2] + row[3])


def run_measure_assembly_completeness(args, logger: Logger):
    context = ProgramContext(
        threads=args.threads,
        working_dir=os.curdir,
        logger=logger,
    )

    os.chdir(context.working_dir)

    results = []

    if args.output_format == "csv":
        output_csv_header()

    for assembly in args.assemblies:
        try:
            mapper = Mapper(
                Genome.from_filename(assembly),
                strategy=MinimapMapper(long_reads=True),
                mapping_params=MappingParameters(
                    mismatches_allowed=0,
                    secondary_mapping_limit=None,
                    sort=False,
                    mapping_stage="completeness",
                    filetype=MappingOutputType.SAM,
                ),
                context=context,
            )

            all_uces = set([s.id for s in read_fasta_sequences(args.uces)])

            context.logger.info(f"Indexing assembly {assembly}...")
            mapper.build_index()

            context.logger.info(f"Mapping UCEs to assembly...")
            sam_filename = mapper.map(os.path.abspath(args.uces))

            uce_mismatches = {}
            uce_scores = {}

            context.logger.info(f"Reading {sam_filename}...")
            alignments = pysam.AlignmentFile(sam_filename, "r")
            for alignment in alignments.fetch():
                if alignment.is_unmapped:
                    continue

                if alignment.query_sequence is None:
                    continue

                nm = count_mismatches(
                    alignment.cigartuples, len(alignment.query_sequence)
                )
                uce_mismatches[alignment.query_name] = nm
                uce_scores[alignment.query_name] = alignment.get_tag("AS")

            n_hits = len(uce_mismatches.keys())
            n_exact = len([v for v in uce_mismatches.values() if v == 0])
            n_miss = len(all_uces - set(uce_mismatches.keys()))

            row = (
                assembly,
                (n_hits, n_exact, n_miss),
                get_stats(np.array(list(uce_mismatches.values()))),
                get_stats(np.array(list(uce_scores.values()))),
            )

            if args.output_format == "csv":
                output_csv_row(row)

            results.append(row)

            if args.output_uces:
                with open(f"{assembly}.hits.csv", "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(
                        ["uce_id", "mismatches", "alignment_score", "is_exact"]
                    )
                    writer.writerows(
                        [
                            [
                                uce_id,
                                mismatches,
                                uce_scores[uce_id],
                                1 if mismatches == 0 else 0,
                            ]
                            for uce_id, mismatches in uce_mismatches.items()
                        ]
                    )

        except Exception as e:
            logger.warning(f"Failed to assess genome {assembly} ({e})!")

    if args.output_format == "human":
        output_human(results)
