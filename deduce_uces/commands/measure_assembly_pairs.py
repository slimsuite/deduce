import csv
import os
from random import random
import sys
from deduce_uces.Genome import Genome
from deduce_uces.Logger import Logger
from deduce_uces.commands.gen_synteny_map import (
    assess_synteny_by_pairs,
    read_unique_uces_from_sam,
)
from deduce_uces.io.fasta import read_fasta_sequences, write_fasta
from deduce_uces.mappers.mapper import Mapper, MappingOutputType, MappingParameters
from deduce_uces.mappers.minimap_mapper import MinimapMapper
from deduce_uces.run import ProgramContext
import networkx as nx


def run_measure_assembly_pairs(args, logger: Logger):
    context = ProgramContext(
        threads=args.threads,
        working_dir=os.curdir,
        logger=logger,
    )

    os.chdir(context.working_dir)

    context.logger.info(f"Reading synteny map...")
    map = nx.read_gexf(args.synteny_map)

    if args.output_format == "csv":
        out = csv.writer(sys.stdout)
        out.writerow(
            [
                "assembly",
                "total",
                "valid",
                "score",
                "average_consensus",
                "unique_uces",
                "nonunique_uces",
            ]
        )

    jobs = []
    for assembly in args.assemblies:
        try:
            mapper = Mapper(
                Genome.from_raw_filename(assembly),
                strategy=MinimapMapper(),
                mapping_params=MappingParameters(
                    mismatches_allowed=args.max_mapping_mismatches,
                    secondary_mapping_limit=0,
                    sort=False,
                    mapping_stage="synteny",
                    filetype=MappingOutputType.SAM,
                ),
                context=context,
            )

            context.logger.info(f"Indexing assembly {assembly}...")
            mapper.build_index()

            context.logger.info(f"Mapping UCEs to assembly...")
            sam_filename = mapper.map(os.path.abspath(args.uces))
            jobs.append((assembly, sam_filename))
        except Exception as e:
            logger.warning(f"Failed to map genome {assembly} ({e})!")

    out = csv.writer(sys.stdout)
    for assembly, sam in jobs:
        unique_uces, nonunique_uces = read_unique_uces_from_sam(
            sam, context.logger, max_mismatches=args.max_mapping_mismatches
        )

        (
            valid,
            score,
            avg_consensus,
            all_valid_pairs,
            total_pairs,
        ) = assess_synteny_by_pairs(
            unique_uces,
            map,
            # This is really the path length, so add 2 for the first and last nodes
            max_dist=2 + args.max_gap,
        )
        out.writerow(
            [
                assembly,
                total_pairs,
                valid,
                score,
                avg_consensus,
                len(unique_uces),
                len(nonunique_uces),
            ]
        )

        if args.output_pairs:
            with open(f"{assembly}.pairs.csv", "w") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        "id_1",
                        "ctg_1",
                        "pos_1",
                        "id_2",
                        "ctg_2",
                        "pos_2",
                        "score",
                        "length",
                        "min_consensus",
                        "avg_consensus",
                    ]
                )
                writer.writerows(
                    [
                        p[0][0].uce_id,
                        p[0][0].contig,
                        p[0][0].position,
                        p[0][1].uce_id,
                        p[0][1].contig,
                        p[0][1].position,
                        p[1],
                        p[2],
                        p[3],
                        p[4],
                    ]
                    for p in all_valid_pairs
                )
