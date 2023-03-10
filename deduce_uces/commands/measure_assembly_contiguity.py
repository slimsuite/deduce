from collections import defaultdict
import csv
from functools import partial
import itertools
import multiprocessing
import os
import random
import sys
from typing import Any, Dict, List, Tuple


from deduce_uces.Genome import Genome
from deduce_uces.io.fasta import read_fasta_sequences, write_fasta
from deduce_uces.io.sam import parse_sam_lazy
from deduce_uces.Logger import ConsoleLogger, Logger
from deduce_uces.mappers.mapper import Mapper, MappingOutputType, MappingParameters
from deduce_uces.mappers.minimap_mapper import MinimapMapper
from deduce_uces.run import ProgramContext


def is_block_valid(
    block: List[str], uce_positions: Dict[str, List[Tuple[str, int]]]
) -> bool:
    assert len(block) >= 2

    if not all(u in uce_positions for u in block):
        return False

    positions_run = [uce_positions[uce_id] for uce_id in block]

    for path in itertools.product(*positions_run):
        if not len(set(p[0] for p in path)) == 1:
            continue

        last_position = path[0]

        valid = True
        for uce in path[1:]:
            if uce[0] != last_position[0] or uce[1] <= last_position[1]:
                valid = False
                break

            last_position = uce
        if valid:
            return True

    return False


def process_job(
    job: Tuple[str, str],
    blocks: List[List[str]],
    pickled_logger: Any,
    max_mapping_mismatches: int,
    working_dir: str,
):
    os.chdir(working_dir)
    assembly, sam_filename = job
    logger = ConsoleLogger.from_pickleable(pickled_logger)
    logger.debug(f"Processing alignments for assembly: {assembly}")

    uce_positions = defaultdict(list)
    for alignment in parse_sam_lazy(
        sam_filename, max_mapping_mismatches, logger, False
    ):
        uce_positions[alignment.query_name].append(
            (
                alignment.reference_name,
                alignment.position,
            )
        )

    logger.debug(f"Read alignments for assembly: {assembly}")

    valid_blocks = len([b for b in blocks if is_block_valid(b, uce_positions)])
    logger.debug(f"Found valid blocks for assembly: {assembly}")
    invalid_blocks = len(blocks) - valid_blocks

    result = [os.path.basename(assembly), valid_blocks, invalid_blocks]

    out = csv.writer(sys.stdout)
    out.writerow(result)
    return os.path.basename(assembly), valid_blocks, invalid_blocks


def run_measure_assembly_contiguity(args, logger: Logger):
    context = ProgramContext(
        threads=args.threads,
        working_dir=os.curdir,
        logger=logger,
    )

    os.chdir(context.working_dir)

    context.logger.info(f"Reading UCEs...")
    sequences = {s.id: str(s.seq) for s in read_fasta_sequences(args.uces)}

    context.logger.info(f"Reading blocks...")
    with open(args.blocks) as f:
        blocks = [line.strip().split(",") for line in f]

    relevant_uce_ids = set(u for b in blocks for u in b)

    block_seqs = [(id, seq) for id, seq in sequences.items() if id in relevant_uce_ids]

    context.logger.info(f"Preparing query...")
    temp_fasta_filename = f"contiguity_query_{random.randint(0,10000)}.deduce.fa"
    write_fasta(temp_fasta_filename, block_seqs)

    if args.output_format == "csv":
        out = csv.writer(sys.stdout)
        out.writerow(["assembly", "valid", "invalid"])

    jobs = []
    for assembly in args.assemblies:
        try:
            mapper = Mapper(
                Genome.from_raw_filename(assembly),
                strategy=MinimapMapper(),
                mapping_params=MappingParameters(
                    mismatches_allowed=0,
                    secondary_mapping_limit=args.max_secondary_alignments,
                    sort=False,
                    mapping_stage="completeness",
                    filetype=MappingOutputType.SAM,
                ),
                context=context,
            )

            context.logger.info(f"Indexing assembly {assembly}...")
            mapper.build_index()

            context.logger.info(f"Mapping UCEs to assembly...")
            sam_filename = mapper.map(os.path.abspath(temp_fasta_filename))
            jobs.append((assembly, sam_filename))
        except Exception as e:
            logger.warning(f"Failed to map genome {assembly} ({e})!")

    process = partial(
        process_job,
        pickled_logger=context.logger.to_pickleable(),
        max_mapping_mismatches=args.max_mapping_mismatches,
        blocks=blocks,
        working_dir=os.curdir,
    )

    with multiprocessing.Pool(args.threads) as p:
        for result in p.imap(process, jobs):
            if args.output_format != "csv":
                print("SYNTENIC BLOCKS\n==============")
                print(f"Aln\t{result[0]}")
                print(f"Valid\t{result[1]}")
                print(f"Invalid\t{result[2]}")

    os.remove(temp_fasta_filename)
