from deduce_uces.commands.gen_syntenic_blocks import run_gen_syntenic_blocks
from deduce_uces.commands.gen_synteny_map import run_gen_synteny_map
from deduce_uces.commands.measure_assembly_blocks import run_measure_assembly_blocks
from deduce_uces.commands.measure_assembly_pairs import run_measure_assembly_pairs
from deduce_uces.commands.validate import run_validate
from deduce_uces.algs.merge_doa import GenomeInfo
from deduce_uces.commands.measure_assembly_contiguity import (
    run_measure_assembly_contiguity,
)
from deduce_uces.commands.measure_assembly_completeness import (
    run_measure_assembly_completeness,
)
import os
import re
import sys
from typing import Set

from deduce_uces.cli import make_argparser, parse_args
from deduce_uces.Genome import Genome
from deduce_uces.io.output import output_uces_bed, output_uces_fasta
from deduce_uces.Logger import ConsoleLogger, Logger, SplotActivity
from deduce_uces.mappers.available_mappers import MAPPING_STRATEGIES
from deduce_uces.mappers.mapper import Mapper, MappingOutputType, MappingParameters
from deduce_uces.run import ProgramContext
from deduce_uces.stages.candidates import find_candidate_uces
from deduce_uces.stages.core import (
    calculate_core_hash_parameters,
    find_core_kmers,
)
from deduce_uces.stages.extend import extend_uces
from deduce_uces.utils import get_canonical, get_files_in_directory


def identity(x):
    return x


def discover_genomes(context: ProgramContext) -> Set[str]:
    return set(
        re.sub(".fa$", "", file)
        for file in get_files_in_directory(context.working_dir)
        if file.endswith(".fa") and not file.endswith(".deduce.fa")
    )


def find_file_for_genome(genome: Genome, extension: str, context: ProgramContext):
    matching = [
        file
        for file in get_files_in_directory(context.working_dir)
        if file.endswith(extension) and file.startswith(genome.name)
    ]

    if len(matching) != 1:
        raise Exception(f"No unambiguous {extension} file for genome: {genome.name}")

    return matching[0]


def run_find(args, logger: Logger):
    context = ProgramContext(
        threads=args.threads,
        working_dir=os.path.abspath(args.dir),
        logger=logger,
    )

    os.chdir(context.working_dir)

    logger.splot_start(SplotActivity.UNKNOWN)

    genomes = set(Genome.from_name(name, context) for name in discover_genomes(context))

    if not args.skip_core_kmer_identification:
        core_kmers = find_core_kmers(genomes, args, context)
    else:
        core_kmers = [
            GenomeInfo(
                name=genome.name,
                bam=find_file_for_genome(genome, ".bam", context),
                fa=find_file_for_genome(genome, ".fa", context),
            )
            for genome in genomes
        ]

    uces = find_candidate_uces(core_kmers, args, context)

    unique_uces = set(get_canonical(u) for u in uces)

    if args.output_format == "fasta" or args.output_format == "all":
        context.logger.info(
            f"Output: save {len(unique_uces)} UCEs in {args.output_format} format"
        )
        logger.splot_start(SplotActivity.WRITING)
        output_uces_fasta(args.output, unique_uces, context)

    if args.output_format == "bed" or args.output_format == "all":
        logger.splot_start(SplotActivity.PROCESSING)
        extended_uces = extend_uces(unique_uces, genomes, args, context)

        if args.output_flanks != 0:
            # uces = find_flanks(uces, args, context)
            # TODO
            pass

        logger.splot_start(SplotActivity.WRITING)
        output_uces_bed(args.output, extended_uces, context)

    logger.splot_end()


def run_preprocess(args, logger: Logger):
    context = ProgramContext(
        threads=args.threads,
        working_dir=os.path.split(os.path.abspath(args.genome))[0],
        logger=logger,
    )

    os.chdir(context.working_dir)
    genome = Genome.from_filename(args.genome)

    context.logger.info("Hashing genome...")
    hash_params = calculate_core_hash_parameters({genome}, args, context)

    genome.hash(hash_params=hash_params, context=context)

    context.logger.info("Indexing genome...")

    mapping_strategy = MAPPING_STRATEGIES[args.mapper]()

    mapper = Mapper(
        genome,
        strategy=mapping_strategy,
        mapping_params=MappingParameters(
            mismatches_allowed=args.core_kmer_mismatches,
            secondary_mapping_limit=None,
            sort=True,
            mapping_stage="core",
            filetype=MappingOutputType.BAM,
        ),
        context=context,
    )
    mapper.build_index()


def main():
    parser = make_argparser()
    args = parse_args(parser)

    if hasattr(args, "debug") and hasattr(args, "no_color"):
        logger = ConsoleLogger(args.debug, args.no_color)
    else:
        logger = ConsoleLogger(False, False)

    if args.command == "find":
        run_find(args, logger)
    elif args.command == "validate":
        run_validate(args, logger)
    elif args.command == "preprocess":
        run_preprocess(args, logger)
    elif args.command == "gen_syntenic_blocks":
        run_gen_syntenic_blocks(args, logger)
    elif args.command == "gen_synteny_map":
        run_gen_synteny_map(args, logger)
    elif args.command == "measure_assembly_contiguity":
        run_measure_assembly_contiguity(args, logger)
    elif args.command == "measure_assembly_completeness":
        run_measure_assembly_completeness(args, logger)
    elif args.command == "measure_assembly_pairs":
        run_measure_assembly_pairs(args, logger)
    elif args.command == "measure_assembly_blocks":
        run_measure_assembly_blocks(args, logger)
    else:
        parser.print_help(sys.stderr)


if __name__ == "__main__":
    main()
