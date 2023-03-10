import math
import os
from deduce_uces.Logger import Logger
from deduce_uces.algs import find_syntenic_blocks
from deduce_uces.run import ProgramContext
import networkx as nx


def run_gen_syntenic_blocks(args, logger: Logger):
    context = ProgramContext(
        threads=args.threads,
        working_dir=os.curdir,
        logger=logger,
    )

    os.chdir(context.working_dir)

    context.logger.info("Loading UCEs from BED files...")
    bed_files = [os.path.abspath(f) for f in args.beds]
    anchors, genome_names = find_syntenic_blocks.load_anchors(bed_files)

    uces = find_syntenic_blocks.group_anchors_by_uce(anchors, genome_names)

    context.logger.info("Identifying syntenic blocks...")

    min_genomes = math.ceil(args.syntenic_block_consensus_support * len(bed_files))

    blocks, synteny_map = find_syntenic_blocks.find_syntenic_blocks(
        uces,
        genome_names,
        min_genomes=min_genomes,
        min_length=args.syntenic_block_min_length,
        include_distance=args.include_distance,
        recursion_limit=args.max_recursion_depth,
    )

    nx.write_gexf(synteny_map, "synteny.gexf")

    if args.minimise_blocks:
        context.logger.info("Minimising syntenic blocks...")
        blocks = find_syntenic_blocks.minimise_blocks(blocks)

    for block in blocks:
        print(",".join(block))
