from deduce_uces.algs.merge_doa import GenomeInfo
from deduce_uces.Logger import SplotActivity
from math import ceil
from typing import NamedTuple, Dict, List, Set, Collection

from deduce_uces.Genome import Genome
from deduce_uces.cli import HashParams, MergeParams, MergeMethod
from deduce_uces.ext.jellyfish import (
    compute_counter_size,
    choose_hash_table_size,
    merge_hashes,
    dump_hash,
)
from deduce_uces.mappers.mapper import Mapper, MappingOutputType, MappingParameters
from deduce_uces.run import ProgramContext
from deduce_uces.types import Strand
from deduce_uces.mappers.available_mappers import MAPPING_STRATEGIES

SingleGenomeCoreKmer = NamedTuple(
    "SingleGenomeCoreKmer",
    [
        ("position", int),
        ("chromosome", str),
        ("sequence", str),
        ("strand", Strand),
        ("mismatches", int),
    ],
)

CoreKmer = NamedTuple(
    "CoreKmer",
    [("name", str), ("instances", Dict[str, SingleGenomeCoreKmer])],
)

GenomeWithCoreKmers = NamedTuple(
    "GenomeWithCoreKmers", [("genome", Genome), ("core_kmers_bam_file", str)]
)


def find_core_kmers(
    genomes: Set[Genome], args, context: ProgramContext
) -> List[GenomeInfo]:

    context.logger.info("Core: starting search for core kmers...")

    # Hash each genome, or load existing hash

    context.logger.info("Core: hashing genomes...")
    for genome in genomes:
        context.logger.debug(f"Core: \t{genome.name}")

        # Recalculate hash parameters for each genome - in case the available virtual memory has changed
        hash_params = calculate_core_hash_parameters(genomes, args, context)

        context.logger.splot_start(
            SplotActivity.PROCESSING, fake_threads=(context.threads, "jellyfish")
        )

        genome.hash(hash_params=hash_params, context=context)

        context.logger.splot_end(fake_threads=(context.threads, "jellyfish"))

    context.logger.info("Core: reducing hash counts...")
    context.logger.splot_start(SplotActivity.PROCESSING)
    for genome in genomes:
        # TODO run in parallel
        genome.reduce_hash(context)
    context.logger.splot_end()

    # Intersect hashes to find core kmers
    context.logger.info("Core: intersecting hashes...")

    core_genomes_required = ceil(args.core_kmer_threshold * len(genomes))

    merge_params = MergeParams(
        min_occurrences=core_genomes_required, method=MergeMethod.SUM
    )

    candidates_hash = merge_hashes(
        "candidates",
        [genome.hash_file for genome in genomes],
        hash_params,
        merge_params,
        context,
    )

    context.logger.info("Core: dumping candidates...")
    candidate_sequence_file = dump_hash(candidates_hash, context)

    context.logger.info("Core: aligning candidates...")

    mapping_strategy = MAPPING_STRATEGIES[args.mapper]()

    instances: List[GenomeInfo] = []
    for genome in genomes:
        mapper = Mapper(
            genome,
            strategy=mapping_strategy,
            mapping_params=MappingParameters(
                mismatches_allowed=args.core_kmer_mismatches,
                secondary_mapping_limit=args.core_kmer_mapping_limit,
                sort=True,
                mapping_stage="core",
                filetype=MappingOutputType.BAM,
            ),
            context=context,
        )

        context.logger.debug(f"Core: \tindexing {genome.name}")
        mapper.build_index()

        context.logger.debug(f"Core: \tmapping {genome.name}")
        bam_file = mapper.map(candidate_sequence_file)

        instances.append(
            GenomeInfo(name=genome.name, bam=bam_file, fa=genome.source_files[0])
        )

    context.logger.info("Core: finished identifying core kmers")

    return instances


def calculate_core_hash_parameters(
    genomes: Set[Genome], args, context: ProgramContext
) -> HashParams:
    counter_len, out_counter_len = compute_counter_size(args.core_kmer_mapping_limit)

    if args.jf_hash_size is not None:
        hash_size = args.jf_hash_size
    else:
        hash_size = max(
            choose_hash_table_size(
                genome.estimate_genome_size(), args.core_kmer_size, context, counter_len
            )
            for genome in genomes
        )

    return HashParams(
        k=args.core_kmer_size,
        min_occurrences=1,
        max_occurrences=None,
        size=hash_size,
        counter_len=counter_len,
        out_counter_len=out_counter_len,
        disk_mode=args.jf_disk,
    )
