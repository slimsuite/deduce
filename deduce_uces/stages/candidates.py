import os
from functools import partial
from multiprocessing import Pool
from random import randint
from typing import Any, Collection, Dict, List, Tuple

import pysam

from deduce_uces.algs.indexed_seq import IndexedSeq
from deduce_uces.algs.merge_doa import (
    GenomeInfo,
    core_kmers_to_uces,
    get_candidate_regions,
)
from deduce_uces.algs.name_mapping import NameMapping, create_name_mapping
from deduce_uces.algs.uce_params import UCEParams, calculate_uce_params
from deduce_uces.ext.samtools import index_fasta
from deduce_uces.io.checkpoint import checkpointed_file
from deduce_uces.io.sam import count_mismatches, index_bam
from deduce_uces.Logger import ConsoleLogger, SplotActivity
from deduce_uces.run import ProgramContext
from deduce_uces.stages.trie import CoreKmer, create_trie, save_trie
from deduce_uces.utils import flatten


def read_bam_without_seqs(
    genome: GenomeInfo, uce_params: UCEParams, name_mapping: NameMapping
):
    f = pysam.AlignmentFile(genome.bam, "rb")  # type: ignore

    for alignment in f.fetch():
        if alignment.is_unmapped or alignment.is_supplementary:
            continue

        if (
            uce_params.core_kmer_mismatches == 0
            and alignment.cigarstring != f"{uce_params.core_kmer_len}="
            and alignment.cigarstring != f"{uce_params.core_kmer_len}M"
        ):
            continue
        elif (
            # TODO use MD count for bowtie
            count_mismatches(alignment.cigartuples, uce_params.core_kmer_len)
            > uce_params.core_kmer_mismatches
        ):
            continue

        yield CoreKmer(
            id=alignment.query_name,
            genome=name_mapping.sequence_to_uchar[genome.name],
            contig=name_mapping.contig_to_uchar[genome.name][alignment.reference_name],
            position=alignment.reference_start,
        )


def build_core_kmer_trie(
    job: Tuple[GenomeInfo, str],
    uce_params: UCEParams,
    name_mapping: NameMapping,
    pickled_logger: Any,
    working_dir: str,
):
    os.chdir(working_dir)
    logger = ConsoleLogger.from_pickleable(pickled_logger)

    logger.debug(f"Candidates: creating core kmer trie for {job[0].name}")

    logger.splot_start(SplotActivity.READING)

    core_kmers_in_regions = get_candidate_regions(
        read_bam_without_seqs(job[0], uce_params, name_mapping), uce_params
    )

    t = create_trie(flatten(core_kmers_in_regions))

    logger.splot_start(SplotActivity.WRITING)
    save_trie(t, job[1])
    logger.splot_end()
    logger.debug(f"Candidates: created core kmer trie for {job[0].name} ({job[1]})")


def read_core_kmers_in_regions(
    genomes: Collection[GenomeInfo],
    uce_params: UCEParams,
    name_mapping: NameMapping,
    context: ProgramContext,
) -> Dict[int, str]:
    core_kmers = {}

    core_kmer_jobs = []
    for genome in genomes:

        with checkpointed_file(
            genome.name,
            "trie",
            {
                "in": genome.bam,
                "uce_params": uce_params,
                "name_mapping": name_mapping,
                "random": randint(1, 100000000),
            },
            context,
        ) as checkpoint:
            core_kmer_jobs.append((genome, checkpoint.filename))

            core_kmers[
                name_mapping.sequence_to_uchar[genome.name]
            ] = checkpoint.filename

    process = partial(
        build_core_kmer_trie,
        uce_params=uce_params,
        name_mapping=name_mapping,
        pickled_logger=context.logger.to_pickleable(),
        working_dir=context.working_dir,
    )

    with Pool(context.threads) as p:
        p.map(process, core_kmer_jobs)

    return core_kmers


def find_candidate_uces(
    genomes: List[GenomeInfo],
    args,
    context: ProgramContext,
) -> List[str]:
    context.logger.info("Candidates: identifying candidate UCEs from core kmers...")

    uce_params = calculate_uce_params(len(genomes), args)

    genomes = sorted(
        genomes,
        key=lambda x: x.name == args.reference,
        reverse=True,
    )

    # Index BAM files if needed
    # This will allow accessing the files per-contig, so that the merging
    # can be done in multiple processes
    context.logger.debug("Merge: indexing BAM files")
    with Pool(min(context.threads, len(genomes))) as p:
        p.map(index_bam, [genome.bam for genome in genomes])

    context.logger.debug("Merge: indexing FASTA files")
    with Pool(min(context.threads, len(genomes))) as p:
        p.map(index_fasta, [genome.fa for genome in genomes])

    # Create name mapping, used to reduce memory required by core kmers
    name_mapping = create_name_mapping(
        {g.name: IndexedSeq(g.fa).get_contigs() for g in genomes}
    )

    # h = hpy()
    # h.setrelheap()
    core_kmers = read_core_kmers_in_regions(
        genomes,
        uce_params,
        name_mapping,
        context,
    )
    # x = h.heap()
    # breakpoint()
    reference_genomes = (
        genomes
        if args.reduced_support_exhaustive
        else genomes[: len(genomes) - uce_params.uce_min_support + 1]
    )

    uces = []
    try:
        for ref in reference_genomes:
            context.logger.debug(f"Candidates: finding UCEs for: {ref.name}")
            uces.extend(
                core_kmers_to_uces(
                    core_kmers,
                    ref,
                    [g for g in genomes if g.name != ref.name],
                    name_mapping,
                    uce_params,
                    context,
                )
            )

    finally:
        # Remove the temporary core kmer tries, even if an exception occurred
        for trie_filename in core_kmers.values():
            os.remove(trie_filename)

    return uces
