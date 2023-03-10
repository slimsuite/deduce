import os
from deduce_uces.commands.gen_synteny_map import sort_and_group
from deduce_uces.io.output import UCE, UCEInstance
import math
from typing import Set, List, Dict


from deduce_uces.Genome import Genome
from deduce_uces.io.checkpoint import checkpointed_file
from deduce_uces.io.fasta import write_fasta
from deduce_uces.io.sam import index_bam, parse_sam_lazy
from deduce_uces.mappers.available_mappers import MAPPING_STRATEGIES
from deduce_uces.mappers.mapper import Mapper, MappingOutputType, MappingParameters
from deduce_uces.run import ProgramContext


def should_include_uce(xs: List[UCEInstance], uce_max_occurrences: int) -> int:
    if len(xs) == 0:
        return False

    return (
        max(len(instances) for _, instances in sort_and_group(xs, lambda x: x.genome))
        <= uce_max_occurrences
    )


def extend_uces(
    uces: Set[str],
    genomes: Set[Genome],
    args,
    context: ProgramContext,
) -> List[UCE]:
    context.logger.info("Extension: searching for UCE appearances...")

    indexed_uces = {i: uce for i, uce in enumerate(uces)}

    with checkpointed_file(f"extension_query", "fa", {}, context) as checkpoint:
        write_fasta(checkpoint.filename, indexed_uces.items())
        query_filename = checkpoint.filename

    # Map the reference sequences to the genomes
    uce_extension_min_homology = (
        args.uce_extension_min_homology or args.uce_core_min_homology
    )

    n_allowable_mismatches = math.floor(
        args.uce_min_length * (1 - uce_extension_min_homology)
    )

    uces_by_id: Dict[int, List[UCEInstance]] = {id: [] for id in indexed_uces.keys()}

    mapping_strategy = MAPPING_STRATEGIES[args.mapper]()

    try:
        for genome in genomes:
            mapper = Mapper(
                genome,
                strategy=mapping_strategy,
                mapping_params=MappingParameters(
                    mismatches_allowed=n_allowable_mismatches,
                    secondary_mapping_limit=args.uce_max_occurrences
                    + 1,  # Used to detect if the limit is hit
                    sort=True,
                    mapping_stage="extension",
                    filetype=MappingOutputType.BAM,
                ),
                context=context,
            )

            mapper.build_index()
            bam_file = mapper.map(query_filename)

            index_bam(bam_file)

            for match in parse_sam_lazy(
                bam_file,
                n_allowable_mismatches,
                context.logger,
                args.mapper == "bowtie",
            ):
                match_length = len(match.query_sequence) - match.number_of_mismatches
                query_length = len(match.query_sequence)

                homology = float(match_length) / float(query_length)

                if homology < args.uce_core_min_homology:
                    continue

                uces_by_id[int(match.query_name)].append(
                    UCEInstance(
                        genome.name,
                        match.reference_name,
                        match.position,
                        match.position + len(match.query_sequence),
                        homology,
                    )
                )
    finally:
        os.remove(query_filename)

    return [
        UCE(id, uce, uces_by_id[id])
        for id, uce in indexed_uces.items()
        if should_include_uce(uces_by_id[id], args.uce_max_occurrences)
    ]
