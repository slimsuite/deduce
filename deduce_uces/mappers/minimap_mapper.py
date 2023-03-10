import os

from deduce_uces.io.checkpoint import checkpointed_file
from deduce_uces.mappers.mapper import (
    MappingOutputType,
    MappingStrategy,
    MappingParameters,
)
from deduce_uces.run import ProgramContext, run_subprocess

from deduce_uces.utils import with_default


class MinimapMapper(MappingStrategy):
    # Most of the time minimap is used for short read mapping to identify core kmer
    # positions.
    # But for mapping whole UCEs as part of assembly assessment, it should use
    # different settings - set `long_reads=True` to enable this.
    def __init__(self, long_reads: bool = False):
        self.long_reads = long_reads

    def build_index(
        self, reference_name: str, reference_file: str, context: ProgramContext
    ) -> str:
        with checkpointed_file(
            reference_name,
            "mmi",
            {},
            context,
        ) as checkpoint:
            if not checkpoint.checkpointed:
                minimap_options = {
                    "t": context.threads,
                    "ax": "sr",
                    "d": checkpoint.filename,
                }

                run_subprocess(
                    ["minimap2"],
                    [reference_file],
                    minimap_options,
                    context=context,
                    short_flags=True,
                )

        return checkpoint.filename

    # Run minimap to align candidates with the target genome
    def map(
        self,
        reads_file: str,
        reference_name: str,
        reference_file: str,
        index: str,
        params: MappingParameters,
        context: ProgramContext,
    ) -> str:
        with checkpointed_file(
            reference_name,
            "bam" if params.filetype == MappingOutputType.BAM else "sam",
            {
                "sm": with_default(params.secondary_mapping_limit, 0),
                "mismatches_allowed": params.mismatches_allowed,
                "sort": params.sort,
                "reads_file": reads_file,
                "reference_file": reference_file,
                "mapping_stage": params.mapping_stage,
                "long_read": self.long_reads,
            },
            context,
        ) as map_checkpoint:
            if not map_checkpoint.checkpointed:
                sam_output_filename = map_checkpoint.filename.replace(".bam", ".sam")

                minimap_options = {
                    "t": context.threads,
                    "a": index,
                    "ax": "sr",
                    "o": sam_output_filename,
                }

                if self.long_reads:
                    minimap_options["f"] = "50000,100000"

                minimap_args = [
                    "--MD",
                    "--eqx",  # Output =/X for match/mismatch in CIGAR
                    "--sam-hit-only",
                    reads_file,
                ]

                if params.secondary_mapping_limit is not None:
                    minimap_args.append("--secondary=yes")
                    minimap_options["N"] = params.secondary_mapping_limit

                run_subprocess(
                    ["minimap2"],
                    minimap_args,
                    minimap_options,
                    context=context,
                    short_flags=True,
                )

                if params.filetype == MappingOutputType.SAM:
                    return sam_output_filename

                # minimap has the "fun" behaviour of not outputting @SQ headers if the index is larger than 8GB
                # if these headers are missing, we need to convert from SAM to BAM with a reference supplied, before
                # doing the standard sort. samtools will repair the file and add the correct headers during the
                # conversion
                run_subprocess(
                    [
                        "samtools",
                        "view",
                    ],
                    [
                        "-b",
                        sam_output_filename,
                    ],
                    {"o": map_checkpoint.filename, "T": reference_file},
                    context=context,
                    short_flags=True,
                )

                os.remove(sam_output_filename)

                run_subprocess(
                    [
                        "samtools",
                        "sort",
                    ],
                    [
                        map_checkpoint.filename,
                    ],
                    {"o": map_checkpoint.filename},
                    context=context,
                    short_flags=True,
                )

            return map_checkpoint.filename
