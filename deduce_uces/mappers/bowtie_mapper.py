import os
from typing import Union, Dict

from deduce_uces.io.checkpoint import checkpointed_file
from deduce_uces.mappers.mapper import MappingStrategy, MappingParameters
from deduce_uces.run import ProgramContext, run_subprocess

# Check if a named bowtie index exists in the working directory
# Indices are multiple files; this just checks the .bt2 file
from deduce_uces.utils import with_default


def does_bowtie_index_exist(name: str, context: ProgramContext) -> bool:
    return os.path.exists(os.path.join(context.working_dir, f"{name}.1.bt2"))


class BowtieMapper(MappingStrategy):
    def build_index(
        self, reference_name: str, reference_file: str, context: ProgramContext
    ) -> str:
        if not does_bowtie_index_exist(reference_name, context=context):

            bowtie_options = {
                "threads": context.threads,
            }

            run_subprocess(
                ["bowtie2-build"],
                [
                    reference_file,
                    reference_name,
                ],
                bowtie_options,
                context=context,
            )

        return reference_name

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
            "bam",
            {
                "sm": with_default(params.secondary_mapping_limit, 0),
                "mismatches_allowed": params.mismatches_allowed,
                "sort": params.sort,
                "reads_file": reads_file,
                "reference_file": reference_file,
                "mapping_stage": params.mapping_stage,
            },
            context,
        ) as map_checkpoint:
            if not map_checkpoint.checkpointed:
                sam_output_filename = map_checkpoint.filename.replace(".bam", ".sam")

                bowtie_options: Dict[str, Union[int, str, None]] = {
                    "x": index,
                    "S": sam_output_filename,
                    "U": reads_file,
                }

                bowtie_args = [
                    # very-fast preset works for end-to-end mapping of 50-mers
                    # (might need to be tested with different values of k)
                    "--very-fast",
                    "--end-to-end",
                    "--quiet",
                    "--no-unal",
                    "--threads",
                    str(context.threads),
                    "-f",  # FASTA input
                ]

                if params.secondary_mapping_limit is not None:
                    bowtie_options["k"] = params.secondary_mapping_limit

                run_subprocess(
                    ["bowtie2"],
                    bowtie_args,
                    bowtie_options,
                    context=context,
                    short_flags=True,
                )

                if params.sort:
                    # Sort in-place and convert to BAM
                    run_subprocess(
                        [
                            "samtools",
                            "sort",
                        ],
                        [
                            sam_output_filename,
                        ],
                        {"o": map_checkpoint.filename},
                        context=context,
                        short_flags=True,
                    )

                else:
                    # Convert to BAM
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
        return map_checkpoint.filename
