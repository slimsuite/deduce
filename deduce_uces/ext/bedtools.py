from typing import List

from deduce_uces.io.checkpoint import checkpointed_file
from deduce_uces.run import run_subprocess, ProgramContext


def get_fasta(
    genome_name: str,
    bed_filename: str,
    reference_filename: str,
    context: ProgramContext,
) -> str:
    with checkpointed_file(
        f"{genome_name}_extension",
        "deduce.fa",
        {"in": bed_filename, "ref": reference_filename},
        context,
    ) as checkpoint:
        if not checkpoint.checkpointed:
            with open(checkpoint.filename, "wb") as f:
                run_subprocess(
                    ["bedtools", "getfasta"],
                    ["-name"],  # Use `name` field rather than coordinates
                    {
                        "fi": reference_filename,
                        "bed": bed_filename,
                    },
                    context,
                    short_flags=True,
                    stdout_file=f,
                )
        return checkpoint.filename


def intersect(
    bed_a: str,
    bed_b: str,
    flags: List[str],
    context: ProgramContext,
) -> str:
    with checkpointed_file(
        f"intersection",
        "bed",
        {"a": bed_a, "b": bed_b, "flags": flags},
        context,
    ) as checkpoint:
        if not checkpoint.checkpointed:
            with open(checkpoint.filename, "wb") as f:
                run_subprocess(
                    ["bedtools", "intersect"],
                    flags,
                    {
                        "a": bed_a,
                        "b": bed_b,
                    },
                    context,
                    short_flags=True,
                    stdout_file=f,
                )

        return checkpoint.filename


def sam_to_bed(sam_filename, context: ProgramContext) -> str:
    bed_filename = sam_filename.replace(".sam", ".bed")
    with open(bed_filename, "wb") as f:
        run_subprocess(
            ["bedtools", "bamtobed"],
            [],
            {
                "i": sam_filename,
            },
            context,
            short_flags=True,
            stdout_file=f,
        )

    return bed_filename
