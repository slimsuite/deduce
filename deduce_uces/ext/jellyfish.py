import math
import os
from math import log2, ceil
from typing import List, Tuple, Mapping, Union, Dict

from psutil import virtual_memory

from deduce_uces.cli import HashParams, MergeParams, MergeMethod
from deduce_uces.filenames import get_hash_filename, get_merged_hash_filename
from deduce_uces.io.checkpoint import checkpointed_file
from deduce_uces.io.fasta import rename_fasta_sequences
from deduce_uces.run import ProgramContext, run_subprocess

# When calculating available virtual memory, leave a bit of room for overhead so that
# Jellyfish doesn't run out of memory
VIRTUAL_MEM_SAFETY_FACTOR = 0.95


# Calculate the Jellyfish counter size to use when counting kmers
# This should be as small as possible while still fitting the maximum count we're interested in, to save
# disk space.
def compute_counter_size(max_occurrences: int) -> Tuple[int, int]:
    # counter_len controls the in-memory counter size in bits.
    counter_len = ceil(log2(max_occurrences + 1)) + 1

    # out_counter_len controls the disk counter size in bytes. We only need one byte to represent a count of 0, 1,
    # or more.
    out_counter_len = ceil(counter_len / 8)

    return counter_len, out_counter_len


def hash_fasta(
    name: str, fasta: str, hash_params: HashParams, context: ProgramContext
) -> str:
    output_filename = get_hash_filename(name, hash_params, context)

    jellyfish_options: Dict[str, Union[int, str, None]] = {
        "output": output_filename,
        "threads": context.threads,
        "mer-len": hash_params.k,
        "size": hash_params.size,
        "counter-len": hash_params.counter_len,
        "out-counter-len": hash_params.out_counter_len,
        "lower-count": hash_params.min_occurrences,
        "upper-count": hash_params.max_occurrences,
    }

    if hash_params.max_occurrences is not None:
        jellyfish_options["upper-count"] = hash_params.max_occurrences

    jellyfish_args: List[str] = (
        ["--canonical"] + ["--disk"]
        if hash_params.disk_mode
        else ["--canonical"] + [os.path.abspath(fasta)]
    )

    run_subprocess(
        ["jellyfish", "count"],
        jellyfish_args,
        jellyfish_options,
        context,
    )

    return output_filename


def merge_hashes(
    name: str,
    hashes: List[str],
    hash_params: HashParams,
    merge_params: MergeParams,
    context: ProgramContext,
) -> str:
    output_filename = get_merged_hash_filename(name, hash_params, merge_params, context)

    if os.path.exists(output_filename):
        context.logger.info(f"Found existing candidates file: {output_filename}")
        return output_filename

    jellyfish_options: Mapping[str, Union[int, str]] = {
        "output": output_filename,
        "lower-count": merge_params.min_occurrences,
    }

    jellyfish_args = ["jellyfish", "merge"]

    if merge_params.method == MergeMethod.MIN:
        jellyfish_args.append("--min")

    run_subprocess(
        ["jellyfish", "merge"],
        [os.path.abspath(hash) for hash in hashes],
        jellyfish_options,
        context,
    )

    return output_filename


def hash_fastas(
    name: str, fastas: List[str], hash_params: HashParams, context: ProgramContext
) -> str:
    if len(fastas) == 1:
        hash = hash_fasta(name, fastas[0], hash_params, context)
        return hash

    intermediate_hashes = []
    for i, fasta in enumerate(fastas):
        intermediate_hashes.append(
            hash_fasta(name + str(i), fasta, hash_params, context)
        )

    merged_hash = merge_hashes(
        name,
        intermediate_hashes,
        hash_params,
        MergeParams(
            min_occurrences=hash_params.min_occurrences, method=MergeMethod.SUM
        ),
        context,
    )

    for intermediate_hash in intermediate_hashes:
        os.remove(intermediate_hash)

    return merged_hash


def dump_hash(hash_filename: str, context: ProgramContext) -> str:
    with checkpointed_file(
        hash_filename + "_unlabelled",
        "deduce.fa",
        depends_on={"hash_filename": hash_filename},
        context=context,
    ) as unlabelled_checkpoint:

        with checkpointed_file(
            hash_filename,
            "deduce.fa",
            depends_on={"hash_filename": hash_filename},
            context=context,
        ) as labelled_checkpoint:
            if not unlabelled_checkpoint.checkpointed:
                jellyfish_options = {
                    "output": unlabelled_checkpoint.filename,
                }

                run_subprocess(
                    ["jellyfish", "dump"],
                    [hash_filename],
                    jellyfish_options,
                    context,
                )

            if not labelled_checkpoint.checkpointed:
                rename_fasta_sequences(
                    unlabelled_checkpoint.filename,
                    labelled_checkpoint.filename,
                    context,
                )

    return labelled_checkpoint.filename


def reduce_hash_to_binary_counts(filename: str, context: ProgramContext) -> str:
    with checkpointed_file(
        filename + "_reduced",
        "jf",
        depends_on={"filename": filename},
        context=context,
    ) as checkpoint:
        if not checkpoint.checkpointed:
            run_subprocess(
                ["jf_bit_counts"],
                [filename, checkpoint.filename],
                {},
                context,
            )

        os.remove(filename)

        return checkpoint.filename


# Estimate the size in bytes of a Jellyfish hash table for the given genome size
# This is the optimal size of the hash table, because if it all fits in memory there are no
# intermediate writes to disk
# Uses the formula from the Jellyfish man page
def get_optimal_hash_table_size(
    genome_size: int, kmer_size: int, max_reprobe: int
) -> int:
    # Hash table is rounded up to the next power of 2 above genome_size
    # https://stackoverflow.com/a/466242
    l = math.floor(math.log2(genome_size)) + 1
    r = math.floor(math.log2(max_reprobe + 1) + 1)

    return math.ceil(pow(2, l) * (2 * kmer_size - l + r + 1) / 8)


def get_available_hash_table_size(
    kmer_size: int, context: ProgramContext, counter_len: int, max_reprobe: int
) -> int:
    available_vmem = math.floor(virtual_memory().available * VIRTUAL_MEM_SAFETY_FACTOR)

    jellyfish_options = {
        "mer-len": kmer_size,
        "mem": available_vmem,
        "counter-len": counter_len,
        "reprobes": max_reprobe,
    }

    jf_mem_output = run_subprocess(
        ["jellyfish", "mem"],
        [],
        jellyfish_options,
        context,
    )

    return int(jf_mem_output.stdout.decode("utf-8").split(" ")[0])


def choose_hash_table_size(
    genome_size: int,
    kmer_size: int,
    context: ProgramContext,
    counter_len: int,
    max_reprobe: int = 62,
) -> int:
    optimal_size = get_optimal_hash_table_size(genome_size, kmer_size, max_reprobe)
    available_size = get_available_hash_table_size(
        kmer_size, context, counter_len, max_reprobe
    )

    context.logger.debug(
        f"Selected hash table size: {min(optimal_size, available_size)} (optimal={optimal_size}, available={available_size})"
    )
    return min(optimal_size, available_size)
