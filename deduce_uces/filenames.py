import hashlib
import os
from pathlib import Path

from deduce_uces.run import ProgramContext
from deduce_uces.cli import HashParams, MergeParams, MergeMethod


def hash_hash_params(params: HashParams) -> str:
    hash = hashlib.md5()
    stringified_params = f"s={params.size},k={params.k},max={params.max_occurrences},min={params.min_occurrences}"
    hash.update(stringified_params.encode("ascii"))
    return hash.hexdigest()


def get_hash_filename(name: str, params: HashParams, context: ProgramContext) -> str:
    return os.path.join(context.working_dir, f"{name}.{hash_hash_params(params)}.jf")


def get_merged_hash_filename(
    name: str,
    hash_params: HashParams,
    merge_params: MergeParams,
    context: ProgramContext,
) -> str:
    hash = hashlib.md5()

    stringified_params = f"min={merge_params.min_occurrences},hash={hash_hash_params(params=hash_params)}"

    # Only add the method to the hashed data if it's not `sum`
    # This maintains backwards compatibility with existing merge names for
    # checkpointing
    if merge_params.method != MergeMethod.SUM:
        stringified_params += f",method={MergeMethod.MIN}"

    hash.update(stringified_params.encode("ascii"))

    return os.path.join(context.working_dir, f"{name}.{hash.hexdigest()}.jf")


def get_genome_sequence_filename(
    name: str, context: ProgramContext, chromosome=None
) -> str:
    if chromosome is not None:
        os.path.join(context.working_dir, f"{name}_chr{chromosome}.fa")

    return os.path.join(context.working_dir, f"{name}.fa")


FILE_EXTENSIONS_TO_STRIP = [".fa", ".jf", ".deduce"]


def strip_filetype(filename: str):
    if any(filename.endswith(x) for x in FILE_EXTENSIONS_TO_STRIP):
        return strip_filetype(Path(filename).stem)

    return filename
