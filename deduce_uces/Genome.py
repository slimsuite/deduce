from dataclasses import dataclass
import os
from pathlib import Path
from typing import Callable, List, Optional, Union

from Bio.Seq import Seq
from deduce_uces.cli import HashParams
from deduce_uces.ext.jellyfish import hash_fastas, reduce_hash_to_binary_counts
from deduce_uces.filenames import get_genome_sequence_filename, get_hash_filename
from deduce_uces.io.fasta import read_fasta_sequences
from deduce_uces.run import ProgramContext
from deduce_uces.utils import find_in_list, get_files_in_directory


@dataclass
class Genome:
    name: str
    source_files: List[str]
    _hash_file: Union[str, None] = None

    @classmethod
    def from_name(cls, name: str, context: ProgramContext):
        filenames = [
            os.path.join(context.working_dir, f)
            for f in get_files_in_directory(context.working_dir)
        ]

        source_files = []
        if get_genome_sequence_filename(name, context) in filenames:
            source_files.append(get_genome_sequence_filename(name, context))
        else:
            i = 1
            while True:
                if get_genome_sequence_filename(name, context, i) in filenames:
                    source_files.append(get_genome_sequence_filename(name, context, i))
                    i += 1
                else:
                    break

        if len(source_files) == 0:
            raise FileNotFoundError(
                f"Could not find sequences for genome {name}: expected ${get_genome_sequence_filename(name, context)} or chromosomes labelled like ${get_genome_sequence_filename(name, context, 1)}"
            )

        return cls(name=name, source_files=source_files, _hash_file=None)

    @classmethod
    def from_filename(cls, filename: str):
        source_file = os.path.abspath(os.path.basename(filename))
        return cls(
            name=Path(filename).stem,
            source_files=[source_file],
            _hash_file=None,
        )

    @classmethod
    def from_raw_filename(cls, filename: str):
        return cls(
            name=Path(filename).stem,
            source_files=[filename],
            _hash_file=None,
        )

    def hash(self, hash_params: HashParams, context: ProgramContext) -> str:
        if self._hash_file is None:
            potential_hash_filename = get_hash_filename(self.name, hash_params, context)

            if os.path.exists(potential_hash_filename):
                self._hash_file = potential_hash_filename
            else:
                self._hash_file = hash_fastas(
                    self.name, self.source_files, hash_params, context
                )

        return self._hash_file

    def reduce_hash(self, context):
        assert self._hash_file is not None

        self._hash_file = reduce_hash_to_binary_counts(self._hash_file, context)

    # Estimate the genome length in nucleotides
    def estimate_genome_size(self) -> int:
        size_in_bytes = sum(os.path.getsize(filename) for filename in self.source_files)

        # Assuming the file is encoded in ASCII or UTF-8, the sequence characters
        # should be one byte each. Any Unicode characters taking more than one byte
        # will be confined to the header, so this function may overestimate very
        # slightly.
        return size_in_bytes

    def n_sequences(self) -> int:
        return len(self.source_files)

    @property
    def hash_file(self) -> str:
        if self._hash_file is None:
            raise Exception("No hash file associated with genome")

        return self._hash_file

    def __hash__(self):
        return hash(repr(self))
