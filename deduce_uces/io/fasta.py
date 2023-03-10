from typing import Generator, Iterable, Callable, Tuple, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from deduce_uces.run import run_subprocess, ProgramContext


def read_fasta_sequences(filename: str) -> Generator[SeqRecord, None, None]:
    return SeqIO.parse(filename, "fasta")


def write_fasta(filename: str, seqs: Iterable[Tuple[Union[str, int], str]]):
    with open(filename, "w") as f:
        SeqIO.write(
            [SeqRecord(Seq(s), id=str(name), description="") for name, s in seqs],
            f,
            "fasta",
        )


def rename_fasta_sequences(
    in_filename: str, out_filename: str, context: ProgramContext
):
    with open(in_filename, "rb") as inf:
        with open(out_filename, "wb") as outf:
            # awk trick from https://www.biostars.org/p/53212/#53219
            run_subprocess(
                ["awk", '/^>/{print ">c" ++i; next}{print}'],
                [],
                {},
                context,
                stdout_file=outf,
                stdin_file=inf,
            )
