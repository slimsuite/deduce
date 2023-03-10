from deduce_uces.run import run_subprocess


def index_fasta(fasta_filename: str):
    run_subprocess(["samtools", "faidx", fasta_filename], [], {}, None)
