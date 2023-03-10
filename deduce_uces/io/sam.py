from deduce_uces.Logger import Logger
from deduce_uces.utils import reverse_complement
import re
from deduce_uces.run import run_subprocess
from typing import Dict, NamedTuple, Iterator, Optional, Tuple, List
import pysam
import os

SAMAlignment = NamedTuple(
    "SAMAlignment",
    [
        ("query_name", str),
        ("flag", int),
        ("reference_name", str),
        ("position", int),  # 0-indexed position
        ("query_sequence", str),
        ("number_of_mismatches", int),
        ("is_mapped", bool),
        ("is_reverse", bool),
    ],
)

# Parses a whole SAM file in one go and yields alignments
# Does not support multiprocessing
def parse_sam_lazy(
    sam_filename: str,
    max_mismatches: int,
    logger: Logger,
    use_md_for_mismatches: Optional[bool] = False,
) -> Iterator[SAMAlignment]:
    filetype = get_filetype(sam_filename)

    # Do a first pass through the file and record all query sequences, from primary
    # mappings
    # This is necessary because secondary mapping records don't include the query
    # sequence
    primary_mapping_sequences: Dict[str, Tuple[str, bool]] = {}

    f1 = pysam.AlignmentFile(  # type: ignore
        sam_filename,
        filetype,
    )

    for alignment in f1.fetch():
        if not alignment.is_secondary and not alignment.is_supplementary:
            primary_mapping_sequences[alignment.query_name] = (
                alignment.query_sequence.upper(),
                alignment.is_reverse,
            )

    # In the second pass, apply filtering constraints
    f2 = pysam.AlignmentFile(sam_filename, filetype)  # type: ignore

    for alignment in f2.fetch():
        if alignment.is_unmapped:
            continue

        primary = primary_mapping_sequences.get(alignment.query_name)

        if primary is None and alignment.query_sequence is None:
            logger.warning(f"no query sequence for query {alignment.query_name}")
            continue

        alignment_query_sequence, query_is_reverse = primary or (
            alignment.query_sequence.upper(),
            alignment.is_reverse,
        )

        if query_is_reverse != alignment.is_reverse:
            alignment_query_sequence = reverse_complement(alignment_query_sequence)

        alignment.query_sequence = alignment_query_sequence

        a = from_pysam_alignment(alignment, use_md_for_mismatches=use_md_for_mismatches)

        if a.number_of_mismatches > max_mismatches:
            continue

        yield a


def get_filetype(sam_filename: str) -> str:
    return "rb" if sam_filename.endswith(".bam") else "r"


def get_index_filename(sam_filename: str) -> str:
    return sam_filename + ".csi"


def get_reference_contigs(sam_filename: str) -> List[str]:
    filetype = get_filetype(sam_filename)

    return pysam.AlignmentFile(sam_filename, filetype).references  # type: ignore


def from_pysam_alignment(
    seg: pysam.AlignedSegment,  # type: ignore
    use_md_for_mismatches: Optional[bool] = False,
) -> SAMAlignment:
    if use_md_for_mismatches:
        number_of_mismatches = count_mismatches_md(
            seg.get_tag("MD"), len(seg.query_sequence)
        )
    else:
        number_of_mismatches = count_mismatches(
            seg.cigartuples, len(seg.query_sequence)
        )

    return SAMAlignment(
        query_name=seg.query_name,
        flag=seg.flag,
        reference_name=seg.reference_name,
        position=seg.reference_start,
        query_sequence=seg.query_sequence,
        number_of_mismatches=number_of_mismatches,
        is_mapped=not seg.is_unmapped,
        is_reverse=seg.is_reverse,
    )


# CIGAR operation mapping
# (From the pysam source code)
# +-----+--------------+-----+
# | M | BAM_CMATCH | 0 |
# +-----+--------------+-----+
# | I | BAM_CINS | 1 |
# +-----+--------------+-----+
# | D | BAM_CDEL | 2 |
# +-----+--------------+-----+
# | N | BAM_CREF_SKIP | 3 |
# +-----+--------------+-----+
# | S | BAM_CSOFT_CLIP | 4 |
# +-----+--------------+-----+
# | H | BAM_CHARD_CLIP | 5 |
# +-----+--------------+-----+
# | P | BAM_CPAD | 6 |
# +-----+--------------+-----+
# |= | BAM_CEQUAL | 7 |
# +-----+--------------+-----+
# | X | BAM_CDIFF | 8 |
# +-----+--------------+-----+
# | B | BAM_CBACK | 9 |
# +-----+--------------+-----+
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 4
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CBACK = 9


def count_mismatches(
    ops: List[Tuple[int, int]],
    sequence_length: int,
) -> int:
    matches = 0

    for op, count in ops:
        if op == BAM_CMATCH or op == BAM_CEQUAL:
            matches += count

    return sequence_length - matches


def count_mismatches_md(
    md: str,
    sequence_length: int,
) -> int:
    matches = sum([int(m) for m in re.findall("[0-9]+", md)])

    return sequence_length - matches


# Create an index in the CSI format, which better supports long chromosomes (e.g. in plants)
def index_bam(filename: str) -> str:
    index_filename = filename.replace(".bam", ".bam.csi")

    if not os.path.exists(index_filename):
        run_subprocess(
            [
                "samtools",
                "index",
            ],
            [filename, "-c", index_filename],
            {},
            context=None,
            short_flags=True,
        )

    return index_filename
