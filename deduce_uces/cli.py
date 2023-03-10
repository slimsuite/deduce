import argparse
import multiprocessing
import os
from enum import Enum
from typing import NamedTuple, Union


def type_uce_min_genomes(arg):
    try:
        genome = int(arg)
    except ValueError:

        raise argparse.ArgumentTypeError(
            "Minimum UCE genomes (uce-min-genomes) must be an integer"
        )
    if genome < 2:
        raise argparse.ArgumentTypeError(
            "Minimum UCE genomes (uce-min-genomes) must be between 2 and the number of genomes provided"
        )
    return genome


def type_uce_min_length(arg):
    try:
        uce_min_length = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Minimum UCE length (uce-min-length) must be an integer"
        )

    if uce_min_length <= 0:
        raise argparse.ArgumentTypeError(
            "Minimum UCE length (uce-min-length) must be greater than 0"
        )

    return uce_min_length


def make_argparser():
    parser = argparse.ArgumentParser(
        description="Find ultraconserved elements across multiple genomes"
    )

    subparsers = parser.add_subparsers(help="Sub-command help")

    find = subparsers.add_parser("find", help="find UCEs in input genomes")
    preprocess = subparsers.add_parser(
        "preprocess", help="preprocess a genome for use with dedUCE"
    )
    validate = subparsers.add_parser(
        "validate", help="compare files to verify dedUCE is working correctly"
    )
    gen_syntenic_blocks = subparsers.add_parser(
        "gen-syntenic-blocks",
        help="generate syntenic blocks of a certain length from UCE output",
    )
    gen_synteny_map = subparsers.add_parser(
        "gen-synteny-map",
        help="generate a consensus synteny map from UCE output",
    )
    measure_assembly_pairs = subparsers.add_parser(
        "measure-assembly-pairs",
        help="use a consensus synteny map and UCE pairs to evaluate the completeness of an assembly",
    )
    measure_assembly_blocks = subparsers.add_parser(
        "measure-assembly-blocks",
        help="use a consensus synteny map to evaluate the ordering of an assembly",
    )
    measure_assembly_contiguity = subparsers.add_parser(
        "map-syntenic-blocks",
        help="map syntenic blocks to evaluate the ordering of an assembly",
    )
    measure_assembly_completeness = subparsers.add_parser(
        "measure-assembly-completeness",
        help="use UCEs to evaluate completeness of an assembly",
    )

    all_parsers = [
        find,
        preprocess,
        validate,
        gen_syntenic_blocks,
        gen_synteny_map,
        measure_assembly_contiguity,
        measure_assembly_pairs,
        measure_assembly_completeness,
        measure_assembly_blocks,
    ]

    # Global configuration
    # Relates to the operation of the program, not domain-specific parameters
    for p in all_parsers:
        p.add_argument(
            "--debug",
            action="store_true",
            help="enable verbose debug logs",
        )

        p.add_argument(
            "--threads",
            type=int,
            default=multiprocessing.cpu_count(),
            help="the number of CPU threads to use",
        )

        p.add_argument(
            "--version",
            action="store_true",
            help="print the version of dedUCE",
        )

        p.add_argument(
            "--no-color",
            action="store_true",
            help="disable color in program output",
        )

    # PREPROCESS subcommand
    preprocess.set_defaults(command="preprocess")
    preprocess.add_argument(
        "genome",
        metavar="GENOME",
        type=str,
        help="the genome FASTA file to preprocess",
    )

    # FIND subcommand
    find.set_defaults(command="find")
    find.add_argument(
        "dir",
        metavar="DIR",
        type=str,
        nargs="?",
        default=os.getcwd(),
        help="the directory containing genome sequences and previously-generated intermediate files",
    )

    find.add_argument(
        "--skip-core-kmer-identification",
        action="store_true",
        help="use existing BAM files",
    )

    find.add_argument(
        "--reduced-support-exhaustive",
        action="store_true",
        help="if finding UCEs in a subset of genomes, check ALL reference genomes to enable comparison with some other methods",
    )

    find.add_argument(
        "--jf-disk",
        action="store_true",
        help="use Jellyfish's --disk mode",
    )

    # GEN_SYNTENY_MAP subcommand
    gen_synteny_map.set_defaults(command="gen_synteny_map")
    gen_synteny_map.add_argument(
        "beds",
        metavar="BED",
        type=str,
        nargs="+",
        help="a list of BED files containing UCE position information",
    )
    gen_synteny_map.add_argument(
        "--output",
        type=str,
        default="consensus.deduce.gexf",
        help="the name of the output graph file",
    )
    gen_synteny_map.add_argument(
        "--dag",
        action="store_true",
        help="remove edges until acyclic",
    )

    # GEN_SYNTENIC_BLOCKS subcommand
    gen_syntenic_blocks.set_defaults(command="gen_syntenic_blocks")
    gen_syntenic_blocks.add_argument(
        "beds",
        metavar="BED",
        type=str,
        nargs="+",
        help="a list of BED files containing UCE position information",
    )
    gen_syntenic_blocks.add_argument(
        "--syntenic-block-min-length",
        type=int,
        default=5,
        help="the minimum number of UCEs which must appear collinearly to form a syntenic block",
    )
    gen_syntenic_blocks.add_argument(
        "--syntenic-block-consensus-support",
        type=float,
        default=1,
        help="the percentage of genomes for which a syntenic relationship must hold to appear in the synteny map",
    )
    gen_syntenic_blocks.add_argument(
        "--max-recursion-depth",
        type=int,
        default=1000,
        help="tune the Python interpreter's max recursion depth, if there are really long blocks",
    )
    gen_syntenic_blocks.add_argument(
        "--include-distance",
        action="store_true",
        help="when removing edges until the synteny map is acyclic, weight closer UCEs higher (quite slow)",
    )
    gen_syntenic_blocks.add_argument(
        "--minimise-blocks",
        action="store_true",
        help="only allow each UCE to appear in a single block (prevents overlap)",
    )

    measure_assembly_blocks.set_defaults(command="measure_assembly_blocks")
    measure_assembly_blocks.add_argument(
        "synteny_map",
        metavar="MAP",
        type=str,
        help="a file containing a consensus synteny map output by gen-synteny-map",
    )
    measure_assembly_blocks.add_argument(
        "assemblies",
        metavar="ASSEMBLY",
        type=str,
        nargs="+",
        help="files containing each target assembly",
    )
    measure_assembly_blocks.add_argument(
        "--input-format",
        type=str,
        choices=["gexf", "bed", "sam", "fasta"],
        default="gexf",
        help="the input format for target assemblies (gexf = synteny map, bed / sam = result of mapping UCEs to assemblies, fasta = assembly sequences)",
    )
    measure_assembly_blocks.add_argument(
        "--output-format",
        type=str,
        choices=["csv"],
        default="csv",
        help="the output format",
    )
    measure_assembly_blocks.add_argument(
        "--output-blocks",
        action="store_true",
        help="output all blocks found in the assembly in BED format",
    )
    measure_assembly_blocks.add_argument(
        "--match-score",
        type=int,
        default=7,
        help="the score for matching nodes in the alignment",
    )
    measure_assembly_blocks.add_argument(
        "--mismatch-score",
        type=int,
        default=-3,
        help="the score for mismatching nodes in the alignment (deletions/mismatches in the target assembly)",
    )
    measure_assembly_blocks.add_argument(
        "--insertion-score",
        type=int,
        default=-1,
        help="the score for inserted nodes (insertions in the target assembly)",
    )
    measure_assembly_blocks.add_argument(
        "--uces",
        type=str,
        required=False,
        help="a FASTA file containing UCE sequences (required if input format is FASTA)",
    )

    measure_assembly_pairs.set_defaults(command="measure_assembly_pairs")
    measure_assembly_pairs.add_argument(
        "uces",
        metavar="UCES",
        type=str,
        help="a FASTA file containing all UCE sequences",
    )
    measure_assembly_pairs.add_argument(
        "synteny_map",
        metavar="MAP",
        type=str,
        help="a file containing a synteny map output by gen-synteny-map",
    )
    measure_assembly_pairs.add_argument(
        "assemblies",
        metavar="ASSEMBLY",
        type=str,
        nargs="+",
        help="FASTA files containing each assembly to measure",
    )
    measure_assembly_pairs.add_argument(
        "--output-format",
        type=str,
        choices=["csv"],
        default="csv",
        help="the output format",
    )
    measure_assembly_pairs.add_argument(
        "--output-pairs",
        action="store_true",
        help="output all pairs found in the assembly",
    )
    measure_assembly_pairs.add_argument(
        "--max-mapping-mismatches",
        type=int,
        default=0,
        help="the maximum mismatches to tolerate when finding UCEs in the assembly",
    )
    measure_assembly_pairs.add_argument(
        "--max-gap",
        type=int,
        default=3,
        help="the maximum gap in the path when checking a pair in the synteny map",
    )

    # MEASURE_ASSEMBLY_CONTIGUITY subcommand
    measure_assembly_contiguity.set_defaults(command="measure_assembly_contiguity")
    measure_assembly_contiguity.add_argument(
        "uces",
        metavar="UCES",
        type=str,
        help="a FASTA file containing all UCE sequences",
    )
    measure_assembly_contiguity.add_argument(
        "blocks",
        metavar="blocks",
        type=str,
        help="a file containing syntenic blocks output by gen-syntenic-blocks",
    )
    measure_assembly_contiguity.add_argument(
        "assemblies",
        metavar="ASSEMBLY",
        type=str,
        nargs="+",
        help="FASTA files containing each assembly to measure",
    )
    measure_assembly_contiguity.add_argument(
        "--output-format",
        type=str,
        choices=["csv", "human"],
        default="human",
        help="whether to output in human-readable or CSV format",
    )
    measure_assembly_contiguity.add_argument(
        "--max-mapping-mismatches",
        type=int,
        default=3,
        help="the maximum mismatches to tolerate when finding blocks in the assembly",
    )
    measure_assembly_contiguity.add_argument(
        "--max-secondary-alignments",
        type=int,
        default=100,
        help="the maximum number of secondary alignments to report during mapping",
    )

    # MEASURE_ASSEMBLY_COMPLETENESS subcommand
    measure_assembly_completeness.set_defaults(command="measure_assembly_completeness")
    measure_assembly_completeness.add_argument(
        "uces",
        metavar="UCES",
        type=str,
        help="a FASTA file containing all UCE sequences",
    )
    measure_assembly_completeness.add_argument(
        "assemblies",
        metavar="ASSEMBLY",
        type=str,
        nargs="+",
        help="FASTA files containing each assembly to measure",
    )
    measure_assembly_completeness.add_argument(
        "--output-format",
        type=str,
        choices=["csv", "human"],
        default="human",
        help="whether to output in human-readable or CSV format",
    )

    measure_assembly_completeness.add_argument(
        "--output-uces",
        action="store_true",
        help="output the IDs of all UCEs found in the target assembly",
    )

    parsers_with_standard_flags = [find, preprocess]

    for p in parsers_with_standard_flags:
        # UCE configuration
        # Relates to the definition of ultraconserved elements used by the program
        p.add_argument(
            "--uce-min-length",
            type=type_uce_min_length,
            default=200,
            help="the minimum number of base pairs required to constitute a UCE",
        )

        p.add_argument(
            "--uce-core-min-homology",
            type=float,
            default=1,
            help="the minimum percent identity required to constitute a UCE during the core stage",
        )

        p.add_argument(
            "--uce-extension-min-homology",
            type=float,
            help="the minimum percent identity required to constitute a UCE during the extension stage (defaults to uce-core-min-homology)",
        )

        p.add_argument(
            "--core-kmer-mapping-limit",
            type=int,
            default=1000,
            help="the maximum number of hits recorded for each core kmer",
        )

        p.add_argument(
            "--uce-max-occurrences",
            type=int,
            default=1,
            help="the maximum number of times a sequence can appear in a single genome and be considered a UCE",
        )

        # Algorithm
        p.add_argument(
            "--core-kmer-size",
            type=int,
            default=50,
            help="the size of kmers to use when finding UCE candidates",
        )

        p.add_argument(
            "--core-kmer-mismatches",
            type=int,
            default=0,
            help="the number of permitted mismatches when mapping core kmers on genomes",
        )

        p.add_argument(
            "--core-kmer-max-gap",
            type=int,
            default=0,
            help="the permitted extra gap length between core kmers during merging (increase this if minimap is missing core kmers)",
        )

        p.add_argument(
            "--core-kmer-threshold",
            type=float,
            default=1,
            help="the % of genomes that must include a core kmer",
        )

        p.add_argument(
            "--allopolyploids",
            action="store_true",
            help="input genomes are allopolyploids (chromosomes duplicated)",
        )

        # Input configuration
        p.add_argument(
            "--genomes",
            type=str,
            nargs="*",
            default=None,
            help="the names of the genomes to process",
        )

        p.add_argument(
            "--reference",
            type=str,
            default=None,
            help="the reference genome to use for FASTA output and %identity calculation. If not set, will choose randomly.",
        )

        # Output configuration
        p.add_argument(
            "--output-flanks",
            type=int,
            default=0,
            help="the size of the flanks of each UCE to output, in nucleotides",
        )

        p.add_argument(
            "--output",
            type=str,
            default="uces.deduce",
            help="the prefix of the output FASTA file containing identified UCEs",
        )

        p.add_argument(
            "--output-format",
            type=str,
            choices=["fasta", "bed", "all"],
            default="fasta",
            help="the output format",
        )

        # Tool configuration
        # Relates to parameters used by external tools called by the program
        p.add_argument(
            "--jf-hash-size",
            type=str,
            default=None,
            help="the size of the Jellyfish hashes to create. Higher = more memory but faster",
        )

        p.add_argument(
            "--mapper",
            default="minimap",
            choices=["minimap", "bowtie"],
            help="the tool to use when aligning candidate UCEs to the input genomes",
        )

    # VALIDATE subcommand
    validate.set_defaults(command="validate")

    validate.add_argument(
        "--format",
        type=str,
        choices=["fasta", "bed"],
        default="fasta",
        help="the type of files to compare",
    )

    validate.add_argument(
        "--minimise",
        action="store_true",
        help="combine UCEs which are substrings of another UCE",
    )

    validate.add_argument(
        "--show-missing",
        action="store_true",
        help="output the exact UCEs which are expected but not found",
    )

    validate.add_argument(
        "deduce_uces",
        metavar="DEDUCE_UCES",
        type=str,
        help="a BED or FASTA file containing the UCE output from dedUCE",
    )

    validate.add_argument(
        "published_uces",
        metavar="PUBLISHED_UCES",
        type=str,
        help="a BED or FASTA file containing published UCEs",
    )

    return parser


def parse_args(parser, args=None):
    if args:
        return parser.parse_args(args)
    else:
        return parser.parse_args()


HashParams = NamedTuple(
    "HashParams",
    [
        ("size", str),
        ("k", int),
        ("max_occurrences", Union[int, None]),
        ("min_occurrences", int),
        ("counter_len", int),
        ("out_counter_len", int),
        ("disk_mode", bool),
    ],
)


class MergeMethod(Enum):
    SUM = 1
    MIN = 2


MergeParams = NamedTuple(
    "MergeParams",
    [
        ("method", MergeMethod),
        ("min_occurrences", int),
    ],
)
