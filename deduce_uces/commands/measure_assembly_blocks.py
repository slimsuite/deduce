import os
import sys
from typing import Dict, List, Literal, NamedTuple, Optional, Tuple, TypeVar, Union

import pysam
from deduce_uces.Genome import Genome
from deduce_uces.Logger import Logger
from deduce_uces.commands.gen_synteny_map import build_simple_synteny_map
from deduce_uces.ext.bedtools import sam_to_bed
from deduce_uces.io.bed import BedItem, read_bed, write_bed
from deduce_uces.mappers.mapper import Mapper, MappingOutputType, MappingParameters
from deduce_uces.mappers.minimap_mapper import MinimapMapper
from deduce_uces.run import ProgramContext
import networkx as nx
import pandas as pd

InputFormat = Union[Literal["gexf"], Literal["bed"], Literal["sam"], Literal["fasta"]]


def read_synteny_map(
    filename: str,
    format: InputFormat,
    uces: Optional[str],
    context: ProgramContext,
) -> nx.DiGraph:
    if format == "gexf":
        return nx.read_gexf(filename)
    elif format == "bed":
        return build_simple_synteny_map(filename)
    elif format == "sam":
        return build_simple_synteny_map(sam_to_bed(filename, context))
    else:
        # FASTA
        if uces is None:
            context.logger.error("FASTA input format requires a UCE file")
            sys.exit(1)

        mapper = Mapper(
            Genome.from_filename(filename),
            strategy=MinimapMapper(long_reads=True),
            mapping_params=MappingParameters(
                mismatches_allowed=0,
                secondary_mapping_limit=None,
                sort=False,
                mapping_stage="completeness",
                filetype=MappingOutputType.SAM,
            ),
            context=context,
        )
        context.logger.info(f"Indexing assembly {filename}...")
        mapper.build_index()

        context.logger.info(f"Mapping UCEs to assembly...")
        sam_filename = mapper.map(os.path.abspath(uces))

        return build_simple_synteny_map(sam_to_bed(sam_filename, context))


# ID -> contig, start, end
UCEPositions = Dict[str, Tuple[str, int, int]]


def read_uce_positions(
    filename: str,
    synteny_map: nx.DiGraph,
    format: InputFormat,
    context: ProgramContext,
) -> UCEPositions:
    positions = {}
    if format == "gexf":
        context.logger.error("Block output not available for GEXF inputs")
        sys.exit(1)
    elif format == "bed":
        with open(filename, "r") as f:
            for b in read_bed(f):
                positions[b.tags["name"]] = (b.chromosome, b.start, b.end)
        return positions
    elif format == "sam":
        alignments = pysam.AlignmentFile(filename, "r")
        for alignment in alignments.fetch():
            if alignment.is_unmapped:
                continue

            positions[alignment.query_name] = (
                alignment.reference_name,
                alignment.reference_start,
                alignment.reference_end,
            )

        return positions
    else:
        context.logger.error("Block output not available for FASTA inputs")
        sys.exit(1)


def get_Ts(genome_map: nx.DiGraph) -> List[List[str]]:
    sources = [n for n, d in genome_map.in_degree() if d == 0]

    Ts = []
    for source in sources:
        T = [source]
        neighbors = list(genome_map.neighbors(source))
        while len(neighbors) > 0:
            T.append(neighbors[0])
            neighbors = list(genome_map.neighbors(neighbors[-1]))
        Ts.append(T)

    return Ts


# Each cell is SCORE, BACKLINK
Cell = Tuple[int, Optional[str]]
Solution = Dict[str, Dict[str, Cell]]

Scores = NamedTuple("Scores", [("match", int), ("mismatch", int), ("insertion", int)])


def dbg(O: Solution):
    table = []
    for x, ys in O.items():
        row = []
        for y in ys:
            row.append(O[x][y])
        table.append(row)

    print("SCORES")

    print("   " + " ".join([f"{x:>3}" for x in list(O.values())[0].keys()]))
    for row, t in zip(table, O.keys()):
        print(t + " |" + "|".join(f"{x[0]:>3}" for x in row))

    print("BACKLINKS")
    print("   " + " ".join([f"{x:>3}" for x in list(O.values())[0].keys()]))
    for row, t in zip(table, O.keys()):
        print(
            t
            + " |"
            + "|".join(f"{x[1]:>3}" if x[1] is not None else " - " for x in row)
        )


def solve(S: nx.DiGraph, T: List[str], scores: Scores) -> Solution:
    O: Solution = {t: {s: (0, None) for s in S} for t in T}

    for i in range(0, len(T)):
        for s in S:
            pred = [x[0] for x in S.in_edges(s)]
            if len(pred) == 0:
                best_s_prime = None
                best_score = None
            else:
                best_s_prime = pred[0]
                best_score = O[T[i - 1]][pred[0]][0]
                for s_prime in pred[1:]:
                    if O[T[i - 1]][s_prime][0] > best_score:
                        best_s_prime = s_prime
                        best_score = O[T[i - 1]][s_prime][0]

            match_score = (
                scores.match + (best_score or 0)
                if T[i] == s
                else scores.mismatch + (best_score or 0)
            )

            if i != 0:
                insertion_score = scores.insertion + O[T[i - 1]][s][0]
                path_score = max(match_score, insertion_score)
                if path_score == match_score:
                    O[T[i]][s] = (match_score, best_s_prime)
                else:
                    O[T[i]][s] = (insertion_score, s)
            else:
                # Has to be the match score
                O[T[i]][s] = (match_score, best_s_prime)

    return O


X = TypeVar("X")


def remove_duplicates(xs: List[X]) -> List[X]:
    seen = set()
    filtered = []
    for x in xs:
        if x in seen:
            continue
        filtered.append(x)
        seen.add(x)
    return filtered


def backtrack(O: Solution) -> Tuple[int, List[str], List[int]]:
    xs = list(O.keys())
    ys = list(list(O.values())[0].keys())

    best_score = O[xs[0]][ys[0]][0]
    best_end = (0, ys[0])

    for i in range(0, len(xs)):
        for y in ys:
            if (
                O[xs[i]][y][0] > best_score
                or O[xs[i]][y][0] == best_score
                and i > best_end[0]
            ):
                best_score = O[xs[i]][y][0]
                best_end = (i, y)

    backtrace = [best_end]

    while O[xs[best_end[0]]][best_end[1]][1] is not None and best_end[0] >= 0:
        best_end = (best_end[0] - 1, O[xs[best_end[0]]][best_end[1]][1])
        backtrace.append(best_end)

    path = remove_duplicates([x[1] for x in reversed(backtrace)])
    t_path = [x[0] for x in reversed(backtrace)][1:]
    return best_score, path, t_path


def run_measure_assembly_blocks(args, logger: Logger):
    context = ProgramContext(
        threads=args.threads,
        working_dir=os.curdir,
        logger=logger,
    )

    os.chdir(context.working_dir)

    scores = Scores(
        match=args.match_score,
        mismatch=args.mismatch_score,
        insertion=args.insertion_score,
    )

    context.logger.info(f"Reading consensus synteny map...")
    consensus_map: nx.DiGraph = nx.read_gexf(args.synteny_map)

    records = []
    for assembly in args.assemblies:
        context.logger.info(f"Building synteny map for assembly: {assembly}")
        target_synteny_map = read_synteny_map(
            assembly, args.input_format, args.uces, context
        )

        # To output blocks, we need their position
        if args.output_blocks:
            positions = read_uce_positions(
                assembly, target_synteny_map, args.input_format, context
            )
        else:
            positions = {}

        Ts = get_Ts(target_synteny_map)
        paths = []

        context.logger.info(f"Evaluating {len(Ts)} paths in assembly: {assembly}")
        for T in Ts:
            best_score, path, t_path = backtrack(solve(consensus_map, T, scores))
            paths.append(([T[i] for i in t_path], best_score))
            records.append([assembly, T[0], best_score, len(path)])

        if args.output_blocks:
            bed_blocks = []

            for path, score in paths:
                if len(path) == 0:
                    continue

                path_start_in_target = None
                path_end_in_target = None
                for p in path:
                    if p in target_synteny_map:
                        path_start_in_target = p
                        break

                if path_start_in_target is None:
                    continue

                for p in reversed(path):
                    if (
                        p in target_synteny_map
                        and target_synteny_map.nodes[p]["contig"]
                        == target_synteny_map.nodes[p]["contig"]
                    ):
                        path_end_in_target = p
                        break

                if path_end_in_target is None:
                    continue

                start_id = (
                    path_start_in_target
                    if positions[path_start_in_target][1]
                    < positions[path_end_in_target][1]
                    else path[-1]
                )
                end_id = (
                    path_end_in_target
                    if positions[path_start_in_target][1]
                    < positions[path_end_in_target][1]
                    else path_start_in_target
                )

                bed_blocks.append(
                    BedItem(
                        chromosome=positions[start_id][0],
                        start=positions[start_id][1],
                        end=positions[end_id][2],
                        tags={
                            "name": len(path),
                            "score": score,
                            "from": start_id,
                            "to": end_id,
                        },
                    )
                )

            with open(f"{assembly}.blocks.bed", "w") as f:
                write_bed(f, bed_blocks)

    df = (
        pd.DataFrame.from_records(
            records, columns=["assembly", "block_start_uce", "score", "path_len"]
        )
        .groupby("assembly")
        .agg(
            n_paths=pd.NamedAgg("path_len", "count"),
            avg_len=pd.NamedAgg("path_len", "mean"),
            median_len=pd.NamedAgg("path_len", "median"),
            min_len=pd.NamedAgg("path_len", "min"),
            max_len=pd.NamedAgg("path_len", "max"),
            avg_score=pd.NamedAgg("score", "mean"),
            median_score=pd.NamedAgg("score", "median"),
            min_score=pd.NamedAgg("score", "min"),
            max_score=pd.NamedAgg("score", "max"),
        )
        .reset_index()
    )

    df.to_csv(sys.stdout)
