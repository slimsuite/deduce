import sys
from typing import Any, Callable, Generator, List, NamedTuple, Optional, Tuple, TypeVar
from deduce_uces.Logger import Logger
from deduce_uces.algs.find_syntenic_blocks import reduce_until_acyclic
from deduce_uces.io.sam import parse_sam_lazy
from deduce_uces.run import ProgramContext
import itertools
import os
import csv
from operator import itemgetter
import networkx as nx

SyntenyUCE = NamedTuple(
    "SyntenyUCE", [("uce_id", str), ("contig", str), ("position", int)]
)

T = TypeVar("T")
K = TypeVar("K")


def sort_and_group(
    xs: List[T], by: Callable[[T], K]
) -> Generator[Tuple[K, List[T]], Any, Any]:
    for g, gxs in itertools.groupby(sorted(xs, key=by), key=by):
        yield g, list(gxs)


def build_simple_synteny_map(bed_file: str) -> nx.DiGraph:
    with open(bed_file) as f:
        r = csv.reader(f, delimiter="\t")

        uces = [SyntenyUCE(line[3], line[0], int(line[1])) for line in r]

    unique_uces, _ = get_unique_uces(uces)
    assert len(unique_uces) == len(
        set(uce.uce_id for uce in unique_uces)
    ), "Uniqueness did not work"

    m = nx.DiGraph()
    for contig, contig_uces in sort_and_group(unique_uces, itemgetter(1)):
        sorted_contig_uces = sorted(contig_uces, key=itemgetter(2))

        for i in range(len(sorted_contig_uces) - 1):
            m.add_node(sorted_contig_uces[i][0], contig=contig)
            m.add_node(sorted_contig_uces[i + 1][0], contig=contig)
            m.add_edge(
                sorted_contig_uces[i][0],
                sorted_contig_uces[i + 1][0],
                weight=1,
                dist=sorted_contig_uces[i + 1][2] - sorted_contig_uces[i][2],
                contig=contig,
            )

    return m


def get_unique_uces(
    uces: List[SyntenyUCE],
) -> Tuple[List[SyntenyUCE], List[SyntenyUCE]]:
    unique = []
    nonunique = []
    for _, instances in sort_and_group(uces, itemgetter(0)):
        if len(instances) > 1:
            nonunique.append(instances[0])
        else:
            unique.append(instances[0])

    return unique, nonunique


def read_unique_uces_from_sam(
    sam_filename: str, logger: Logger, max_mismatches: int = 0
) -> Tuple[List[SyntenyUCE], List[SyntenyUCE]]:
    uces = [
        SyntenyUCE(a.query_name, a.reference_name, a.position)
        for a in parse_sam_lazy(sam_filename, max_mismatches, logger, False)
    ]
    return get_unique_uces(uces)


def avg(xs):
    return sum(xs) / len(xs)


def scored_shortest_path_or_none(
    g: nx.DiGraph, x: str, y: str, max_len: int
) -> Tuple[Optional[float], Optional[int], Optional[float], Optional[float]]:
    if x not in g:
        return None, None, None, None
    if y not in g:
        return None, None, None, None

    try:
        p = nx.shortest_path(g, source=x, target=y)
        if len(p) > max_len:
            return None, None, None, None

        consensus = []
        for i in range(len(p) - 1):
            consensus.append(g[p[i]][p[i + 1]]["weight"])
        return min(consensus) * (1 / len(p)), len(p), min(consensus), avg(consensus)

    except nx.NetworkXNoPath:
        return None, None, None, None


def assess_synteny_by_pairs(
    uces: List[SyntenyUCE], synteny_map: nx.DiGraph, max_dist=1
) -> Tuple[int, float, float, List[Any], int]:
    syntenic = 0
    score = 0
    consensus_l = []
    total_pairs = 0
    valid_pairs = []
    for _, contig_uces in sort_and_group(uces, itemgetter(1)):
        sorted_uces = list(sorted(contig_uces, key=itemgetter(2)))

        for i in range(len(sorted_uces) - 1):
            pair = (sorted_uces[i], sorted_uces[i + 1])
            total_pairs += 1
            s, d, mc, ac = scored_shortest_path_or_none(
                synteny_map, pair[0][0], pair[1][0], max_dist
            )
            if s and d:
                score += s
                syntenic += 1
                consensus_l.append(mc)
                valid_pairs.append((pair, s, d, mc, ac))

    return syntenic, score, avg(consensus_l), valid_pairs, total_pairs


def merge_synteny_maps(maps: List[nx.DiGraph], dag: bool) -> nx.DiGraph:
    merged = nx.DiGraph()

    for m in maps:
        for n in m.nodes:
            merged.add_node(n)

    for m in maps:
        for e in m.edges(data="weight"):
            if merged.has_edge(e[0], e[1]):
                merged[e[0]][e[1]]["weight"] = (
                    merged[e[0]][e[1]]["weight"] + m[e[0]][e[1]]["weight"]
                )
            else:
                merged.add_edge(e[0], e[1], weight=m[e[0]][e[1]]["weight"])

    if dag and not nx.algorithms.dag.is_directed_acyclic_graph(merged):
        merged = reduce_until_acyclic(
            merged, len(maps), False, consensus_label="weight"
        )

    return merged


def run_gen_synteny_map(args, logger: Logger):
    context = ProgramContext(
        threads=args.threads,
        working_dir=os.curdir,
        logger=logger,
    )

    os.chdir(context.working_dir)

    context.logger.info("Loading UCEs from BED files...")
    bed_files = [os.path.abspath(f) for f in args.beds]

    maps = []
    for b in bed_files:
        context.logger.info(f"Building synteny map for genome: {b}")
        maps.append(build_simple_synteny_map(b))

    context.logger.info(f"Creating consensus synteny map...")
    consensus = merge_synteny_maps(maps, args.dag)

    nx.write_gexf(consensus, args.output)
