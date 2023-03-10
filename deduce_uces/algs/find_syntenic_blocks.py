import itertools
from typing import Iterable, List, NamedTuple, Optional, Tuple
import numpy as np
import csv
import networkx as nx
import sys

# A block contains a list of UCEs, in order
SyntenicBlock = List[str]

Anchor = NamedTuple(
    "Anchor",
    [
        ("anchor_id", int),
        ("genome", str),
        ("contig", str),
        ("start", int),
        ("end", int),
        ("uce_id", str),
        ("identity", float),
    ],
)


def load_anchors(bed_files: List[str]) -> Tuple[List[Anchor], List[str]]:
    anchors = []
    genome_names = []
    for fname in bed_files:
        name = fname.replace(".bed", "")

        with open(fname, "r") as f:
            records = [
                Anchor(
                    i,
                    name,
                    line[0],
                    int(line[1]),
                    int(line[2]),
                    line[3],
                    float(line[4]),
                )
                for i, line in enumerate(csv.reader(f, delimiter="\t"))
            ]

        anchors.extend(records)
        genome_names.append(name)

    return anchors, genome_names


UCEAnchor = NamedTuple(
    "UCEAnchor",
    [("uce_id", str), ("anchors", List[Anchor]), ("contig_filter", List[str])],
)


def apply_contigs_filter(
    uce_anchors: List[Anchor], contigs_filter: List[str], genome_names: List[str]
) -> List[Anchor]:
    contig_by_genome = {k: v for k, v in zip(genome_names, contigs_filter)}
    return [a for a in uce_anchors if contig_by_genome[a.genome] == a.contig]


def group_anchors_by_uce(
    anchors: List[Anchor], genome_names: List[str]
) -> List[UCEAnchor]:
    uces = []
    for uce_id, uce_anchors_i in itertools.groupby(
        sorted(anchors, key=lambda x: x.uce_id), key=lambda x: x.uce_id
    ):
        uce_anchors = list(uce_anchors_i)
        relevant_contigs = [set([None])] * len(genome_names)
        for genome_sort, genome_anchors in itertools.groupby(
            sorted(uce_anchors, key=lambda x: (genome_names.index(x.genome), x.genome)),
            key=lambda x: (genome_names.index(x.genome), x.genome),
        ):
            genome_i, _ = genome_sort
            relevant_contigs[genome_i] = set(x.contig for x in genome_anchors)

        for contig_filter in itertools.product(*relevant_contigs):
            if None in contig_filter:
                print("WARNING: contig_filter contains None", uce_id, relevant_contigs)

            else:
                uces.append(
                    UCEAnchor(
                        uce_id,
                        apply_contigs_filter(uce_anchors, contig_filter, genome_names),
                        contig_filter,
                    )
                )

    return uces


def anchors_to_tuple(
    uce: UCEAnchor, genome_names: List[str]
) -> Tuple[str, List[Optional[int]]]:
    positions = [None] * len(genome_names)
    for genome, genome_anchors_i in itertools.groupby(
        sorted(uce.anchors, key=lambda x: x.genome), key=lambda x: x.genome
    ):
        genome_anchors = list(genome_anchors_i)
        positions[genome_names.index(genome)] = (
            genome_anchors[0].start if len(genome_anchors) > 0 else None
        )

    return uce.uce_id, positions


def reduce_until_acyclic(
    synteny_map: nx.DiGraph,
    max_genomes: int,
    include_distance: bool,
    consensus_label: str = "support",
) -> nx.DiGraph:
    edges = synteny_map.edges(data=True)
    if include_distance:
        edges_by_support = list(
            sorted(
                [e for e in edges if e[2].get(consensus_label) < max_genomes],
                key=lambda x: (
                    x[2].get(consensus_label),
                    -x[2].get("average_distance"),
                ),
            )
        )
    else:
        edges_by_support = list(
            sorted(
                [e for e in edges if e[2].get(consensus_label) < max_genomes],
                key=lambda x: x[2].get(consensus_label),
            )
        )

    while synteny_map.number_of_edges() > 0:
        if len(edges_by_support) > 0:
            edge_to_delete = edges_by_support.pop()
        else:
            break

        synteny_map.remove_edge(edge_to_delete[0], edge_to_delete[1])

        if nx.algorithms.dag.is_directed_acyclic_graph(synteny_map):
            return synteny_map

    raise Exception("Could not reduce any more")


def find_paths_dfs(
    g_unred: nx.DiGraph, max_genomes: int, include_distance: bool
) -> List[SyntenicBlock]:
    if not nx.algorithms.dag.is_directed_acyclic_graph(g_unred):
        g_unred = reduce_until_acyclic(g_unred, max_genomes, include_distance)

    g = nx.algorithms.transitive_reduction(g_unred)

    sources = [n for n, d in g.in_degree() if d == 0]
    sinks = [n for n, d in g.out_degree() if d == 0]

    def dfs(node):
        if g.out_degree(node) == 0:
            # Leaf
            return [[node]]
        else:
            return [
                [node] + path
                for neighbor in g.neighbors(node)
                for path in dfs(neighbor)
            ]

    paths = []

    for source in sources:
        paths.extend(dfs(source))

    return paths


def compute_score(edges: Iterable[Tuple[str, str, int, float]]):
    # Right now this is just binary!
    return 1


def build_combined_synteny_map(edges: List[Tuple[str, str, int, float]]):
    g = nx.DiGraph()

    sorted_edges = sorted(edges, key=lambda x: (x[0], x[1]))

    for nodes, es in itertools.groupby(sorted_edges, key=lambda x: (x[0], x[1])):
        score = compute_score(es)

        g.add_node(nodes[0])
        g.add_node(nodes[1])
        g.add_edge(nodes[0], nodes[1], weight=score)

    return g


# A syntenic block contains >=min_size UCEs such that:
#
# * the UCEs appear on the same contig in each genome
# * the UCEs appear in the same order in each genome
def find_syntenic_blocks(
    uces: List[UCEAnchor],
    genome_names: List[str],
    min_genomes: int,
    min_length: int,
    include_distance: bool,
    recursion_limit: int,
):
    combined_synteny_edges = []

    blocks = []

    uces_by_contigs = sorted(uces, key=lambda x: x.contig_filter)

    for _, uce_group_i in itertools.groupby(
        uces_by_contigs, key=lambda x: x.contig_filter
    ):
        g = nx.DiGraph()

        nodes = [anchors_to_tuple(uce, genome_names) for uce in uce_group_i]

        for ref_i in range(len(genome_names) - min_genomes + 1):
            sorted_nodes = sorted(nodes, key=lambda x: x[1][ref_i])
            for i in range(len(sorted_nodes)):
                j = i + 1

                while j < len(sorted_nodes):
                    relationship_support = len(
                        [
                            1
                            for n, p in zip(sorted_nodes[j][1], sorted_nodes[i][1])
                            if n > p
                        ]
                    )

                    if relationship_support >= min_genomes:
                        average_distance = np.mean(
                            [
                                n - p
                                for n, p in zip(sorted_nodes[j][1], sorted_nodes[i][1])
                                if n > p
                            ]
                        )

                        g.add_nodes_from(
                            [
                                sorted_nodes[i][0],
                                sorted_nodes[j][0],
                            ]
                        )
                        g.add_edge(
                            sorted_nodes[i][0],
                            sorted_nodes[j][0],
                            support=relationship_support,
                            average_distance=average_distance,
                        )

                        combined_synteny_edges.append(
                            (
                                sorted_nodes[i][0],
                                sorted_nodes[j][0],
                                relationship_support,
                                average_distance,
                            )
                        )

                    j += 1

        sys.setrecursionlimit(recursion_limit)

        blocks.extend(
            p
            for p in find_paths_dfs(g, len(genome_names), include_distance)
            if len(p) >= min_length
        )

    return blocks, build_combined_synteny_map(combined_synteny_edges)


def minimise_blocks(blocks: List[SyntenicBlock]) -> List[SyntenicBlock]:
    minimal = []
    seen = set()

    for block in sorted(blocks, key=lambda x: len(x), reverse=True):
        accept = True
        for uce_id in block:
            if uce_id in seen:
                accept = False
        if accept:
            minimal.append(block)
            for uce_id in block:
                seen.add(uce_id)
    return minimal
