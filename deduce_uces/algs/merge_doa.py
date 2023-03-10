from deduce_uces.stages.trie import (
    CoreKmer,
    iter_trie,
    load_trie_mem,
    load_trie_mmap,
)
import math

# from fast_bio.distance import n_matching_positions
from functools import partial
from itertools import groupby
from multiprocessing import Pool
from operator import attrgetter
from typing import (
    Any,
    Counter,
    Dict,
    Generator,
    Iterable,
    Iterator,
    List,
    NamedTuple,
    Optional,
    Set,
    Tuple,
    TypeVar,
)


from deduce_uces.algs.indexed_seq import IndexedSeq
from deduce_uces.algs.name_mapping import NameMapping
from deduce_uces.algs.uce_params import UCEParams
from deduce_uces.io.sam import SAMAlignment
from deduce_uces.Logger import ConsoleLogger, SplotActivity
from deduce_uces.run import ProgramContext
from deduce_uces.utils import (
    chunk,
    find_in_list,
    flatten,
    hamming,
    reverse_complement,
    unique_by,
)


def sam_alignment_to_core_kmer(
    genome_name: str, name_mapping: NameMapping, a: SAMAlignment
) -> CoreKmer:
    return CoreKmer(
        id=int(a.query_name),
        genome=name_mapping.sequence_to_uchar[genome_name],
        contig=name_mapping.contig_to_uchar[genome_name][a.reference_name],
        position=a.position,
    )


X = TypeVar("X")


def count_differences(a: List[X], b: List[X]) -> int:
    assert len(a) == len(b), "count_differences called on lists of different sizes"

    return sum([1 for ax, bx in zip(a, b) if ax != bx])


def generate_contig_sequences_indexed(
    core_kmers: Iterator[CoreKmer],
    min_len: int,
    uce_params: UCEParams,
    index: IndexedSeq,
) -> Generator[str, None, None]:
    run = None

    ck = safe_next_in_generator(core_kmers)
    while ck is not None:
        if run is None:
            run = [ck]
            ck = safe_next_in_generator(core_kmers)
        else:
            if (
                ck.position
                <= run[-1].position + uce_params.core_kmer_len + uce_params.uce_max_gap
            ):
                run.append(ck)
                ck = safe_next_in_generator(core_kmers)
            else:
                if (
                    run[-1].position
                    - run[0].position
                    + uce_params.core_kmer_len
                    + uce_params.uce_max_gap
                    >= min_len
                ):
                    seq = index.fetch(
                        ck.contig,
                        run[0].position,
                        run[-1].position,
                        uce_params.core_kmer_len,
                    )
                    for subseq in chunk(seq, min_len):
                        yield subseq

                run = None

    if (
        run
        and run[-1].position
        - run[0].position
        + uce_params.core_kmer_len
        + uce_params.uce_max_gap
        >= min_len
    ):
        seq = index.fetch(
            run[0].contig,
            run[0].position,
            run[-1].position,
            uce_params.core_kmer_len + uce_params.uce_max_gap,
        )
        for subseq in chunk(seq, min_len):
            yield subseq


G = TypeVar("G")


def safe_next_in_generator(gen: Iterator[G]) -> Optional[G]:
    try:
        return next(gen)
    except StopIteration:
        return None


# Precondition: core kmers are for a single sequence and contig
# Precondition: core kmers are sorted by start position, ascending
# TODO: make a generator?
def get_candidate_regions(
    core_kmers: Iterable[CoreKmer], uce_params: UCEParams
) -> List[List[CoreKmer]]:
    regions: List[List[CoreKmer]] = []
    in_region = False

    cur_position = 1
    cur_region_len = 0
    cur_region: List[CoreKmer] = []

    # core_kmers may be a regular iterator or a generator
    # if the former, convert to a generator
    core_kmer_generator = (x for x in core_kmers)
    a = safe_next_in_generator(core_kmer_generator)
    while a is not None:
        if not in_region:
            # Start of region
            in_region = True
            cur_position = a.position
            cur_region_len = 1
            cur_region = [a]
            a = safe_next_in_generator(core_kmer_generator)

        else:
            if a.position == cur_position:
                # Multiple core kmers may map to the same position
                # This occurs when allowable mismatches > 0
                cur_region.append(a)
                a = safe_next_in_generator(core_kmer_generator)

            elif (
                a.position
                <= cur_position + uce_params.core_kmer_len + uce_params.uce_max_gap
            ):
                # Still in run, adding to the length
                cur_region_len += a.position - cur_position
                cur_position = a.position
                cur_region.append(a)
                a = safe_next_in_generator(core_kmer_generator)

            else:
                # Out of run
                if (
                    cur_region_len
                    + uce_params.core_kmer_len
                    - 1
                    + uce_params.uce_max_gap
                    >= uce_params.uce_min_length
                ):
                    regions.append(cur_region)

                in_region = False

    if in_region:
        if (
            cur_region_len + uce_params.core_kmer_len - 1 + uce_params.uce_max_gap
            >= uce_params.uce_min_length
        ):
            regions.append(cur_region)
    return regions


# Core Kmer ID -> Genome -> (CoreKmerStart, CoreKmerSeq)
CoreKmerPosition = NamedTuple("CoreKmerPosition", [("contig", int), ("position", int)])
CoreKmerPositions = Dict[int, Dict[int, List[CoreKmerPosition]]]


def get_reference_sequence(
    ref_core_kmers: List[CoreKmer], core_kmer_len: int, index: IndexedSeq
) -> str:
    return index.fetch(
        ref_core_kmers[0].contig,
        ref_core_kmers[0].position,
        ref_core_kmers[-1].position,
        core_kmer_len,
    )


def calculate_homology(ref: str, seq: str) -> float:
    assert len(ref) == len(seq)

    return (
        max(
            hamming(ref, seq),
            hamming(ref, reverse_complement(seq)),
        )
        / float(len(seq))
    )


def is_window_a_uce(
    window_instances: List[CoreKmer],
    other_sequence_names: Set[int],
    core_kmer_positions: CoreKmerPositions,
    uce_params: UCEParams,
    ref_index: IndexedSeq,
    indices: Dict[int, IndexedSeq],
):

    window_size = window_instances[-1].position - window_instances[0].position

    if window_size + uce_params.core_kmer_len + 1 < uce_params.uce_min_length:
        return None

    valid_sequences = other_sequence_names.copy()

    # 0. Map each core kmer in the reference to all matching core kmers in other sequences
    ref_position_to_ref_cks: Dict[int, List[CoreKmer]] = {}
    ref_position_to_matching_cks: Dict[int, Dict[int, List[CoreKmerPosition]]] = {}

    for ck in window_instances:
        if ck.position not in ref_position_to_ref_cks:
            ref_position_to_ref_cks[ck.position] = []

        if ck.position not in ref_position_to_matching_cks:
            ref_position_to_matching_cks[ck.position] = {}

        ref_position_to_ref_cks[ck.position].append(ck)

        for seq_name, seq_cks in core_kmer_positions.get(ck.id, {}).items():
            if seq_name not in ref_position_to_matching_cks[ck.position]:
                ref_position_to_matching_cks[ck.position][seq_name] = []

            ref_position_to_matching_cks[ck.position][seq_name].extend(seq_cks)

    # 1. Find shared contigs for the other sequences
    other_sequence_contigs: Dict[int, Counter[int]] = {
        s: Counter() for s in other_sequence_names
    }

    for other_cks in ref_position_to_matching_cks.values():
        for seq_name, cks_at_position in other_cks.items():
            # Record that a core kmer has appeared at this position
            other_sequence_contigs[seq_name].update(
                set(ck.contig for ck in cks_at_position)
            )

    # 1a. Reject sequences which don't have at least one shared contig for enough positions
    valid_sequences -= set(
        seq_name
        for seq_name, seq_contigs in other_sequence_contigs.items()
        if len(seq_contigs.most_common()) > 0
        and window_size - uce_params.core_kmer_len - 1 - seq_contigs.most_common()[0][1]
        >= uce_params.uce_max_gap
    )

    # 2. Find the homology of the window, reference vs. each other genome
    reference_sequence = get_reference_sequence(
        unique_by(window_instances, key=attrgetter("position")),
        uce_params.core_kmer_len,
        ref_index,
    )

    # If the reference sequence isn't long enough, it's been cut off at the
    # end of the genome. This window can't be a UCE.
    if len(reference_sequence) < uce_params.uce_min_length:
        return None

    invalid_sequences = set()

    for seq_name in valid_sequences:
        seq_core_kmers = []
        for matching in ref_position_to_matching_cks.values():
            seq_core_kmers.extend(matching.get(seq_name, []))

        def does_contig_contain_uce(contig_core_kmers) -> bool:
            for contig_sequence in generate_contig_sequences_indexed(
                (x for x in contig_core_kmers),
                min_len=len(reference_sequence),
                uce_params=uce_params,
                index=indices[seq_name],
            ):
                if (
                    calculate_homology(reference_sequence, contig_sequence)
                    >= uce_params.uce_min_homology
                ):
                    # Found a suitable UCE, no need to keep looking
                    return True

            return False

        for contig in other_sequence_contigs[seq_name]:
            # Get a sorted list of core kmers in the contig
            # There will only be one core kmer per position; the rest are thrown away
            contig_core_kmers = sorted(
                unique_by(
                    [ck for ck in seq_core_kmers if ck.contig == contig],
                    key=attrgetter("position"),
                ),
                key=attrgetter("position"),
            )
            if does_contig_contain_uce(contig_core_kmers):
                break
        else:
            # No suitable UCE found, remove this from the set of genomes
            invalid_sequences.add(seq_name)

    valid_sequences -= invalid_sequences

    # 3. Accept if >= genome threshold
    if len(valid_sequences) >= uce_params.uce_min_support - 1:
        # Got a UCE
        return reference_sequence
    else:
        return None


def build_core_kmer_positions(
    reference_genome_id: int,
    region_core_kmers: Set[int],
    core_kmer_trie_filenames: Dict[int, str],
) -> CoreKmerPositions:
    other_tries = [
        load_trie_mmap(f)
        for i, f in core_kmer_trie_filenames.items()
        if i != reference_genome_id
    ]

    positions: CoreKmerPositions = {}
    for ck_id in region_core_kmers:
        if not ck_id in positions:
            positions[ck_id] = {}

        for t in other_tries:
            for appearance in t.get(ck_id, []):
                sequence, contig, position = appearance
                if not sequence in positions[ck_id]:
                    positions[ck_id][sequence] = []

                positions[ck_id][sequence].append(CoreKmerPosition(contig, position))

    return positions


# Preconditions:
# - region is on a single contig only
# - region is sorted in ascending position
def get_uces_from_region(
    region_info: Tuple[str, List[CoreKmer]],
    ref_genome_id: int,
    name_mapping: NameMapping,
    uce_params: UCEParams,
    core_kmer_trie_filenames: Dict[int, str],
    genome_fasta_filenames: Dict[int, str],
    pickled_logger: Any,
) -> List[str]:
    logger = ConsoleLogger.from_pickleable(pickled_logger)

    # This comes in as a tuple for boring reasons involving pickling
    # in the multiprocessing library
    region_name, region = region_info

    logger.debug(
        f"Identify: finding UCEs for region {region_name} (core_kmers={len(region)})"
    )
    logger.splot_start(SplotActivity.READING)

    region_core_kmers = set(ck.id for ck in region)
    core_kmer_positions = build_core_kmer_positions(
        ref_genome_id, region_core_kmers, core_kmer_trie_filenames
    )

    ref_index = IndexedSeq(
        genome_fasta_filenames[ref_genome_id],
        name_mapping.uchar_to_contig[ref_genome_id],  # type: ignore
    )

    indices = {
        i: IndexedSeq(filename, name_mapping.uchar_to_contig[i])  # type: ignore
        for i, filename in genome_fasta_filenames.items()
    }

    other_sequence_names = set(name_mapping.uchar_to_sequence.keys())

    logger.splot_start(SplotActivity.PROCESSING)

    uces = []
    is_uce = partial(
        is_window_a_uce,
        other_sequence_names=other_sequence_names,
        core_kmer_positions=core_kmer_positions,
        uce_params=uce_params,
        ref_index=ref_index,
        indices=indices,
    )

    # In the unlikely case there is a UCE right on the end of the sequence,
    # it might not be picked up due to the window extension logic
    # To fix this, generate some "dummy" core kmers for the overhanging end
    overhang_core_kmers = []
    for pos in range(uce_params.core_kmer_len - 1):
        overhang_core_kmers.append(
            CoreKmer(
                -pos,  # Negative to mark overhang
                region[-1].genome,
                region[-1].contig,
                region[-1].position + (pos + 1),
            )
        )

    region.extend(overhang_core_kmers)

    region_start = region[0].position
    region_end = region[-1].position

    window_start = region_start
    while window_start <= region_end:
        # Start a new window, of the minimum length
        logger.splot_start(SplotActivity.WINDOW_SEARCHING)
        window_end = window_start + uce_params.uce_min_length - uce_params.core_kmer_len
        if window_end > region_end:
            break
        window = [
            row
            for row in region
            if row.position >= window_start and row.position <= window_end
        ]

        maybe_uce = is_uce(window) if len(window) > 0 else None
        if maybe_uce is not None:
            # This window is a valid UCE
            # Start extending the window until it is no longer valid
            logger.splot_start(SplotActivity.WINDOW_EXTENDING)
            while window_end <= region_end:
                extended_window = window.copy()
                while extended_window == window and window_end <= region_end:
                    window_end += 1
                    extended_window = window + [
                        row for row in region if row.position == window_end - 1
                    ]

                if not is_uce(extended_window):
                    # We've flown too close to the sun
                    break

                window = extended_window

            uces.append(is_uce(window))

            window_start = window_end
        else:
            # Otherwise, all done and no UCE found
            # Advance the window one place forward and try again
            # Gee it's hard work to be a computer
            window_start += 1

    logger.splot_end()
    return uces


GenomeInfo = NamedTuple("GenomeInfo", [("name", str), ("bam", str), ("fa", str)])


def core_kmers_to_uces(
    core_kmer_tries: Dict[int, str],
    ref: GenomeInfo,
    others: List[GenomeInfo],
    name_mapping: NameMapping,
    uce_params: UCEParams,
    context: ProgramContext,
) -> List[str]:
    context.logger.debug("Identify: mapping core kmer positions")
    context.logger.splot_start(SplotActivity.PROCESSING)
    context.logger.debug("Identify: sorting reference sequence core kmers")

    ref_genome_id = name_mapping.sequence_to_uchar[ref.name]
    ref_core_kmers = load_trie_mem(core_kmer_tries[ref_genome_id])
    ref_seq_core_kmers = sorted(
        (ck for ck in iter_trie(ref_core_kmers) if ck.genome == ref_genome_id),
        key=lambda ck: (ck.contig, ck.position),
    )

    context.logger.debug("Identify: filtering other sequence core kmers")

    genome_fastas = {
        name_mapping.sequence_to_uchar[name]: find_in_list(
            others + [ref], lambda x: x.name == name
        ).fa  # type: ignore
        for name in name_mapping.sequence_to_uchar.keys()
    }

    def get_region_id(contig_name: str, idx: int, n_regions: int) -> str:
        digits = math.floor(math.log10(n_regions) + 1)
        return f"{contig_name}_{str(idx).rjust(digits, '0')}"

    # Reference regions are identified by key contig_IDX
    ref_regions = []
    ref_region_core_kmers = {}
    for contig, contig_core_kmers in groupby(ref_seq_core_kmers, attrgetter("contig")):
        context.logger.debug(f"Identify: finding regions for contig {contig}")

        contig_regions = get_candidate_regions(contig_core_kmers, uce_params)

        for i, region in enumerate(contig_regions):
            region_id = get_region_id(contig, i, len(contig_regions))
            ref_region_core_kmers[region_id] = set(ck.id for ck in region)
            ref_regions.append((region_id, region))

    context.logger.debug("Identify: finished identifying regions")

    context.logger.splot_start(SplotActivity.PROCESSING)
    get_uces = partial(
        get_uces_from_region,
        ref_genome_id=ref_genome_id,
        uce_params=uce_params,
        name_mapping=name_mapping,
        core_kmer_trie_filenames=core_kmer_tries,
        genome_fasta_filenames=genome_fastas,
        pickled_logger=context.logger.to_pickleable(),
    )

    context.logger.debug("Identify: sending to thread pool")

    if context.threads == 1:
        # When running single-threaded, avoiding the pool is good
        # for debugging / profiling
        uces = [get_uces(r) for r in ref_regions]
    else:
        with Pool(context.threads) as p:
            uces = list(p.imap_unordered(get_uces, ref_regions))

    context.logger.splot_end()

    return flatten(uces)
