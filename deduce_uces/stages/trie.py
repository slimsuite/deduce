from typing import Generator, Iterable, NamedTuple
import marisa_trie

CoreKmer = NamedTuple(
    "CoreKmer",
    [("id", int), ("genome", int), ("contig", int), ("position", int)],
)

# Struct specification for packing core kmers
CORE_KMER_TRIE_FORMAT = f"<BII"


def create_trie(core_kmers: Iterable[CoreKmer]) -> marisa_trie.RecordTrie:
    trie_pairs = [
        (
            core_kmer.id,
            (
                core_kmer.genome,
                core_kmer.contig,
                core_kmer.position,
            ),
        )
        for core_kmer in core_kmers
    ]

    return marisa_trie.RecordTrie(CORE_KMER_TRIE_FORMAT, trie_pairs)


def save_trie(t: marisa_trie.RecordTrie, filename: str):
    t.save(filename)


def load_trie_mem(filename: str) -> marisa_trie.RecordTrie:
    t = marisa_trie.RecordTrie(CORE_KMER_TRIE_FORMAT)
    t.load(filename)
    return t


def load_trie_mmap(filename: str) -> marisa_trie.RecordTrie:
    return marisa_trie.RecordTrie(CORE_KMER_TRIE_FORMAT).mmap(filename)


def iter_trie(t: marisa_trie.RecordTrie) -> Generator[CoreKmer, None, None]:
    for id, args in t.items():
        yield CoreKmer(id=id, genome=args[0], contig=args[1], position=args[2])
