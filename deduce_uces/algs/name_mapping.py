from deduce_uces.utils import invert_dict
from typing import Dict, List, NamedTuple


NameMapping = NamedTuple(
    "NameMapping",
    [
        ("sequence_to_uchar", Dict[str, int]),
        ("contig_to_uchar", Dict[str, Dict[str, int]]),
        # And the inverse too, for efficient lookup
        ("uchar_to_sequence", Dict[int, str]),
        ("uchar_to_contig", Dict[int, Dict[int, str]]),
    ],
)


def create_name_mapping(genomes: Dict[str, List[str]]) -> NameMapping:
    # To save space in memory and core kmer tries, and account for different lengths in the sequence names,
    # convert genome and contig names to integers

    assert len(genomes) < 255  # Max number of sequences

    sequence_to_uchar = {name: i for i, name in enumerate(genomes.keys())}
    contig_to_uchar = {
        genome_name: {ctg: i for i, ctg in enumerate(contigs)}
        for genome_name, contigs in genomes.items()
    }

    return NameMapping(
        sequence_to_uchar=sequence_to_uchar,
        contig_to_uchar=contig_to_uchar,
        uchar_to_sequence=invert_dict(sequence_to_uchar),
        uchar_to_contig={
            sequence_to_uchar[sequence_name]: {
                id: contig_name for contig_name, id in contig.items()
            }
            for sequence_name, contig in contig_to_uchar.items()
        },
    )
