from typing import Dict, List, Optional, Union

from pyfaidx import Fasta


class IndexedSeq:
    def __init__(
        self,
        fa: str,
        contig_name_mapping: Optional[Dict[Union[int, str], str]] = None,
        build_index: Optional[bool] = False,
    ) -> None:
        self.index = Fasta(fa, rebuild=False, build_index=build_index)
        self.contig_name_mapping = contig_name_mapping

    def get_contigs(self) -> List[str]:
        return [c.name for c in self.index]

    def fetch(
        self, contig: Union[str, int], start_0: int, end_0: int, core_kmer_len: int
    ):
        resolved_contig = (
            self.contig_name_mapping[contig] if self.contig_name_mapping else contig
        )

        # print(f"\tFETCH: {resolved_contig}:{start_0 + 1}-{end_0 + core_kmer_len}")

        return str(
            self.index.faidx.fetch(resolved_contig, start_0 + 1, end_0 + core_kmer_len)
        ).upper()
