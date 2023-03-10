from enum import Enum
from typing import Protocol, Union, NamedTuple, Set

from deduce_uces.Genome import Genome
from deduce_uces.run import ProgramContext


class MappingOutputType(Enum):
    SAM = 1
    BAM = 2


MatchPosition = NamedTuple(
    "MatchPosition",
    [
        ("start", int),
        ("end", int),
        ("chromosome", str),
        ("core_kmers", Set[str]),
    ],
)

MappingParameters = NamedTuple(
    "MappingParameters",
    [
        ("secondary_mapping_limit", Union[int, None]),
        ("mismatches_allowed", int),
        ("sort", bool),
        ("mapping_stage", str),
        ("filetype", MappingOutputType),
    ],
)


class MappingStrategy(Protocol):
    def __call__(self, *args, **kwargs):
        pass

    def build_index(
        self, reference_name: str, reference_file: str, context: ProgramContext
    ) -> str:
        return ""

    def map(
        self,
        reads_file: str,
        reference_name: str,
        reference_file: str,
        index: str,
        params: MappingParameters,
        context: ProgramContext,
    ) -> str:
        return ""


class Mapper:
    def __init__(
        self,
        genome: Genome,
        strategy: MappingStrategy,
        mapping_params: MappingParameters,
        context: ProgramContext,
    ):
        self.genome = genome
        self.index = None
        self.mapped = None
        self.strategy = strategy
        self.params = mapping_params
        self.context = context

    def build_index(self):
        self.index = self.strategy.build_index(
            self.genome.name, self.genome.source_files[0], self.context
        )

    def map(self, reads_file: str) -> str:
        assert self.index is not None

        self.mapped = self.strategy.map(
            reads_file,
            self.genome.name,
            self.genome.source_files[0],
            self.index,
            self.params,
            self.context,
        )

        return self.mapped

    @property
    def mapped_reads_file(self) -> str:
        if self.mapped is None:
            raise FileNotFoundError("Mapped reads file does not exist")

        return self.mapped
