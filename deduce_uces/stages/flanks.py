from deduce_uces.io.output import UCE, UCEInstance
from typing import List

from deduce_uces.run import ProgramContext


def find_flanks(
    uces: List[UCE],
    args,
    context: ProgramContext,
) -> List[UCE]:
    context.logger.info("Flanks: extending UCE flanks...")

    def extend_flanks(instance: UCEInstance) -> UCEInstance:
        genome, chromosome, start, end, homology = instance
        new_start = max(0, start - args.output_flanks)
        new_end = end + args.output_flanks

        return UCEInstance(genome, chromosome, new_start, new_end, homology)

    uces = [
        UCE(u.id, u.consensus_sequence, [extend_flanks(i) for i in u.instances])
        for u in uces
    ]

    context.logger.info("Flanks: extended UCE flanks")

    return uces
