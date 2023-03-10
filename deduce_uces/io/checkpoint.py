import hashlib
import os
from typing import Mapping, NamedTuple, Any
from contextlib import contextmanager

from deduce_uces.run import ProgramContext

CheckpointedFile = NamedTuple(
    "CheckpointedFile", [("checkpointed", bool), ("filename", str)]
)

# When a breaking change is made to file contents between dedUCE versions, such that old
# files created by dedUCE are no longer compatible, increment this counter.
# All checkpointed files will then be recreated.
MAGIC_DEPENDS_VERSION = 1


@contextmanager
def checkpointed_file(
    label: str, extension: str, depends_on: Mapping[str, Any], context: ProgramContext
):
    depends_on_stringified = (
        ",".join([f"{k}={v}" for k, v in depends_on.items()])
        + ","
        + str(MAGIC_DEPENDS_VERSION)
    )
    hash = hashlib.md5()
    hash.update(depends_on_stringified.encode("utf-8"))

    filename = os.path.join(
        context.working_dir, f"{label}.{hash.hexdigest()}.{extension}"
    )

    if os.path.exists(filename):
        context.logger.debug(f"Using checkpointed file: {filename}")
        yield CheckpointedFile(True, filename)
    else:
        yield CheckpointedFile(False, filename)
