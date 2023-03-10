import os
import subprocess
from typing import NamedTuple, Union, List, BinaryIO, Mapping, Optional

from deduce_uces.Logger import Logger

ProgramContext = NamedTuple(
    "ProgramContext",
    [
        ("threads", int),
        ("working_dir", str),
        ("logger", Logger),
    ],
)


def run_subprocess(
    cmd: List[str],
    args: List[str],
    flags: Mapping[str, Union[int, str, None]],
    context: Optional[ProgramContext],
    short_flags: bool = False,
    stdout_file: BinaryIO = None,
    stdin_file: BinaryIO = None,
) -> subprocess.CompletedProcess:
    cmd_with_args = cmd

    for arg, val in flags.items():
        if val is None:
            continue

        if short_flags:
            cmd_with_args.append(f"-{arg}")
        else:
            cmd_with_args.append(f"--{arg}")
        cmd_with_args.append(str(val))

    cmd_with_args.extend(args)

    working_dir = context.working_dir if context else os.getcwd()

    try:
        if context:
            context.logger.debug("Command: " + " ".join(cmd_with_args))

        if stdout_file is not None or stdin_file is not None:
            return subprocess.run(
                cmd_with_args,
                cwd=working_dir,
                check=True,
                stdout=stdout_file,
                stdin=stdin_file,
            )
        else:
            return subprocess.run(
                cmd_with_args, cwd=working_dir, capture_output=True, check=True
            )
    except subprocess.CalledProcessError as e:
        err_message = (
            f"Command failed with exit code {e.returncode}: {' '.join(cmd_with_args)}"
        )

        if context:
            context.logger.warning(e.stderr)
            context.logger.error(err_message)
        raise Exception(err_message)
