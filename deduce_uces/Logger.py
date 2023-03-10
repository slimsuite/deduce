import os
import sys
from datetime import datetime
from enum import Enum
from time import localtime, strftime
from typing import Literal, Optional, Protocol, Tuple, Union

LOG_LEVEL_NORMAL = "NORMAL"
LOG_LEVEL_WARNING = "WARNING"
LOG_LEVEL_ERROR = "ERROR"


class SplotActivity(Enum):
    PROCESSING = 1
    READING = 2
    WRITING = 3
    WINDOW_SEARCHING = 4
    WINDOW_EXTENDING = 5
    UNKNOWN = 6


SPLOT_ACTIVITY_TO_COLOR = {
    SplotActivity.PROCESSING: "green",
    SplotActivity.READING: "blue",
    SplotActivity.WRITING: "brown",
    SplotActivity.WINDOW_SEARCHING: "orange",
    SplotActivity.WINDOW_EXTENDING: "red",
    SplotActivity.UNKNOWN: "gray",
}


class Logger(Protocol):
    def debug(self, message: str) -> None:
        """debug"""

    def warning(self, message: str) -> None:
        """warning"""

    def error(self, message: str) -> None:
        """error"""

    def info(self, message: str) -> None:
        """info"""

    def to_pickleable(self):
        """to_pickleable"""

    def splot_start(
        self, activity: SplotActivity, fake_threads: Optional[Tuple[int, str]] = None
    ):
        """splot_start"""

    def splot_end(self, fake_threads: Optional[Tuple[int, str]] = None):
        """splot_end"""

    @classmethod
    def from_pickleable(cls, pickled):
        """from_pickleable"""


class ConsoleLogger:
    def __init__(
        self,
        enable_debug: bool = False,
        show_color: bool = False,
        to_stdout: bool = False,
        logfile_name: Optional[str] = None,
    ):
        self.enable_debug = enable_debug
        self.show_color = show_color
        self.to_stdout = to_stdout
        self.file = sys.stdout if to_stdout else sys.stderr

        if enable_debug:
            self.logfile_name = (
                logfile_name
                if logfile_name is not None
                else os.path.abspath(
                    "deduce_log_" + strftime("%Y-%m-%d_%H-%M-%s", localtime())
                )
            )
            self.logfile = open(self.logfile_name + ".log", "a")
            self.splotfile = open(self.logfile_name + ".splot", "a")

            # self.debug("Initialised logger in debug mode")
            # self.debug(f"Logging to {self.logfile_name}.log")
            # self.debug(f"Splotting to {self.logfile_name}.splot")
        else:
            self.logfile_name = None

    def log(self, message, level):
        tag = ""
        if level == LOG_LEVEL_ERROR:
            tag = "ERROR: "
        elif level == LOG_LEVEL_WARNING:
            tag = "WARNING: "

        if self.enable_debug:
            out = f"[{strftime('%Y-%m-%d %H:%M:%s', localtime())}] [{os.getpid()}] {tag}{message}"

            print(
                out,
                file=self.file,
            )
            self.logfile.write(out + "\n")
            self.logfile.flush()

        else:
            print(f"{tag}{message}", file=self.file)

    def splot_start(
        self, activity: SplotActivity, fake_threads: Optional[Tuple[int, str]] = None
    ):
        if self.enable_debug:
            self._splot(">", activity, fake_threads)

    def splot_end(self, fake_threads: Optional[Tuple[int, str]] = None):
        if self.enable_debug:
            self._splot("<", activity=None, fake_threads=fake_threads)

    def _splot(
        self,
        modifier: str,
        activity: Optional[SplotActivity] = None,
        fake_threads: Optional[Tuple[int, str]] = None,
    ):
        # Splot uses a weird format: 2010-10-21 16:45:09,431
        ts = (
            datetime.utcnow()
            .isoformat(sep=" ", timespec="milliseconds")
            .replace(".", ",")
        )

        if fake_threads is not None:
            thread_ids = [f"{fake_threads[1]}_{i}" for i in range(fake_threads[0])]
        else:
            thread_ids = [os.getpid()]

        out = ""
        for thread_id in thread_ids:
            if activity:
                out = (
                    f"{ts} {modifier}{thread_id} {SPLOT_ACTIVITY_TO_COLOR[activity]}\n"
                )
            else:
                out = f"{ts} {modifier}{thread_id}\n"

        self.splotfile.write(out)
        self.splotfile.flush()

    def debug(self, message):
        if self.enable_debug:
            self.log(message, LOG_LEVEL_NORMAL)

    def warning(self, message):
        self.log(message, LOG_LEVEL_WARNING)

    def error(self, message):
        self.log(message, LOG_LEVEL_ERROR)

    def info(self, message):
        self.log(message, LOG_LEVEL_NORMAL)

    def to_pickleable(self):
        return (
            self.enable_debug,
            self.show_color,
            self.to_stdout,
            # None isn't directly pickleable
            self.logfile_name if self.logfile_name is not None else "NONE",
        )

    @classmethod
    def from_pickleable(cls, pickled):
        enable_debug, show_color, to_stdout, logfile_name = pickled
        return cls(
            enable_debug,
            show_color,
            to_stdout,
            logfile_name if logfile_name != "NONE" else None,
        )


class NullLogger:
    def debug(self, message: str) -> None:
        pass

    def warning(self, message: str) -> None:
        pass

    def error(self, message: str) -> None:
        pass

    def info(self, message: str) -> None:
        pass

    def to_pickleable(self):
        return 0, 0, 0, 0

    def splot_start(
        self, activity: SplotActivity, fake_threads: Optional[Tuple[int, str]] = None
    ):
        pass

    def splot_end(self, fake_threads: Optional[Tuple[int, str]] = None):
        pass

    @classmethod
    def from_pickleable(cls, pickled):
        return cls()
