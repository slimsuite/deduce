from os import walk
from typing import (
    Dict,
    OrderedDict,
    TypeVar,
    Callable,
    List,
    Optional,
    Any,
)


def get_files_in_directory(dir: str) -> List[str]:
    # Cool! https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
    _, _, filenames = next(walk(dir))
    return filenames


T = TypeVar("T")


def find_in_list(l: List[T], f: Callable[[T], bool]) -> Optional[T]:
    for x in l:
        if f(x):
            return x

    return None


X = TypeVar("X")


# Will get the last occurrence of each item, bucketed by key
def unique_by(xs: List[X], key: Callable[[X], Any]) -> List[X]:
    # https://stackoverflow.com/a/10024750
    seen = OrderedDict()

    for x in xs:
        seen[key(x)] = x

    return list(seen.values())


def with_default(x: Optional[T], default: T) -> T:
    if x is None:
        return default
    return x


def flatten(nested: List[List[T]]) -> List[T]:
    return [x for l in nested for x in l]


def reverse_complement(seq: str) -> str:
    m = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return seq.translate(seq.maketrans(m))[::-1]  # type: ignore


def chunk(xs: str, size: int) -> List[str]:
    # https://stackoverflow.com/a/38163917
    return [xs[i : i + size] for i in range(0, len(xs) - size + 1)]


A = TypeVar("A")
B = TypeVar("B")


def invert_dict(d: Dict[A, B]) -> Dict[B, A]:
    return {b: a for a, b in d.items()}


def hamming(x: str, y: str) -> int:
    assert len(x) == len(y)

    same = 0
    for a, b in zip(x, y):
        if a == b:
            same += 1

    return same


def get_canonical(s):
    return min(s, reverse_complement(s))
