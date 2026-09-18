"""
Reader for IPG CarMaker "Info File" format.

Format (see any file under <CM_INSTALL>/Data/Sensor, e.g. Radar_Default,
RCS_Car, RadarRSI_Default, MaterialLib):

    #INFOFILE1.1 - Do not remove this line!
    Key = value value value
    Key.Sub = value
    BlockKey:
    <TAB>row of values
    <TAB>row of values

Rules implemented here (derived from inspection of the shipped files):
  * ``#`` starts a comment, ``#INFOFILE`` header line is ignored.
  * ``Key = ...`` is a scalar/vector entry; the value is the rest of the line.
  * ``Key:`` starts a block; every following line that starts with a TAB (or
    other leading whitespace) belongs to the block. Blocks are returned as a
    list of raw lines (whitespace stripped).
  * Keys are case sensitive.  Duplicate keys: last one wins.

This module is intentionally dependency-free so that the same parsing logic can
be ported to C++ with a straight line-by-line translation.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence


class InfoFile:
    """Parsed IPG Info File."""

    def __init__(self) -> None:
        self.scalars: Dict[str, str] = {}
        self.blocks: Dict[str, List[str]] = {}

    # ------------------------------------------------------------------ load
    @classmethod
    def load(cls, path: str) -> "InfoFile":
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return cls.loads(fh.read())

    @classmethod
    def loads(cls, text: str) -> "InfoFile":
        inf = cls()
        current_block: Optional[str] = None
        for raw in text.splitlines():
            if raw.startswith("#"):
                continue
            if not raw.strip():
                current_block = None
                continue

            indented = raw[0] in " \t"
            if indented and current_block is not None:
                inf.blocks[current_block].append(raw.strip())
                continue

            current_block = None
            line = raw.strip()
            # strip trailing comments (keep '#' inside quoted strings out of scope)
            hash_pos = line.find(" #")
            if hash_pos >= 0:
                line = line[:hash_pos].strip()

            if "=" in line:
                key, _, value = line.partition("=")
                inf.scalars[key.strip()] = value.strip()
            elif line.endswith(":"):
                current_block = line[:-1].strip()
                inf.blocks[current_block] = []
        return inf

    # ---------------------------------------------------------------- access
    def has(self, key: str) -> bool:
        return key in self.scalars or key in self.blocks

    def str(self, key: str, default: Optional[str] = None) -> str:
        if key in self.scalars:
            return self.scalars[key]
        if default is None:
            raise KeyError(f"InfoFile: missing key '{key}'")
        return default

    def int(self, key: str, default: Optional[int] = None) -> int:
        if key in self.scalars:
            return int(float(self.scalars[key].split()[0]))
        if default is None:
            raise KeyError(f"InfoFile: missing key '{key}'")
        return default

    def float(self, key: str, default: Optional[float] = None) -> float:
        if key in self.scalars:
            return float(self.scalars[key].split()[0])
        if default is None:
            raise KeyError(f"InfoFile: missing key '{key}'")
        return default

    def floats(self, key: str, default: Optional[Sequence[float]] = None) -> List[float]:
        """Scalar line *or* block, flattened into one list of floats."""
        if key in self.scalars:
            return [float(t) for t in self.scalars[key].split()]
        if key in self.blocks:
            out: List[float] = []
            for line in self.blocks[key]:
                out.extend(float(t) for t in line.split())
            return out
        if default is None:
            raise KeyError(f"InfoFile: missing key '{key}'")
        return list(default)

    def matrix(self, key: str) -> List[List[float]]:
        """Block interpreted as a matrix (one row per line)."""
        if key not in self.blocks:
            raise KeyError(f"InfoFile: missing block '{key}'")
        return [[float(t) for t in line.split()] for line in self.blocks[key]]


def find_data_file(name: str, search_dirs: Sequence[str]) -> str:
    """Resolve an IPG data-file reference (e.g. ``RCS_Car``) against a list of
    directories, mimicking CarMaker's project/data-pool/installation lookup
    order."""
    for d in search_dirs:
        cand = os.path.join(d, name)
        if os.path.isfile(cand):
            return cand
    raise FileNotFoundError(f"data file '{name}' not found in {list(search_dirs)}")
