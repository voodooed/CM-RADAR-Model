"""
Locating a CarMaker installation.

Nothing in this package requires CarMaker: the models run standalone. A
CarMaker installation is only needed to load IPG's own data files
(``Data/Sensor/Radar_Default``, ``RCS_*``, ``RadarRSI_Default``,
``MaterialLib``) and to run the validation scripts against them.

Resolution order:
  1. an explicit ``--ipg <dir>`` command-line argument;
  2. the ``CARMAKER_DIR`` environment variable;
  3. nothing - the caller falls back to the internal models / synthetic data.

``CARMAKER_DIR`` must point at the versioned installation directory, i.e. the
one that contains ``Data/``, ``include/`` and ``doc/``:

    export CARMAKER_DIR=/opt/ipg/carmaker/linux64-15.1     # Linux
    set CARMAKER_DIR=C:\\IPG\\carmaker\\win64-15.1          # Windows
"""

from __future__ import annotations

import os
from typing import Optional

ENV_VAR = "CARMAKER_DIR"


def default_carmaker_dir() -> Optional[str]:
    """Installation directory from ``CARMAKER_DIR``, or ``None`` if unset."""
    d = os.environ.get(ENV_VAR, "").strip()
    return d or None


def sensor_data_dir(carmaker_dir: Optional[str]) -> Optional[str]:
    """``<carmaker_dir>/Data/Sensor`` if it exists, else ``None``."""
    if not carmaker_dir:
        return None
    d = os.path.join(carmaker_dir, "Data", "Sensor")
    return d if os.path.isdir(d) else None


def short(path: str, keep: int = 3) -> str:
    """Last ``keep`` components of a path, for log output that should not
    depend on where the installation happens to live."""
    parts = os.path.normpath(path).split(os.sep)
    return ("..." + os.sep if len(parts) > keep else "") + os.sep.join(parts[-keep:])


def missing_message(what: str = "CarMaker installation") -> str:
    return (f"[skip] no {what} found. Pass --ipg <dir> or set "
            f"{ENV_VAR}=<...>/carmaker/<arch>-<version>")
