"""
Antenna gain map and RCS map for the CarMaker "Radar Sensor".

Two data sets are needed by the model:

1. **Antenna gain map** - ``FileIdent = CarMaker-AntennaGainMap``
   (Reference Manual -> Sensors -> Radar Sensor -> ... -> Antenna Gain Map File).
   A 2-D map ``Gain[nAzimuth][nElevation]`` in dB, plus the generator parameters
   ``ScanRange``, ``BeamWidth`` and ``AntennaEff``.
   CarMaker ships ``Data/Sensor/Radar_Default``.

2. **RCS map** - ``FileIdent = CarMaker-RCS Map``
   (Reference Manual -> Sensors -> Radar Sensor -> ... -> RCS map files).
   ``Azim.data`` / ``AzimRCS.data`` give RCS over the azimuth of incidence in the
   *object* body frame, plus ``ProbExist`` and ``OcclusionFactor``.
   CarMaker ships ``RCS_Car``, ``RCS_Truck``, ``RCS_Pedestrian``, ``RCS_Bicycle``,
   ``RCS_GuardRailPost``.

Both readers accept the original IPG files unchanged, so the model can be run
against the exact data CarMaker uses.
"""

from __future__ import annotations

import math
from typing import List, Optional, Sequence

from common.infofile import InfoFile

# Half-power beamwidth constant of a uniform aperture:
#   sinc(pi*nu) drops 3 dB at nu = 0.442946  ->  BW = 2*0.442946*lambda/a
APERTURE_BW_CONST = 0.8858929413


def _sinc(x: float) -> float:
    """sin(x)/x with the removable singularity handled."""
    if abs(x) < 1e-12:
        return 1.0
    return math.sin(x) / x


# ===========================================================================
#  Antenna gain map
# ===========================================================================
class AntennaGainMap:
    """One-way antenna power gain [dB] over (azimuth, elevation).

    Storage layout matches the IPG file: ``gain_db[i_az][i_el]`` with equidistant
    azimuth and elevation sample points (Reference Manual: "Specify equidistant
    sample points ... from -90 deg to 90 deg").
    """

    __slots__ = ("az", "el", "gain_db", "beam_width", "scan_range", "antenna_eff")

    def __init__(self, az: Sequence[float], el: Sequence[float],
                 gain_db: Sequence[Sequence[float]],
                 beam_width: Sequence[float] = (20.0, 15.0),
                 scan_range: Sequence[float] = (0.0, 0.0),
                 antenna_eff: float = 1.0):
        self.az = list(az)
        self.el = list(el)
        self.gain_db = [list(r) for r in gain_db]
        self.beam_width = (float(beam_width[0]), float(beam_width[1]))
        self.scan_range = (float(scan_range[0]), float(scan_range[1]))
        self.antenna_eff = float(antenna_eff)

    # ----------------------------------------------------------------- load
    @classmethod
    def load(cls, path: str) -> "AntennaGainMap":
        inf = InfoFile.load(path)
        ident = inf.str("FileIdent", "")
        if not ident.startswith("CarMaker-AntennaGainMap"):
            raise ValueError(f"{path}: not a CarMaker-AntennaGainMap (FileIdent='{ident}')")
        unit = inf.str("Type", "rad").strip().lower()
        scale = 1.0 if unit.startswith("rad") else math.pi / 180.0
        az = [v * scale for v in inf.floats("Azimuth")]
        el = [v * scale for v in inf.floats("Elevation")]
        gain = inf.matrix("Gain")
        if len(gain) != len(az) or len(gain[0]) != len(el):
            raise ValueError(f"{path}: Gain must be nAzimuth x nElevation "
                             f"({len(az)}x{len(el)}), got {len(gain)}x{len(gain[0])}")
        bw = inf.floats("BeamWidth", (20.0, 15.0))
        sr = inf.floats("ScanRange", (0.0, 0.0))
        eff = inf.float("AntennaEff", 1.0)
        return cls(az, el, gain, bw, sr, eff)

    # ------------------------------------------------------------- generate
    @classmethod
    def generate(cls, beam_width_deg: Sequence[float],
                 scan_range_deg: Sequence[float] = (0.0, 0.0),
                 antenna_eff: float = 1.0,
                 n_samples: int = 181,
                 obliquity: bool = True) -> "AntennaGainMap":
        """Internal uniform-rectangular-aperture model.

        Reference Manual -> Sensors -> Radar Sensor -> Antenna Gain:

            f(Theta, Phi) = sin(pi nu_y)/(pi nu_y) * sin(pi nu_z)/(pi nu_z)   (Eq. 507)
            nu_y = a/lambda * sin(Theta) cos(Phi)
            nu_z = b/lambda * sin(Theta) sin(Phi)                             (Eq. 508)

        with Theta the polar angle away from boresight and Phi the roll angle
        about it, so that ``sin(Theta)cos(Phi) = cos(el) sin(az)`` and
        ``sin(Theta)sin(Phi) = sin(el)`` in the azimuth/elevation frame used
        everywhere else.

        The aperture size follows from the 3 dB one-way main-lobe width, and the
        peak gain from the standard uniform-aperture directivity
        ``G0 = 4*pi*a*b/lambda^2``, scaled by AntennaEff.

        ``obliquity`` adds a Huygens/Kirchhoff obliquity factor
        ``(1 + cos Theta) / 2`` to the field.  This factor is **[INFERRED]**: it
        is not in the documentation, but including it reduces the mismatch
        against the shipped ``Data/Sensor/Radar_Default`` map from 1.17 dB to
        0.12 dB mean (see validation/validate_antenna_map.py).
        """
        bw_h = math.radians(beam_width_deg[0])
        bw_v = math.radians(beam_width_deg[1])
        if bw_h <= 0.0 or bw_v <= 0.0:
            raise ValueError("BeamWidth must have a positive non-zero value")
        a_lam = APERTURE_BW_CONST / bw_h          # a / lambda
        b_lam = APERTURE_BW_CONST / bw_v          # b / lambda
        g0_db = 10.0 * math.log10(4.0 * math.pi * a_lam * b_lam * max(antenna_eff, 1e-12))

        # beam steering: shift of the pattern in direction-cosine space
        u0 = math.sin(math.radians(scan_range_deg[0]))
        v0 = math.sin(math.radians(scan_range_deg[1]))

        half = math.pi / 2.0
        az = [-half + 2.0 * half * i / (n_samples - 1) for i in range(n_samples)]
        el = list(az)
        gain = []
        for a in az:
            row = []
            for e in el:
                uy = math.cos(e) * math.sin(a)
                uz = math.sin(e)
                f = _sinc(math.pi * a_lam * (uy - u0)) * _sinc(math.pi * b_lam * (uz - v0))
                if obliquity:
                    cos_theta = math.cos(a) * math.cos(e)
                    f *= 0.5 * (1.0 + cos_theta)
                row.append(g0_db + 20.0 * math.log10(max(abs(f), 1e-9)))
            gain.append(row)
        return cls(az, el, gain, beam_width_deg, scan_range_deg, antenna_eff)

    # --------------------------------------------------------------- lookup
    def gain(self, az: float, el: float) -> float:
        """Bilinearly interpolated one-way gain [dB] (clamped outside the map)."""
        return _bilinear(self.az, self.el, self.gain_db, az, el)

    @property
    def peak_db(self) -> float:
        return max(max(r) for r in self.gain_db)


def _bilinear(xs: List[float], ys: List[float],
              z: List[List[float]], x: float, y: float) -> float:
    """Bilinear interpolation on an equidistant grid; clamps at the edges."""
    nx, ny = len(xs), len(ys)
    dx = (xs[-1] - xs[0]) / (nx - 1)
    dy = (ys[-1] - ys[0]) / (ny - 1)
    fx = (x - xs[0]) / dx
    fy = (y - ys[0]) / dy
    i = int(math.floor(fx))
    j = int(math.floor(fy))
    i = max(0, min(nx - 2, i))
    j = max(0, min(ny - 2, j))
    tx = max(0.0, min(1.0, fx - i))
    ty = max(0.0, min(1.0, fy - j))
    z00, z10 = z[i][j], z[i + 1][j]
    z01, z11 = z[i][j + 1], z[i + 1][j + 1]
    return (z00 * (1 - tx) * (1 - ty) + z10 * tx * (1 - ty)
            + z01 * (1 - tx) * ty + z11 * tx * ty)


# ===========================================================================
#  RCS map
# ===========================================================================
class RcsMap:
    """Azimuth-dependent RCS look-up table of one object class.

    ``rcs_lin[i]`` is the RCS in m^2 at incidence azimuth ``azim[i]`` (rad),
    measured in the *object* body frame (x along the object's longitudinal axis,
    0 deg = illuminated from the front, +/-180 deg = from the rear).
    """

    __slots__ = ("name", "azim", "rcs_lin", "prob_exist", "occlusion_factor")

    def __init__(self, name: str, azim: Sequence[float], rcs_lin: Sequence[float],
                 prob_exist: int = 0, occlusion_factor: float = 0.0):
        self.name = name
        self.azim = list(azim)
        self.rcs_lin = list(rcs_lin)
        self.prob_exist = int(prob_exist)
        self.occlusion_factor = float(occlusion_factor)

    # ----------------------------------------------------------------- load
    @classmethod
    def load(cls, path: str, name: Optional[str] = None) -> "RcsMap":
        inf = InfoFile.load(path)
        ident = inf.str("FileIdent", "")
        if not ident.startswith("CarMaker-RCS Map"):
            raise ValueError(f"{path}: not a CarMaker-RCS Map (FileIdent='{ident}')")
        atype = inf.str("Azim.type", "deg").strip().lower()
        scale = 1.0 if atype.startswith("rad") else math.pi / 180.0
        azim = [v * scale for v in inf.floats("Azim.data")]
        rtype = inf.str("AzimRCS.type", "10log10").strip().lower()
        raw = inf.floats("AzimRCS.data")
        if len(raw) != len(azim):
            raise ValueError(f"{path}: AzimRCS.data ({len(raw)}) != Azim.data ({len(azim)})")
        if rtype == "linear":
            lin = [max(v, 0.0) for v in raw]
        elif rtype == "10log10":
            lin = [10.0 ** (v / 10.0) for v in raw]
        elif rtype == "20log10":
            lin = [10.0 ** (v / 20.0) for v in raw]
        else:
            raise ValueError(f"{path}: unknown AzimRCS.type '{rtype}'")
        return cls(name or path.rsplit("/", 1)[-1], azim, lin,
                   inf.int("ProbExist", 0), inf.float("OcclusionFactor", 0.0))

    # --------------------------------------------------------------- lookup
    def rcs(self, azimuth: float) -> float:
        """Linear RCS [m^2] at one incidence azimuth (wrapped to the table range)."""
        return _interp_wrapped(self.azim, self.rcs_lin, azimuth)

    def rcs_mean(self, az_min: float, az_max: float, n: int = 9) -> float:
        """Mean linear RCS over an azimuth interval.

        Used for the documented "extended objects (object size and sensor
        resolution)" effect: "The RCS map is evaluated over a range of angles
        depending on the distance, sensor resolution and orientation of the
        object".  The averaging rule itself is **[INFERRED]** - we average the
        *linear* RCS (i.e. power), which is the physically meaningful average
        and is consistent with the documented merging rule (Eq. 511).
        """
        if n < 2 or az_max <= az_min:
            return self.rcs(0.5 * (az_min + az_max))
        acc = 0.0
        for i in range(n):
            a = az_min + (az_max - az_min) * i / (n - 1)
            acc += self.rcs(a)
        return acc / n


def _interp_wrapped(xs: List[float], ys: List[float], x: float) -> float:
    """Linear interpolation on a monotone table covering a full 2*pi period."""
    lo, hi = xs[0], xs[-1]
    period = hi - lo
    if period <= 0.0:
        return ys[0]
    # wrap x into [lo, hi)
    t = (x - lo) % period
    x = lo + t
    n = len(xs)
    # binary search
    a, b = 0, n - 1
    while b - a > 1:
        m = (a + b) // 2
        if xs[m] <= x:
            a = m
        else:
            b = m
    span = xs[b] - xs[a]
    if span <= 0.0:
        return ys[a]
    w = (x - xs[a]) / span
    return ys[a] * (1.0 - w) + ys[b] * w
