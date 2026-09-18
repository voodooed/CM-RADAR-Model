"""
Antenna / VRx configuration of the Radar RSI.

File format (Reference Manual -> Sensors -> Radar RSI -> ... -> Antenna / VRx
File), ``FileIdent = CarMaker-RadarRSI_TranceiverConfig``:

    TransmitGain.FoV      = <az deg> <el deg>
    TransmitGain.NSamples = <nAzimuth> <mElevation>
    TransmitGain.Map:     <nAzimuth x mElevation matrix, dB>
    ReceiveGain.FoV       = ...
    ReceiveGain.NSamples  = ...
    ReceiveGain.Map:      ...
    VRxCalibration:       <n x 4: y [m], z [m], weight amplitude, weight phase [rad]>

CarMaker ships ``Data/Sensor/RadarRSI_Default`` (180x180 deg FoV, 181x181
samples, 16 virtual receivers spaced ~1.95 mm = lambda/2 at 77 GHz along y).

Interpolation rules are documented and implemented as stated:
  * transmit  : "computed by linearly interpolating the parameterized gain map.
                 For angles not covered within the field of view, the ray is not
                 sent into the scene."
  * receive   : "evaluates the gain map at the next smaller neighboring sample
                 point."
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np

from common.infofile import InfoFile


class GainMap2D:
    """(azimuth, elevation) -> gain [dB], equidistant samples inside a FoV."""

    def __init__(self, fov_deg: Sequence[float], n_samples: Sequence[int],
                 gain_db: np.ndarray, floor_db: float = -300.0):
        self.fov = (math.radians(fov_deg[0]), math.radians(fov_deg[1]))
        self.n_az = int(n_samples[0])
        self.n_el = int(n_samples[1])
        self.gain_db = np.asarray(gain_db, dtype=np.float64).reshape(self.n_az, self.n_el)
        self.floor_db = floor_db
        self.az = np.linspace(-0.5 * self.fov[0], 0.5 * self.fov[0], self.n_az)
        self.el = np.linspace(-0.5 * self.fov[1], 0.5 * self.fov[1], self.n_el)

    # --------------------------------------------------------------- lookup
    def _indices(self, az: np.ndarray, el: np.ndarray):
        daz = (self.az[-1] - self.az[0]) / max(self.n_az - 1, 1)
        de = (self.el[-1] - self.el[0]) / max(self.n_el - 1, 1)
        fi = (az - self.az[0]) / daz
        fj = (el - self.el[0]) / de
        return fi, fj

    def linear(self, az, el) -> np.ndarray:
        """Bilinear interpolation (transmit path)."""
        az = np.atleast_1d(np.asarray(az, dtype=np.float64))
        el = np.atleast_1d(np.asarray(el, dtype=np.float64))
        fi, fj = self._indices(az, el)
        inside = ((fi >= 0) & (fi <= self.n_az - 1) & (fj >= 0) & (fj <= self.n_el - 1))
        i = np.clip(np.floor(fi).astype(int), 0, self.n_az - 2)
        j = np.clip(np.floor(fj).astype(int), 0, self.n_el - 2)
        tx = np.clip(fi - i, 0.0, 1.0)
        ty = np.clip(fj - j, 0.0, 1.0)
        g = (self.gain_db[i, j] * (1 - tx) * (1 - ty)
             + self.gain_db[i + 1, j] * tx * (1 - ty)
             + self.gain_db[i, j + 1] * (1 - tx) * ty
             + self.gain_db[i + 1, j + 1] * tx * ty)
        return np.where(inside, g, self.floor_db)

    def nearest_lower(self, az, el) -> np.ndarray:
        """Value at the next smaller neighbouring sample point (receive path)."""
        az = np.atleast_1d(np.asarray(az, dtype=np.float64))
        el = np.atleast_1d(np.asarray(el, dtype=np.float64))
        fi, fj = self._indices(az, el)
        inside = ((fi >= 0) & (fi <= self.n_az - 1) & (fj >= 0) & (fj <= self.n_el - 1))
        i = np.clip(np.floor(fi).astype(int), 0, self.n_az - 1)
        j = np.clip(np.floor(fj).astype(int), 0, self.n_el - 1)
        return np.where(inside, self.gain_db[i, j], self.floor_db)

    @property
    def peak_db(self) -> float:
        return float(self.gain_db.max())


class VRxArray:
    """Virtual-receiver configuration: positions [m] and complex weights."""

    def __init__(self, rows: Sequence[Sequence[float]]):
        a = np.asarray(rows, dtype=np.float64).reshape(-1, 4)
        self.y = a[:, 0]
        self.z = a[:, 1]
        self.weight = a[:, 2] * np.exp(1j * a[:, 3])

    def __len__(self) -> int:
        return int(self.y.size)


class TransceiverConfig:
    def __init__(self, tx: GainMap2D, rx: GainMap2D, vrx: Optional[VRxArray] = None):
        self.tx = tx
        self.rx = rx
        self.vrx = vrx

    # ----------------------------------------------------------------- load
    @classmethod
    def load(cls, path: str) -> "TransceiverConfig":
        inf = InfoFile.load(path)
        ident = inf.str("FileIdent", "")
        if "RadarRSI_TranceiverConfig" not in ident:
            raise ValueError(f"{path}: FileIdent '{ident}' is not a "
                             f"CarMaker-RadarRSI_TranceiverConfig")

        def gmap(prefix: str) -> GainMap2D:
            fov = inf.floats(f"{prefix}.FoV")
            ns = inf.floats(f"{prefix}.NSamples")
            data = np.array(inf.floats(f"{prefix}.Map"), dtype=np.float64)
            return GainMap2D(fov, (int(ns[0]), int(ns[1])), data)

        tx = gmap("TransmitGain")
        rx = gmap("ReceiveGain")
        vrx = VRxArray(inf.matrix("VRxCalibration")) if "VRxCalibration" in inf.blocks else None
        return cls(tx, rx, vrx)

    # ------------------------------------------------------------- generate
    @classmethod
    def generate(cls, fov_deg: Sequence[float], beam_width_deg: Sequence[float],
                 n_samples: int = 181, n_vrx: int = 16,
                 wavelength: float = 299792458.0 / 77e9) -> "TransceiverConfig":
        """Uniform-rectangular-aperture gain maps + a lambda/2 ULA.

        Uses the same internal antenna model as the object-list Radar Sensor
        (Reference Manual Eq. 507/508; "An example gain map for an uniform
        rectangular aperture antenna is available.  Details are described in
        Detection Based on Signal-to-Noise Ratio via Equation 507 and
        Equation 508").
        """
        from radar_object_list.maps import APERTURE_BW_CONST
        bw_h = math.radians(beam_width_deg[0])
        bw_v = math.radians(beam_width_deg[1])
        a_lam = APERTURE_BW_CONST / bw_h
        b_lam = APERTURE_BW_CONST / bw_v
        g0 = 10.0 * math.log10(4.0 * math.pi * a_lam * b_lam)
        az = np.linspace(-0.5 * math.radians(fov_deg[0]), 0.5 * math.radians(fov_deg[0]),
                         n_samples)
        el = np.linspace(-0.5 * math.radians(fov_deg[1]), 0.5 * math.radians(fov_deg[1]),
                         n_samples)
        AZ, EL = np.meshgrid(az, el, indexing="ij")
        uy = np.cos(EL) * np.sin(AZ)
        uz = np.sin(EL)
        f = np.sinc(a_lam * uy) * np.sinc(b_lam * uz)     # np.sinc(x)=sin(pi x)/(pi x)
        f = f * 0.5 * (1.0 + np.cos(AZ) * np.cos(EL))     # obliquity factor (inferred)
        gain = g0 + 20.0 * np.log10(np.maximum(np.abs(f), 1e-9))
        gm_tx = GainMap2D(fov_deg, (n_samples, n_samples), gain)
        gm_rx = GainMap2D(fov_deg, (n_samples, n_samples), gain.copy())
        d = 0.5 * wavelength
        rows = [[(i - (n_vrx - 1) / 2.0) * d, 0.0, 1.0, 0.0] for i in range(n_vrx)]
        return cls(gm_tx, gm_rx, VRxArray(rows))

    # --------------------------------------------------------------- access
    def tx_gain(self, az, el) -> np.ndarray:
        return self.tx.linear(az, el)

    def rx_gain(self, az, el) -> np.ndarray:
        return self.rx.nearest_lower(az, el)
