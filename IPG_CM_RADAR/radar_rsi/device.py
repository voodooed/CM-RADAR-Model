"""
Device Model of the Radar RSI: Front-End Hardware + Back-End Software.

Reference Manual -> Sensors -> Radar RSI -> Device Model / Front-End Hardware /
Back-End Software.  The chain implemented here is exactly the documented one:

    interaction points
      -> Range-Doppler-Processing   (2-D DFT, windows, zero padding, aliasing)
      -> [+ Front-End noise]        (Eq. 525, scaled by NoiseScaling LUTs)
      -> Range-Doppler-Filtering    (OS-CFAR, Eq. 526)  [+ 2D Peak-Finder]
      -> either  VRx output (complex amplitudes per virtual receiver)
         or      Angular-Processing (spatial DFT over a lambda/2 ULA)
                 -> Angular-Filtering (CA-CFAR) [+ Peak-Finder]
                 -> Peak-Interpolation (Eq. 527-529)
      -> point cloud (Cartesian or Spherical)

Data-cube layout matches ``tUserRadarCube`` of the shipped GPU Coding Interface
sample (``Examples/GPUCodingInterface/.../GPUCodingInterface.h``):

    index = index_range + index_velocity * n_r + index_vrx * n_r * n_v

Normalisation: "All detections are normalized to compensate for any signal
processing gain or loss.  The power of the output therefore corresponds to the
signal power at the end of the Front-End Hardware."  Every DFT here is therefore
divided by the coherent gain of its window, so cube values stay
input-referred complex voltage amplitudes in sqrt(W).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np

from radar_rsi.config import (K_BOLTZMANN, CfarConfig, RadarRSIConfig,
                              WindowFunction)
from radar_rsi.propagation import InteractionPoints

REFERENCE_IMPEDANCE = 50.0      # [Ohm], for the VRx output in millivolt [INFERRED]


# ===========================================================================
#  outputs
# ===========================================================================
@dataclass
class DetectionPoint:
    """``tDetPoint`` of ``Sensor_RadarRSI.h``."""
    coordinates: Tuple[float, float, float]   # cartesian [m] or (r, az, el)
    power_dbm: float
    velocity: float                           # relative radial velocity [m/s]


@dataclass
class VRxDetection:
    """``tDetVRx`` of ``Sensor_RadarRSI.h``."""
    range_m: float
    velocity: float
    amp_vrx: np.ndarray                       # complex, millivolt


@dataclass
class RadarRSIOutput:
    """``tRadarRSI`` of ``Sensor_RadarRSI.h``."""
    n_detections: int = 0
    time_fired: float = 0.0
    det_points: List[DetectionPoint] = None      # type: ignore[assignment]
    det_vrx: List[VRxDetection] = None           # type: ignore[assignment]

    def __post_init__(self):
        if self.det_points is None:
            self.det_points = []
        if self.det_vrx is None:
            self.det_vrx = []


# ===========================================================================
#  window functions
# ===========================================================================
def window(kind: WindowFunction, n: int) -> np.ndarray:
    """Reference Manual: "Either a Von-Hann-window (Hann) or a rectangular
    window (Rectangular) can be applied." """
    if kind == WindowFunction.HANN:
        return 0.5 - 0.5 * np.cos(2.0 * math.pi * np.arange(n) / max(n - 1, 1))
    return np.ones(n)


def dft_kernel(n_samples: int, n_bins: int, w: np.ndarray,
               delta_bins: np.ndarray) -> np.ndarray:
    """Normalised DTFT of a windowed complex tone, evaluated ``delta_bins`` away
    from its true (fractional) bin position.

    ``X(delta) = (1 / sum w) * sum_{n<N_s} w[n] exp(j 2 pi n delta / N_bins)``

    This is the analytic result of the DFT the device model performs, so using
    it to synthesise the Range-Doppler map reproduces the documented effects -
    finite resolution, spectral leakage, windowing and (through the periodicity
    in ``delta``) aliasing - without simulating the time-domain samples.
    """
    n = np.arange(n_samples)
    phase = np.exp(2j * math.pi * np.outer(delta_bins, n) / n_bins)
    return (phase @ w) / np.sum(w)


# ===========================================================================
#  Front-End: noise
# ===========================================================================
def noise_power(cfg: RadarRSIConfig) -> float:
    """Reference Manual Equation 525:
        P_Noise = k_B * T_Noise * B_Noise * 10^(NF/10)
    """
    return (K_BOLTZMANN * cfg.noise_temperature_k * cfg.noise_bandwidth_hz
            * 10.0 ** (cfg.noise_figure_db / 10.0))


def _lut_db(lut: Sequence[Tuple[float, float]], x: np.ndarray) -> np.ndarray:
    """Linear interpolation of a (sample point, dB) look-up table, clamped."""
    if not lut:
        return np.zeros_like(x)
    xs = np.array([p[0] for p in lut], dtype=np.float64)
    ys = np.array([p[1] for p in lut], dtype=np.float64)
    return np.interp(x, xs, ys)


# ===========================================================================
#  Back-End
# ===========================================================================
class DeviceModel:
    def __init__(self, cfg: RadarRSIConfig, vrx_y: np.ndarray,
                 vrx_weight: np.ndarray, seed: int = 1):
        cfg.validate()
        self.cfg = cfg
        self.rng = np.random.default_rng(seed)
        self.vrx_y = np.asarray(vrx_y, dtype=np.float64)
        self.vrx_w = np.asarray(vrx_weight, dtype=np.complex128)
        self.n_vrx = self.vrx_y.size

        self.n_r = cfg.n_range_bins
        self.n_v = cfg.n_doppler_bins
        self.w_r = window(cfg.window_function, cfg.range_samples)
        self.w_v = window(cfg.window_function, cfg.doppler_samples)
        self.w_a = window(cfg.window_function, min(cfg.azimuth_samples, self.n_vrx))

        self.range_axis = np.linspace(0.0, cfg.range_max, self.n_r)
        self.vel_axis = np.linspace(cfg.doppler_min_max[0], cfg.doppler_min_max[1], self.n_v)

        self.p_noise = noise_power(cfg)

    # ------------------------------------------------- 1. signal generation
    def build_cube(self, ip: InteractionPoints) -> np.ndarray:
        """Range-Doppler-VRx data cube from the interaction points.

        Two synthesis modes (``cfg.rd_synthesis``):
          * ``"spectral"`` - each IA point is added as a windowed DFT kernel
            around its exact fractional range/Doppler bin.  Models resolution,
            leakage and aliasing.  (default; [INFERRED but equivalent to what a
            real 2-D DFT of the beat signal produces])
          * ``"binning"``  - the amplitude is added to the nearest bin only,
            exactly as the shipped GPU-Coding-Interface sample
            ``Samples/4_Radar/.../Radar.cu`` does.
        """
        cfg = self.cfg
        cube = np.zeros((self.n_r, self.n_v, self.n_vrx), dtype=np.complex128)
        if len(ip) == 0:
            return cube

        # ---- axis mapping (same convention as the shipped Radar.cu sample) --
        r = np.clip(ip.range_m, 0.0, cfg.range_max)
        v_lo, v_hi = cfg.doppler_min_max
        if cfg.doppler_aliasing:
            # v_min / v_max are the *unambiguous* velocities: anything outside
            # folds back, exactly as a real FMCW Doppler DFT does.
            v = v_lo + np.mod(ip.velocity - v_lo, v_hi - v_lo)
        else:
            v = np.clip(ip.velocity, v_lo, v_hi)
        p_r = r / cfg.range_max * (self.n_r - 1)
        p_v = ((v - cfg.doppler_min_max[0])
               / (cfg.doppler_min_max[1] - cfg.doppler_min_max[0])) * (self.n_v - 1)

        # ---- ULA steering vector for every IA point -----------------------
        # d = lambda/2 spacing -> phase = 2 pi y sin(az) / lambda
        steer = np.exp(2j * math.pi * np.outer(np.sin(ip.azimuth), self.vrx_y)
                       / cfg.wavelength) * self.vrx_w[None, :]

        if cfg.rd_synthesis == "binning":
            ir = np.rint(p_r).astype(int)
            iv = np.rint(p_v).astype(int)
            contrib = ip.amplitude[:, None] * steer
            np.add.at(cube, (ir, iv), contrib)
            return cube

        # ---- spectral synthesis -------------------------------------------
        h = cfg.rd_kernel_halfwidth
        offs = np.arange(-h, h + 1)
        ir0 = np.rint(p_r).astype(int)
        iv0 = np.rint(p_v).astype(int)
        # kernels evaluated at (bin - exact position)
        d_r = (ir0[:, None] + offs[None, :]) - p_r[:, None]
        d_v = (iv0[:, None] + offs[None, :]) - p_v[:, None]
        k_r = dft_kernel(cfg.range_samples, self.n_r, self.w_r, d_r.ravel()).reshape(d_r.shape)
        k_v = dft_kernel(cfg.doppler_samples, self.n_v, self.w_v, d_v.ravel()).reshape(d_v.shape)

        idx_r = np.clip(ir0[:, None] + offs[None, :], 0, self.n_r - 1)
        idx_v = np.clip(iv0[:, None] + offs[None, :], 0, self.n_v - 1)

        m = len(ip)
        chunk = max(1, 2_000_000 // (offs.size * offs.size * max(self.n_vrx, 1)))
        for s in range(0, m, chunk):
            e = min(m, s + chunk)
            amp = ip.amplitude[s:e]
            # (chunk, nr_w, nv_w, n_vrx)
            block = (amp[:, None, None, None]
                     * k_r[s:e, :, None, None] * k_v[s:e, None, :, None]
                     * steer[s:e, None, None, :])
            ri = np.repeat(idx_r[s:e, :, None], offs.size, axis=2)
            vi = np.repeat(idx_v[s:e, None, :], offs.size, axis=1)
            ri = np.repeat(ri[:, :, :, None], self.n_vrx, axis=3)
            vi = np.repeat(vi[:, :, :, None], self.n_vrx, axis=3)
            ai = np.broadcast_to(np.arange(self.n_vrx), ri.shape)
            np.add.at(cube, (ri, vi, ai), block)
        return cube

    # ------------------------------------------------------------- 2. noise
    def add_noise(self, cube: np.ndarray) -> np.ndarray:
        """White Gaussian noise floor with the documented range/velocity scaling.

        Reference Manual -> Front-End Hardware: "The noise floor is modeled as
        white Gaussian noise (WGN) ... To account for sensor specific
        characteristics, a scaling factor can be parameterized with the help of
        a look-up table.  The logarithmic noise factors NF_R (depending on
        range) and NF_Vel (depending on velocity) are used to scale the noise
        power P_Noise."
        """
        cfg = self.cfg
        p = np.full((self.n_r, self.n_v), self.p_noise)
        if cfg.noise_scaling:
            p = p * 10.0 ** (_lut_db(cfg.noise_scaling_range, self.range_axis)[:, None] / 10.0)
            p = p * 10.0 ** (_lut_db(cfg.noise_scaling_doppler, self.vel_axis)[None, :] / 10.0)
        sigma = np.sqrt(0.5 * p)[:, :, None]
        noise = (self.rng.normal(0.0, 1.0, cube.shape)
                 + 1j * self.rng.normal(0.0, 1.0, cube.shape)) * sigma
        return cube + noise

    def noise_map(self) -> np.ndarray:
        cfg = self.cfg
        p = np.full((self.n_r, self.n_v), self.p_noise)
        if cfg.noise_scaling:
            p = p * 10.0 ** (_lut_db(cfg.noise_scaling_range, self.range_axis)[:, None] / 10.0)
            p = p * 10.0 ** (_lut_db(cfg.noise_scaling_doppler, self.vel_axis)[None, :] / 10.0)
        return p

    # --------------------------------------------- 3. Range-Doppler filtering
    def os_cfar(self, power: np.ndarray, c: CfarConfig) -> np.ndarray:
        """Ordered-Statistic CFAR on the 2-D Range-Doppler map (Eq. 526).

        "Its parameters include the number of surrounding layers n_lay around the
         cell-under-test (CUT), as well as the guarded layers n_guard which are
         excluded from evaluation.  Furthermore, a threshold T_OS-CFAR is used to
         distinguish detection from noise.  To change the false alarm rate ... a
         scaling factor S_OS-CFAR is applied to the CUT.
         ... the CUT has to be larger than the m-th value, given as
             m = T_OS-CFAR * N_lay^tot"
        """
        nr, nv = power.shape
        lay, gd = c.layers, c.guards
        scale = 10.0 ** (c.snr_db / 10.0)
        det = np.zeros_like(power, dtype=bool)
        # training-cell offsets: square ring between the guard band and n_lay
        offs = [(i, j)
                for i in range(-lay, lay + 1)
                for j in range(-lay, lay + 1)
                if max(abs(i), abs(j)) > gd]
        n_tot = len(offs)
        if n_tot == 0:
            return det
        m = int(min(max(round(c.threshold_pct / 100.0 * n_tot), 1), n_tot)) - 1
        # gather all training cells as a stack, then take the m-th order statistic
        stack = np.empty((n_tot, nr, nv))
        for k, (di, dj) in enumerate(offs):
            stack[k] = _shift2d(power, di, dj)
        ref = np.partition(stack, m, axis=0)[m]
        det = power > ref * scale
        return det

    def ca_cfar_1d(self, power: np.ndarray, c: CfarConfig) -> np.ndarray:
        """Cell-averaging CFAR on the angular spectrum (Angular-Filtering)."""
        n = power.size
        lay, gd = c.layers, c.guards
        scale = 10.0 ** (c.snr_db / 10.0)
        det = np.zeros(n, dtype=bool)
        for i in range(n):
            lo = max(0, i - lay)
            hi = min(n, i + lay + 1)
            acc, cnt = 0.0, 0
            for j in range(lo, hi):
                if abs(j - i) <= gd:
                    continue
                acc += power[j]
                cnt += 1
            if cnt == 0:
                continue
            det[i] = power[i] > (acc / cnt) * scale
        return det

    # ------------------------------------------------------ 4. peak finding
    @staticmethod
    def peak_finder_2d(power: np.ndarray, mask: np.ndarray, layers: int) -> np.ndarray:
        """"The Peak-Finder ensures that for each detection after filtering only
        a local maximum is reported as a detection.  ... the surrounding is
        spanned up via a square defined by layers." """
        if layers <= 0:
            return mask
        local_max = np.ones_like(mask)
        for di in range(-layers, layers + 1):
            for dj in range(-layers, layers + 1):
                if di == 0 and dj == 0:
                    continue
                local_max &= power >= _shift2d(power, di, dj, fill=-np.inf)
        return mask & local_max

    @staticmethod
    def peak_finder_1d(power: np.ndarray, mask: np.ndarray, layers: int) -> np.ndarray:
        if layers <= 0:
            return mask
        n = power.size
        out = mask.copy()
        for i in range(n):
            if not out[i]:
                continue
            lo = max(0, i - layers)
            hi = min(n, i + layers + 1)
            if power[i] < power[lo:hi].max():
                out[i] = False
        return out

    # ------------------------------------------------- 5. peak interpolation
    @staticmethod
    def parabolic_peak(alpha: float, beta: float, gamma: float) -> Tuple[float, float]:
        """Reference Manual Equations 527/528.

            p    = 1/2 * (alpha - gamma) / (alpha - 2 beta + gamma)
            y(p) = beta - 1/4 * p * (alpha - gamma)

        with ``beta`` the peak bin and ``alpha``/``gamma`` its neighbours.
        Returns ``(offset_in_bins, interpolated_magnitude)``.
        """
        den = alpha - 2.0 * beta + gamma
        if abs(den) < 1e-30:
            return 0.0, beta
        p = 0.5 * (alpha - gamma) / den
        p = max(-0.5, min(0.5, p))
        y = beta - 0.25 * p * (alpha - gamma)
        return p, y


def _shift2d(a: np.ndarray, di: int, dj: int, fill: float = 0.0) -> np.ndarray:
    """Shift a 2-D array by (di, dj) filling the border with ``fill``."""
    out = np.full_like(a, fill)
    nr, nv = a.shape
    si0, si1 = max(0, di), min(nr, nr + di)
    di0, di1 = max(0, -di), min(nr, nr - di)
    sj0, sj1 = max(0, dj), min(nv, nv + dj)
    dj0, dj1 = max(0, -dj), min(nv, nv - dj)
    if si1 > si0 and sj1 > sj0:
        out[di0:di1, dj0:dj1] = a[si0:si1, sj0:sj1]
    return out
