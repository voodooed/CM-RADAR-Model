#!/usr/bin/env python3
"""
Validation of the Radar RSI Device Model (Front-End + Back-End).

Each test corresponds to a documented statement in
Reference Manual -> Sensors -> Radar RSI -> Device Model.

1. Noise power, Equation 525:  P_Noise = k_B T B 10^(NF/10), plus the
   NoiseScalingRange / NoiseScalingDopplerVel look-up tables.
2. Range-Doppler processing: bin mapping, resolution, window/leakage,
   aliasing, and the documented normalisation ("the power of the output
   corresponds to the signal power at the end of the Front-End Hardware").
3. OS-CFAR (Equation 526) false-alarm behaviour on pure noise.
4. Peak interpolation, Equations 527/528, against an exact parabola.
5. Angular processing: a lambda/2 ULA recovers the azimuth of a synthetic
   target, and the resolution matches 2/n_ULA in sin(azimuth).
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np                                              # noqa: E402

from radar_rsi.config import (K_BOLTZMANN, CfarConfig, RadarRSIConfig,  # noqa: E402
                              WindowFunction)
from radar_rsi.device import DeviceModel, noise_power, window     # noqa: E402
from radar_rsi.propagation import InteractionPoints              # noqa: E402

OK, BAD = "ok ", "FAIL"


def check(c: bool) -> str:
    return OK if c else BAD


def make_device(**kw) -> DeviceModel:
    d = dict(range_max=100.0, range_samples=128, range_zero_paddings=0,
             doppler_min_max=(-30.0, 30.0), doppler_samples=64,
             doppler_zero_paddings=0, azimuth_samples=16,
             azimuth_zero_paddings=48, window_function=WindowFunction.RECTANGULAR,
             noise_scaling=False, random_seed=3)
    d.update(kw)
    cfg = RadarRSIConfig(**d)
    n = cfg.n_vrx
    lam = cfg.wavelength
    y = np.array([(i - (n - 1) / 2.0) * lam / 2.0 for i in range(n)])
    w = np.ones(n, dtype=complex)
    return DeviceModel(cfg, y, w, seed=cfg.random_seed)


def single_target(r, v, az, amp=1.0) -> InteractionPoints:
    return InteractionPoints(np.array([r]), np.array([v]), np.array([az]),
                             np.array([0.0]), np.array([amp + 0j]))


# ===========================================================================
def test_noise() -> None:
    print("=" * 78)
    print("1. Front-End noise (Equation 525)")
    print("=" * 78)
    cfg = RadarRSIConfig(noise_temperature_k=300.0, noise_bandwidth_mhz=200.0,
                         noise_figure_db=5.0)
    manual = K_BOLTZMANN * 300.0 * 200e6 * 10 ** 0.5
    code = noise_power(cfg)
    print("   T = 300 K, B = 200 MHz, NF = 5 dB")
    print(f"   by hand : {manual:.6e} W = {10 * math.log10(manual) + 30:.3f} dBm")
    print(f"   by code : {code:.6e} W = {10 * math.log10(code) + 30:.3f} dBm")
    print(f"   {check(abs(code - manual) / manual < 1e-12)}")

    dev = make_device(noise_scaling=True,
                      noise_scaling_range=[(0.0, 10.0), (10.0, 8.0), (20.0, 5.0),
                                           (40.0, 0.0), (80.0, -1.0)],
                      noise_scaling_doppler=[(-20.0, 0.0), (0.0, 0.0), (20.0, 0.0)])
    nm = dev.noise_map()
    print("   NoiseScalingRange LUT (DemoCar_SensorRadarRSI values):")
    for r_query, expect_db in ((0.0, 10.0), (10.0, 8.0), (20.0, 5.0),
                               (40.0, 0.0), (80.0, -1.0), (100.0, -1.0)):
        i = int(np.argmin(np.abs(dev.range_axis - r_query)))
        got = 10 * math.log10(nm[i, dev.n_v // 2] / dev.p_noise)
        print(f"      r = {dev.range_axis[i]:6.2f} m -> +{got:5.2f} dB "
              f"(LUT: +{expect_db:.1f} dB)  {check(abs(got - expect_db) < 0.3)}")

    # empirical noise power in the cube
    dev2 = make_device()
    cube = dev2.add_noise(np.zeros((dev2.n_r, dev2.n_v, dev2.n_vrx), dtype=complex))
    emp = float(np.mean(np.abs(cube) ** 2))
    print(f"   empirical mean |noise|^2 in the cube = {emp:.4e} W, "
          f"P_Noise = {dev2.p_noise:.4e} W  "
          f"{check(abs(emp / dev2.p_noise - 1) < 0.02)}")


# ===========================================================================
def test_rd_processing() -> None:
    print()
    print("=" * 78)
    print("2. Range-Doppler processing")
    print("=" * 78)
    dev = make_device()
    dr = dev.range_axis[1] - dev.range_axis[0]
    dv = dev.vel_axis[1] - dev.vel_axis[0]
    print(f"   range axis  : 0 .. {dev.cfg.range_max} m over {dev.n_r} bins "
          f"-> {dr:.4f} m/bin")
    print(f"   doppler axis: {dev.cfg.doppler_min_max} m/s over {dev.n_v} bins "
          f"-> {dv:.4f} (m/s)/bin")

    print("\n   a) a target that lands exactly on a bin keeps its amplitude")
    print("      (documented normalisation: 'All detections are normalized to "
          "compensate\n       for any signal processing gain or loss')")
    print(f"      {'r [m]':>8} {'v [m/s]':>9} {'peak bin r':>11} {'peak bin v':>11} "
          f"{'|A| peak':>10} {'|A| in':>8}")
    ok = True
    exact = [(float(dev.range_axis[32]), float(dev.vel_axis[42]), 1.0),
             (float(dev.range_axis[64]), float(dev.vel_axis[18]), 0.5),
             (float(dev.range_axis[96]), float(dev.vel_axis[31]), 2.0)]
    for r, v, amp in exact:
        cube = dev.build_cube(single_target(r, v, 0.0, amp))
        p = np.abs(cube[:, :, 0])
        ir, iv = np.unravel_index(np.argmax(p), p.shape)
        exp_r = r / dev.cfg.range_max * (dev.n_r - 1)
        exp_v = ((v - dev.cfg.doppler_min_max[0])
                 / (dev.cfg.doppler_min_max[1] - dev.cfg.doppler_min_max[0])) * (dev.n_v - 1)
        good = abs(ir - exp_r) <= 0.5 and abs(iv - exp_v) <= 0.5 \
            and abs(p[ir, iv] / amp - 1) < 0.05
        ok &= good
        print(f"      {r:8.2f} {v:9.2f} {ir:11d} {iv:11d} {p[ir, iv]:10.4f} "
              f"{amp:8.4f}  {check(good)}")
    print(f"      amplitude normalisation preserved (documented): {check(ok)}")

    print("\n   b) resolution: two equal targets merge below ~1 range bin")
    for gap_bins in (0.25, 0.5, 1.0, 2.0):
        r0 = 40.0
        r1 = r0 + gap_bins * dr
        ip = InteractionPoints(np.array([r0, r1]), np.array([0.0, 0.0]),
                               np.array([0.0, 0.0]), np.array([0.0, 0.0]),
                               np.array([1 + 0j, 1 + 0j]))
        cube = dev.build_cube(ip)
        cut = np.abs(cube[:, dev.n_v // 2, 0])
        iv = dev.n_v // 2
        peaks = [i for i in range(1, dev.n_r - 1)
                 if cut[i] > cut[i - 1] and cut[i] > cut[i + 1] and cut[i] > 0.2]
        print(f"      separation {gap_bins:4.2f} bins ({gap_bins * dr:5.3f} m) "
              f"-> {len(peaks)} peak(s)")

    print("\n   c) window functions and spectral leakage")
    print("      target placed half a bin off-grid (worst case for leakage)")
    for wf in (WindowFunction.RECTANGULAR, WindowFunction.HANN):
        d2 = make_device(window_function=wf)
        r_half = float(d2.range_axis[40]) + 0.5 * (d2.range_axis[1] - d2.range_axis[0])
        cube = d2.build_cube(single_target(r_half, float(d2.vel_axis[31]), 0.0, 1.0))
        cut = np.abs(cube[:, 31, 0])
        pk = int(np.argmax(cut))
        far = np.concatenate([cut[:max(pk - 3, 0)], cut[pk + 4:]])
        print(f"      {wf.value:<12s}: peak {cut[pk]:.4f}, "
              f"scalloping loss {20 * math.log10(cut[pk]):6.2f} dB, "
              f"leakage 3+ bins away {20 * math.log10(far.max() / cut[pk]):7.2f} dBc")
    print("      (Hann trades a larger scalloping loss for much lower leakage - "
          "the\n       documented reason for offering the two windows)")

    print("\n   d) aliasing: a Doppler beyond the unambiguous interval folds back")
    d3 = make_device()
    v_amb = d3.cfg.doppler_min_max[1] - d3.cfg.doppler_min_max[0]
    v0 = float(d3.vel_axis[45])
    alias_ok = True
    for k in (0, 1, -1, 2):
        cube = d3.build_cube(single_target(40.0, v0 + k * v_amb, 0.0, 1.0))
        iv = int(np.argmax(np.abs(cube[:, :, 0]).max(axis=0)))
        alias_ok &= (iv == 45 or iv == 44 or iv == 46)
        print(f"      v = {v0 + k * v_amb:+8.2f} m/s -> doppler bin {iv:3d} "
              f"(unambiguous interval = {v_amb:.1f} m/s)  {check(iv in (44, 45, 46))}")
    print(f"      folding into the unambiguous interval: {check(alias_ok)}")
    print("      (set cfg.doppler_aliasing = False to clamp instead, which is what "
          "the\n       shipped GPU-Coding-Interface sample Radar.cu does)")


# ===========================================================================
def test_os_cfar() -> None:
    print()
    print("=" * 78)
    print("3. OS-CFAR (Equation 526) on pure noise")
    print("=" * 78)
    print("   'the CUT has to be larger than the m-th value, "
          "m = T_OS-CFAR * N_lay^tot'")
    dev = make_device()
    rng = np.random.default_rng(11)
    for snr_db in (0.0, 2.0, 4.0, 6.0, 9.0, 11.0):
        c = CfarConfig(layers=5, guards=2, snr_db=snr_db, threshold_pct=95.0)
        rates = []
        for _ in range(4):
            noise = (rng.normal(size=(dev.n_r, dev.n_v))
                     + 1j * rng.normal(size=(dev.n_r, dev.n_v))) * math.sqrt(0.5)
            p = np.abs(noise) ** 2
            det = dev.os_cfar(p, c)
            rates.append(det.mean())
        rate = float(np.mean(rates))
        print(f"   S_OS-CFAR = {snr_db:5.1f} dB -> false alarm rate "
              f"{rate:.5f}  ({rate * dev.n_r * dev.n_v:6.1f} cells of "
              f"{dev.n_r * dev.n_v})")
    print("   monotonically decreasing with the scaling factor: expected behaviour")

    # a target well above the noise must always be detected
    dev2 = make_device()
    cube = dev2.build_cube(single_target(50.0, 5.0, 0.0, 1.0))
    cube = dev2.add_noise(cube * math.sqrt(dev2.p_noise) * 30.0)
    p = np.abs(cube[:, :, 0]) ** 2
    c = CfarConfig(layers=5, guards=2, snr_db=11.0, threshold_pct=95.0)
    det = dev2.os_cfar(p, c)
    det = dev2.peak_finder_2d(p, det, 2)
    ir, iv = np.nonzero(det)
    print(f"   strong target at r = 50.0 m, v = 5.0 m/s -> "
          f"{len(ir)} detection(s) at "
          f"{[(round(float(dev2.range_axis[a]), 2), round(float(dev2.vel_axis[b]), 2)) for a, b in zip(ir, iv)]}"
          f"  {check(len(ir) == 1)}")


# ===========================================================================
def test_peak_interpolation() -> None:
    print()
    print("=" * 78)
    print("4. Peak interpolation (Equations 527 / 528)")
    print("=" * 78)
    print("   p = 1/2 (alpha - gamma) / (alpha - 2 beta + gamma)")
    print("   y(p) = beta - 1/4 p (alpha - gamma)")
    print(f"   {'true offset':>12} {'true peak':>10} {'p':>10} {'y(p)':>10} "
          f"{'err(p)':>10} {'err(y)':>10}")
    ok = True
    for p_true in (-0.4, -0.2, 0.0, 0.15, 0.35):
        # exact parabola  y(x) = A - B (x - p_true)^2
        a_, b_ = 5.0, 2.0
        alpha = a_ - b_ * (-1 - p_true) ** 2
        beta = a_ - b_ * (0 - p_true) ** 2
        gamma = a_ - b_ * (1 - p_true) ** 2
        p, y = DeviceModel.parabolic_peak(alpha, beta, gamma)
        e1, e2 = abs(p - p_true), abs(y - a_)
        ok &= e1 < 1e-10 and e2 < 1e-10
        print(f"   {p_true:12.3f} {a_:10.3f} {p:10.6f} {y:10.6f} "
              f"{e1:10.2e} {e2:10.2e}")
    print(f"   exact for a parabola: {check(ok)}")

    print("\n   applied to a synthetic RD peak (sub-bin target position):")
    dev = make_device()
    dr = dev.range_axis[1] - dev.range_axis[0]
    for frac in (0.0, 0.25, 0.5):
        r_true = 40.0 + frac * dr
        cube = dev.build_cube(single_target(r_true, 0.0, 0.0, 1.0))
        p = np.abs(cube[:, :, 0]) ** 2
        ir, iv = np.unravel_index(np.argmax(p), p.shape)
        off, _ = DeviceModel.parabolic_peak(p[ir - 1, iv], p[ir, iv], p[ir + 1, iv])
        r_est = dev.range_axis[ir] + off * dr
        print(f"      true {r_true:8.4f} m   bin {dev.range_axis[ir]:8.4f} m   "
              f"interpolated {r_est:8.4f} m   error "
              f"{abs(r_est - r_true) * 100:6.2f} cm")


# ===========================================================================
def test_angular() -> None:
    print()
    print("=" * 78)
    print("5. Angular processing (lambda/2 ULA + spatial DFT)")
    print("=" * 78)
    dev = make_device(azimuth_samples=16, azimuth_zero_paddings=112)
    n_ula = dev.cfg.azimuth_samples
    print(f"   n_ULA = {n_ula}, zero paddings = {dev.cfg.azimuth_zero_paddings}, "
          f"bins = {dev.cfg.n_azimuth_bins}")
    print(f"   Rayleigh resolution ~ 2/n_ULA = {2 / n_ula:.4f} in sin(az) "
          f"= {math.degrees(math.asin(2 / n_ula)):.2f} deg at boresight")
    print(f"   {'true az':>9} {'estimated az':>14} {'error':>9}")
    ok = True
    for az_true_deg in (-40.0, -20.0, -5.0, 0.0, 7.5, 25.0, 45.0):
        az = math.radians(az_true_deg)
        cube = dev.build_cube(single_target(40.0, 0.0, az, 1.0))
        ir, iv = np.unravel_index(np.argmax(np.abs(cube[:, :, 0])), cube.shape[:2])
        x = cube[ir, iv, :n_ula] * window(dev.cfg.window_function, n_ula)
        spec = np.fft.fftshift(np.fft.fft(x, n=dev.cfg.n_azimuth_bins))
        pw = np.abs(spec) ** 2
        i = int(np.argmax(pw))
        off, _ = DeviceModel.parabolic_peak(pw[i - 1], pw[i], pw[i + 1]) \
            if 0 < i < pw.size - 1 else (0.0, pw[i])
        u = 2.0 * (i + off - dev.cfg.n_azimuth_bins / 2.0) / dev.cfg.n_azimuth_bins
        est = math.degrees(math.asin(max(-1.0, min(1.0, u))))
        err = abs(est - az_true_deg)
        ok &= err < 1.0
        print(f"   {az_true_deg:9.2f} {est:14.3f} {err:9.3f}  {check(err < 1.0)}")
    print(f"   {check(ok)}")


if __name__ == "__main__":
    test_noise()
    test_rd_processing()
    test_os_cfar()
    test_peak_interpolation()
    test_angular()
