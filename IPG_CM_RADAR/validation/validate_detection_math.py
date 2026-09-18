#!/usr/bin/env python3
"""
Validation of the detection mathematics of the object-list Radar Sensor.

Checks that can be made purely against the Reference Manual:

1. ``erfcinv`` implementation and the ten values of ``ProbFalseAlarmIdx``
   (1..10 -> P_FA = 1e-1 .. 1e-10) that enter Equation 504.
2. Round-trip consistency of Equation 504 and its inversion ``prob_detect()``:
   at ``SNR = SNR_min`` the probability of detection must equal
   ``ProbDetectMin`` exactly - i.e. the documented statements
   "detection <=> SNR > SNR_min" (Eq. 502) and "ProbDetect ... resulting from
   Equation 504" describe the same threshold.
3. Radar equation (Eq. 505) against a hand-computed link budget.
4. Thermal noise (Eq. 506) against kTBF.
5. Monotonicity / sanity of the resulting detection ranges.
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common.paths import default_carmaker_dir, sensor_data_dir   # noqa: E402
from radar_object_list import detection as det                   # noqa: E402
from radar_object_list.config import (K_BOLTZMANN, RadarSensorConfig)  # noqa: E402
from radar_object_list.maps import AntennaGainMap                # noqa: E402
from radar_object_list.model import RadarSensorModel             # noqa: E402

OK = "ok "
BAD = "FAIL"


def check(cond: bool) -> str:
    return OK if cond else BAD


def main() -> None:
    print("=" * 78)
    print("1. erfc^-1 and the ProbFalseAlarmIdx table")
    print("=" * 78)
    print("   ProbFalseAlarmIdx maps 1..10 to P_FA = 1e-1 .. 1e-10 "
          "(Reference Manual,\n   Parameterization of Radar Sensor).")
    print(f"   {'idx':>4} {'P_FA':>10} {'erfc^-1(2 P_FA)':>18} "
          f"{'erfc(x) - 2 P_FA':>20}  ")
    worst = 0.0
    for idx in range(1, 11):
        pfa = 10.0 ** (-idx)
        x = det.erfcinv(2.0 * pfa)
        resid = math.erfc(x) - 2.0 * pfa
        worst = max(worst, abs(resid) / (2.0 * pfa))
        print(f"   {idx:4d} {pfa:10.0e} {x:18.10f} {resid:20.3e}  "
              f"{check(abs(resid) < 1e-12 * max(1.0, 2 * pfa))}")
    print(f"   worst relative residual: {worst:.2e}")
    print(f"   table constant ERFCINV_PFA matches: "
          f"{check(all(abs(det.ERFCINV_PFA[i] - det.erfcinv(2.0 * 10.0 ** (-i))) < 1e-12 for i in range(1, 11)))}")

    print()
    print("=" * 78)
    print("2. Equation 504 and its inversion")
    print("=" * 78)
    print("   SNR_min = 2 (erfc^-1(2 P_FA) - erfc^-1(2 P_Dmin))^2")
    print("   P_D(SNR_min) must equal P_Dmin")
    print(f"   {'P_FA':>8} {'P_Dmin':>8} {'SNR_min':>10} {'SNR_min[dB]':>12} "
          f"{'P_D(SNR_min)':>14} {'err':>10}")
    max_err = 0.0
    for pfa_idx in (3, 6, 9):
        for pdmin in (0.5, 0.9, 0.99):
            pfa = 10.0 ** (-pfa_idx)
            s = det.snr_min(pfa, pdmin)
            pd = det.prob_detect(s, pfa)
            err = abs(pd - pdmin)
            max_err = max(max_err, err)
            print(f"   {pfa:8.0e} {pdmin:8.2f} {s:10.4f} "
                  f"{10 * math.log10(s):12.3f} {pd:14.10f} {err:10.2e}")
    print(f"   max round-trip error: {max_err:.2e}  {check(max_err < 1e-9)}")

    print()
    print("=" * 78)
    print("3. Radar equation (Eq. 505) - hand check")
    print("=" * 78)
    p_dbm, g_db, rcs, r, la, latm = 14.0, 20.365, 1.0, 100.0, 0.0, 0.0
    lam = 299792458.0 / 77e9
    p_w = 10.0 ** ((p_dbm - 30.0) / 10.0)
    g = 10.0 ** (g_db / 10.0)
    manual = p_w * g * g * lam ** 2 * rcs / ((4 * math.pi) ** 3 * r ** 4)
    code = det.signal_strength(p_w, g_db, lam, rcs, r, la, latm)
    print(f"   P = {p_dbm} dBm, G = {g_db} dB, RCS = {rcs} m^2, r = {r} m")
    print(f"   by hand : {manual:.6e} W  ({10 * math.log10(manual):.3f} dBW)")
    print(f"   by code : {code:.6e} W  ({10 * math.log10(code):.3f} dBW)")
    print(f"   relative error {abs(code - manual) / manual:.2e}  "
          f"{check(abs(code - manual) / manual < 1e-12)}")
    print("   1/r^4 law:")
    prev = None
    for rr in (50.0, 100.0, 200.0):
        s = det.signal_strength(p_w, g_db, lam, rcs, rr, la, latm)
        db = 10 * math.log10(s)
        note = "" if prev is None else f"  delta = {db - prev:+.3f} dB (expect -12.041)"
        print(f"      r = {rr:6.1f} m -> {db:9.3f} dBW{note}")
        prev = db

    print()
    print("=" * 78)
    print("4. Thermal noise (Eq. 506)")
    print("=" * 78)
    t0, nf, bn = 293.15, 4.8, 25000.0
    manual = t0 * (10 ** (nf / 10)) * K_BOLTZMANN * bn
    code = det.thermal_noise(t0, nf, bn)
    print(f"   T0 = {t0} K, F = {nf} dB, B = {bn} Hz")
    print(f"   by hand : {manual:.6e} W ({10 * math.log10(manual):.3f} dBW)")
    print(f"   by code : {code:.6e} W ({10 * math.log10(code):.3f} dBW)")
    print(f"   {check(abs(code - manual) / manual < 1e-12)}")

    print()
    print("=" * 78)
    print("5. Detection ranges with the DemoCar_SensorRadar parameter set")
    print("=" * 78)
    cfg = RadarSensorConfig(random_seed=1)      # defaults = DemoCar_SensorRadar
    sensor_dir = sensor_data_dir(default_carmaker_dir())
    ipg = os.path.join(sensor_dir, "Radar_Default") if sensor_dir else ""
    antenna = AntennaGainMap.load(ipg) if ipg and os.path.isfile(ipg) else None
    model = RadarSensorModel(cfg, antenna=antenna)
    print(f"   antenna: {'Radar_Default (IPG file)' if antenna else 'internal model'}"
          f", peak {model.antenna.peak_db:.2f} dB")
    print(f"   SNR_min          = {10 * math.log10(model.snr_min):.2f} dB")
    print(f"   thermal noise    = {10 * math.log10(model.noise_thermal):.2f} dBW")
    print(f"   {'RCS [dBm^2]':>12} {'range (thermal only)':>22} {'range (with clutter)':>22}")
    cfg_nc = RadarSensorConfig(random_seed=1, clutter_enabled=False)
    model_nc = RadarSensorModel(cfg_nc, antenna=antenna)
    mono = True
    prev = 1e9
    for rcs_db in (20.0, 10.0, 0.0, -10.0, -20.0):
        s = 10.0 ** (rcs_db / 10.0)
        r_nc = model_nc.detection_range(s)
        r_c = model.detection_range(s)
        mono &= r_nc < prev
        prev = r_nc
        print(f"   {rcs_db:12.1f} {r_nc:22.1f} {r_c:22.1f}")
    print(f"   monotonically decreasing with RCS: {check(mono)}")
    print("   note: Sensor.Param.<n>.Range_max = 200 m only limits the "
          "candidate search;\n         detection itself depends only on SNR "
          "(Reference Manual, Range_max).")


if __name__ == "__main__":
    main()
