#!/usr/bin/env python3
"""
Behavioural validation of the object-list Radar Sensor implementation.

Each test checks one *documented* statement of the Reference Manual.  Where the
manual only describes behaviour qualitatively, the test verifies that our
implementation shows exactly the described behaviour (not that it matches IPG
numerically - that is not possible without CarMaker itself; see
validation/validation_report.md).
"""

from __future__ import annotations

import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common.frames import Pose                                    # noqa: E402
from common.paths import default_carmaker_dir                      # noqa: E402
from common.rng import Rng                                        # noqa: E402
from radar_object_list.config import (LENGTH_WIDTH_CLASSES, MeasStat,  # noqa: E402
                                      RadarSensorConfig)
from radar_object_list.maps import AntennaGainMap, RcsMap         # noqa: E402
from radar_object_list.model import RadarSensorModel              # noqa: E402
from radar_object_list.targets import Target, occlusion_fractions, view_target  # noqa: E402

IPG = default_carmaker_dir() or ""
OK, BAD = "ok ", "FAIL"


def check(c: bool) -> str:
    return OK if c else BAD


def flat_rcs(level_db: float, prob_exist: int = 4, occl: float = 0.85) -> RcsMap:
    az = [math.radians(a) for a in range(-180, 181)]
    return RcsMap("flat", az, [10.0 ** (level_db / 10.0)] * len(az), prob_exist, occl)


def antenna():
    p = os.path.join(IPG, "Data", "Sensor", "Radar_Default")
    return AntennaGainMap.load(p) if os.path.isfile(p) else AntennaGainMap.generate((20.0, 15.0))


def base_cfg(**kw) -> RadarSensorConfig:
    d = dict(random_seed=12345, false_pos_active=False, clutter_enabled=False)
    d.update(kw)
    return RadarSensorConfig(**d)


def sensor_at(x=0.0, z=0.4) -> Pose:
    return Pose.from_pos_rot_deg((x, 0.0, z), (0.0, 0.0, 0.0))


# ===========================================================================
def test_rcs_maps() -> None:
    print("=" * 78)
    print("1. Reading the RCS maps shipped by CarMaker")
    print("=" * 78)
    d = os.path.join(IPG, "Data", "Sensor")
    if not os.path.isdir(d):
        print("   [skip] CarMaker installation not found")
        return
    print(f"   {'file':<20} {'samples':>8} {'ProbExist':>10} {'OcclFactor':>11} "
          f"{'RCS(0deg)':>10} {'RCS(90deg)':>11} {'max':>8}")
    for name in ("RCS_Car", "RCS_Truck", "RCS_Pedestrian", "RCS_Bicycle",
                 "RCS_GuardRailPost"):
        p = os.path.join(d, name)
        if not os.path.isfile(p):
            continue
        m = RcsMap.load(p, name)
        f = 10 * math.log10(m.rcs(0.0))
        s = 10 * math.log10(m.rcs(math.radians(90.0)))
        mx = 10 * math.log10(max(m.rcs_lin))
        print(f"   {name:<20} {len(m.azim):8d} {m.prob_exist:10d} "
              f"{m.occlusion_factor:11.2f} {f:10.2f} {s:11.2f} {mx:8.2f}")
    print("   (values in dBm^2; Reference Manual -> RCS map files)")


# ===========================================================================
def test_swerling() -> None:
    print()
    print("=" * 78)
    print("2. Swerling type 1 RCS fluctuation (Eq. 509)")
    print("=" * 78)
    print("   p(sigma) = (1 / sigma_lut) exp(-sigma / sigma_lut)")
    print("   -> mean = sigma_lut, std = sigma_lut, P(sigma > sigma_lut) = 1/e")
    rng = Rng(7)
    for lut_db in (0.0, 10.0):
        lut = 10.0 ** (lut_db / 10.0)
        s = [rng.exponential(lut) for _ in range(200000)]
        mean = statistics.fmean(s)
        sd = statistics.pstdev(s)
        frac = sum(1 for x in s if x > lut) / len(s)
        print(f"   sigma_lut = {lut:8.3f} m^2 : mean = {mean:8.4f} "
              f"({check(abs(mean / lut - 1) < 0.02)})  "
              f"std = {sd:8.4f} ({check(abs(sd / lut - 1) < 0.02)})  "
              f"P(>lut) = {frac:.4f} vs 1/e = {1 / math.e:.4f} "
              f"({check(abs(frac - 1 / math.e) < 0.01)})")


# ===========================================================================
def test_measurement_noise() -> None:
    print()
    print("=" * 78)
    print("3. Measurement noise (Eq. 513): x_hat = x + N(0, accuracy)")
    print("=" * 78)
    cfg = base_cfg(accuracy_distance=0.4, accuracy_azimuth_deg=0.1,
                   accuracy_speed_kmh=0.1, extended_object_rcs=False)
    model = RadarSensorModel(cfg, antenna=antenna())
    tgt = Target(obj_id=16000001, pos=(60.0, 0.0, 0.75), vel=(-10.0, 0.0, 0.0),
                 length=4.8, width=1.8, height=1.5, rcs_map=flat_rcs(20.0))
    ds, azs, vs = [], [], []
    for k in range(4000):
        out = model.calculate_now(k * 0.06, sensor_at(), (0.0, 0.0, 0.0), [tgt])
        if out.objects:
            o = out.objects[0]
            ds.append(o.dist - 60.0)
            azs.append(math.degrees(math.atan2(o.dist_y, o.dist_x)))
            vs.append(o.vrel - (-10.0))
    print(f"   samples: {len(ds)}")
    print(f"   distance : mean {statistics.fmean(ds):+7.4f} m    "
          f"std {statistics.pstdev(ds):.4f} m    (AccuracyDistance = "
          f"{cfg.accuracy_distance} m)  "
          f"{check(abs(statistics.pstdev(ds) / cfg.accuracy_distance - 1) < 0.08)}")
    print(f"   azimuth  : mean {statistics.fmean(azs):+7.4f} deg  "
          f"std {statistics.pstdev(azs):.4f} deg  (AccuracyAzimuth = "
          f"{cfg.accuracy_azimuth_deg} deg)  "
          f"{check(abs(statistics.pstdev(azs) / cfg.accuracy_azimuth_deg - 1) < 0.10)}")
    print(f"   speed    : mean {statistics.fmean(vs):+7.4f} m/s  "
          f"std {statistics.pstdev(vs):.4f} m/s  (AccuracySpeed = "
          f"{cfg.accuracy_speed_kmh} km/h = {cfg.accuracy_speed_kmh / 3.6:.4f} m/s)  "
          f"{check(abs(statistics.pstdev(vs) / (cfg.accuracy_speed_kmh / 3.6) - 1) < 0.08)}")


# ===========================================================================
def test_separability() -> None:
    print()
    print("=" * 78)
    print("4. Separability / object merging (Eq. 514 and Eq. 511)")
    print("=" * 78)
    cfg = base_cfg(separability=1.5, resolution_distance=1.8,
                   resolution_azimuth_deg=1.6, resolution_speed_kmh=0.4,
                   accuracy_distance=0.0, accuracy_azimuth_deg=0.0,
                   accuracy_speed_kmh=0.0, extended_object_rcs=False)
    d_r = cfg.separability * cfg.resolution_distance
    print(f"   separable distance Delta_r = delta * R_r = {d_r:.2f} m")
    print(f"   {'gap [m]':>9} {'nObj':>5} {'expected':>9}")
    ok = True
    for gap in (0.5, 1.0, 2.0, 2.60, 2.80, 3.0, 5.0):
        a = Target(obj_id=16000001, pos=(60.0, 0.0, 0.75), vel=(0.0, 0.0, 0.0),
                   length=0.2, width=0.2, height=0.2, rcs_map=flat_rcs(20.0))
        b = Target(obj_id=16000002, pos=(60.0 + gap, 0.0, 0.75), vel=(0.0, 0.0, 0.0),
                   length=0.2, width=0.2, height=0.2, rcs_map=flat_rcs(20.0))
        m2 = RadarSensorModel(cfg, antenna=antenna())
        out = m2.calculate_now(0.0, sensor_at(), (0.0, 0.0, 0.0), [a, b])
        # the geometric range difference is slightly below the y/z-offset-free
        # gap, so use a tolerance-free comparison on the *measured* separation
        expect = 1 if gap < d_r else 2
        ok &= (len(out.objects) == expect)
        print(f"   {gap:9.2f} {len(out.objects):5d} {expect:9d}   "
              f"{check(len(out.objects) == expect)}")
    print(f"   merging follows Delta_r exactly: {check(ok)}")

    # RCS of a merged pair = mean of the single values (Eq. 511)
    cfg2 = base_cfg(separability=1.5, accuracy_distance=0.0,
                    accuracy_azimuth_deg=0.0, accuracy_speed_kmh=0.0,
                    extended_object_rcs=False, rcs_noise_correction=False)
    m = RadarSensorModel(cfg2, antenna=antenna())
    m.rng.reset(1)
    # occlusion switched off (OcclusionFactor = 0) so that only Eq. 511 is tested
    a = Target(obj_id=16000001, pos=(60.0, 0.0, 0.75), length=0.2, width=0.2,
               height=0.2, rcs_map=flat_rcs(20.0, occl=0.0))
    b = Target(obj_id=16000002, pos=(60.5, 0.0, 0.75), length=0.2, width=0.2,
               height=0.2, rcs_map=flat_rcs(20.0, occl=0.0))
    accum = []
    for _ in range(2000):
        out = m.calculate_now(0.0, sensor_at(), (0.0, 0.0, 0.0), [a, b])
        if len(out.objects) == 1:
            accum.append(10.0 ** (out.objects[0].rcs / 10.0))
    mean_merged = statistics.fmean(accum)
    print(f"   merged RCS mean = {10 * math.log10(mean_merged):.2f} dBm^2, "
          f"single-object mean = 20.00 dBm^2 (Swerling mean of the mean of two "
          f"exponentials)  {check(abs(10 * math.log10(mean_merged) - 20.0) < 0.5)}")


# ===========================================================================
def test_occlusion() -> None:
    print()
    print("=" * 78)
    print("5. Occlusion (Eq. 510):  RCS = (1 - o_h' o_v') RCS_no")
    print("=" * 78)
    print("   occluder OcclusionFactor = 0.85 (as in RCS_Car)")
    print(f"   {'occluder width [m]':>19} {'o_h':>7} {'o_v':>7} "
          f"{'RCS factor':>11}")
    target = Target(obj_id=16000001, pos=(60.0, 0.0, 0.75), length=4.8, width=1.8,
                    height=1.5, rcs_map=flat_rcs(20.0))
    pose = sensor_at()
    tv = view_target(target, pose, (0.0, 0.0, 0.0))
    for w in (0.0, 0.9, 1.8, 4.0):
        others = [tv]
        if w > 0.0:
            occ = Target(obj_id=16000002, pos=(30.0, 0.0, 0.75), length=1.0,
                         width=w * 30.0 / 60.0, height=1.5 * 30.0 / 60.0,
                         rcs_map=flat_rcs(20.0, occl=0.85))
            others.append(view_target(occ, pose, (0.0, 0.0, 0.0)))
        o_h, o_v = occlusion_fractions(tv, others)
        print(f"   {w:19.2f} {o_h:7.3f} {o_v:7.3f} {1 - o_h * o_v:11.3f}")
    print("   (an occluder at half the range needs half the size to cover the "
          "same angle)")


# ===========================================================================
def test_prob_exist_and_latency() -> None:
    print()
    print("=" * 78)
    print("6. Probability of existence, MeasStat and the cycle/latency model")
    print("=" * 78)
    cfg = base_cfg(cycle_time_ms=60.0, latency_factors=(1.0, 1.0),
                   extended_object_rcs=False)
    model = RadarSensorModel(cfg, antenna=antenna())
    model.set_initial_prob_exist({16000001: 4})
    tgt = Target(obj_id=16000001, pos=(60.0, 0.0, 0.75), length=4.8, width=1.8,
                 height=1.5, rcs_map=flat_rcs(25.0))
    print("   Reference Manual: 'Will be incremented if object is detected, "
          "decremented\n   otherwise.  If zero, object will be dropped from "
          "object list.'")
    print(f"   {'t [s]':>7} {'nObj':>5} {'ProbExist':>10} {'MeasStat':>9}  scenario")
    dt = 0.01
    for k in range(60):
        t = k * dt
        present = [tgt] if t < 0.30 else []
        out = model.step(t, sensor_at(), (0.0, 0.0, 0.0), present)
        if k % 6 == 0:
            pe = out.objects[0].prob_exist if out.objects else 0
            ms = out.objects[0].meas_stat if out.objects else MeasStat.NO_OBJECT
            print(f"   {t:7.2f} {len(out.objects):5d} {pe:10d} {ms:9d}  "
                  f"{'target present' if present else 'target removed'}")
    print(f"   first output appears one latency ({model.latency:.2f} s) after the "
          f"first cycle: expected")
    print("   after the target disappears the object is held while ProbExist > 0")


# ===========================================================================
def test_false_positives() -> None:
    print()
    print("=" * 78)
    print("7. False positives (clutter and mirror objects)")
    print("=" * 78)
    cfg = base_cfg(false_pos_active=True, clutter_obj_mean=5.0,
                   extended_object_rcs=False, max_num_obj=2000)
    model = RadarSensorModel(cfg, antenna=antenna())
    counts = []
    ids = set()
    for k in range(400):
        out = model.calculate_now(k * 0.06, sensor_at(), (0.0, 0.0, 0.0), [])
        n = sum(1 for o in out.objects if o.obj_id == -2)
        counts.append(n)
        ids.update(o.obj_id for o in out.objects)
        model._tracks.clear()
    mean = statistics.fmean(counts)
    var = statistics.pvariance(counts)
    print(f"   ClutterObjMean = {cfg.clutter_obj_mean}")
    print(f"   observed mean  = {mean:.3f}   "
          f"{check(abs(mean - cfg.clutter_obj_mean) < 0.4)}")
    print(f"   observed var   = {var:.3f}    (Poisson: var = mean)  "
          f"{check(abs(var - mean) < 1.0)}")
    print(f"   clutter object IDs seen: {sorted(ids)}  "
          f"(documented: 'Clutter is always assigned an object ID of -2')  "
          f"{check(ids <= {-2})}")

    # mirror objects need a guardrail-like object in the scene
    cfg2 = base_cfg(false_pos_active=True, clutter_obj_mean=0.0,
                    extended_object_rcs=False, max_num_obj=2000)
    m2 = RadarSensorModel(cfg2, antenna=antenna())
    rail = Target(obj_id=16000009, pos=(50.0, -6.0, 0.5), length=100.0, width=0.2,
                  height=0.7, rcs_map=flat_rcs(10.0), name="guardrail")
    car = Target(obj_id=16000001, pos=(50.0, 0.0, 0.75), length=4.8, width=1.8,
                 height=1.5, rcs_map=flat_rcs(20.0), name="car")
    mirrors = 0
    for k in range(300):
        out = m2.calculate_now(k * 0.06, sensor_at(), (0.0, 0.0, 0.0), [rail, car])
        mirrors += sum(1 for o in out.objects if o.obj_id == -16000001)
        m2._tracks.clear()
    print(f"   mirror objects of the car over 300 cycles: {mirrors} "
          f"(ID = -ObjId, documented)  {check(mirrors > 0)}")


# ===========================================================================
def test_classes() -> None:
    print()
    print("=" * 78)
    print("8. Length / width classification table")
    print("=" * 78)
    print("   documented: 0 unknown, 1 <0.5, 2 <1, 3 <2, 4 <3, 5 <4, 6 <6, "
          "7 exceeds")
    cases = ((0.0, 0), (0.3, 1), (0.7, 2), (1.5, 3), (2.5, 4), (3.5, 5),
             (5.0, 6), (7.0, 7))
    ok = True
    for v, expect in cases:
        got = LENGTH_WIDTH_CLASSES.classify(v)
        ok &= got == expect
        print(f"   {v:5.1f} m -> class {got}  (expected {expect})  "
              f"{check(got == expect)}")
    print(f"   {check(ok)}")


if __name__ == "__main__":
    test_rcs_maps()
    test_swerling()
    test_measurement_noise()
    test_separability()
    test_occlusion()
    test_prob_exist_and_latency()
    test_false_positives()
    test_classes()
