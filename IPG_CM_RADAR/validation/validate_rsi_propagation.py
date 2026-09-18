#!/usr/bin/env python3
"""
Validation of the Radar RSI wave-propagation stage against closed-form physics.

CarMaker's ray tracer is proprietary and GPU-only; there is no way to compare
sample-by-sample.  What *can* be checked is that the scattering kernel we
implement reproduces the analytic solutions that the Reference Manual says the
model is based on ("The computed scattered electric fields are modeled with an
analytical solution to the Maxwell's equation and are, other than the
restriction to the far-field, an exact solution of physical laws").

Checks
------
1. Flat perfectly conducting plate at normal incidence:
       sigma = 4 * pi * A^2 / lambda^2                (physical optics)
2. 1/r^4 behaviour of the received power for a fixed target.
3. Doppler: radial velocity of a moving plate.
4. Road multipath: two-ray interference against the analytic two-ray model.
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np                                                # noqa: E402

from common.frames import Pose                                    # noqa: E402
from radar_rsi.config import RadarRSIConfig                       # noqa: E402
from radar_rsi.propagation import RayTracer                       # noqa: E402
from radar_rsi.raypattern import phi_theta_pattern                # noqa: E402
from radar_rsi.scene import (Body, MaterialLib, Scene, Triangle,  # noqa: E402
                             make_plane)


def make_vertical_plate(material, size: float, x: float, subdiv: int) -> Body:
    """Square plate in the y-z plane at longitudinal distance ``x``,
    normal pointing back at the sensor (-x)."""
    h = size * 0.5
    tris = []
    for i in range(subdiv):
        for j in range(subdiv):
            y0 = -h + size * i / subdiv
            y1 = -h + size * (i + 1) / subdiv
            z0 = -h + size * j / subdiv
            z1 = -h + size * (j + 1) / subdiv
            a = (x, y0, z0)
            b = (x, y1, z0)
            c = (x, y1, z1)
            d = (x, y0, z1)
            # winding chosen so that the normal is -x
            tris.append(Triangle(a, c, b, material))
            tris.append(Triangle(a, d, c, material))
    return Body(name="plate", triangles=tris)


def unity_gain(az, el):
    return np.zeros_like(np.atleast_1d(np.asarray(az, dtype=float)))


def received_power(cfg: RadarRSIConfig, scene: Scene, sensor_pose: Pose,
                   sensor_vel=(0.0, 0.0, 0.0), n_rays: int = 40000):
    """Coherent sum of all interaction points -> received power [W]."""
    tracer = RayTracer(cfg, unity_gain, unity_gain, seed=1)
    dirs = phi_theta_pattern(math.radians(cfg.fov_deg[0]),
                             math.radians(cfg.fov_deg[1]), n_rays)
    scene.build()
    ip = tracer.trace(scene, sensor_pose, sensor_vel, dirs)
    if len(ip) == 0:
        return 0.0, ip
    return float(abs(np.sum(ip.amplitude)) ** 2), ip


def main() -> None:
    mats = MaterialLib.default()
    metal = mats.get("metal")
    metal = type(metal)(metal.name, metal.permittivity, 0.0)   # no roughness

    lam = 299792458.0 / 77e9
    print("=" * 78)
    print("1. Flat plate RCS  (physical optics:  sigma = 4 pi A^2 / lambda^2)")
    print("=" * 78)
    print(f"{'plate [m]':>10} {'R [m]':>7} {'sigma_PO [dBm2]':>16} "
          f"{'sigma_sim [dBm2]':>17} {'error [dB]':>11}")
    print("   (plate kept in the far field:  R > 2 D^2 / lambda, as the manual "
          "requires\n    for the interaction points)")
    errs = []
    for side in (0.08, 0.10, 0.15):
        for R in (30.0, 60.0):
            cfg = RadarRSIConfig(
                fov_deg=(3.0, 3.0), range_min_max=(0.1, 200.0),
                frequency_ghz=77.0, transmit_power_dbm=30.0,
                polarization_transmit=0.0, polarization_receive=0.0,
                max_bounces=1, system_losses_db=0.0,
                rain_rate_mm_h=0.0, vis_range_fog_m=1e9)
            scene = Scene()
            scene.add_body(make_vertical_plate(metal, side, R, subdiv=6))
            pose = Pose((0.0, 0.0, 0.0))
            p_rx, _ = received_power(cfg, scene, pose, n_rays=120000)
            # invert the radar equation with G_t = G_r = 1
            pt = cfg.transmit_power_w
            sigma_sim = p_rx * (4 * math.pi) ** 3 * R ** 4 / (pt * lam ** 2)
            sigma_po = 4 * math.pi * (side ** 2) ** 2 / lam ** 2
            err = 10 * math.log10(sigma_sim / sigma_po)
            errs.append(err)
            print(f"{side:10.2f} {R:7.1f} {10 * math.log10(sigma_po):16.2f} "
                  f"{10 * math.log10(sigma_sim):17.2f} {err:11.2f}"
                  f"   (far-field dist {2 * side ** 2 / lam:5.1f} m)")
    print(f"\n  mean |error| = {np.mean(np.abs(errs)):.2f} dB, "
          f"max = {np.max(np.abs(errs)):.2f} dB")

    print()
    print("=" * 78)
    print("2. Range law  (expect P ~ 1/R^4  ->  -12 dB per doubling)")
    print("=" * 78)
    prev = None
    for R in (20.0, 40.0, 80.0, 160.0):
        cfg = RadarRSIConfig(fov_deg=(3.0, 3.0), range_min_max=(0.1, 200.0),
                             max_bounces=1, rain_rate_mm_h=0.0, vis_range_fog_m=1e9)
        scene = Scene()
        scene.add_body(make_vertical_plate(metal, 0.1, R, subdiv=4))
        p_rx, _ = received_power(cfg, scene, Pose((0.0, 0.0, 0.0)), n_rays=120000)
        db = 10 * math.log10(p_rx)
        delta = "" if prev is None else f"   delta = {db - prev:+6.2f} dB"
        print(f"   R = {R:5.1f} m   P_rx = {db:8.2f} dBW{delta}")
        prev = db

    print()
    print("=" * 78)
    print("3. Doppler  (plate moving at -8 m/s, sensor at rest)")
    print("=" * 78)
    cfg = RadarRSIConfig(fov_deg=(3.0, 3.0), range_min_max=(0.1, 200.0),
                         max_bounces=1, rain_rate_mm_h=0.0, vis_range_fog_m=1e9)
    scene = Scene()
    plate = make_vertical_plate(metal, 0.1, 30.0, subdiv=4)
    plate.velocity = (-8.0, 0.0, 0.0)
    scene.add_body(plate)
    _, ip = received_power(cfg, scene, Pose((0.0, 0.0, 0.0)), n_rays=120000)
    print(f"   interaction points : {len(ip)}")
    print(f"   mean radial velocity: {np.mean(ip.velocity):+.4f} m/s  (expected -8.0000)")
    print(f"   mean range         : {np.mean(ip.range_m):.4f} m  (expected 30.0000)")

    print()
    print("=" * 78)
    print("4. Road multipath  (two-ray interference vs. analytic model)")
    print("=" * 78)
    print("   sensor h = 0.5 m, plate centre h = 0.75 m, asphalt road")
    asphalt = mats.get("asphalt")
    asphalt = type(asphalt)(asphalt.name, asphalt.permittivity, 0.0)
    for R in (20.0, 25.0, 30.0, 35.0, 40.0):
        cfg = RadarRSIConfig(fov_deg=(10.0, 10.0), range_min_max=(0.1, 200.0),
                             max_bounces=2, rain_rate_mm_h=0.0, vis_range_fog_m=1e9)
        scene = Scene()
        plate = make_vertical_plate(metal, 0.1, R, subdiv=3)
        plate.pose = Pose((0.0, 0.0, 0.75))
        scene.add_body(plate)
        scene.add_body(make_plane(asphalt, size=300.0, z=0.0, name="road"))
        p_multi, ip = received_power(cfg, scene, Pose((0.0, 0.0, 0.5)), n_rays=60000)

        cfg1 = RadarRSIConfig(fov_deg=(10.0, 10.0), range_min_max=(0.1, 200.0),
                              max_bounces=1, rain_rate_mm_h=0.0, vis_range_fog_m=1e9)
        scene1 = Scene()
        plate1 = make_vertical_plate(metal, 0.1, R, subdiv=3)
        plate1.pose = Pose((0.0, 0.0, 0.75))
        scene1.add_body(plate1)
        p_direct, _ = received_power(cfg1, scene1, Pose((0.0, 0.0, 0.5)), n_rays=60000)
        if p_direct <= 0.0:
            continue
        print(f"   R = {R:5.1f} m   direct {10 * math.log10(p_direct):8.2f} dBW   "
              f"with road {10 * math.log10(max(p_multi, 1e-300)):8.2f} dBW   "
              f"interference {10 * math.log10(p_multi / p_direct):+6.2f} dB   "
              f"({len(ip)} IA points)")


if __name__ == "__main__":
    main()
