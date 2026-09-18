#!/usr/bin/env python3
"""
Example / demo for the CarMaker Radar RSI reference implementation.

Scenario:
    * flat asphalt road plane (gives the documented road multipath)
    * one metal car body 40 m ahead, closing at 5 m/s
    * one metal guardrail 6 m to the right (gives mirror / ghost targets)

The sensor is parameterised like CarMaker's shipped example vehicle
``Data/Vehicle/Examples/DemoCar_SensorRadarRSI``, with a smaller data cube and
ray count so that the demo runs in a few seconds on a CPU.

Run:
    python3 radar_rsi/example.py [--ipg <CarMaker install dir>] [--rays N]
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np                                                # noqa: E402

from common.frames import Pose                                    # noqa: E402
from radar_rsi.config import (OutputType, RadarRSIConfig, RayPattern,          # noqa: E402
                              WindowFunction)
from radar_rsi.model import RadarRSIModel                          # noqa: E402
from radar_rsi.scene import (MaterialLib, Scene, make_box,        # noqa: E402
                             make_plane)
from radar_rsi.transceiver import TransceiverConfig                # noqa: E402


def build_scene(mats: MaterialLib, car_x: float, car_v: float,
                with_guardrail: bool) -> Scene:
    scene = Scene()

    road = make_plane(mats.get("asphalt"), size=400.0, z=0.0, name="road")
    scene.add_body(road)

    car = make_box(mats.get("metal"), length=4.5, width=1.8, height=1.4,
                   name="car", subdiv=3)
    car.pose = Pose((car_x, 0.0, 0.7))
    car.velocity = (car_v, 0.0, 0.0)
    scene.add_body(car)

    if with_guardrail:
        rail = make_box(mats.get("metal"), length=120.0, width=0.1, height=0.7,
                        name="guardrail", subdiv=2)
        rail.pose = Pose((60.0, -6.0, 0.6))
        scene.add_body(rail)

    return scene


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ipg", help="CarMaker installation dir to load "
                                  "Data/Sensor/RadarRSI_Default and MaterialLib")
    ap.add_argument("--rays", type=int, default=6000)
    ap.add_argument("--no-guardrail", action="store_true")
    ap.add_argument("--vrx", action="store_true", help="use VRx output mode")
    args = ap.parse_args()

    # ---------------------------------------------------------- materials
    if args.ipg:
        mats = MaterialLib.load(os.path.join(args.ipg, "Data", "Sensor", "MaterialLib"))
        print(f"[info] MaterialLib: {len(mats.materials)} materials loaded")
    else:
        mats = MaterialLib.default()
        print("[info] using the built-in material subset")

    # ------------------------------------------------------------- sensor
    cfg = RadarRSIConfig(
        name="RARS00", pos=(4.4, 0.0, 0.5), rot_deg=(0.0, 0.0, 0.0),
        fov_deg=(40.0, 10.0),                 # narrowed for a fast demo
        range_min_max=(0.1, 120.0),
        n_rays=args.rays, ray_pattern=RayPattern.FIBONACCI,
        frequency_ghz=77.0, transmit_power_dbm=30.0,
        polarization_transmit=0.5, polarization_receive=0.5,
        noise_bandwidth_mhz=200.0, noise_figure_db=5.0, noise_temperature_k=300.0,
        noise_scaling=True,
        range_max=120.0, range_samples=128, range_zero_paddings=0,
        doppler_min_max=(-40.0, 40.0), doppler_samples=64, doppler_zero_paddings=0,
        azimuth_samples=16, azimuth_zero_paddings=48,
        process_elevation_angle=False,
        window_function=WindowFunction.HANN,
        output_type=OutputType.VRX if args.vrx else OutputType.SPHERICAL,
        n_max_detections=200,
        max_bounces=3,
        random_seed=7)
    cfg.cfar_range_doppler.layers = 5
    cfg.cfar_range_doppler.guards = 2
    cfg.cfar_range_doppler.snr_db = 11.0
    cfg.cfar_range_doppler.threshold_pct = 95.0
    cfg.cfar_azimuth.layers = 8
    cfg.cfar_azimuth.guards = 2
    cfg.cfar_azimuth.snr_db = 12.0

    if args.ipg:
        tc = TransceiverConfig.load(os.path.join(args.ipg, "Data", "Sensor",
                                                 "RadarRSI_Default"))
        print(f"[info] transceiver: Tx peak {tc.tx.peak_db:.2f} dB, "
              f"Rx peak {tc.rx.peak_db:.2f} dB, "
              f"{len(tc.vrx) if tc.vrx else 0} virtual receivers")
    else:
        tc = TransceiverConfig.generate(cfg.fov_deg, (20.0, 10.0),
                                        n_vrx=cfg.n_vrx, wavelength=cfg.wavelength)
        print("[info] transceiver generated from the internal aperture model")

    model = RadarRSIModel(cfg, transceiver=tc)

    print(f"wavelength                : {cfg.wavelength * 1000:.3f} mm")
    print(f"noise power (Eq. 525)     : {10 * math.log10(model.device.p_noise) + 30:.2f} dBm")
    print(f"range   resolution (bin)  : {model.range_resolution():.3f} m")
    print(f"doppler resolution (bin)  : {model.velocity_resolution():.3f} m/s")
    print(f"angular resolution (ULA)  : {math.degrees(model.angular_resolution()):.2f} deg")
    print(f"cube                      : {cfg.n_range_bins} x {cfg.n_doppler_bins} "
          f"x {model.device.n_vrx}")
    print()

    # ----------------------------------------------------------- scenario
    sensor_pose = Pose.from_pos_rot_deg(cfg.pos, cfg.rot_deg)
    sensor_vel = (0.0, 0.0, 0.0)          # ego at rest, target closing

    for step, (car_x, car_v) in enumerate([(45.0, -5.0), (30.0, -5.0)]):
        scene = build_scene(mats, car_x, car_v, not args.no_guardrail)
        t0 = time.time()
        out = model.calculate(step * 0.06, scene, sensor_pose, sensor_vel)
        dt = time.time() - t0
        ip = model.last_interaction_points
        r_centre = car_x - cfg.pos[0]
        r_front = (car_x - 2.25) - cfg.pos[0]     # nearest face of the 4.5 m box
        rd = np.abs(model.last_cube[:, :, 0]) ** 2
        n_cfar = int(model.device.os_cfar(rd, cfg.cfar_range_doppler).sum())
        print(f"--- step {step}: car centre at x={car_x:.1f} m "
              f"(range: front face {r_front:.2f} m, centre {r_centre:.2f} m), "
              f"v={car_v:+.1f} m/s")
        print(f"    ray tracing      : {len(ip) if ip else 0} interaction points "
              f"[{dt:.2f} s]")
        print(f"    RD map peak      : "
              f"{10 * math.log10(rd.max() / model.device.p_noise):.1f} dB over the "
              f"noise floor")
        print(f"    OS-CFAR          : {n_cfar} cells")
        print(f"    nDetections      : {out.n_detections}")
        if cfg.output_type == OutputType.VRX:
            for d in out.det_vrx[:8]:
                mag = np.abs(d.amp_vrx)
                print(f"     r={d.range_m:7.2f} m  v={d.velocity:+6.2f} m/s  "
                      f"|A| = [{mag[0]:.3e} ... {mag[-1]:.3e}] mV")
        else:
            pts = sorted(out.det_points, key=lambda p: -p.power_dbm)[:10]
            print("       range[m]   az[deg]   el[deg]     v[m/s]   power[dBm]")
            for p in pts:
                r, az, el = p.coordinates
                print(f"      {r:9.2f} {math.degrees(az):9.2f} {math.degrees(el):9.2f} "
                      f"{p.velocity:10.2f} {p.power_dbm:12.2f}")
        print()


if __name__ == "__main__":
    main()
