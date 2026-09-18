#!/usr/bin/env python3
"""
Example / demo for the CarMaker Radar Sensor reference implementation.

Scenario (a CarMaker-style "ACC" test run):

    ego  : 25 m/s, front-bumper radar at x = 4.2 m, z = 0.4 m
    lead : car, same lane, 60 m ahead, 22 m/s
    truck: 90 m ahead, right lane (+/- 3.5 m), 20 m/s
    ped  : pedestrian on the shoulder, 40 m ahead, 4 m lateral

The sensor is parameterised exactly like CarMaker's shipped example vehicle
``Data/Vehicle/Examples/DemoCar_SensorRadar``.

Run:
    python3 radar_object_list/example.py [--ipg <CarMaker install dir>]
"""

from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common.frames import Pose                                   # noqa: E402
from radar_object_list.config import RadarSensorConfig           # noqa: E402
from radar_object_list.maps import AntennaGainMap, RcsMap        # noqa: E402
from radar_object_list.model import RadarSensorModel             # noqa: E402
from radar_object_list.targets import Target                     # noqa: E402


def synthetic_rcs_map(name: str, level_db: float, prob_exist: int,
                      occlusion: float) -> RcsMap:
    """Fallback RCS map when the CarMaker installation is not available:
    isotropic in azimuth, with front/rear lobes like the shipped maps."""
    azim = [math.radians(a) for a in range(-180, 181)]
    lin = []
    for a in azim:
        # +8 dB towards the broadside faces, as in RCS_Car
        lobe = 8.0 * abs(math.sin(a))
        lin.append(10.0 ** ((level_db + lobe) / 10.0))
    return RcsMap(name, azim, lin, prob_exist, occlusion)


def load_maps(ipg_dir: str | None):
    if ipg_dir:
        d = os.path.join(ipg_dir, "Data", "Sensor")
        return (AntennaGainMap.load(os.path.join(d, "Radar_Default")),
                RcsMap.load(os.path.join(d, "RCS_Car"), "RCS_Car"),
                RcsMap.load(os.path.join(d, "RCS_Truck"), "RCS_Truck"),
                RcsMap.load(os.path.join(d, "RCS_Pedestrian"), "RCS_Pedestrian"))
    print("[info] no --ipg given: using the internal antenna model and "
          "synthetic RCS maps")
    return (AntennaGainMap.generate((20.0, 15.0)),
            synthetic_rcs_map("RCS_Car", 10.0, 4, 0.85),
            synthetic_rcs_map("RCS_Truck", 20.0, 5, 0.95),
            synthetic_rcs_map("RCS_Pedestrian", -8.0, 2, 0.30))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ipg", help="CarMaker installation dir "
                                  "(…/carmaker/win64-15.1) to load the original "
                                  "Data/Sensor maps")
    ap.add_argument("--steps", type=int, default=40)
    args = ap.parse_args()

    antenna, rcs_car, rcs_truck, rcs_ped = load_maps(args.ipg)

    # --- sensor: DemoCar_SensorRadar parameters -------------------------
    cfg = RadarSensorConfig(
        name="RA00", pos=(4.2, 0.0, 0.4), rot_deg=(0.0, 0.0, 0.0),
        cycle_time_ms=60.0, latency_factors=(1.0, 1.0),
        fov_deg=(40.0, 30.0), range_min=0.2, range_max=200.0, max_num_obj=200,
        frequency_ghz=77.0, transmit_power_dbm=14.0, system_losses_db=0.0,
        noise_bandwidth_hz=25000.0, noise_figure_db=4.8,
        prob_detect_min=0.5, prob_false_alarm_idx=6,
        accuracy_distance=0.4, accuracy_azimuth_deg=0.1, accuracy_speed_kmh=0.1,
        resolution_distance=1.8, resolution_azimuth_deg=1.6,
        resolution_speed_kmh=0.4, separability=1.5,
        false_pos_active=False, clutter_obj_mean=5.0,
        random_seed=42)

    model = RadarSensorModel(cfg, antenna=antenna, road_z=0.0)
    model.set_initial_prob_exist({16000001: rcs_car.prob_exist,
                                  16000002: rcs_truck.prob_exist,
                                  16000003: rcs_ped.prob_exist})

    print(f"wavelength                 : {cfg.wavelength * 1000:.3f} mm")
    print(f"antenna peak gain          : {antenna.peak_db:.3f} dB")
    print(f"P_FA                       : {cfg.prob_false_alarm:.0e}")
    print(f"SNR_min (Eq. 504)          : {model.snr_min:.4f} "
          f"({10 * math.log10(model.snr_min):.2f} dB)")
    print(f"thermal noise (Eq. 506)    : {10 * math.log10(model.noise_thermal):.2f} dBW")
    print(f"two-way damping            : {model.damping_db_per_km:.3f} dB/km")
    print(f"detection range, 10 dBm^2  : {model.detection_range(10.0):.1f} m")
    print(f"detection range,  0 dBm^2  : {model.detection_range(1.0):.1f} m")
    print(f"detection range,-10 dBm^2  : {model.detection_range(0.1):.1f} m")
    print()

    # --- scenario -------------------------------------------------------
    ego_speed = 25.0
    dt = 0.01
    for k in range(args.steps):
        t = k * dt
        ego_x = ego_speed * t
        sensor_pose = Pose.from_pos_rot_deg((ego_x + cfg.pos[0], cfg.pos[1], cfg.pos[2]),
                                            cfg.rot_deg)
        sensor_vel = (ego_speed, 0.0, 0.0)

        targets = [
            Target(obj_id=16000001, pos=(ego_x + 64.2, 0.0, 0.75),
                   vel=(22.0, 0.0, 0.0), yaw=0.0,
                   length=4.8, width=1.8, height=1.5,
                   rcs_map=rcs_car, name="lead_car"),
            Target(obj_id=16000002, pos=(ego_x + 94.2, -3.5, 1.6),
                   vel=(20.0, 0.0, 0.0), yaw=0.0,
                   length=12.0, width=2.5, height=3.2,
                   rcs_map=rcs_truck, name="truck"),
            Target(obj_id=16000003, pos=(ego_x + 44.2, 4.0, 0.9),
                   vel=(0.0, 1.2, 0.0), yaw=math.pi / 2,
                   length=0.5, width=0.6, height=1.8,
                   rcs_map=rcs_ped, name="pedestrian"),
        ]

        out = model.step(t, sensor_pose, sensor_vel, targets,
                         ego_speed=ego_speed, ego_yaw_rate=0.0, ego_width=1.8)

        if k in (0, 6, 12, 39):
            print(f"t = {t:5.2f} s   RolCount={out.info.rol_count}  "
                  f"nObj={out.info.n_obj}  RelvTgt={out.info.relv_tgt}")
            print("   ObjId      Dist   DistX   DistY    Vrel     RCS      SNR   "
                  "Pd    PoE  Cls(L/W) MeasStat  ProbObst")
            for o in out.objects:
                print(f"  {o.obj_id:>9d} {o.dist:8.2f} {o.dist_x:7.2f} {o.dist_y:7.2f} "
                      f"{o.vrel:7.2f} {o.rcs:7.2f} {o.snr:8.2f} {o.prob_detect:5.3f} "
                      f"{o.prob_exist:4d}   {o.length_class}/{o.width_class}    "
                      f"{o.meas_stat:5d}   {o.prob_obst:8.3f}")
            print()


if __name__ == "__main__":
    main()
