"""
Independent reference implementation of the CarMaker **Radar RSI**
(raw-signal-interface radar, ``Sensor.Param.<n>.Type = RadarRSI``).

``RadarRSIModel.calculate()`` runs one sensor cycle:

    ray pattern  ->  wave propagation (SBR)  ->  interaction points
                 ->  radar data cube (range x doppler x VRx)
                 ->  Front-End noise
                 ->  OS-CFAR + Peak-Finder + Peak-Interpolation
                 ->  VRx output   *or*   Angular processing + CA-CFAR
                 ->  point cloud (Cartesian / Spherical)

Output mirrors ``tRadarRSI`` / ``tDetPoint`` / ``tDetVRx`` of
``include/Vehicle/Sensor_RadarRSI.h``.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np

from common.frames import Pose
from radar_rsi.config import OutputType, RadarRSIConfig, RayPattern
from radar_rsi.device import (REFERENCE_IMPEDANCE, DetectionPoint, DeviceModel,
                              RadarRSIOutput, VRxDetection, window)
from radar_rsi.propagation import InteractionPoints, RayTracer, intersect_batch
from radar_rsi.raypattern import make_pattern
from radar_rsi.scene import Scene
from radar_rsi.transceiver import TransceiverConfig


class RadarRSIModel:
    def __init__(self, cfg: RadarRSIConfig,
                 transceiver: Optional[TransceiverConfig] = None):
        cfg.validate()
        self.cfg = cfg
        if transceiver is not None:
            self.tc = transceiver
        elif cfg.antenna_fname:
            self.tc = TransceiverConfig.load(cfg.antenna_fname)
        else:
            self.tc = TransceiverConfig.generate(cfg.fov_deg, (20.0, 10.0),
                                                 n_vrx=cfg.n_vrx,
                                                 wavelength=cfg.wavelength)
        vrx = self.tc.vrx
        if vrx is None or len(vrx) == 0:
            d = 0.5 * cfg.wavelength
            n = cfg.n_vrx
            y = np.array([(i - (n - 1) / 2.0) * d for i in range(n)])
            w = np.ones(n, dtype=complex)
        else:
            y, w = vrx.y, vrx.weight
        self.device = DeviceModel(cfg, y, w, seed=cfg.random_seed)
        self.tracer = RayTracer(cfg, self.tc.tx_gain, self.tc.rx_gain,
                                seed=cfg.random_seed)
        self.last_cube: Optional[np.ndarray] = None
        self.last_interaction_points: Optional[InteractionPoints] = None

    # ==================================================================== run
    def calculate(self, t: float, scene: Scene, sensor_pose: Pose,
                  sensor_vel: Sequence[float] = (0.0, 0.0, 0.0)) -> RadarRSIOutput:
        cfg = self.cfg
        scene.build()

        # ---- 1. ray pattern ------------------------------------------------
        dirs = self._ray_pattern(scene, sensor_pose)

        # ---- 2. wave propagation -------------------------------------------
        ip = self.tracer.trace(scene, sensor_pose, sensor_vel, dirs)
        self.last_interaction_points = ip

        # ---- 3. signal generation: radar cube ------------------------------
        cube = self.device.build_cube(ip)
        cube = self.device.add_noise(cube)
        self.last_cube = cube

        # ---- 4. Range-Doppler filtering ------------------------------------
        # "the Range-Doppler-Filtering is applied to the Range-Doppler-Map of a
        #  single virtual receiver"
        rd = cube[:, :, 0]
        power = np.abs(rd) ** 2
        mask = self.device.os_cfar(power, cfg.cfar_range_doppler)
        mask = self.device.peak_finder_2d(power, mask,
                                          cfg.peak_finder_range_doppler_layers)
        det_r, det_v = np.nonzero(mask)

        out = RadarRSIOutput(time_fired=t)
        if det_r.size == 0:
            return out

        # ---- 5. per-detection processing ------------------------------------
        for ir, iv in zip(det_r.tolist(), det_v.tolist()):
            rng, vel, p_rd = self._interpolate_rd(power, ir, iv)
            if cfg.output_type == OutputType.VRX:
                amp_mv = cube[ir, iv, :] * math.sqrt(2.0 * REFERENCE_IMPEDANCE) * 1e3
                out.det_vrx.append(VRxDetection(float(rng), float(vel), amp_mv))
                if len(out.det_vrx) >= cfg.n_max_detections:
                    break
                continue

            for az, el, p_pt in self._angular_processing(cube[ir, iv, :], p_rd):
                if len(out.det_points) >= cfg.n_max_detections:
                    break
                p_dbm = 10.0 * math.log10(max(p_pt, 1e-300)) + 30.0
                if cfg.output_type == OutputType.SPHERICAL:
                    coords = (float(rng), float(az), float(el))
                else:
                    coords = (float(rng * math.cos(el) * math.cos(az)),
                              float(rng * math.cos(el) * math.sin(az)),
                              float(rng * math.sin(el)))
                out.det_points.append(DetectionPoint(coords, float(p_dbm),
                                                     float(vel)))

        out.n_detections = len(out.det_points) if cfg.output_type != OutputType.VRX \
            else len(out.det_vrx)
        return out

    # ------------------------------------------------------------ ray pattern
    def _ray_pattern(self, scene: Scene, sensor_pose: Pose):
        cfg = self.cfg
        if cfg.ray_pattern != RayPattern.DYNAMIC:
            return make_pattern(cfg)

        # scene scan: one ray through each bucket centre (documented)
        from radar_rsi.propagation import SceneArrays
        sa = SceneArrays(scene)
        S = np.array(sensor_pose.t)
        R = np.array(sensor_pose.R)
        r_max = cfg.range_min_max[1]

        def scan(az: float, el: float) -> bool:
            d_local = np.array([[math.cos(el) * math.cos(az),
                                 math.cos(el) * math.sin(az),
                                 math.sin(el)]])
            d = d_local @ R.T
            hit, _, _ = intersect_batch(sa, S[None, :], d, cfg.range_min_max[0],
                                        np.array([r_max]))
            return bool(hit[0])

        boxes = self._instance_boxes(scene, sensor_pose)
        return make_pattern(cfg, scene_scan=scan, instance_boxes=boxes)

    def _instance_boxes(self, scene: Scene, sensor_pose: Pose):
        """Angular bounding boxes of the scene instances, from their bounding
        spheres (documented: "Each instance has a bounding sphere which encloses
        all its vertices ... its minimum and maximum horizontal and vertical
        angles within the field of view are used as a bounding box")."""
        out = []
        for body in scene.bodies:
            if not body.triangles:
                continue
            pts = np.array([body.pose.to_parent(v)
                            for tri in body.triangles for v in (tri.v0, tri.v1, tri.v2)])
            centre = pts.mean(axis=0)
            radius = float(np.linalg.norm(pts - centre, axis=1).max())
            if radius > 25.0:
                continue                        # documented Movie NX limit
            c_local = np.array(sensor_pose.to_child(centre))
            dist = float(np.linalg.norm(c_local))
            if dist < 1e-6 or dist > self.cfg.range_min_max[1]:
                continue
            az = math.atan2(c_local[1], c_local[0])
            el = math.atan2(c_local[2], math.hypot(c_local[0], c_local[1]))
            half = math.asin(min(radius / dist, 1.0))
            out.append((az - half, az + half, el - half, el + half))
        return out

    # --------------------------------------------------- RD peak interpolation
    def _interpolate_rd(self, power: np.ndarray, ir: int, iv: int):
        """Two-fold parabolic interpolation of the RD peak (Eq. 527-529)."""
        cfg = self.cfg
        dev = self.device
        r_axis, v_axis = dev.range_axis, dev.vel_axis
        rng = float(r_axis[ir])
        vel = float(v_axis[iv])
        p_cfar = float(power[ir, iv])
        if not cfg.peak_interpolation:
            return rng, vel, p_cfar

        p_max_r, p_max_v = p_cfar, p_cfar
        if 0 < ir < power.shape[0] - 1:
            d, p_max_r = dev.parabolic_peak(power[ir - 1, iv], power[ir, iv],
                                            power[ir + 1, iv])
            rng += d * (r_axis[1] - r_axis[0])
        if 0 < iv < power.shape[1] - 1:
            d, p_max_v = dev.parabolic_peak(power[ir, iv - 1], power[ir, iv],
                                            power[ir, iv + 1])
            vel += d * (v_axis[1] - v_axis[0])
        # Equation 529: P_peak = P_maxV * P_maxR / P_CFAR
        p_peak = p_max_v * p_max_r / max(p_cfar, 1e-300)
        return rng, vel, p_peak

    # ---------------------------------------------------- angular processing
    def _angular_processing(self, vrx_slice: np.ndarray, p_rd: float):
        """Spatial DFT over the ULA -> CA-CFAR -> Peak-Finder -> interpolation.

        Returns a list of ``(azimuth, elevation, power)`` detections.
        Elevation is 0 unless ``ProcessElevationAngle`` is enabled; the manual
        performs "the same steps as for the azimuthal angle" in that case, which
        requires a 2-D virtual array (not part of the default parameterisation).
        """
        cfg = self.cfg
        dev = self.device
        n_ula = min(cfg.azimuth_samples, vrx_slice.size)
        n_bins = cfg.n_azimuth_bins
        w = window(cfg.window_function, n_ula)
        x = vrx_slice[:n_ula] * w
        spec = np.fft.fftshift(np.fft.fft(x, n=n_bins)) / np.sum(w)
        power = np.abs(spec) ** 2

        mask = dev.ca_cfar_1d(power, cfg.cfar_azimuth)
        mask = dev.peak_finder_1d(power, mask, cfg.peak_finder_azimuth_layers)
        idx = np.nonzero(mask)[0]
        if idx.size == 0:
            # no angular detection survived - report the strongest bin so that
            # the RD detection is not silently lost
            idx = np.array([int(np.argmax(power))])

        # bin -> sin(azimuth):  d = lambda/2 ULA is unambiguous over +/-90 deg
        results = []
        total = float(power.sum()) or 1.0
        for i in idx.tolist():
            off, p_int = 0.0, float(power[i])
            if cfg.peak_interpolation and 0 < i < n_bins - 1:
                off, p_int = dev.parabolic_peak(power[i - 1], power[i], power[i + 1])
            u = 2.0 * (i + off - n_bins / 2.0) / n_bins       # sin(az)
            u = max(-1.0, min(1.0, u))
            az = math.asin(u)
            # distribute the (already normalised) RD power over the angular peaks
            results.append((float(az), 0.0, float(p_rd * p_int / total)))
        return results

    # ------------------------------------------------------------- utilities
    def range_resolution(self) -> float:
        """c / (2 B) equivalent expressed through the RD-map discretisation."""
        return self.cfg.range_max / max(self.cfg.range_samples - 1, 1)

    def velocity_resolution(self) -> float:
        lo, hi = self.cfg.doppler_min_max
        return (hi - lo) / max(self.cfg.doppler_samples - 1, 1)

    def angular_resolution(self) -> float:
        """Rayleigh resolution of a uniform lambda/2 ULA: ~2/N rad at boresight."""
        return 2.0 / max(self.cfg.azimuth_samples, 1)
