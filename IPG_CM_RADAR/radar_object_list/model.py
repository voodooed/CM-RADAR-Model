"""
Independent reference implementation of the CarMaker **Radar Sensor**
(object-list "HiFi" radar, ``Sensor.Param.<n>.Type = Radar``).

The class ``RadarSensorModel`` reproduces the pipeline documented in
Reference Manual -> Sensors -> Radar Sensor:

    candidates -> antenna gain -> RCS -> radar equation -> SNR -> threshold
               -> grouping/merging -> measurement noise -> object management
               -> false positives -> object list (+ latency)

The output mirrors ``tOutQuants`` / ``tOutQuantsGlob`` of
``include/Vehicle/Sensor_Radar.h``.

Design notes for a later C++/CARLA port
---------------------------------------
* No Python-only constructs in the algorithm: plain loops, plain floats, small
  fixed-size vectors.  ``dataclass`` is used for records only (-> ``struct``).
* All state lives in ``RadarSensorModel`` (-> class members); ``step()`` is a
  pure function of (time, ego pose, target list, state).
* Randomness goes through ``common.rng.Rng`` (-> ``std::mt19937_64``).
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from common.frames import Pose, azimuth_elevation, vnorm
from common.rng import Rng
from radar_object_list import atmosphere as atm
from radar_object_list import detection as det
from radar_object_list.config import (LENGTH_WIDTH_CLASSES, OBJ_ID_CLUTTER,
                                      DynProp, MeasStat, RadarSensorConfig)
from radar_object_list.maps import AntennaGainMap
from radar_object_list.targets import Target, TargetView, occlusion_fractions, view_target

KMH = 1.0 / 3.6


# ===========================================================================
#  output records  (mirror tOutQuants / tOutQuantsGlob of Sensor_Radar.h)
# ===========================================================================
@dataclass
class RadarObject:
    obj_id: int = -1
    dist: float = 0.0
    dist_x: float = 0.0
    dist_y: float = 0.0
    dist_z: float = 0.0
    vrel: float = 0.0
    vrel_x: float = 0.0
    vrel_y: float = 0.0
    arel_x: float = 0.0
    length: float = 0.0
    width: float = 0.0
    length_class: int = 0
    width_class: int = 0
    theta: float = 0.0                 # elevation [rad]
    rcs: float = 0.0                   # [dBm^2]
    signal_strength: float = 0.0       # [dBW]
    snr: float = 0.0                   # [dB]
    prob_detect: float = 0.0
    dyn_prop: int = DynProp.NOT_CLASSIFIED
    prob_exist: int = 0
    prob_obst: float = 0.0
    meas_stat: int = MeasStat.NO_OBJECT
    group_id: int = -1
    left_in_group: int = -1
    rel_course_angle: float = 0.0


@dataclass
class RadarGlobalInfo:
    rol_count: int = 0
    n_obj: int = 0
    relv_tgt: int = -1
    n_lanes_l: int = 0
    n_lanes_r: int = 0
    dist_to_left_border: float = 0.0
    dist_to_right_border: float = 0.0
    transmit_power: float = 0.0
    noise_bandwidth: float = 0.0
    noise_figure: float = 0.0
    separability: float = 0.0
    resolution_speed: float = 0.0
    accuracy_speed: float = 0.0
    resolution_azimuth: float = 0.0
    accuracy_azimuth: float = 0.0
    resolution_distance: float = 0.0
    accuracy_distance: float = 0.0
    range_min: float = 0.0
    range_max: float = 0.0


@dataclass
class RadarOutput:
    time: float = 0.0
    objects: List[RadarObject] = field(default_factory=list)
    info: RadarGlobalInfo = field(default_factory=RadarGlobalInfo)


# ===========================================================================
#  internal per-object tracking state
# ===========================================================================
@dataclass
class _TrackState:
    prob_exist: int = 0
    seen_before: bool = False
    last_output: Optional[RadarObject] = None


# ===========================================================================
#  merged detection candidate
# ===========================================================================
@dataclass
class _Candidate:
    ids: List[int]                       # object ids merged into this candidate
    views: List[TargetView]
    range_m: float
    azimuth: float
    elevation: float
    range_rate: float
    rcs_lin: float                       # linear RCS [m^2]
    signal_w: float
    noise_w: float
    snr: float
    prob_detect: float
    ds: Tuple[float, float, float]
    dv: Tuple[float, float, float]
    da: Tuple[float, float, float]
    length: float
    width: float
    rel_course_angle: float
    group_id: int = -1
    left_in_group: int = -1
    is_false_positive: bool = False
    fp_obj_id: int = 0


# ===========================================================================
#  the model
# ===========================================================================
class RadarSensorModel:
    """CarMaker Radar Sensor - independent reference implementation."""

    def __init__(self, cfg: RadarSensorConfig,
                 antenna: Optional[AntennaGainMap] = None,
                 road_z: float = 0.0):
        cfg.validate()
        self.cfg = cfg
        self.road_z = road_z
        if antenna is not None:
            self.antenna = antenna
        elif cfg.gain_fname:
            self.antenna = AntennaGainMap.load(cfg.gain_fname)
        else:
            self.antenna = AntennaGainMap.generate(cfg.beam_width_deg,
                                                   cfg.scan_range_deg,
                                                   cfg.antenna_eff)
        self.rng = Rng(cfg.random_seed)

        # constant parts of the link budget
        self.snr_min = det.snr_min(cfg.prob_false_alarm, cfg.prob_detect_min)
        self.noise_thermal = det.thermal_noise(cfg.temperature_k,
                                               cfg.noise_figure_db,
                                               cfg.noise_bandwidth_hz)
        self.damping_db_per_km = atm.total_damping_db_per_km(
            cfg.frequency_ghz, cfg.rain_rate_mm_h, cfg.vis_range_fog_m,
            cfg.rain_k, cfg.rain_alpha, cfg.fog_k_db_km_per_gm3,
            cfg.atm_gamma_db_km)

        # tracking / timing state
        self._tracks: Dict[int, _TrackState] = {}
        self._rol_count = 0
        self._t_last_calc = -1e30
        self._pending: List[Tuple[float, RadarOutput]] = []   # (release time, output)
        self._current = RadarOutput(info=self._make_info())

    # ------------------------------------------------------------------ API
    @property
    def cycle_time(self) -> float:
        return self.cfg.cycle_time_ms * 1e-3

    @property
    def latency(self) -> float:
        """Documented as Latency[0..1] * CycleTime; we use the lower bound
        (``[1.0 1.0]`` = "quantity update without jitter")."""
        return self.cfg.latency_factors[0] * self.cycle_time

    def step(self, t: float,
             sensor_pose_world: Pose,
             sensor_vel_world: Sequence[float],
             targets: Sequence[Target],
             ego_speed: float = 0.0,
             ego_yaw_rate: float = 0.0,
             ego_width: float = 1.8,
             road_info: Optional[Dict[str, float]] = None) -> RadarOutput:
        """Advance the sensor to simulation time ``t`` and return the currently
        visible object list.

        A new radar cycle is computed every ``CycleTime``; its result becomes
        visible ``Latency`` later (Reference Manual -> Time Delay: "The time
        delay results in an age of output data between dt_latency and
        dt_latency + T_radar").
        """
        offset = self.cfg.cycle_offset_ms * 1e-3
        if t + 1e-12 >= self._t_last_calc + self.cycle_time or self._t_last_calc < -1e29:
            if t + 1e-12 >= offset:
                self._t_last_calc = t
                out = self._calculate(t, sensor_pose_world, sensor_vel_world,
                                      targets, ego_speed, ego_yaw_rate,
                                      ego_width, road_info)
                self._pending.append((t + self.latency, out))

        while self._pending and self._pending[0][0] <= t + 1e-12:
            _, out = self._pending.pop(0)
            self._current = out
        return self._current

    def calculate_now(self, t: float, sensor_pose_world: Pose,
                      sensor_vel_world: Sequence[float],
                      targets: Sequence[Target], **kw) -> RadarOutput:
        """Run one radar cycle immediately, bypassing cycle time and latency.
        Useful for unit tests and for validating single scenes."""
        return self._calculate(t, sensor_pose_world, sensor_vel_world, targets,
                               kw.get("ego_speed", 0.0), kw.get("ego_yaw_rate", 0.0),
                               kw.get("ego_width", 1.8), kw.get("road_info"))

    # =================================================================== core
    def _calculate(self, t: float, sensor_pose_world: Pose,
                   sensor_vel_world: Sequence[float],
                   targets: Sequence[Target],
                   ego_speed: float, ego_yaw_rate: float, ego_width: float,
                   road_info: Optional[Dict[str, float]]) -> RadarOutput:
        cfg = self.cfg
        sensor_height = max(sensor_pose_world.t[2] - self.road_z, 0.0)

        # ---- 1. observation area: range gate + field of view ---------------
        views: List[TargetView] = []
        for tgt in targets:
            if not tgt.detect_mask or tgt.rcs_map is None:
                continue                       # cannot be detected (documented)
            v = view_target(tgt, sensor_pose_world, sensor_vel_world)
            if v.range_m < cfg.range_min or v.range_m > cfg.range_max:
                continue
            if not self._in_fov(v):
                continue
            views.append(v)
        views.sort(key=lambda v: v.range_m)

        # ---- 2. per-object link budget ------------------------------------
        candidates: List[_Candidate] = []
        for v in views:
            cand = self._link_budget(v, views, sensor_height)
            if cand.snr > self.snr_min:
                candidates.append(cand)

        # ---- 3. grouping / merging (separability, Eq. 514) -----------------
        candidates = self._merge_candidates(candidates)

        # ---- 4. false positives (mirror objects + clutter) -----------------
        if cfg.false_pos_active:
            candidates.extend(self._false_positives(views, sensor_height))

        # ---- 5. measurement noise + output records ------------------------
        objects: List[RadarObject] = []
        for cand in candidates:
            objects.append(self._to_output(cand, ego_speed, ego_yaw_rate, ego_width))

        # ---- 6. object management: probability of existence, list size -----
        objects = self._object_management(objects)

        # ---- 7. global information ----------------------------------------
        self._rol_count = (self._rol_count + 1) & 0xFFFFFFFF
        info = self._make_info()
        info.rol_count = self._rol_count
        info.n_obj = len(objects)
        info.relv_tgt = self._relevant_target(objects)
        if road_info:
            info.n_lanes_l = int(road_info.get("nLanesL", 0))
            info.n_lanes_r = int(road_info.get("nLanesR", 0))
            info.dist_to_left_border = road_info.get("DistToLeftBorder", 0.0)
            info.dist_to_right_border = road_info.get("DistToRightBorder", 0.0)
        return RadarOutput(time=t, objects=objects, info=info)

    # ------------------------------------------------------------ FoV gate
    def _in_fov(self, v: TargetView) -> bool:
        """Observation area test.

        The Reference Manual states that detection itself is *not* geometric
        ("The detection is explicitly not based on geometrical analyses"), and
        that FoV/range only select detection candidates.  We therefore accept an
        object as soon as *any* part of its bounding box falls inside the
        FoV wedge.
        """
        h = math.radians(self.cfg.fov_deg[0]) * 0.5
        v_ = math.radians(self.cfg.fov_deg[1]) * 0.5
        return (v.az_min <= h and v.az_max >= -h
                and v.el_min <= v_ and v.el_max >= -v_)

    # ------------------------------------------------------- link budget
    def _link_budget(self, v: TargetView, all_views: Sequence[TargetView],
                     sensor_height: float) -> _Candidate:
        cfg = self.cfg
        lam = cfg.wavelength

        gain_db = self.antenna.gain(v.azimuth, v.elevation)
        damping_db = atm.two_way_damping_db(self.damping_db_per_km, v.range_m)

        # ---- RCS -------------------------------------------------------
        rcs_map = v.target.rcs_map
        assert rcs_map is not None
        if cfg.extended_object_rcs:
            # extended-object handling: widen the incidence-azimuth interval to
            # at least the angular resolution projected onto the object
            res = math.radians(cfg.resolution_azimuth_deg)
            span = max(v.rcs_az_max - v.rcs_az_min, res)
            mid = 0.5 * (v.rcs_az_min + v.rcs_az_max)
            rcs_lut = rcs_map.rcs_mean(mid - 0.5 * span, mid + 0.5 * span)
        else:
            rcs_lut = rcs_map.rcs(v.rcs_azimuth)

        # Swerling type 1 fluctuation (Eq. 509; the manual prints the pdf
        # without the minus sign - an exponential pdf requires it)
        rcs_fluct = self.rng.exponential(rcs_lut)

        # occlusion (Eq. 510)
        o_h, o_v = occlusion_fractions(v, all_views)
        rcs_lin = (1.0 - o_h * o_v) * rcs_fluct

        # ---- signal and noise -----------------------------------------
        signal_w = det.signal_strength(cfg.transmit_power_w, gain_db, lam,
                                       rcs_lin, v.range_m,
                                       cfg.system_losses_db, damping_db)
        noise_w = self.noise_thermal
        if cfg.clutter_enabled:
            noise_w += det.clutter_noise(
                v.range_m, cfg.transmit_power_w, gain_db, lam,
                cfg.clutter_sigma0_db, math.radians(cfg.resolution_azimuth_deg),
                cfg.resolution_distance, sensor_height, cfg.system_losses_db,
                damping_db, math.radians(cfg.clutter_grazing_min_deg))

        snr = signal_w / noise_w if noise_w > 0.0 else 0.0
        pd = det.prob_detect(snr, cfg.prob_false_alarm)

        # reported RCS "corrected for the total noise": the receiver measures
        # S + N, so the RCS estimate that follows from the radar equation is
        # inflated by (1 + 1/SNR).  [INFERRED]
        rcs_reported = rcs_lin * (1.0 + 1.0 / snr) if (cfg.rcs_noise_correction and snr > 0.0) \
            else rcs_lin

        return _Candidate(ids=[v.target.obj_id], views=[v], range_m=v.range_m,
                          azimuth=v.azimuth, elevation=v.elevation,
                          range_rate=v.range_rate, rcs_lin=rcs_reported,
                          signal_w=signal_w, noise_w=noise_w, snr=snr,
                          prob_detect=pd, ds=v.ds, dv=v.dv, da=v.da,
                          length=v.target.length, width=v.target.width,
                          rel_course_angle=v.rel_course_angle)

    # ------------------------------------------------------------- merging
    def _merge_candidates(self, cands: List[_Candidate]) -> List[_Candidate]:
        """Separability (Reference Manual Eq. 514).

        "If two objects differ in all three states less than the respective
         separable value  Delta_i = delta * R_i  ... they are merged into a
         single object with a new bounding box enveloping the original ones.
         The output quantities are calculated for the new emerged object, the
         RCS corresponds to the mean of the original values."  (Eq. 511)
        """
        cfg = self.cfg
        d_r = cfg.separability * cfg.resolution_distance
        d_v = cfg.separability * cfg.resolution_speed_kmh * KMH
        d_a = cfg.separability * math.radians(cfg.resolution_azimuth_deg)

        merged: List[_Candidate] = []
        used = [False] * len(cands)
        for i, ci in enumerate(cands):
            if used[i]:
                continue
            group = [ci]
            used[i] = True
            # single linkage against the growing group
            changed = True
            while changed:
                changed = False
                for j, cj in enumerate(cands):
                    if used[j]:
                        continue
                    if any(abs(cj.range_m - g.range_m) < d_r
                           and abs(cj.range_rate - g.range_rate) < d_v
                           and abs(cj.azimuth - g.azimuth) < d_a for g in group):
                        group.append(cj)
                        used[j] = True
                        changed = True
            merged.append(self._combine(group))
        return merged

    def _combine(self, group: List[_Candidate]) -> _Candidate:
        if len(group) == 1:
            return group[0]
        n = float(len(group))
        # envelope of the merged bounding boxes in the sensor frame
        xs, ys, zs = [], [], []
        for g in group:
            for v in g.views:
                for c in v.corners_s:
                    xs.append(c[0])
                    ys.append(c[1])
                    zs.append(c[2])
        cx = 0.5 * (min(xs) + max(xs))
        cy = 0.5 * (min(ys) + max(ys))
        cz = 0.5 * (min(zs) + max(zs))
        ds = (cx, cy, cz)
        rng = math.sqrt(cx * cx + cy * cy + cz * cz)
        az, el = azimuth_elevation(ds)

        # weight the kinematics by signal power (strongest scatterer dominates)
        wsum = sum(max(g.signal_w, 0.0) for g in group) or 1.0
        dv = tuple(sum(g.dv[k] * g.signal_w for g in group) / wsum for k in range(3))
        da = tuple(sum(g.da[k] * g.signal_w for g in group) / wsum for k in range(3))
        rr = sum(g.range_rate * g.signal_w for g in group) / wsum

        rep = max(group, key=lambda g: g.signal_w)      # representative object
        # left-most object in sensor view (documented: LeftInGroup)
        left = max(group, key=lambda g: g.azimuth)

        return _Candidate(
            ids=[i for g in group for i in g.ids],
            views=[v for g in group for v in g.views],
            range_m=rng, azimuth=az, elevation=el, range_rate=rr,
            rcs_lin=sum(g.rcs_lin for g in group) / n,          # Eq. 511
            signal_w=sum(g.signal_w for g in group),
            noise_w=rep.noise_w,
            snr=sum(g.signal_w for g in group) / rep.noise_w if rep.noise_w > 0 else 0.0,
            prob_detect=det.prob_detect(
                sum(g.signal_w for g in group) / rep.noise_w if rep.noise_w > 0 else 0.0,
                self.cfg.prob_false_alarm),
            ds=ds, dv=dv, da=da,
            length=max(xs) - min(xs), width=max(ys) - min(ys),
            rel_course_angle=rep.rel_course_angle,
            group_id=rep.ids[0], left_in_group=left.ids[0])

    # ------------------------------------------------------- false positives
    def _false_positives(self, views: Sequence[TargetView],
                         sensor_height: float) -> List[_Candidate]:
        """Reference Manual -> Processing Effects -> False positives.

        Documented facts implemented here:
          * mirror objects (reflection via guardrail / tunnel wall); their
            probability of existence grows with truck-like objects and depends on
            the sensor-object and object-wall distances; object id = -ObjId;
          * clutter: number ~ Poisson(ClutterObjMean), distance and object
            dimensions ~ uniform; object id = -2.

        Everything else (the exact probability law, the wall geometry lookup) is
        **[INFERRED]**; the wall position must be supplied by the caller through
        ``Target.name == 'guardrail'`` style markers or is skipped.
        """
        cfg = self.cfg
        out: List[_Candidate] = []

        # ---- mirror objects -------------------------------------------
        walls = [v for v in views if v.target.name.lower() in ("guardrail", "wall", "tunnel")]
        for v in views:
            if v.target.name.lower() in ("guardrail", "wall", "tunnel"):
                continue
            for w in walls:
                d_wall = abs(v.ds[1] - w.ds[1])            # lateral distance object <-> wall
                if d_wall <= 0.0:
                    continue
                # documented: probability grows for trucks, falls with distances
                truck_boost = 2.0 if v.target.length > 8.0 else 1.0
                p = truck_boost * math.exp(-d_wall / 5.0) * math.exp(-v.range_m / cfg.range_max)
                if self.rng.uniform() >= min(p, 1.0):
                    continue
                # mirrored position: reflect the object across the wall plane
                mirror_y = 2.0 * w.ds[1] - v.ds[1]
                ds = (v.ds[0], mirror_y, v.ds[2])
                rng = vnorm(ds)
                az, el = azimuth_elevation(ds)
                out.append(_Candidate(
                    ids=[-v.target.obj_id], views=[], range_m=rng, azimuth=az,
                    elevation=el, range_rate=v.range_rate,
                    rcs_lin=0.5 * v.target.length * v.target.width,
                    signal_w=0.0, noise_w=self.noise_thermal, snr=self.snr_min * 1.1,
                    prob_detect=cfg.prob_detect_min,
                    ds=ds, dv=(v.dv[0], -v.dv[1], v.dv[2]), da=v.da,
                    length=v.target.length, width=v.target.width,
                    rel_course_angle=v.rel_course_angle,
                    is_false_positive=True, fp_obj_id=-v.target.obj_id))

        # ---- clutter ---------------------------------------------------
        n_clutter = self.rng.poisson(cfg.clutter_obj_mean)
        h = math.radians(cfg.fov_deg[0]) * 0.5
        vfov = math.radians(cfg.fov_deg[1]) * 0.5
        for _ in range(n_clutter):
            rng = self.rng.uniform(cfg.range_min, cfg.range_max)
            az = self.rng.uniform(-h, h)
            el = self.rng.uniform(-vfov, vfov)
            length = self.rng.uniform(0.3, 6.0)
            width = self.rng.uniform(0.3, 3.0)
            ds = (rng * math.cos(el) * math.cos(az),
                  rng * math.cos(el) * math.sin(az),
                  rng * math.sin(el))
            out.append(_Candidate(
                ids=[OBJ_ID_CLUTTER], views=[], range_m=rng, azimuth=az,
                elevation=el, range_rate=self.rng.uniform(-1.0, 1.0),
                rcs_lin=10.0 ** (self.rng.uniform(-2.0, 0.5)),
                signal_w=0.0, noise_w=self.noise_thermal, snr=self.snr_min * 1.05,
                prob_detect=cfg.prob_detect_min, ds=ds,
                dv=(0.0, 0.0, 0.0), da=(0.0, 0.0, 0.0),
                length=length, width=width, rel_course_angle=0.0,
                is_false_positive=True, fp_obj_id=OBJ_ID_CLUTTER))
        return out

    # ----------------------------------------------------- output + noise
    def _to_output(self, c: _Candidate, ego_speed: float,
                   ego_yaw_rate: float, ego_width: float) -> RadarObject:
        """Apply measurement noise (Eq. 513) and fill one ``tOutQuants`` record."""
        cfg = self.cfg
        r = c.range_m + self.rng.normal(0.0, cfg.accuracy_distance)
        r = max(r, 0.0)
        az = c.azimuth + self.rng.normal(0.0, math.radians(cfg.accuracy_azimuth_deg))
        el = c.elevation + self.rng.normal(0.0, math.radians(cfg.accuracy_azimuth_deg))
        v_r = c.range_rate + self.rng.normal(0.0, cfg.accuracy_speed_kmh * KMH)

        # measured cartesian position derived from the noisy polar measurement
        dx = r * math.cos(el) * math.cos(az)
        dy = r * math.cos(el) * math.sin(az)
        dz = r * math.sin(el)

        # measured lateral/longitudinal velocity: keep the ground-truth
        # tangential part, replace the radial part by the measured one
        ux, uy = math.cos(az), math.sin(az)
        v_true_radial = c.dv[0] * ux + c.dv[1] * uy
        dvx = c.dv[0] + (v_r - v_true_radial) * ux
        dvy = c.dv[1] + (v_r - v_true_radial) * uy

        # size measurement (documented: "enhanced by a resolution and distance
        # based clustering model with noise")  [INFERRED noise law]
        size_sigma = 0.5 * math.radians(cfg.resolution_azimuth_deg) * max(r, 1.0)
        length = max(c.length + self.rng.normal(0.0, size_sigma), 0.0)
        width = max(c.width + self.rng.normal(0.0, size_sigma), 0.0)

        obj = RadarObject()
        obj.obj_id = c.fp_obj_id if c.is_false_positive else c.ids[0]
        obj.dist = r
        obj.dist_x, obj.dist_y, obj.dist_z = dx, dy, dz
        obj.vrel = v_r
        obj.vrel_x, obj.vrel_y = dvx, dvy
        obj.arel_x = c.da[0]
        obj.length, obj.width = length, width
        obj.length_class = LENGTH_WIDTH_CLASSES.classify(length)
        obj.width_class = LENGTH_WIDTH_CLASSES.classify(width)
        obj.theta = el
        obj.rcs = det.to_db(c.rcs_lin)                 # dBm^2
        obj.signal_strength = det.to_dbw(c.signal_w)   # dBW
        obj.snr = det.to_db(c.snr)                     # dB
        obj.prob_detect = c.prob_detect
        obj.dyn_prop = self._dyn_prop(c, ego_speed)
        obj.prob_obst = self._obstacle_probability(c, ego_yaw_rate, ego_speed, ego_width)
        obj.group_id = c.group_id
        obj.left_in_group = c.left_in_group
        obj.rel_course_angle = c.rel_course_angle
        return obj

    def _dyn_prop(self, c: _Candidate, ego_speed: float) -> int:
        """Documented: "Based on velocity direction of target and sensor.  The
        target is classified as oncoming, if the velocity vectors have opposite
        directions." (stationary / stopped / moving / oncoming)."""
        if not c.views:
            return DynProp.NOT_CLASSIFIED
        v_abs = vnorm(c.views[0].target.vel)
        if v_abs < 0.1:
            # never moved vs. temporarily stopped cannot be distinguished from a
            # single frame; CarMaker keeps history.  [INFERRED simplification]
            return DynProp.STATIONARY
        # oncoming: target velocity has a component against the sensor x axis
        if c.dv[0] < 0.0 and c.ds[0] > 0.0 and c.views[0].target.vel[0] * ego_speed < 0.0:
            return DynProp.ONCOMING
        return DynProp.MOVING

    def _obstacle_probability(self, c: _Candidate, ego_yaw_rate: float,
                              ego_speed: float, ego_width: float) -> float:
        """Reference Manual -> Object Quantities -> Obstacle probability:

            P_obst = 1 - g / (1 m)

        with ``g`` the gap between the object's bounding box and the driving
        path.  The driving path is "defined by the width of the ego vehicle, an
        additional clearance of 0.3 m and the current curvature radius of the
        trajectory".
        """
        cfg = self.cfg
        half_path = 0.5 * ego_width + cfg.obstacle_extra_clearance_m
        x, y = c.ds[0], c.ds[1]
        # lateral offset of the curved driving path at longitudinal distance x
        if abs(ego_yaw_rate) > 1e-6 and abs(ego_speed) > 1e-3:
            radius = ego_speed / ego_yaw_rate
            y_path = math.copysign(1.0, radius) * (abs(radius)
                                                   - math.sqrt(max(radius * radius - x * x, 0.0)))
        else:
            y_path = 0.0
        gap = abs(y - y_path) - half_path - 0.5 * c.width
        gap = max(gap, 0.0)
        if gap >= cfg.obstacle_safety_gap_m:
            return 0.0
        return 1.0 - gap / cfg.obstacle_safety_gap_m

    # ------------------------------------------------------ object management
    def _object_management(self, objects: List[RadarObject]) -> List[RadarObject]:
        """Probability of existence, measurement status and list-size limit.

        Reference Manual -> Object Quantities -> Probability of existence:
          "Object specific initialization after first detection.  Will be
           incremented if object is detected, decremented otherwise.  If zero,
           object will be dropped from object list.  Therefore, an object can
           remain on the object list for several cycles, even if it is no longer
           detected."
        """
        cfg = self.cfg
        seen = set()
        for obj in objects:
            seen.add(obj.obj_id)
            st = self._tracks.get(obj.obj_id)
            if st is None:
                st = _TrackState()
                # initialised from the RCS map's ProbExist (documented)
                st.prob_exist = self._initial_prob_exist(obj.obj_id)
                self._tracks[obj.obj_id] = st
                obj.meas_stat = MeasStat.NEW_OBJECT
            else:
                obj.meas_stat = MeasStat.MEASURED
            if obj.prob_detect >= cfg.prob_detect_min:
                st.prob_exist = min(cfg.prob_exist_max, st.prob_exist + cfg.prob_exist_step)
            else:
                st.prob_exist = max(0, st.prob_exist - cfg.prob_exist_step)
            st.seen_before = True
            obj.prob_exist = st.prob_exist
            st.last_output = copy.copy(obj)      # snapshot; never alias the output

        # objects that were not detected in this cycle: decrement, keep while > 0
        stale: List[RadarObject] = []
        for oid, st in list(self._tracks.items()):
            if oid in seen:
                continue
            st.prob_exist = max(0, st.prob_exist - cfg.prob_exist_step)
            if st.prob_exist <= 0 or st.last_output is None:
                del self._tracks[oid]
                continue
            held = copy.copy(st.last_output)
            held.prob_exist = st.prob_exist
            held.meas_stat = MeasStat.NOT_MEASURED
            st.last_output = held
            stale.append(copy.copy(held))

        objects = objects + stale
        objects.sort(key=lambda o: o.dist)
        if len(objects) > cfg.max_num_obj:
            objects = objects[:cfg.max_num_obj]
        return objects

    def _initial_prob_exist(self, obj_id: int) -> int:
        st = self._rcs_prob_exist.get(obj_id, 0) if hasattr(self, "_rcs_prob_exist") else 0
        return st

    def set_initial_prob_exist(self, mapping: Dict[int, int]) -> None:
        """Register ``ProbExist`` values from the objects' RCS maps."""
        self._rcs_prob_exist = dict(mapping)

    def _relevant_target(self, objects: Sequence[RadarObject]) -> int:
        """``RelvTgt`` - "relevant target, dep. on Dist and ObstacleProbability"
        (``Sensor_Radar.h``).  We take the closest object with non-zero obstacle
        probability.  [INFERRED tie-breaking]"""
        best = -1
        best_dist = float("inf")
        for i, o in enumerate(objects):
            if o.prob_obst > 0.0 and o.dist < best_dist:
                best, best_dist = i, o.dist
        return best

    # ------------------------------------------------------------- helpers
    def _make_info(self) -> RadarGlobalInfo:
        cfg = self.cfg
        return RadarGlobalInfo(
            transmit_power=cfg.transmit_power_dbm,
            noise_bandwidth=cfg.noise_bandwidth_hz,
            noise_figure=cfg.noise_figure_db,
            separability=cfg.separability,
            resolution_speed=cfg.resolution_speed_kmh,
            accuracy_speed=cfg.accuracy_speed_kmh,
            resolution_azimuth=cfg.resolution_azimuth_deg,
            accuracy_azimuth=cfg.accuracy_azimuth_deg,
            resolution_distance=cfg.resolution_distance,
            accuracy_distance=cfg.accuracy_distance,
            range_min=cfg.range_min,
            range_max=cfg.range_max)

    # --------------------------------------------------------- diagnostics
    def detection_range(self, rcs_m2: float, azimuth: float = 0.0,
                        elevation: float = 0.0, tol: float = 1e-3) -> float:
        """Range at which a point target of the given RCS reaches ``SNR_min``.

        Solves the radar equation for r; useful for validating the link budget
        and for comparing against CarMaker's ``Range_max`` parameterisation.
        """
        cfg = self.cfg
        lo, hi = 0.1, 10000.0
        gain_db = self.antenna.gain(azimuth, elevation)

        def snr_at(r: float) -> float:
            damp = atm.two_way_damping_db(self.damping_db_per_km, r)
            s = det.signal_strength(cfg.transmit_power_w, gain_db, cfg.wavelength,
                                    rcs_m2, r, cfg.system_losses_db, damp)
            n = self.noise_thermal
            if cfg.clutter_enabled:
                n += det.clutter_noise(r, cfg.transmit_power_w, gain_db, cfg.wavelength,
                                       cfg.clutter_sigma0_db,
                                       math.radians(cfg.resolution_azimuth_deg),
                                       cfg.resolution_distance, 0.4,
                                       cfg.system_losses_db, damp,
                                       math.radians(cfg.clutter_grazing_min_deg))
            return s / n

        if snr_at(lo) < self.snr_min:
            return 0.0
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            if snr_at(mid) >= self.snr_min:
                lo = mid
            else:
                hi = mid
        return lo
