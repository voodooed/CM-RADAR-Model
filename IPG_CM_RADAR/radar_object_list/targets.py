"""
Target (traffic object) representation and the geometric quantities the
CarMaker Radar Sensor derives from it.

Object model
------------
CarMaker represents every detectable object by an oriented **bounding box**
(``Basics.Dimension = l w h`` of the Traffic Object Template Info File) plus an
**RCS map** (``RCSMap.FName``).  The radar model uses:

  * the bounding-box **centre** (``BBC``) as the reference point for range,
    velocity and angles - Reference Manual -> Object Distance:
    "The object's bounding box center BBC is selected as reference point";
  * the eight bounding-box **corners** for the azimuth/elevation extent and for
    occlusion - ``Sensor_Radar_protected.h`` documents ``tPTraffic`` with the 8
    box vertices and the "most left, right, top, bottom vertex" pointers plus
    ``azimuth_min`` / ``azimuth_max`` ("is needed for Range resolution") and
    ``OccPhi`` / ``OccTheta`` ("part of occlusion horizontal / vertical");
  * the **azimuth of incidence** in the object body frame as the RCS-map
    argument - Reference Manual -> Radar Cross Section -> Direction of Incidence.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

from common.frames import (Pose, Vec3, azimuth_elevation, bbox_corners_local,
                           vdot, vnorm, vsub)
from radar_object_list.maps import RcsMap


@dataclass
class Target:
    """One detectable object in the sensor surrounding.

    All kinematics are given in the *world* frame (Fr0); the model transforms
    them into the sensor frame itself.
    """
    obj_id: int
    pos: Tuple[float, float, float]              # bounding-box centre, world [m]
    vel: Tuple[float, float, float] = (0.0, 0.0, 0.0)   # world [m/s]
    acc: Tuple[float, float, float] = (0.0, 0.0, 0.0)   # world [m/s^2]
    yaw: float = 0.0                             # heading about z [rad]
    pitch: float = 0.0
    roll: float = 0.0
    length: float = 4.8                          # Basics.Dimension l [m]
    width: float = 1.8                           # Basics.Dimension w [m]
    height: float = 1.5                          # Basics.Dimension h [m]
    rcs_map: Optional[RcsMap] = None             # RCSMap.FName; None -> undetectable
    detect_mask: bool = True                     # Traffic.<Id>.DetectMask[Sensor]
    name: str = ""

    # ------------------------------------------------------------------ pose
    def pose(self) -> Pose:
        from common.frames import rot_zyx
        return Pose(self.pos, rot_zyx(self.roll, self.pitch, self.yaw))

    def corners_world(self) -> List[Vec3]:
        p = self.pose()
        return [p.to_parent(c) for c in
                bbox_corners_local(self.length, self.width, self.height)]


@dataclass
class TargetView:
    """Geometry of one target as seen from the sensor (sensor frame)."""
    target: Target
    ds: Vec3                 # BBC position in sensor frame [m]
    dv: Vec3                 # relative velocity in sensor frame [m/s]
    da: Vec3                 # relative acceleration in sensor frame [m/s^2]
    range_m: float
    range_rate: float        # radial (approach) velocity, + = receding
    azimuth: float           # [rad]
    elevation: float         # [rad]
    az_min: float            # angular extent of the box, sensor frame [rad]
    az_max: float
    el_min: float
    el_max: float
    rcs_azimuth: float       # incidence azimuth in the *object* frame [rad]
    rcs_az_min: float        # incidence-azimuth span covered by the box [rad]
    rcs_az_max: float
    rel_course_angle: float  # direction of the relative velocity in sensor frame
    corners_s: List[Vec3] = field(default_factory=list)


def view_target(target: Target, sensor_pose_world: Pose,
                sensor_vel_world: Sequence[float]) -> TargetView:
    """Express one target in the sensor frame and derive all angular quantities."""
    ds = sensor_pose_world.to_child(target.pos)
    dv = sensor_pose_world.dir_to_child(vsub(target.vel, sensor_vel_world))
    da = sensor_pose_world.dir_to_child(target.acc)

    rng = vnorm(ds)
    az, el = azimuth_elevation(ds)
    # radial velocity: projection of the relative velocity on the line of sight
    rr = vdot(dv, ds) / rng if rng > 1e-9 else 0.0

    corners_s = [sensor_pose_world.to_child(c) for c in target.corners_world()]
    azs, els = [], []
    for c in corners_s:
        a, e = azimuth_elevation(c)
        azs.append(a)
        els.append(e)
    # unwrap the azimuths around the centre azimuth so that a box straddling
    # +/-pi does not produce a bogus 2*pi extent
    azs = [_unwrap_to(a, az) for a in azs]

    # ---- incidence azimuth in the object body frame (Reference Manual:
    #      "measured with respect to the object's body fixed reference frame in
    #      the bounding box center BBC with the x-axis along the object's
    #      longitudinal axis and is parametrized by azimuth phi")
    obj_pose = target.pose()
    sensor_in_obj = obj_pose.to_child(sensor_pose_world.t)
    rcs_az, _ = azimuth_elevation(sensor_in_obj)
    # span of incidence azimuth covered by the box corners (used for the
    # extended-object RCS averaging)
    rcs_az_min, rcs_az_max = _incidence_azimuth_span(target, sensor_pose_world.t, rcs_az)

    # relative course angle: direction of the (relative) velocity in sensor frame
    rca = math.atan2(dv[1], dv[0]) if (abs(dv[0]) + abs(dv[1])) > 1e-9 else 0.0

    return TargetView(target=target, ds=ds, dv=dv, da=da, range_m=rng,
                      range_rate=rr, azimuth=az, elevation=el,
                      az_min=min(azs), az_max=max(azs),
                      el_min=min(els), el_max=max(els),
                      rcs_azimuth=rcs_az, rcs_az_min=rcs_az_min,
                      rcs_az_max=rcs_az_max, rel_course_angle=rca,
                      corners_s=corners_s)


def _unwrap_to(a: float, ref: float) -> float:
    while a - ref > math.pi:
        a -= 2.0 * math.pi
    while a - ref < -math.pi:
        a += 2.0 * math.pi
    return a


def _incidence_azimuth_span(target: Target, sensor_pos_world: Sequence[float],
                            centre_az: float) -> Tuple[float, float]:
    """Range of RCS-map arguments spanned by the object's own extent.

    The sensor sees the object over a finite angular sector; expressed in the
    object frame this corresponds to a range of incidence azimuths.  Computed
    from the four bottom bounding-box corners.  [INFERRED - the documentation
    states that such an interval is used but not how it is constructed.]
    """
    obj_pose = target.pose()
    s_obj = obj_pose.to_child(sensor_pos_world)
    lo = hi = None
    for sx, sy in ((-0.5, -0.5), (-0.5, 0.5), (0.5, 0.5), (0.5, -0.5)):
        cx, cy = sx * target.length, sy * target.width
        a = math.atan2(s_obj[1] - cy, s_obj[0] - cx)
        a = _unwrap_to(a, centre_az)
        lo = a if lo is None else min(lo, a)
        hi = a if hi is None else max(hi, a)
    return lo, hi   # type: ignore[return-value]


# ===========================================================================
#  Occlusion
# ===========================================================================
def occlusion_fractions(view: TargetView, others: Sequence[TargetView]) -> Tuple[float, float]:
    """Relative horizontal and vertical occlusion ``(o_h, o_v)``, each in [0, 1].

    Reference Manual -> Radar Cross Section -> Occlusion:
      "Partial occlusion is considered based on the objects' bounding boxes.
       Using the bounding box corners relative horizontal o_h as well as
       vertical occlusion o_v (ranging from 0 to 1) is evaluated conservatively.
       Effects of multiple objects occluding the target and object-specific
       transparency properties are also considered."

    Implementation [INFERRED in its details]:
      * an occluder counts only if it is *closer* to the sensor than the target;
      * the angular intervals [az_min, az_max] and [el_min, el_max] of occluder
        and target are intersected;
      * the covered fraction is weighted with the occluder's ``OcclusionFactor``
        (RCS map file parameter: "how much the object in the bounding box
        obstructs the view of objects behind", 0 = transparent, 1 = opaque);
      * multiple occluders are combined "conservatively" as the maximum, not the
        sum, so that o_h, o_v stay in [0, 1].
    """
    o_h = 0.0
    o_v = 0.0
    span_az = max(view.az_max - view.az_min, 1e-9)
    span_el = max(view.el_max - view.el_min, 1e-9)
    for other in others:
        if other.target.obj_id == view.target.obj_id:
            continue
        if other.range_m >= view.range_m:
            continue
        factor = other.target.rcs_map.occlusion_factor if other.target.rcs_map else 1.0
        if factor <= 0.0:
            continue
        oaz_lo = _unwrap_to(other.az_min, view.azimuth)
        oaz_hi = _unwrap_to(other.az_max, view.azimuth)
        cov_az = min(view.az_max, oaz_hi) - max(view.az_min, oaz_lo)
        cov_el = min(view.el_max, other.el_max) - max(view.el_min, other.el_min)
        if cov_az <= 0.0 or cov_el <= 0.0:
            continue
        o_h = max(o_h, factor * min(1.0, cov_az / span_az))
        o_v = max(o_v, factor * min(1.0, cov_el / span_el))
    return o_h, o_v
