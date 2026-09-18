"""
Coordinate frames, rigid-body transforms and bounding boxes.

Conventions follow CarMaker (Reference Manual -> Sensors -> Sensor Assembly and
the figure "Direction by azimuth and elevation"):

  * Vehicle / sensor frames are right handed:  x forward, y left, z up.
  * ``Sensor.<n>.rot = rx ry rz`` [deg] is applied in **z-y-x order**, i.e.
    R = Rz(rz) * Ry(ry) * Rx(rx), and describes the orientation of the sensor
    (transceiver) frame with respect to the mounting frame.
  * The sensor boresight is the **+x axis** of the sensor frame.
  * Azimuth  phi   = atan2(d_y, d_x)                (0 on boresight, + to the left)
  * Elevation theta = asin(d_z / |d|)               (0 in the x-y plane, + up)
  * Direction cosines used by the aperture-antenna model (Reference Manual
    Equation 508) are
        u_y = cos(theta) * sin(phi)   (== sin(Theta)cos(Phi) in IPG notation)
        u_z = sin(theta)              (== sin(Theta)sin(Phi) in IPG notation)
    where IPG's Theta is the polar angle away from boresight and Phi the roll
    angle about the boresight.  Both parameterisations describe the same unit
    vector; see docs in radar_object_list/explanation.md.

Everything here is written with plain tuples/lists of floats so it maps 1:1 to
C++ (``std::array<double,3>``, ``double[3][3]``).
"""

from __future__ import annotations

import math
from typing import List, Sequence, Tuple

Vec3 = Tuple[float, float, float]
Mat3 = Tuple[Vec3, Vec3, Vec3]

# --------------------------------------------------------------------- basics


def vadd(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vsub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vscale(a: Sequence[float], s: float) -> Vec3:
    return (a[0] * s, a[1] * s, a[2] * s)


def vdot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vcross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def vnorm(a: Sequence[float]) -> float:
    return math.sqrt(vdot(a, a))


def vunit(a: Sequence[float]) -> Vec3:
    n = vnorm(a)
    if n <= 0.0:
        return (0.0, 0.0, 0.0)
    return (a[0] / n, a[1] / n, a[2] / n)


# ----------------------------------------------------------------- rotations


def rot_x(a: float) -> Mat3:
    c, s = math.cos(a), math.sin(a)
    return ((1.0, 0.0, 0.0), (0.0, c, -s), (0.0, s, c))


def rot_y(a: float) -> Mat3:
    c, s = math.cos(a), math.sin(a)
    return ((c, 0.0, s), (0.0, 1.0, 0.0), (-s, 0.0, c))


def rot_z(a: float) -> Mat3:
    c, s = math.cos(a), math.sin(a)
    return ((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0))


def mat_mul(a: Mat3, b: Mat3) -> Mat3:
    return tuple(  # type: ignore[return-value]
        tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )


def mat_vec(m: Mat3, v: Sequence[float]) -> Vec3:
    return (m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2])


def mat_t_vec(m: Mat3, v: Sequence[float]) -> Vec3:
    """m^T * v (i.e. rotate from parent frame into the child frame)."""
    return (m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
            m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
            m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2])


def mat_transpose(m: Mat3) -> Mat3:
    return ((m[0][0], m[1][0], m[2][0]),
            (m[0][1], m[1][1], m[2][1]),
            (m[0][2], m[1][2], m[2][2]))


def rot_zyx(rx: float, ry: float, rz: float) -> Mat3:
    """CarMaker ``rot`` convention: R = Rz(rz) * Ry(ry) * Rx(rx), angles in rad."""
    return mat_mul(rot_z(rz), mat_mul(rot_y(ry), rot_x(rx)))


# ------------------------------------------------------------------- spherical


def azimuth_elevation(d: Sequence[float]) -> Tuple[float, float]:
    """(azimuth, elevation) [rad] of a direction expressed in the sensor frame."""
    r_xy = math.hypot(d[0], d[1])
    az = math.atan2(d[1], d[0])
    el = math.atan2(d[2], r_xy)
    return az, el


def direction_cosines(az: float, el: float) -> Tuple[float, float]:
    """Aperture direction cosines (u_y, u_z) used by Reference Manual Eq. 508."""
    return math.cos(el) * math.sin(az), math.sin(el)


def unit_from_az_el(az: float, el: float) -> Vec3:
    return (math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el))


# ------------------------------------------------------------------ rigid pose


class Pose:
    """Rigid transform child -> parent:  x_parent = R * x_child + t."""

    __slots__ = ("t", "R")

    def __init__(self, t: Sequence[float] = (0.0, 0.0, 0.0),
                 R: Mat3 = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))):
        self.t: Vec3 = (float(t[0]), float(t[1]), float(t[2]))
        self.R: Mat3 = R

    @classmethod
    def from_pos_rot_deg(cls, pos: Sequence[float], rot_deg: Sequence[float]) -> "Pose":
        """CarMaker ``Sensor.<n>.pos`` / ``Sensor.<n>.rot`` (degrees, z-y-x)."""
        return cls(pos, rot_zyx(math.radians(rot_deg[0]),
                               math.radians(rot_deg[1]),
                               math.radians(rot_deg[2])))

    def to_parent(self, p: Sequence[float]) -> Vec3:
        return vadd(mat_vec(self.R, p), self.t)

    def to_child(self, p: Sequence[float]) -> Vec3:
        return mat_t_vec(self.R, vsub(p, self.t))

    def dir_to_child(self, d: Sequence[float]) -> Vec3:
        return mat_t_vec(self.R, d)

    def dir_to_parent(self, d: Sequence[float]) -> Vec3:
        return mat_vec(self.R, d)

    def compose(self, child: "Pose") -> "Pose":
        """self ∘ child  (child expressed in self's child frame)."""
        return Pose(self.to_parent(child.t), mat_mul(self.R, child.R))


# --------------------------------------------------------------- bounding box

# Vertex order matches Sensor_Radar_protected.h (ePntType):
#   0 BBR back-bottom-right, 1 BBL, 2 FBL, 3 FBR,
#   4 BTR back-top-right,    5 BTL, 6 FTL, 7 FTR
BBOX_SIGNS: Tuple[Tuple[float, float, float], ...] = (
    (-0.5, -0.5, -0.5),   # 0 back  bottom right
    (-0.5, +0.5, -0.5),   # 1 back  bottom left
    (+0.5, +0.5, -0.5),   # 2 front bottom left
    (+0.5, -0.5, -0.5),   # 3 front bottom right
    (-0.5, -0.5, +0.5),   # 4 back  top   right
    (-0.5, +0.5, +0.5),   # 5 back  top   left
    (+0.5, +0.5, +0.5),   # 6 front top   left
    (+0.5, -0.5, +0.5),   # 7 front top   right
)


def bbox_corners_local(length: float, width: float, height: float) -> List[Vec3]:
    """8 corners of an axis-aligned box centred at the bounding-box centre,
    expressed in the object body frame (x forward, y left, z up)."""
    return [(sx * length, sy * width, sz * height) for (sx, sy, sz) in BBOX_SIGNS]
