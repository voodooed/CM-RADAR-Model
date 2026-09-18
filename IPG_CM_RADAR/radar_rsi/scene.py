"""
Scene representation for the Radar RSI ray tracer.

CarMaker feeds the RSI ray tracer with the *3D geometry* of the scene (Movie NX
/ IPGMovie instances), not with bounding boxes.  Reference Manual -> Sensors ->
Radar RSI -> Object Model Requirements:

    "The faces have to be small enough to represent the curvature and to contain
     relevant corner reflectors ... The normals have to point perpendicularly
     outward from a closed surface ... The faces have to be tagged with a
     suitable material (cf. Material Parameters)."

Materials come from ``Data/Sensor/MaterialLib`` (``FileIdent = MaterialLib 15``),
which defines per material:

    Material.<i>.Radar.Permittivity = <relative permittivity>
    Material.<i>.Radar.Scattering   = <surface roughness, degrees>

(e.g. asphalt: 4.4 / 20.0 deg, metal: 1e9 / 0.01 deg, concrete-and-steel:
1e9 / 0.5 deg).  The permittivity drives the Fresnel reflection coefficient,
the scattering value the "material dependent stochastic component on the normal
of the surface" (Reference Manual -> Material-Dependent Reflection).

Geometry here is a plain indexed triangle soup with per-triangle material and
per-body rigid motion, so that Doppler (including micro-Doppler from bodies that
move relative to each other) falls out naturally.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

from common.frames import (Pose, Vec3, vadd, vcross, vdot, vnorm, vscale,
                           vsub, vunit)
from common.infofile import InfoFile


# ===========================================================================
#  materials
# ===========================================================================
@dataclass
class RadarMaterial:
    name: str
    permittivity: float = 1e9        # Material.<i>.Radar.Permittivity
    scattering_deg: float = 0.0      # Material.<i>.Radar.Scattering


class MaterialLib:
    """Loader for CarMaker's ``Data/Sensor/MaterialLib``."""

    def __init__(self, materials: Optional[Sequence[RadarMaterial]] = None):
        self.materials: List[RadarMaterial] = list(materials or [])
        self._by_name = {m.name: m for m in self.materials}

    @classmethod
    def load(cls, path: str) -> "MaterialLib":
        inf = InfoFile.load(path)
        n = inf.int("Materials.N", 0)
        mats = []
        for i in range(1, n + 1):
            mats.append(RadarMaterial(
                name=inf.str(f"Material.{i}.Name", f"material{i}"),
                permittivity=inf.float(f"Material.{i}.Radar.Permittivity", 1e9),
                scattering_deg=inf.float(f"Material.{i}.Radar.Scattering", 0.0)))
        return cls(mats)

    @classmethod
    def default(cls) -> "MaterialLib":
        """Subset of the shipped MaterialLib, for use without a CarMaker install."""
        return cls([
            RadarMaterial("asphalt", 4.4, 20.0),
            RadarMaterial("concrete-and-steel", 1e9, 0.5),
            RadarMaterial("stone", 5.0, 5.0),
            RadarMaterial("metal", 1e9, 0.01),
            RadarMaterial("metal_scatter", 1e9, 20.0),
            RadarMaterial("glass", 6.0, 0.5),
            RadarMaterial("plastic", 3.0, 5.0),
            RadarMaterial("rubber", 4.0, 10.0),
            RadarMaterial("human", 40.0, 30.0),
        ])

    def get(self, name: str) -> RadarMaterial:
        m = self._by_name.get(name)
        if m is None:
            raise KeyError(f"material '{name}' not in MaterialLib")
        return m


# ===========================================================================
#  geometry
# ===========================================================================
@dataclass
class Triangle:
    """One face, expressed in the *body* frame of its owning ``Body``."""
    v0: Vec3
    v1: Vec3
    v2: Vec3
    material: RadarMaterial

    def normal(self) -> Vec3:
        return vunit(vcross(vsub(self.v1, self.v0), vsub(self.v2, self.v0)))


@dataclass
class Body:
    """A rigid group of triangles with its own pose and velocity.

    Separate bodies with different velocities produce the micro-Doppler the
    Reference Manual mentions ("If the 3D model comprises of multiple bodies
    with relative motion, micro-Doppler will also occur").
    """
    name: str
    pose: Pose = field(default_factory=Pose)
    velocity: Vec3 = (0.0, 0.0, 0.0)              # world frame [m/s]
    angular_velocity: Vec3 = (0.0, 0.0, 0.0)      # world frame [rad/s], about pose.t
    triangles: List[Triangle] = field(default_factory=list)

    def point_velocity(self, p_world: Sequence[float]) -> Vec3:
        """Velocity of the material point at ``p_world`` (rigid-body motion)."""
        r = vsub(p_world, self.pose.t)
        return vadd(self.velocity, vcross(self.angular_velocity, r))


class Scene:
    """Flattened triangle soup with brute-force intersection.

    A real implementation uses a BVH on the GPU; the loop here is written so
    that only ``intersect()`` needs replacing by an accelerated query when the
    algorithm is ported to C++/CARLA.
    """

    def __init__(self) -> None:
        self.bodies: List[Body] = []
        # flattened world-space triangles: (p0, p1, p2, normal, material, body_index)
        self._tris: List[Tuple[Vec3, Vec3, Vec3, Vec3, RadarMaterial, int]] = []

    def add_body(self, body: Body) -> None:
        self.bodies.append(body)

    def build(self) -> None:
        """Transform all triangles into world space (call after moving bodies)."""
        self._tris = []
        for bi, body in enumerate(self.bodies):
            for tri in body.triangles:
                p0 = body.pose.to_parent(tri.v0)
                p1 = body.pose.to_parent(tri.v1)
                p2 = body.pose.to_parent(tri.v2)
                n = vunit(vcross(vsub(p1, p0), vsub(p2, p0)))
                self._tris.append((p0, p1, p2, n, tri.material, bi))

    @property
    def n_triangles(self) -> int:
        return len(self._tris)

    # --------------------------------------------------------- intersection
    def intersect(self, origin: Sequence[float], direction: Sequence[float],
                  t_min: float = 1e-4, t_max: float = 1e9):
        """Closest ray/triangle hit.

        Returns ``(t, point, normal, material, body_index)`` or ``None``.
        Moeller-Trumbore, double sided (the normal is flipped towards the ray so
        that Fresnel always sees a front face).
        """
        best_t = t_max
        best = None
        ox, oy, oz = origin
        dx, dy, dz = direction
        for (p0, p1, p2, n, mat, bi) in self._tris:
            e1 = (p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2])
            e2 = (p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2])
            # pvec = d x e2
            px = dy * e2[2] - dz * e2[1]
            py = dz * e2[0] - dx * e2[2]
            pz = dx * e2[1] - dy * e2[0]
            det = e1[0] * px + e1[1] * py + e1[2] * pz
            if -1e-12 < det < 1e-12:
                continue
            inv = 1.0 / det
            tv = (ox - p0[0], oy - p0[1], oz - p0[2])
            u = (tv[0] * px + tv[1] * py + tv[2] * pz) * inv
            if u < 0.0 or u > 1.0:
                continue
            # qvec = tv x e1
            qx = tv[1] * e1[2] - tv[2] * e1[1]
            qy = tv[2] * e1[0] - tv[0] * e1[2]
            qz = tv[0] * e1[1] - tv[1] * e1[0]
            v = (dx * qx + dy * qy + dz * qz) * inv
            if v < 0.0 or u + v > 1.0:
                continue
            t = (e2[0] * qx + e2[1] * qy + e2[2] * qz) * inv
            if t_min < t < best_t:
                best_t = t
                hp = (ox + dx * t, oy + dy * t, oz + dz * t)
                nn = n if vdot(n, direction) < 0.0 else vscale(n, -1.0)
                best = (t, hp, nn, mat, bi)
        return best

    def is_occluded(self, a: Sequence[float], b: Sequence[float],
                    eps: float = 1e-3) -> bool:
        """True if the segment a->b is blocked (used for the line-of-sight test
        of the back-propagation path to the receiver)."""
        d = vsub(b, a)
        dist = vnorm(d)
        if dist <= eps:
            return False
        hit = self.intersect(a, vscale(d, 1.0 / dist), eps, dist - eps)
        return hit is not None


# ===========================================================================
#  primitive builders (boxes, planes) - enough for a reproducible demo scene
# ===========================================================================
def make_plane(material: RadarMaterial, size: float = 400.0,
               z: float = 0.0, name: str = "road") -> Body:
    """Large horizontal plane (road surface) centred on the origin."""
    h = size * 0.5
    v = [(-h, -h, z), (h, -h, z), (h, h, z), (-h, h, z)]
    tris = [Triangle(v[0], v[1], v[2], material),
            Triangle(v[0], v[2], v[3], material)]
    return Body(name=name, triangles=tris)


def make_box(material: RadarMaterial, length: float, width: float, height: float,
             name: str = "box", subdiv: int = 1) -> Body:
    """Axis-aligned box centred on its own origin, optionally subdivided so that
    the faces are "small enough" (Object Model Requirements)."""
    hx, hy, hz = length * 0.5, width * 0.5, height * 0.5
    corners = {
        "xn": ((-hx, -hy, -hz), (-hx, hy, -hz), (-hx, hy, hz), (-hx, -hy, hz)),
        "xp": ((hx, hy, -hz), (hx, -hy, -hz), (hx, -hy, hz), (hx, hy, hz)),
        "yn": ((-hx, -hy, -hz), (hx, -hy, -hz), (hx, -hy, hz), (-hx, -hy, hz)),
        "yp": ((hx, hy, -hz), (-hx, hy, -hz), (-hx, hy, hz), (hx, hy, hz)),
        "zn": ((-hx, hy, -hz), (hx, hy, -hz), (hx, -hy, -hz), (-hx, -hy, -hz)),
        "zp": ((-hx, -hy, hz), (hx, -hy, hz), (hx, hy, hz), (-hx, hy, hz)),
    }
    tris: List[Triangle] = []
    for quad in corners.values():
        tris.extend(_subdivide_quad(quad, material, subdiv))
    return Body(name=name, triangles=tris)


def _subdivide_quad(q, material: RadarMaterial, n: int) -> List[Triangle]:
    a, b, c, d = q
    tris: List[Triangle] = []
    for i in range(n):
        for j in range(n):
            s0, s1 = i / n, (i + 1) / n
            t0, t1 = j / n, (j + 1) / n
            p00 = _bilerp(a, b, c, d, s0, t0)
            p10 = _bilerp(a, b, c, d, s1, t0)
            p11 = _bilerp(a, b, c, d, s1, t1)
            p01 = _bilerp(a, b, c, d, s0, t1)
            tris.append(Triangle(p00, p10, p11, material))
            tris.append(Triangle(p00, p11, p01, material))
    return tris


def _bilerp(a, b, c, d, s: float, t: float) -> Vec3:
    ab = tuple(a[k] + (b[k] - a[k]) * s for k in range(3))
    dc = tuple(d[k] + (c[k] - d[k]) * s for k in range(3))
    return tuple(ab[k] + (dc[k] - ab[k]) * t for k in range(3))   # type: ignore[return-value]
