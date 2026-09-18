"""
Ray patterns of the Radar RSI.

Reference Manual -> Sensors -> Radar RSI -> Ray Pattern:

  * **Phi-Theta (static)** - "the numbers of used horizontal and vertical rays
    result from the chosen field of view and the condition of equidistant
    angles".
  * **Fibonacci (static)** - "the rays are distributed using the spherical
    fibonacci lattice.  The Fibonacci pattern reduces the regularity of the
    pattern.  This helps to reduce exaggerated interference effects that could
    occur for certain alignments between the sensor and planes in the 3D
    environment."
  * **Dynamic** - the FoV is split into buckets of at most 1 deg side length; a
    scene scan (one ray through each bucket centre) and/or an instance check
    marks buckets as occupied; the requested number of rays is then distributed
    "as evenly as possible" over the occupied buckets only.  "At no time during
    the simulation this number is larger [than nRays]."

All directions are returned as (azimuth, elevation) pairs in the sensor frame.
"""

from __future__ import annotations

import math
from typing import Callable, List, Optional, Sequence, Tuple

from radar_rsi.config import OccupancyCheckMode, RadarRSIConfig, RayPattern

AzEl = Tuple[float, float]

GOLDEN_ANGLE = math.pi * (3.0 - math.sqrt(5.0))


# ---------------------------------------------------------------- static ----
def phi_theta_pattern(fov_h_rad: float, fov_v_rad: float, n_rays: int) -> List[AzEl]:
    """Equidistant angular grid; the actual number of rays can be smaller than
    ``n_rays`` because of the "condition of equidistant angles"."""
    if n_rays <= 0:
        return []
    aspect = fov_h_rad / max(fov_v_rad, 1e-9)
    n_v = max(1, int(math.floor(math.sqrt(n_rays / max(aspect, 1e-9)))))
    n_h = max(1, int(math.floor(n_rays / n_v)))
    out: List[AzEl] = []
    for j in range(n_v):
        el = -0.5 * fov_v_rad + fov_v_rad * (j + 0.5) / n_v
        for i in range(n_h):
            az = -0.5 * fov_h_rad + fov_h_rad * (i + 0.5) / n_h
            out.append((az, el))
    return out


def fibonacci_pattern(fov_h_rad: float, fov_v_rad: float, n_rays: int) -> List[AzEl]:
    """Spherical Fibonacci lattice restricted to the field of view.

    The lattice is generated on the unit sphere cap that covers the FoV and
    mapped onto the (azimuth, elevation) rectangle; this keeps the low-
    discrepancy property that motivates the pattern while respecting the FoV
    rectangle CarMaker uses.  [Mapping detail is INFERRED; the manual only names
    the "spherical fibonacci lattice".]
    """
    if n_rays <= 0:
        return []
    out: List[AzEl] = []
    n = n_rays
    for i in range(n):
        # low-discrepancy pair (u, v) in [0,1)^2
        u = (i + 0.5) / n
        v = (i * GOLDEN_ANGLE / (2.0 * math.pi)) % 1.0
        # equal-solid-angle mapping in elevation, uniform in azimuth
        el = math.asin(math.sin(0.5 * fov_v_rad) * (2.0 * u - 1.0))
        az = fov_h_rad * (v - 0.5)
        out.append((az, el))
    return out


# --------------------------------------------------------------- dynamic ----
def dynamic_pattern(cfg: RadarRSIConfig,
                    scene_scan: Optional[Callable[[float, float], bool]] = None,
                    instance_boxes: Optional[Sequence[Tuple[float, float, float, float]]] = None
                    ) -> List[AzEl]:
    """Dynamic ray pattern (bucket occupancy grid).

    ``scene_scan(az, el)``  -> True if a ray through that bucket centre hits
                               anything within range.
    ``instance_boxes``      -> angular bounding boxes ``(az_min, az_max,
                               el_min, el_max)`` of the instances in the scene.
    """
    fov_h = math.radians(cfg.fov_deg[0])
    fov_v = math.radians(cfg.fov_deg[1])
    n_h = cfg.n_horizontal_buckets if cfg.n_horizontal_buckets > 0 \
        else int(math.ceil(cfg.fov_deg[0]))
    n_v = cfg.n_vertical_buckets if cfg.n_vertical_buckets > 0 \
        else int(math.ceil(cfg.fov_deg[1]))
    n_h, n_v = max(1, n_h), max(1, n_v)
    dh, dv = fov_h / n_h, fov_v / n_v

    occupied = [[False] * n_v for _ in range(n_h)]
    mode = cfg.occupancy_check_mode

    if mode in (OccupancyCheckMode.INSTANCES_AND_SCENE_SCAN,
                OccupancyCheckMode.ONLY_SCENE_SCAN) and scene_scan is not None:
        for i in range(n_h):
            az = -0.5 * fov_h + (i + 0.5) * dh
            for j in range(n_v):
                el = -0.5 * fov_v + (j + 0.5) * dv
                if scene_scan(az, el):
                    occupied[i][j] = True

    if mode in (OccupancyCheckMode.INSTANCES_AND_SCENE_SCAN,
                OccupancyCheckMode.ONLY_INSTANCES) and instance_boxes:
        for (a0, a1, e0, e1) in instance_boxes:
            i0 = max(0, int(math.floor((a0 + 0.5 * fov_h) / dh)))
            i1 = min(n_h - 1, int(math.floor((a1 + 0.5 * fov_h) / dh)))
            j0 = max(0, int(math.floor((e0 + 0.5 * fov_v) / dv)))
            j1 = min(n_v - 1, int(math.floor((e1 + 0.5 * fov_v) / dv)))
            for i in range(i0, i1 + 1):
                for j in range(j0, j1 + 1):
                    occupied[i][j] = True

    cells = [(i, j) for i in range(n_h) for j in range(n_v) if occupied[i][j]]
    if not cells:
        # "If no bucket is classified as occupied, the rays will be distributed
        #  within the total field of view."
        return phi_theta_pattern(fov_h, fov_v, cfg.n_rays)

    per_bucket = cfg.n_rays // len(cells)
    if cfg.max_rays_per_bucket > 0:
        per_bucket = min(per_bucket, cfg.max_rays_per_bucket)
    per_bucket = max(per_bucket, 1)
    k = max(1, int(math.floor(math.sqrt(per_bucket))))    # k x k sub-grid per bucket

    out: List[AzEl] = []
    for (i, j) in cells:
        az0 = -0.5 * fov_h + i * dh
        el0 = -0.5 * fov_v + j * dv
        for a in range(k):
            for b in range(k):
                if len(out) >= cfg.n_rays:      # never exceed nRays (documented)
                    return out
                out.append((az0 + dh * (a + 0.5) / k, el0 + dv * (b + 0.5) / k))
    return out


def make_pattern(cfg: RadarRSIConfig,
                 scene_scan: Optional[Callable[[float, float], bool]] = None,
                 instance_boxes: Optional[Sequence[Tuple[float, float, float, float]]] = None
                 ) -> List[AzEl]:
    fov_h = math.radians(cfg.fov_deg[0])
    fov_v = math.radians(cfg.fov_deg[1])
    if cfg.ray_pattern == RayPattern.FIBONACCI:
        return fibonacci_pattern(fov_h, fov_v, cfg.n_rays)
    if cfg.ray_pattern == RayPattern.PHI_THETA:
        return phi_theta_pattern(fov_h, fov_v, cfg.n_rays)
    return dynamic_pattern(cfg, scene_scan, instance_boxes)
