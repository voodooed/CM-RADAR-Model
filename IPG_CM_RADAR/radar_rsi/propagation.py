"""
Wave-propagation model of the Radar RSI (shooting-and-bouncing rays).

Reference Manual -> Sensors -> Radar RSI -> Wave Propagation Model and its
sub-sections:

  * *3D Geometry Illumination* - "The electromagnetic wave is discretized by
    rays, the interaction with the scenario is determined for each ray
    separately.  Each ray incorporates information about its polarization as
    well as its complex-valued amplitude.  To ensure accurate simulation of
    interference effects a fully complex-valued calculation is used".
  * *Material-Dependent Reflection* - "The reflected ray is modeled via Fresnel
    reflection incorporating the permittivity of the hit material ... Damping of
    the electric field amplitude ... depends on frequency, rain rate and visual
    range in fog.  Damping due to free-space propagation is considered as well.
    In order to account for scattering effects of rough surfaces, a material
    dependent stochastic component on the normal of the surface is introduced."
  * *Multipath Propagation* - "Each individual interaction is tracked in
    outgoing direction, for further tracing, and on paths towards the receiver,
    for reception.  The latter ... consider both Line-of-Sight (LOS), if
    available, and back-propagation of the incoming ray to the receiver."
  * *Doppler Shift* - "The Radar RSI models the Doppler shift incorporating all
    reflections relative to the transmitter."
  * *Scattered Electromagnetic Fields* - "each individual initial ray produces
    one or more interaction points (IA-points).  For each IA-point, the effect
    on the electromagnetic field at receiver is calculated ... Each IA-point has
    to be in the far-field".

What is documented vs. inferred
------------------------------
Documented : ray discretisation, complex amplitudes, polarisation, Fresnel with
             material permittivity, stochastic normal perturbation driven by the
             material, free-space + atmospheric damping, LOS *and*
             back-propagated reception paths, Doppler over the whole path,
             far-field validity, IA-point output (range, radial velocity,
             direction, complex field).
Inferred   : the concrete scattering kernel.  We use the standard
             physical-optics ray-tube formulation

                 E_s = j k (L^2 dOmega) / (2 pi) * Gamma * Lambda(q)
                       * E_inc * exp(-j k r2) / r2

             with ``Lambda(q) = sinc(q_t1 a/2) sinc(q_t2 a/2)`` the PO aperture
             factor of the illuminated patch.  This choice reproduces the exact
             physical-optics RCS of a flat plate, ``sigma = 4 pi A^2 / lambda^2``
             (verified in validation/validate_rsi_propagation.py), which is the
             standard analytic solution the manual refers to.

The output of this module is the list of interaction points that feeds the
Device Model, matching ``tUserRadarInteractionPoint`` of
``Examples/GPUCodingInterface/.../GPUCodingInterface.h``:
``{float range; float vel; float2 direction; float2 electricField;}``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np

from common.frames import Pose
from radar_object_list import atmosphere as atm
from radar_rsi.config import RadarRSIConfig
from radar_rsi.scene import Scene


# ===========================================================================
#  flattened scene (numpy) - the ray loop is data parallel, exactly as it is
#  on the GPU inside CarMaker; a C++ port maps this onto one kernel per bounce.
# ===========================================================================
class SceneArrays:
    def __init__(self, scene: Scene):
        tris = scene._tris                                   # noqa: SLF001
        n = len(tris)
        self.n = n
        self.p0 = np.array([t[0] for t in tris], dtype=np.float64).reshape(n, 3)
        self.p1 = np.array([t[1] for t in tris], dtype=np.float64).reshape(n, 3)
        self.p2 = np.array([t[2] for t in tris], dtype=np.float64).reshape(n, 3)
        self.nrm = np.array([t[3] for t in tris], dtype=np.float64).reshape(n, 3)
        self.permittivity = np.array([t[4].permittivity for t in tris])
        self.scattering = np.radians(np.array([t[4].scattering_deg for t in tris]))
        self.body = np.array([t[5] for t in tris], dtype=np.int64)
        self.e1 = self.p1 - self.p0
        self.e2 = self.p2 - self.p0
        # rigid-body motion of every triangle's owner
        self.body_vel = np.array([b.velocity for b in scene.bodies], dtype=np.float64) \
            if scene.bodies else np.zeros((1, 3))
        self.body_omega = np.array([b.angular_velocity for b in scene.bodies],
                                   dtype=np.float64) if scene.bodies else np.zeros((1, 3))
        self.body_org = np.array([b.pose.t for b in scene.bodies], dtype=np.float64) \
            if scene.bodies else np.zeros((1, 3))


def intersect_batch(sa: SceneArrays, origins: np.ndarray, dirs: np.ndarray,
                    t_min: float, t_max: np.ndarray):
    """Moeller-Trumbore for M rays against all triangles.

    Returns ``(hit_mask, t, tri_index)``.
    """
    m = origins.shape[0]
    best_t = np.full(m, np.inf)
    best_i = np.full(m, -1, dtype=np.int64)
    for i in range(sa.n):
        e1 = sa.e1[i]
        e2 = sa.e2[i]
        pvec = np.cross(dirs, e2)                    # (m,3)
        det = pvec @ e1
        ok = np.abs(det) > 1e-12
        if not ok.any():
            continue
        inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
        tvec = origins - sa.p0[i]
        u = np.einsum("ij,ij->i", tvec, pvec) * inv
        ok &= (u >= 0.0) & (u <= 1.0)
        qvec = np.cross(tvec, e1)
        v = np.einsum("ij,ij->i", dirs, qvec) * inv
        ok &= (v >= 0.0) & (u + v <= 1.0)
        t = (qvec @ e2) * inv
        ok &= (t > t_min) & (t < t_max) & (t < best_t)
        best_t = np.where(ok, t, best_t)
        best_i = np.where(ok, i, best_i)
    return best_i >= 0, best_t, best_i


# ===========================================================================
#  interaction points
# ===========================================================================
@dataclass
class InteractionPoints:
    """Output of the wave-propagation stage - one entry per (ray, reception path).

    Mirrors ``tUserRadarInteractionPoint``:
        range [m], vel [m/s] (radial), direction (azimuth, elevation) [rad],
        electricField (complex voltage amplitude at the receiver, sqrt(W)).
    """
    range_m: np.ndarray
    velocity: np.ndarray
    azimuth: np.ndarray
    elevation: np.ndarray
    amplitude: np.ndarray            # complex

    def __len__(self) -> int:
        return int(self.range_m.size)

    @staticmethod
    def empty() -> "InteractionPoints":
        z = np.zeros(0)
        return InteractionPoints(z, z, z, z, np.zeros(0, dtype=complex))

    def concat(self, other: "InteractionPoints") -> "InteractionPoints":
        return InteractionPoints(
            np.concatenate([self.range_m, other.range_m]),
            np.concatenate([self.velocity, other.velocity]),
            np.concatenate([self.azimuth, other.azimuth]),
            np.concatenate([self.elevation, other.elevation]),
            np.concatenate([self.amplitude, other.amplitude]))


# ===========================================================================
#  polarisation helpers (Reference Manual Eq. 521-524)
# ===========================================================================
def polarization_vector(az: float | np.ndarray, el: float | np.ndarray,
                        weight: float) -> np.ndarray:
    """Linear polarisation vector for a direction (azimuth, elevation).

    Eq. 521/522:  A_TE = sqrt(w_p),  A_TM = sqrt(1 - w_p)
    Eq. 523   :  E_TM = (-sin(theta) cos(phi), -sin(theta) sin(phi), cos(theta))
    Eq. 524   :  E_TE = (-sin(phi), cos(phi), 0)

    NOTE: the manual prints Eq. 524 with ``theta`` instead of ``phi``; with
    ``theta`` the two vectors would not be orthogonal, so we use ``phi``
    (the standard azimuthal unit vector).  This is an inferred correction of an
    apparent typo, flagged in the documentation.
    """
    az = np.asarray(az, dtype=np.float64)
    el = np.asarray(el, dtype=np.float64)
    a_te = math.sqrt(max(weight, 0.0))
    a_tm = math.sqrt(max(1.0 - weight, 0.0))
    e_tm = np.stack([-np.sin(el) * np.cos(az), -np.sin(el) * np.sin(az),
                     np.cos(el) * np.ones_like(az)], axis=-1)
    e_te = np.stack([-np.sin(az), np.cos(az), np.zeros_like(az)], axis=-1)
    v = a_te * e_te + a_tm * e_tm
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.maximum(n, 1e-15)


def fresnel(permittivity: np.ndarray, cos_i: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Fresnel reflection coefficients (r_s, r_p) for an air/dielectric boundary.

    ``permittivity`` is the relative permittivity from ``MaterialLib``
    (``Material.<i>.Radar.Permittivity``); the huge values used for metal
    (1e9) drive both coefficients to +/-1 as they should.
    """
    eps = np.asarray(permittivity, dtype=np.complex128)
    ci = np.asarray(cos_i, dtype=np.complex128)
    st2 = 1.0 - ci * ci
    ct = np.sqrt(np.maximum(np.real(eps - st2), 0.0) + 1j * np.imag(eps - st2)) / np.sqrt(eps)
    r_s = (ci - np.sqrt(eps) * ct) / (ci + np.sqrt(eps) * ct)
    r_p = (np.sqrt(eps) * ci - ct) / (np.sqrt(eps) * ci + ct)
    return r_s, r_p


# ===========================================================================
#  the tracer
# ===========================================================================
class RayTracer:
    def __init__(self, cfg: RadarRSIConfig, tx_gain, rx_gain, seed: int = 1):
        self.cfg = cfg
        self.tx_gain = tx_gain          # callable(az, el) -> dB   (array capable)
        self.rx_gain = rx_gain
        self.rng = np.random.default_rng(seed)
        self.k = 2.0 * math.pi / cfg.wavelength
        # one-way field attenuation coefficient [1/m] from the dB/km damping
        gamma_db_km = atm.total_damping_db_per_km(
            cfg.frequency_ghz, cfg.rain_rate_mm_h, cfg.vis_range_fog_m,
            0.9, 0.74, 4.0, None)
        self.alpha_np_per_m = gamma_db_km / 1000.0 / 8.685889638

    # -------------------------------------------------------------- helpers
    def _atten(self, path_len: np.ndarray) -> np.ndarray:
        return np.exp(-self.alpha_np_per_m * path_len)

    def _point_velocity(self, sa: SceneArrays, tri_idx: np.ndarray,
                        points: np.ndarray) -> np.ndarray:
        b = sa.body[tri_idx]
        r = points - sa.body_org[b]
        return sa.body_vel[b] + np.cross(sa.body_omega[b], r)

    # ---------------------------------------------------------------- trace
    def trace(self, scene: Scene, sensor_pose: Pose,
              sensor_vel: Sequence[float],
              directions_az_el: Sequence[Tuple[float, float]]) -> InteractionPoints:
        cfg = self.cfg
        sa = SceneArrays(scene)
        if sa.n == 0 or not directions_az_el:
            return InteractionPoints.empty()

        m = len(directions_az_el)
        az0 = np.array([d[0] for d in directions_az_el])
        el0 = np.array([d[1] for d in directions_az_el])

        # solid angle represented by one ray (FoV rectangle on the unit sphere)
        fov_h = math.radians(cfg.fov_deg[0])
        fov_v = math.radians(cfg.fov_deg[1])
        omega_fov = 2.0 * fov_h * math.sin(0.5 * fov_v)
        d_omega = omega_fov / m

        # launch directions in world coordinates
        R = np.array(sensor_pose.R, dtype=np.float64)
        S = np.array(sensor_pose.t, dtype=np.float64)
        v_s = np.array(sensor_vel, dtype=np.float64)
        d_local = np.stack([np.cos(el0) * np.cos(az0),
                            np.cos(el0) * np.sin(az0),
                            np.sin(el0)], axis=-1)
        dirs = d_local @ R.T
        launch_dirs = dirs.copy()

        # transmit field: |E| = sqrt(P_t * G_t / (4 pi)) / L, phase -k L
        g_tx = 10.0 ** (self.tx_gain(az0, el0) / 10.0)
        e0 = np.sqrt(cfg.transmit_power_w * g_tx / (4.0 * math.pi))
        pol_tx_local = polarization_vector(az0, el0, cfg.polarization_transmit)
        pol = pol_tx_local @ R.T                      # world frame, real unit vectors
        field = pol.astype(np.complex128) * e0[:, None]

        origins = np.repeat(S[None, :], m, axis=0)
        path_len = np.full(m, 0.0)
        path_rate = np.full(m, 0.0)                   # d(path length)/dt
        prev_vel = np.repeat(v_s[None, :], m, axis=0)  # velocity of the last
        #                                               scattering point (or of
        #                                               the sensor for bounce 1)
        alive = np.ones(m, dtype=bool)
        # ray must start at Range[0] and stop at Range[1]  (documented)
        r_min, r_max = cfg.range_min_max

        out = InteractionPoints.empty()

        for bounce in range(1, cfg.max_bounces + 1):
            if not alive.any():
                break
            idx = np.nonzero(alive)[0]
            o = origins[idx]
            d = dirs[idx]
            t_lim = np.full(idx.size, r_max) - path_len[idx]
            t_lim = np.maximum(t_lim, 0.0)
            hit, t, tri = intersect_batch(sa, o, d, r_min if bounce == 1 else 1e-4, t_lim)
            if not hit.any():
                alive[idx] = False
                continue
            sel = idx[hit]
            t = t[hit]
            tri = tri[hit]
            p = origins[sel] + dirs[sel] * t[:, None]

            # ---- accumulated geometry ---------------------------------
            L = path_len[sel] + t
            # velocity of the illuminated material point (rigid body motion)
            v_p = self._point_velocity(sa, tri, p)
            # rate of change of this segment's length
            seg_dir = dirs[sel]
            v_prev = prev_vel[sel]
            seg_rate = np.einsum("ij,ij->i", v_p - v_prev, seg_dir)
            L_rate = path_rate[sel] + seg_rate

            # ---- surface normal with material roughness ----------------
            n = sa.nrm[tri].copy()
            flip = np.einsum("ij,ij->i", n, seg_dir) > 0.0
            n[flip] = -n[flip]
            sig = sa.scattering[tri]
            if np.any(sig > 0.0):
                n = _perturb_normals(n, sig, self.rng)

            cos_i = -np.einsum("ij,ij->i", n, seg_dir)
            cos_i = np.clip(cos_i, 1e-4, 1.0)

            # illuminated patch: ray tube cross-section / cos(theta_i)
            area_proj = (L ** 2) * d_omega                    # = A * cos(theta_i)
            patch_side = np.sqrt(np.maximum(area_proj / cos_i, 1e-12))

            # ---- reception path 1: direct line of sight ----------------
            to_rx = S[None, :] - p
            r2 = np.linalg.norm(to_rx, axis=1)
            u_rx = to_rx / np.maximum(r2, 1e-12)[:, None]
            los = _visible(sa, p, S, r2)
            if los.any():
                ip = self._reception(sa, tri, p, seg_dir, n, field[sel], cos_i,
                                     area_proj, patch_side, L, L_rate, r2, u_rx,
                                     R, S, v_s, v_p, los, back_prop=False,
                                     launch_dir=launch_dirs[sel])
                out = out.concat(ip)

            # ---- reception path 2: back-propagation (bounce >= 2) ------
            if bounce >= 2:
                u_back = -seg_dir
                ip = self._reception(sa, tri, p, seg_dir, n, field[sel], cos_i,
                                     area_proj, patch_side, L, L_rate,
                                     L.copy(), u_back, R, S, v_s, v_p,
                                     np.ones(sel.size, dtype=bool), back_prop=True,
                                     launch_dir=launch_dirs[sel])
                out = out.concat(ip)

            # ---- continue the ray: specular reflection ----------------
            r_s, r_p = fresnel(sa.permittivity[tri], cos_i)
            e_in = field[sel]
            s_hat = np.cross(seg_dir, n)
            s_norm = np.linalg.norm(s_hat, axis=1, keepdims=True)
            s_hat = np.where(s_norm > 1e-9, s_hat / np.maximum(s_norm, 1e-12),
                             np.array([0.0, 1.0, 0.0]))
            e_s_comp = np.einsum("ij,ij->i", e_in, s_hat.astype(np.complex128))
            e_perp = e_in - e_s_comp[:, None] * s_hat
            d_ref = seg_dir - 2.0 * np.einsum("ij,ij->i", seg_dir, n)[:, None] * n
            p_hat_out = np.cross(d_ref, s_hat)
            p_hat_in = np.cross(seg_dir, s_hat)
            e_p_comp = np.einsum("ij,ij->i", e_perp, p_hat_in.astype(np.complex128))
            e_new = (r_s * e_s_comp)[:, None] * s_hat + (r_p * e_p_comp)[:, None] * p_hat_out
            # spherical spreading is applied through 1/L at reception, so the
            # propagating field keeps only the reflection factors and damping
            e_new = e_new * self._atten(t)[:, None]

            # update state
            alive[:] = False
            alive[sel] = True
            origins[sel] = p + n * 1e-4
            dirs[sel] = d_ref
            field[sel] = e_new
            path_len[sel] = L
            path_rate[sel] = L_rate
            prev_vel[sel] = v_p

            # energy cut-off: stop tracing rays whose field has decayed
            rel_db = 20.0 * np.log10(np.maximum(
                np.linalg.norm(np.abs(field[sel]), axis=1)
                / np.maximum(e0[sel], 1e-30), 1e-12))
            alive[sel[rel_db < cfg.ray_energy_cutoff_db]] = False

        return out

    # ------------------------------------------------------------ reception
    def _reception(self, sa: SceneArrays, tri, p, d_in, n, e_in, cos_i,
                   area_proj, patch_side, L, L_rate, r2, u_out,
                   R, S, v_s, v_p, mask, back_prop: bool, launch_dir):
        """Scattered field at the receiver for one set of IA points."""
        cfg = self.cfg
        k = self.k
        sel = np.nonzero(mask)[0]
        if sel.size == 0:
            return InteractionPoints.empty()

        n_ = n[sel]
        d_in_ = d_in[sel]
        u_out_ = u_out[sel] if u_out.ndim == 2 else u_out
        r2_ = r2[sel] if np.ndim(r2) else np.full(sel.size, float(r2))
        cos_i_ = cos_i[sel]
        L_ = L[sel]
        L_rate_ = L_rate[sel]
        e_in_ = e_in[sel]
        area_ = area_proj[sel]
        side_ = patch_side[sel]

        cos_o = np.einsum("ij,ij->i", n_, u_out_)
        valid = cos_o > 1e-4
        if not valid.any():
            return InteractionPoints.empty()

        # ---- physical-optics aperture factor Lambda(q) -----------------
        # q = k (d_in - u_out); only its tangential part matters
        q = k * (d_in_ - u_out_)
        q_t = q - np.einsum("ij,ij->i", q, n_)[:, None] * n_
        t1 = _tangent(n_)
        t2 = np.cross(n_, t1)
        lam = _sinc(0.5 * np.einsum("ij,ij->i", q_t, t1) * side_) * \
              _sinc(0.5 * np.einsum("ij,ij->i", q_t, t2) * side_)

        # ---- Fresnel for the outgoing direction ------------------------
        r_s, r_p = fresnel(sa.permittivity[tri][sel], cos_i_)
        s_hat = np.cross(d_in_, n_)
        s_norm = np.linalg.norm(s_hat, axis=1, keepdims=True)
        s_hat = np.where(s_norm > 1e-9, s_hat / np.maximum(s_norm, 1e-12),
                         np.array([0.0, 1.0, 0.0]))
        e_s_comp = np.einsum("ij,ij->i", e_in_, s_hat.astype(np.complex128))
        e_perp = e_in_ - e_s_comp[:, None] * s_hat
        p_hat_out = np.cross(u_out_, s_hat)
        p_hat_in = np.cross(d_in_, s_hat)
        e_p_comp = np.einsum("ij,ij->i", e_perp, p_hat_in.astype(np.complex128))
        e_scat_pol = (r_s * e_s_comp)[:, None] * s_hat + (r_p * e_p_comp)[:, None] * p_hat_out

        # ---- PO scattered field at the receiver ------------------------
        total_path = L_ + r2_
        # E_s = j k (L^2 dOmega) / (2 pi) * Gamma * Lambda * E_inc(P) * e^{-jk r2} / r2
        # with E_inc(P) = E_launch / L  (spherical spreading applied here)
        amp = (1j * k * area_ / (2.0 * math.pi) * lam
               / (np.maximum(L_, 1e-6) * np.maximum(r2_, 1e-6)))
        phase = np.exp(-1j * k * total_path)
        e_rx = e_scat_pol * (amp * phase * self._atten(total_path))[:, None]

        # ---- receive antenna: direction, gain, polarisation ------------
        # Look direction = the direction the sensor *believes* the echo comes
        # from.  For a direct (line-of-sight) return that is the direction from
        # the sensor to the interaction point; for a back-propagated return the
        # wave re-enters along the original launch path, which is what produces
        # ghost targets and NLOS detections.
        if back_prop:
            d_look_world = launch_dir[sel]
        else:
            d_look_world = -u_out_
        d_rx_local = d_look_world @ R           # world -> sensor frame (R^T v)
        r_xy = np.hypot(d_rx_local[:, 0], d_rx_local[:, 1])
        az_rx = np.arctan2(d_rx_local[:, 1], d_rx_local[:, 0])
        el_rx = np.arctan2(d_rx_local[:, 2], r_xy)

        g_rx = 10.0 ** (self.rx_gain(az_rx, el_rx) / 10.0)
        pol_rx = polarization_vector(az_rx, el_rx, cfg.polarization_receive) @ R.T
        proj = np.einsum("ij,ij->i", e_rx, pol_rx.astype(np.complex128))
        v_rx = np.sqrt(g_rx * cfg.wavelength ** 2 / (4.0 * math.pi)) * proj
        v_rx = v_rx / (10.0 ** (cfg.system_losses_db / 20.0))

        # ---- range and radial velocity ---------------------------------
        if back_prop:
            rng = L_                                  # (L + L) / 2
            rate = L_rate_                            # d/dt of the one-way path
        else:
            rate_r2 = np.einsum("ij,ij->i", (v_s[None, :] - v_p[sel]), u_out_)
            rng = 0.5 * (L_ + r2_)
            rate = 0.5 * (L_rate_ + rate_r2)

        good = valid & (np.abs(v_rx) > 0.0) & np.isfinite(np.abs(v_rx))
        return InteractionPoints(rng[good], rate[good], az_rx[good], el_rx[good],
                                 v_rx[good])


# ----------------------------------------------------------------- helpers
def _sinc(x: np.ndarray) -> np.ndarray:
    return np.sinc(x / math.pi)          # np.sinc(z) = sin(pi z)/(pi z)


def _tangent(n: np.ndarray) -> np.ndarray:
    """Any unit vector orthogonal to n (numerically stable)."""
    a = np.where(np.abs(n[:, 0:1]) < 0.9, np.array([1.0, 0.0, 0.0]),
                 np.array([0.0, 1.0, 0.0]))
    t = np.cross(n, a)
    return t / np.maximum(np.linalg.norm(t, axis=1, keepdims=True), 1e-12)


def _perturb_normals(n: np.ndarray, sigma: np.ndarray, rng) -> np.ndarray:
    """Material-dependent stochastic component on the surface normal
    (Reference Manual -> Material-Dependent Reflection).  ``sigma`` is the
    material's ``Radar.Scattering`` value converted to radians."""
    t1 = _tangent(n)
    t2 = np.cross(n, t1)
    a = rng.normal(0.0, 1.0, size=n.shape[0]) * sigma
    b = rng.normal(0.0, 1.0, size=n.shape[0]) * sigma
    out = n + t1 * a[:, None] + t2 * b[:, None]
    return out / np.maximum(np.linalg.norm(out, axis=1, keepdims=True), 1e-12)


def _visible(sa: SceneArrays, p: np.ndarray, s: np.ndarray,
             dist: np.ndarray) -> np.ndarray:
    """Line-of-sight test from every IA point to the sensor."""
    d = (s[None, :] - p) / np.maximum(dist, 1e-12)[:, None]
    hit, _, _ = intersect_batch(sa, p + d * 1e-3, d, 1e-4, dist - 1e-2)
    return ~hit
