# Implementation logic — Radar Model 2 (Radar RSI)

Implementation-oriented specification. Reference implementation:
`config.py`, `scene.py`, `raypattern.py`, `transceiver.py`, `propagation.py`,
`device.py`, `model.py`.

---

## 1. Coordinate systems and conventions

Same as Model 1: sensor frame with **boresight +x**, y left, z up;
`azimuth = atan2(d_y, d_x)`, `elevation = atan2(d_z, hypot(d_x, d_y))`;
`Sensor.<n>.rot` applied **z-y-x**.

A world-frame vector is rotated into the sensor frame with `Rᵀ v`
(in row-vector form: `v @ R`). Getting this wrong mirrors every reported
azimuth — the single most likely bug in a port.

Virtual receivers lie along **y** at spacing `d = λ/2` (`VRxCalibration`
columns are `y [m]`, `z [m]`, weight amplitude, weight phase [rad]).

---

## 2. Data structures

```c
struct RadarMaterial { const char* name; double permittivity, scattering_rad; };

struct Triangle { double v0[3], v1[3], v2[3]; const RadarMaterial* mat; };

struct Body { Pose pose; double velocity[3], omega[3]; Triangle* tris; int n; };

struct SceneArrays {              // world space, rebuilt whenever bodies move
    double p0[T][3], e1[T][3], e2[T][3], n[T][3];
    double permittivity[T], scattering[T];  int body[T];
    double body_vel[B][3], body_omega[B][3], body_org[B][3];
};

struct InteractionPoint {         // == tUserRadarInteractionPoint
    float range;                  // one-way range  [m]
    float vel;                    // relative radial velocity [m/s]
    float az, el;                 // look direction in the sensor frame [rad]
    complex<float> amplitude;     // complex voltage at the receiver [sqrt(W)]
};

struct RadarCube { int n_r, n_v, n_vrx; complex<double>* data; };
// index = i_r + i_v*n_r + i_vrx*n_r*n_v          (documented layout)

struct DetectionPoint { double coords[3]; double power_dbm, velocity; };
struct VRxDetection   { double range, velocity; complex<double> amp[n_vrx]; };
```

---

## 3. Ray patterns

```
Omega_FoV = 2 * FoV_h * sin(FoV_v / 2)          # solid angle of the FoV rectangle
dOmega    = Omega_FoV / n_rays_actual           # solid angle per ray
```

**PhiTheta**: `n_v = floor(sqrt(nRays / (FoV_h/FoV_v)))`, `n_h = floor(nRays/n_v)`;
cell centres of an `n_h × n_v` grid.

**Fibonacci**: for `i = 0 … n-1`
```
u  = (i + 0.5) / n
v  = frac( i * golden_angle / (2 pi) ),   golden_angle = pi (3 - sqrt 5)
el = asin( sin(FoV_v/2) * (2u - 1) )      # equal solid angle in elevation
az = FoV_h * (v - 0.5)
```

**Dynamic**:
```
n_h = nHorizontalBuckets  (default ceil(FoV_h_deg)), likewise n_v
occupied[n_h][n_v] = false
if mode != OnlyInstances:      shoot one ray through each bucket centre;
                               mark the bucket if it hits anything within Range[1]
if mode != OnlySceneScan:      for each instance with bounding-sphere radius <= 25 m
                                   and centre within Range[1]:
                                   half = asin(radius / dist)
                                   mark all buckets inside
                                   [az +- half] x [el +- half]
if no bucket occupied:         fall back to PhiTheta over the whole FoV
per_bucket = min(nRays / n_occupied, MaxNumberOfRaysPerBucket)
k          = floor(sqrt(per_bucket))          # k x k sub-grid per bucket
emit k*k directions per occupied bucket, stopping at nRays
```

---

## 4. Wave propagation (shooting and bouncing rays)

### 4.1 State per ray

```
origin[3], direction[3]          world frame
field[3]                         complex E, WITHOUT the 1/L spreading factor
launch_dir[3]                    original launch direction (for back-propagation)
path_len L                       accumulated one-way path length
path_rate dL/dt                  accumulated time derivative of L
prev_vel[3]                      velocity of the previous scattering point
                                 (sensor velocity for the first segment)
```

Initialisation:

```
E0        = sqrt( P_tx * G_tx(az, el) / (4 pi) )      # G from the Tx gain map, linear
pol_tx    = A_TE * E_TE(az, el) + A_TM * E_TM(az, el),  normalised
field     = pol_tx (rotated to world) * E0
L         = 0,  dL/dt = 0,  prev_vel = sensor velocity
```

`A_TE = sqrt(w_p)`, `A_TM = sqrt(1 - w_p)`, and (Eq. 523/524, with the corrected
Eq. 524):

```
E_TM(phi, theta) = (-sin(theta) cos(phi), -sin(theta) sin(phi), cos(theta))
E_TE(phi, theta) = (-sin(phi),             cos(phi),             0)
```

If the launch direction is outside the transmit gain map's FoV, the ray is not
sent (documented).

### 4.2 Per bounce

```
1. intersect the scene (Moeller-Trumbore); t limited to Range[1] - L
   (first bounce additionally starts at t >= Range[0])
2. hit point P, triangle index k, body index b
3. L      += t
   v_P     = body_vel[b] + cross(body_omega[b], P - body_org[b])
   dL/dt  += dot(v_P - prev_vel, direction)          # segment length rate
   prev_vel = v_P
4. normal n = face normal flipped towards the incoming ray, then perturbed:
      n' = normalize( n + t1 * N(0, sigma) + t2 * N(0, sigma) )
      sigma = material.Radar.Scattering [rad], t1, t2 orthonormal tangents
5. cos_i    = -dot(n', direction),  clamped to [1e-4, 1]
   area_proj = L^2 * dOmega                 # ray tube cross-section = A * cos_i
   patch_side a = sqrt(area_proj / cos_i)   # side of the illuminated patch
6. emit reception paths (§4.3)
7. continue: Fresnel-reflect the field, reflect the direction,
      d_ref = direction - 2 dot(direction, n') n'
      apply the atmospheric damping over t
      stop when the relative field level drops below ray_energy_cutoff_db
```

Fresnel (air → dielectric, relative permittivity `eps`):

```
cos_t = sqrt(eps - (1 - cos_i^2)) / sqrt(eps)
r_s   = (cos_i - sqrt(eps) cos_t) / (cos_i + sqrt(eps) cos_t)
r_p   = (sqrt(eps) cos_i - cos_t) / (sqrt(eps) cos_i + cos_t)
```

with the s/p basis
`s_hat = normalize(cross(direction, n'))`,
`p_hat_in = cross(direction, s_hat)`, `p_hat_out = cross(d_out, s_hat)`
(fall back to `s_hat = (0,1,0)` at exact normal incidence).
`Permittivity = 1e9` (metal) drives both coefficients to ±1, as intended.

### 4.3 Reception paths

Two per interaction point:

| path | outgoing direction `u_out` | `r2` | reported look direction | condition |
|---|---|---|---|---|
| direct (LOS) | `(S − P)/|S − P|` | `|S − P|` | `−u_out` | segment `P→S` not blocked |
| back-propagation | `−direction` | `L` | `+launch_dir` | bounce ≥ 2 |

For each:

```
cos_o = dot(n', u_out);  skip if <= 0

# physical-optics aperture factor
q     = k (direction - u_out)                 # k = 2 pi / lambda
q_t   = q - dot(q, n') n'
Lambda = sinc( 0.5 * dot(q_t, t1) * a ) * sinc( 0.5 * dot(q_t, t2) * a )
         # sinc(x) = sin(x)/x ; == 1 in the specular direction

# Fresnel into the outgoing direction (same decomposition as §4.2)
E_scat = r_s (E·s_hat) s_hat + r_p (E·p_hat_in) p_hat_out

# scattered field at the receiver
E_rx = E_scat * j k * area_proj / (2 pi) * Lambda / (L * r2)
              * exp(-j k (L + r2)) * atten(L + r2)

# receive antenna
(az_rx, el_rx) = angles of the look direction in the sensor frame
G_rx           = 10^(rx_gain_map.nearest_lower(az_rx, el_rx) / 10)
pol_rx         = A_TE E_TE(az_rx, el_rx) + A_TM E_TM(az_rx, el_rx)   (receive w_p)
V              = sqrt(G_rx * lambda^2 / (4 pi)) * dot(E_rx, pol_rx)
V             /= 10^(SystemLosses / 20)

# range and radial velocity
direct : range = (L + r2)/2 ,  vel = (dL/dt + dot(v_sensor - v_P, u_out)) / 2
back   : range = L           ,  vel = dL/dt
```

`atten(x) = exp(-alpha x)` with
`alpha [1/m] = gamma[dB/km] / 1000 / 8.6858896`.

### 4.4 Why this kernel

`area_proj = L² dΩ`, `E_inc = E0/L`, so for a flat plate at normal incidence in
the far field, summing over all rays that hit it gives

```
sum E_s = j k A Gamma / (2 pi R^2) * E0        ->   sigma = 4 pi A^2 |Gamma|^2 / lambda^2
```

which is exactly the physical-optics plate RCS. This is the calibration point of
the whole propagation stage (validated: mean error 0.76 dB).

---

## 5. Device model

### 5.1 Radar cube

Dimensions `n_r = Range.Samples + Range.ZeroPaddings`,
`n_v = DopplerVel.Samples + DopplerVel.ZeroPaddings`, `n_vrx = |VRxCalibration|`.

Axes (same convention as the shipped `Radar.cu`):

```
range_axis[i] = i * Range.Max / (n_r - 1)
vel_axis[j]   = v_min + j * (v_max - v_min) / (n_v - 1)
```

Fractional bin positions of an IA point:

```
r  = clamp(range, 0, Range.Max)
v  = v_min + mod(velocity - v_min, v_max - v_min)      # unambiguous folding
p_r = r / Range.Max * (n_r - 1)
p_v = (v - v_min) / (v_max - v_min) * (n_v - 1)
```

ULA steering vector (per IA point, per virtual receiver `m`):

```
steer[m] = exp( j 2 pi * y_m * sin(azimuth) / lambda ) * w_m
```

Accumulation, `"spectral"` mode (default):

```
K(delta, N_s, N_bins, w) = ( sum_{n<N_s} w[n] exp(j 2 pi n delta / N_bins) ) / sum(w)

for each IA point, for each bin (i_r, i_v) within +-halfwidth of (p_r, p_v):
    cube[i_r][i_v][m] += A * K(i_r - p_r, ...) * K(i_v - p_v, ...) * steer[m]
```

`K` is the analytic DTFT of a windowed tone — it reproduces resolution, spectral
leakage, the window trade-off and (through its periodicity) aliasing without
generating time samples. `"binning"` mode instead adds `A * steer[m]` to the
nearest bin only, matching the shipped GPU sample.

### 5.2 Noise (Eq. 525)

```
P_Noise = k_B * T_Noise * B_Noise * 10^(NF/10)
if NoiseScaling:
    P(i_r, i_v) = P_Noise * 10^(NF_R(range_axis[i_r])/10) * 10^(NF_V(vel_axis[i_v])/10)
cube += (N(0,1) + j N(0,1)) * sqrt(P/2)                 # per cell, per VRx
```

`NF_R`, `NF_V` are linear interpolations of `NoiseScalingRange` /
`NoiseScalingDopplerVel`, clamped outside the table.

### 5.3 OS-CFAR (Eq. 526)

```
power = |cube[:, :, 0]|^2
training cells: all (di, dj) with max(|di|, |dj|) in (n_guard, n_lay]
N_tot = number of training cells
m     = clamp( round(Threshold[%] / 100 * N_tot), 1, N_tot ) - 1     # 0-based
Z     = m-th smallest training value
detect if  power > Z * 10^(SNR_dB / 10)
```

### 5.4 Peak finder

2-D: keep a detection only if `power >= power(shifted by (di,dj))` for all
`|di|, |dj| <= PeakFinder.RangeDoppler.Layers`.
1-D (angular): keep only if it is the max in `± PeakFinder.AzimuthAngle.Layers`.

### 5.5 Peak interpolation (Eq. 527–529)

```
p    = 0.5 * (alpha - gamma) / (alpha - 2 beta + gamma)     clamped to [-0.5, 0.5]
y(p) = beta - 0.25 * p * (alpha - gamma)
```

Applied along range and along Doppler independently; the interpolated power is

```
P_peak = P_maxV * P_maxR / P_CFAR
```

and the same scheme is used for the 2-D angular spectrum.

### 5.6 Angular processing

```
x    = cube[i_r][i_v][0 : n_ULA] * window(n_ULA)
spec = fftshift( FFT(x, n = AzimuthAngle.Samples + AzimuthAngle.ZeroPaddings) ) / sum(window)
power = |spec|^2
CA-CFAR:  detect if power[i] > mean(training cells) * 10^(SNR_dB/10)
peak finder + parabolic interpolation -> fractional bin i + off
sin(azimuth) = 2 * (i + off - n_bins/2) / n_bins
azimuth      = asin( clamp(sin(azimuth), -1, 1) )
```

If nothing survives the angular CFAR, report the strongest bin so the RD
detection is not silently lost (implementation choice).

### 5.7 Output

```
Cartesian: (r cos(el) cos(az), r cos(el) sin(az), r sin(el))
Spherical: (r, az, el)
Power [dBm] = 10 log10(P) + 30
VRx mode  : amp[m] = cube[i_r][i_v][m] * sqrt(2 * 50 Ohm) * 1000   [mV]
truncate at nMaxDetections
```

---

## 6. Randomness

| Where | Distribution |
|---|---|
| Surface-roughness normal tilt | `N(0, Radar.Scattering)` per bounce, per ray |
| Front-end noise | complex `N(0, P_Noise/2)` per cube cell |

Both seeded from the Test Run `RandomSeed`.

---

## 7. Edge cases

* No triangles or no rays → zero detections.
* A ray that leaves the scene simply dies (no environment map).
* Exact normal incidence → `cross(direction, n) = 0`; use a fallback `s_hat`.
* `cos_o <= 0` (receiver on the back side of the facet) → no reception path.
* Occluded LOS → only the back-propagated path contributes.
* Near-field IA points (`R < 2D²/λ` for the illuminated patch) violate the
  model's own far-field assumption; increase `nRays` so the patches shrink.
* Velocities beyond `DopplerVel.MinMax` fold back (set `doppler_aliasing = False`
  to clamp like the GPU sample instead).
* Cube memory is `n_r · n_v · n_vrx · 16 bytes` — the shipped example
  (512 × 256 × 16) is 32 MB; with 60 VRx it is 126 MB.
* `nMaxDetections` truncates silently (documented).

---

## 8. Expected outputs

Per cycle: `nDetections` × (`tDetPoint` or `tDetVRx`) plus `TimeFired`.
Deterministic for a given seed, scene and parameter set.

---

## 9. Assumptions required for an independent implementation

1. The PO scattering kernel and its aperture factor (§4.4) — IPG only states
   that an analytic far-field solution is used.
2. The number of bounces and the ray energy cut-off.
3. The RD-map synthesis: `"spectral"` is mathematically equivalent to a real
   2-D DFT of the beat signal but is not literally what CarMaker computes;
   `"binning"` matches the shipped sample.
4. Damping coefficients (shared with Model 1).
5. VRx millivolt reference impedance (50 Ω).
6. Elevation angular processing needs a 2-D virtual array; only the azimuth
   chain is implemented.
7. IPG's ray tracer runs on the GPU against the Movie NX scene graph, including
   skeleton models for pedestrians and cylinder approximations for tree trunks.
   The reference implementation uses an explicit triangle soup instead.

---

## 10. Notes for the C++/CARLA port

* The per-ray algorithm in §4 is embarrassingly parallel — one thread per ray,
  one kernel launch per bounce, exactly as CarMaker does it. The numpy code is
  written in that shape (arrays of rays, one loop over bounces).
* Replace `intersect_batch` with a BVH query (or CARLA's `cast_ray` / an
  Embree/OptiX traversal). Nothing else in `propagation.py` depends on it.
* The device model is a small amount of dense linear algebra: three FFTs and two
  sliding-window filters. `std::complex<double>` + FFTW/cuFFT map directly.
* Geometry in CARLA: each actor's mesh with a per-material permittivity /
  roughness table replaces `MaterialLib`. Static geometry can be baked once;
  only `Body.pose` / `Body.velocity` change per tick.
* Keep the interaction-point struct byte-compatible with
  `tUserRadarInteractionPoint` if you ever want to swap CarMaker's GPU Coding
  Interface for your own kernels.
