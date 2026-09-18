# Implementation logic — Radar Model 1 (object-list Radar Sensor)

Implementation-oriented specification. Follow it top to bottom and you get the
model without re-reading the IPG material. Everything is expressed so that it
maps 1:1 onto C++ (`struct`, `std::vector`, `double`).

Reference implementation: `config.py`, `maps.py`, `atmosphere.py`,
`targets.py`, `detection.py`, `model.py`.

---

## 1. Coordinate systems

| Frame | Definition |
|---|---|
| World `Fr0` | right-handed, z up |
| Mounting frame | vehicle body frame (`Fr1A`/`Fr1B`), x forward, y left, z up |
| Sensor frame | `pos` + `rot` (deg, applied **z-y-x**: `R = Rz(rz)·Ry(ry)·Rx(rx)`) relative to the mounting frame; **boresight = +x** |
| Object body frame | x along the object's longitudinal axis, origin at the bounding-box centre |

Angles in the sensor frame, for a direction `d`:

```
azimuth   phi   = atan2(d_y, d_x)
elevation theta = atan2(d_z, hypot(d_x, d_y))
```

Aperture direction cosines used by the antenna model:

```
u_y = cos(theta) * sin(phi)      # == sin(Theta) cos(Phi) of Eq. 508
u_z = sin(theta)                 # == sin(Theta) sin(Phi) of Eq. 508
```

---

## 2. Data structures

```c
struct RadarConfig {            // all Sensor.Param.<n>.* keys, see config.py
    double fov[2], range_min, range_max;  int max_num_obj;
    double frequency_ghz, transmit_power_dbm, system_losses_db;
    double noise_bandwidth_hz, noise_figure_db;
    double prob_detect_min;     int prob_false_alarm_idx;   // 1..10
    double acc_dist, acc_az_deg, acc_speed_kmh;
    double res_dist, res_az_deg, res_speed_kmh, separability;
    int    false_pos_active;    double clutter_obj_mean;
    double cycle_time_s, cycle_offset_s, latency_factor[2];
    double temperature_k, rain_rate_mm_h, vis_range_fog_m;
    /* inferred, see explanation.md §10 */
    double clutter_sigma0_db, rain_k, rain_alpha, fog_k;
    int    seed;
};

struct AntennaGainMap {  double az[N], el[M], gain_db[N][M];
                         double beam_width[2], scan_range[2], antenna_eff; };

struct RcsMap        {  double azim[K], rcs_lin[K];     // linear m^2
                        int prob_exist; double occlusion_factor; };

struct Target {                       // one entry of the sensor surrounding
    int obj_id;  double pos[3], vel[3], acc[3];  double yaw, pitch, roll;
    double length, width, height;     // Basics.Dimension
    const RcsMap* rcs_map;            // NULL -> undetectable
    int detect_mask;  const char* name;
};

struct TargetView {                   // derived, sensor frame
    double ds[3], dv[3], da[3];
    double range, range_rate, azimuth, elevation;
    double az_min, az_max, el_min, el_max;
    double rcs_azimuth, rcs_az_min, rcs_az_max, rel_course_angle;
    double corners_s[8][3];
};

struct RadarObject { /* = tOutQuants of Sensor_Radar.h */ };
struct TrackState  { int prob_exist; RadarObject last_output; int valid; };
```

Persistent model state: `map<obj_id, TrackState>`, rolling counter,
last-calculation time, queue of `(release_time, output)` pairs, RNG.

---

## 3. Initialisation

1. Validate parameters (`MaxNumObj ∈ (0, 2000]`, `ProbFalseAlarmIdx ∈ [1,10]`,
   `FoV ∈ (0, 180]°`, `ResolutionAzimuth ≥ 1e-5`, `NoiseBandWidth > 0`, …).
2. `lambda = c / (Frequency * 1e9)`, `c = 299 792 458 m/s`.
3. Load or generate the antenna gain map (§4).
4. Load one `RcsMap` per object class.
5. Precompute
   * `SNR_min = 2 (erfcinv(2 P_FA) − erfcinv(2 P_Dmin))²`
   * `N_thermal = T0 · 10^(NF/10) · k_B · B_n`, `k_B = 1.380649e-23`
   * `gamma_damp [dB/km] = gamma_atm(f) + k·R^alpha + K_fog·M(d_vis)`
6. Seed the RNG from the Test Run `RandomSeed`.

`erfcinv(y)` for `y ∈ (0,2)`: start from `−Phi^-1(y/2)/sqrt(2)` (Acklam's
rational approximation) and run 3 Newton steps on `erfc(x) − y`
(derivative `−2/sqrt(pi)·exp(−x²)`). Precomputed table for the ten `P_FA`:

```
idx   1        2        3        4        5        6        7        8        9        10
x  0.906194 1.644976 2.185124 2.629742 3.015733 3.361179 3.676487 3.968284 4.241090 4.498147
```

---

## 4. Antenna gain map

**Loading** (`FileIdent = CarMaker-AntennaGainMap`): `Type` (`rad`/`deg`),
`Azimuth` (1×n, equidistant, −90…90°), `Elevation` (1×m), `Gain` (n×m, dB),
plus `ScanRange`, `BeamWidth`, `AntennaEff`. Lookup = bilinear interpolation,
clamped at the borders.

**Generating** (uniform rectangular aperture, Eq. 507/508):

```
a_lam = 0.8858929 / BW_h_rad          # 3 dB half-power constant of a sinc
b_lam = 0.8858929 / BW_v_rad
G0_dB = 10 log10( 4 pi * a_lam * b_lam * AntennaEff )
u0 = sin(ScanRange_h),  v0 = sin(ScanRange_v)

for each (az, el):
    u_y = cos(el) sin(az);   u_z = sin(el)
    f   = sinc(pi a_lam (u_y - u0)) * sinc(pi b_lam (u_z - v0))
    f  *= 0.5 * (1 + cos(az) cos(el))          # INFERRED obliquity factor
    G_dB(az, el) = G0_dB + 20 log10(|f|)
```

with `sinc(x) = sin(x)/x`, `sinc(0) = 1`.

---

## 5. Per-cycle algorithm

```
function calculate(t, sensor_pose, sensor_vel, targets, ego_state):

  ## (1) candidates
  views = []
  for tgt in targets:
      if !tgt.detect_mask or tgt.rcs_map == NULL: continue
      v = view_target(tgt, sensor_pose, sensor_vel)          # §6
      if v.range < Range_min or v.range > Range_max: continue
      if v.az_max < -FoV_h/2 or v.az_min > +FoV_h/2: continue
      if v.el_max < -FoV_v/2 or v.el_min > +FoV_v/2: continue
      views.push(v)
  sort views by range ascending           # nearer objects occlude farther ones

  ## (2) link budget
  cands = []
  for v in views:
      c = link_budget(v, views)                              # §7
      if c.snr > SNR_min: cands.push(c)

  ## (3) merging
  cands = merge(cands)                                       # §8

  ## (4) false positives
  if FalsePosActive: cands += false_positives(views)         # §9

  ## (5) measurement noise + output record
  objects = [to_output(c) for c in cands]                    # §10

  ## (6) object management
  objects = object_management(objects)                       # §11

  ## (7) global info
  RolCount++;  nObj = len(objects);  RelvTgt = index of nearest ProbObst > 0
  return objects, global_info
```

Cycle/latency wrapper:

```
step(t):
    if t >= t_last_calc + CycleTime:
        t_last_calc = t
        pending.push( (t + Latency[0]*CycleTime, calculate(t, ...)) )
    while pending.front().release <= t: current = pending.pop()
    return current
```

---

## 6. `view_target` — geometry

```
ds  = R_s^T (p_obj - p_sensor)                 # sensor frame
dv  = R_s^T (v_obj - v_sensor)
da  = R_s^T a_obj
range      = |ds|
azimuth    = atan2(ds_y, ds_x)
elevation  = atan2(ds_z, hypot(ds_x, ds_y))
range_rate = dot(dv, ds) / range               # + = receding

corners (object frame, order of ePntType in Sensor_Radar_protected.h):
    (±l/2, ±w/2, ±h/2)                         # centred on the BBox centre
corners_s[i] = R_s^T (R_obj corner_i + p_obj - p_sensor)
az_min/max, el_min/max = extrema over corners_s, azimuths unwrapped
                         around the centre azimuth

# RCS-map argument: sensor position in the object body frame
s_obj      = R_obj^T (p_sensor - p_obj)
rcs_azimuth = atan2(s_obj_y, s_obj_x)
rcs_az_min/max = extrema of atan2(s_obj_y - c_y, s_obj_x - c_x)
                 over the four bottom corners c (unwrapped)

rel_course_angle = atan2(dv_y, dv_x)
```

---

## 7. `link_budget`

```
G_dB      = antenna.gain(azimuth, elevation)
damp_dB   = gamma_damp[dB/km] * 2 * range / 1000        # two-way

# --- RCS
if extended_object_rcs:
    span = max(rcs_az_max - rcs_az_min, ResolutionAzimuth_rad)
    mid  = 0.5*(rcs_az_min + rcs_az_max)
    rcs_lut = mean over 9 samples of rcs_map(mid ± span/2)      # LINEAR average
else:
    rcs_lut = rcs_map(rcs_azimuth)

rcs = exponential_random(mean = rcs_lut)                # Swerling 1, Eq. 509
(o_h, o_v) = occlusion(view, all_views)                 # §7a
rcs *= (1 - o_h * o_v)                                  # Eq. 510

# --- signal, Eq. 505
S = P_w * 10^(G_dB/5) * lambda^2 * rcs
    / ( (4 pi)^3 * range^4 * 10^(L_A/10) * 10^(damp_dB/10) )
        # note 10^(G_dB/5) == (10^(G_dB/10))^2  (G appears squared)

# --- noise, Eq. 506 + clutter
N = N_thermal
if clutter_enabled:
    psi   = max(atan2(h_sensor, range), psi_min)
    A_c   = range * ResolutionAzimuth_rad * ResolutionDistance / cos(psi)
    sig_c = 10^(sigma0_dB/10) * A_c
    N    += radar_equation(P_w, G_dB, lambda, sig_c, range, L_A, damp_dB)

SNR       = S / N
ProbDetect= 0.5 * erfc( erfcinv(2 P_FA) - sqrt(SNR/2) )
rcs_out   = rcs * (1 + 1/SNR)                           # noise correction
```

### 7a. Occlusion

```
o_h = o_v = 0
for other in views with other.range < view.range:
    f = other.rcs_map.occlusion_factor            # 0 = transparent, 1 = opaque
    cov_az = min(view.az_max, other.az_max) - max(view.az_min, other.az_min)
    cov_el = min(view.el_max, other.el_max) - max(view.el_min, other.el_min)
    if cov_az <= 0 or cov_el <= 0: continue
    o_h = max(o_h, f * min(1, cov_az / (view.az_max - view.az_min)))
    o_v = max(o_v, f * min(1, cov_el / (view.el_max - view.el_min)))
```

---

## 8. Merging (separability)

```
d_r = Separability * ResolutionDistance                    # [m]
d_v = Separability * ResolutionSpeed / 3.6                 # [m/s]
d_a = Separability * ResolutionAzimuth_rad                 # [rad]

single-linkage clustering: i and j belong together iff
    |r_i - r_j| < d_r  AND  |vr_i - vr_j| < d_v  AND  |az_i - az_j| < d_a
```

Combining a cluster (Eq. 511 and *Separability*):

* new bounding box = axis-aligned envelope of all member corners in the sensor
  frame → new `ds`, `range`, `azimuth`, `elevation`, `Length`, `Width`;
* `RCS = mean(RCS_i)`;
* `S = sum(S_i)`, `N` from the strongest member, `SNR = S/N`, `ProbDetect` from `SNR`;
* kinematics = signal-power-weighted mean of the members;
* `GroupId` = id of the strongest member, `LeftInGroup` = id of the member with
  the largest azimuth ("most left in sensor view").

---

## 9. False positives

**Clutter** (documented): `n ~ Poisson(ClutterObjMean)`; for each,
`range ~ U(Range_min, Range_max)`, `az ~ U(−FoV_h/2, +FoV_h/2)`,
`el ~ U(−FoV_v/2, +FoV_v/2)`, length/width uniform, `ObjId = −2`.

**Mirror objects** (partly inferred): for every non-wall object and every
wall-like object (guardrail / tunnel wall) in the scene,

```
d_wall = |y_obj - y_wall|                    (sensor frame)
p      = truck_boost * exp(-d_wall / 5 m) * exp(-range / Range_max)
truck_boost = 2 if length > 8 m else 1
if uniform() < p:  emit a candidate mirrored across the wall plane,
                   ObjId = -obj_id, lateral velocity sign flipped
```

Objects not on the road or the RoadMargin are excluded, and false positives are
skipped when the reference point is exactly on a junction link (documented).
False positives are never grouped.

---

## 10. Measurement noise and output record

```
r   = max(range     + normal(0, AccuracyDistance), 0)
az  = azimuth       + normal(0, AccuracyAzimuth_rad)
el  = elevation     + normal(0, AccuracyAzimuth_rad)
v_r = range_rate    + normal(0, AccuracySpeed/3.6)

DistX = r cos(el) cos(az);  DistY = r cos(el) sin(az);  DistZ = r sin(el)
u     = (cos az, sin az)
VrelX = dv_x + (v_r - dot(dv_xy, u)) * u_x        # keep the tangential part
VrelY = dv_y + (v_r - dot(dv_xy, u)) * u_y
Length/Width = ground truth + normal(0, 0.5 * ResolutionAzimuth_rad * r)   # INFERRED
LengthClass/WidthClass: 0 if <=0; else first of (0.5,1,2,3,4,6) that exceeds; else 7
RCS            = 10 log10(rcs_out)        [dBm^2]
SignalStrength = 10 log10(S)              [dBW]
SNR            = 10 log10(SNR)            [dB]
DynProp: 1 stationary / 2 stopped / 3 moving / 4 oncoming
ProbObst: see below
```

`ProbObst` (documented):

```
half_path = ego_width/2 + 0.3
y_path    = 0                                    if |yaw_rate| ~ 0
          = sign(R) (|R| - sqrt(R^2 - x^2)),  R = v_ego / yaw_rate   otherwise
g         = max(|DistY - y_path| - half_path - Width/2, 0)
ProbObst  = 0                if g >= 1 m
          = 1 - g / 1 m      otherwise
```

---

## 11. Object management

```
for each detected object:
    if not tracked: create TrackState, prob_exist = rcs_map.ProbExist,
                    MeasStat = 1 (new object)
    else:           MeasStat = 3 (object measured)
    prob_exist += (ProbDetect >= ProbDetectMin) ? +1 : -1, clamped to [0, 7]

for each tracked object NOT detected this cycle:
    prob_exist -= 1
    if prob_exist <= 0: erase the track
    else: re-emit a COPY of the last output with MeasStat = 2

sort the list by range;  truncate to MaxNumObj
```

The copy matters: the held record must not alias a record already handed to the
user.

---

## 12. Randomness

One RNG stream per sensor instance, seeded by the Test Run `RandomSeed`
(`std::mt19937_64` is a faithful stand-in for CarMaker's generator; the *stream*
will differ from IPG's, only the distributions match).

| Where | Distribution |
|---|---|
| Swerling-1 RCS | `Exp(mean = RCS_lut)` → `-mean · ln(U)` |
| Measurement noise | `N(0, accuracy)` |
| Size noise | `N(0, 0.5·Θ_res·r)` |
| Clutter count | `Poisson(ClutterObjMean)` (Knuth for λ < 30) |
| Clutter placement | uniform |
| Mirror occurrence | Bernoulli(p) |

---

## 13. Edge cases

* `rcs_map == NULL` or `DetectMask[Sensor] == 0` → object is invisible.
* `range < Range_min` → skipped (sensor origin inside the vehicle body).
* Object straddling ±180° azimuth → unwrap corner azimuths around the centre.
* Normal incidence / degenerate cross products in the geometry → guard with 1e-9.
* `SNR = 0` → `ProbDetect = 0`, no RCS noise correction.
* `noise = 0` (impossible with `B_n > 0`) → guard division.
* Merging is single-linkage: a chain of nearby objects becomes one object.
* Held objects (`ProbExist > 0`, not detected) still occupy `MaxNumObj` slots.
* `Range_max` changes at run time (DVA) only affect static objects after the
  surrounding update (documented restriction).

---

## 14. Expected outputs

Per cycle: `nObj` records of `tOutQuants` plus one `tOutQuantsGlob`.
Deterministic given `(RandomSeed, scenario, parameters)`.
Output age between `Latency` and `Latency + CycleTime`.

---

## 15. Assumptions required for an independent implementation

1. Road clutter: choose `sigma0_dB` — the default −35 dB gives a
   thermal-noise-dominated link at short range and a clutter-limited one beyond
   ~150 m for small targets. **Not an IPG value.**
2. Damping coefficients (rain `k`, `alpha`; fog `K_l`; clear air) are public
   literature values, not IPG's.
3. The RCS noise correction, the extended-object averaging rule, the occlusion
   combination rule, the object ordering and `RelvTgt` tie-breaking are
   interpretations (see `explanation.md` §10).
4. `nLanesL/R` and `DistToLeft/RightBorder` require a road model; the reference
   implementation takes them as optional caller-supplied inputs.
5. The "tracking behaviour" mentioned for `Length`/`Width` is modelled only as
   noise; no filter is documented.

---

## 16. Notes for the C++/CARLA port

* `RadarSensorModel::_calculate` is a pure function of
  `(t, sensor pose, sensor velocity, target array, ego state, model state)`.
  Everything else is bookkeeping.
* Replace `dataclass` records by `struct`s, `dict` tracks by
  `std::unordered_map<int, TrackState>`, `Rng` by `std::mt19937_64` +
  `std::normal_distribution` / `std::exponential_distribution` /
  `std::poisson_distribution`.
* In CARLA the `Target` list is built from `carla::Actor` bounding boxes and
  velocities; the RCS map is a per-blueprint asset. Occlusion can keep using
  bounding boxes, or be replaced by CARLA ray casts without touching the rest.
* Cost is `O(n_objects²)` for occlusion and merging — fine for the ≤ 2000
  objects the model allows.
