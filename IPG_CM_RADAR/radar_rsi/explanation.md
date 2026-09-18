# Radar Model 2 — CarMaker "Radar RSI" (Raw Signal Interface radar)

CarMaker identifier: `Sensor.Param.<no>.Type = "RadarRSI"`
C interface: `include/Vehicle/Sensor_RadarRSI.h`
GPU interface + sample source: `Examples/GPUCodingInterface/src.GPUCodingInterface/`
Documentation: *Reference Manual → Sensors → Radar RSI* (chapters 1137–1163 of the TOC)

---

## 1. What the model simulates

A **radar device**, not an object list.

> "The Radar Raw Signal Interface (Radar RSI) is a radar device simulation. It
> includes a physics-based polarimetric electromagnetic wave propagation as well
> as a Device Model. The computed electromagnetic vector fields stimulate the
> parameterizable Device Model, which outputs data in the form of a point cloud."
> — *Radar RSI → Introduction*

It runs on the GPU (Movie NX / IPGMovie back end) and is meant for **component**
development: radar signal processing, MIMO/beamforming, sensor-model validation.
The output is a point cloud (or raw virtual-receiver amplitudes), not tracked
objects.

The two halves are strictly separated:

1. **Wave propagation** — ray tracing over the real 3D geometry, producing
   *interaction points* (IA points) with a complex field amplitude, a range, a
   radial velocity and a direction.
2. **Device Model** — antenna → front-end hardware → back-end software
   (DFTs, CFAR filters, peak finding), producing detections.

---

## 2. Inputs

| Input | Source |
|---|---|
| 3D scene geometry | Movie NX / IPGMovie instances — real triangle meshes, not bounding boxes |
| Per-face material | `Data/Sensor/MaterialLib`: `Material.<i>.Radar.Permittivity`, `.Radar.Scattering` |
| Body motion | per-instance velocity (multiple bodies → micro-Doppler) |
| Sensor pose | `Sensor.<n>.pos` / `.rot`, optionally `ExtMotion` (`t_ext`, `rot_zyx_ext`) |
| Transceiver config | `Sensor.Param.<n>.Antenna.FName` → `Data/Sensor/<file>`: Tx gain map, Rx gain map, VRx calibration |
| Sensor parameters | `Sensor.Param.<n>.*` — see §9 |
| Environment | rain rate, visual range in fog, via the same damping model as the object-list radar |

Object-model requirements are documented explicitly: faces small enough to carry
the curvature and corner reflectors, outward normals on a closed surface, and a
material tag on every face. Textured detail (e.g. a painted radiator grille) does
**not** act as a corner reflector.

---

## 3. Processing pipeline

```
                      +-----------------------------+
                      |  ray pattern                |   PhiTheta / Fibonacci / Dynamic
                      +--------------+--------------+
                                     v
  +------------------------- WAVE PROPAGATION -------------------------+
  |  launch ray: Tx gain, Tx polarisation weight, complex amplitude    |
  |  for each bounce (<= N):                                           |
  |     intersect scene -> interaction point                           |
  |     perturb the normal by the material's Radar.Scattering          |
  |     Fresnel reflection from Radar.Permittivity                     |
  |     free-space + atmospheric/rain/fog damping                      |
  |     emit reception paths:  (a) direct LOS to the receiver          |
  |                            (b) back-propagation along the ray      |
  |     accumulate path length and its time derivative -> Doppler      |
  +--------------------------------+-----------------------------------+
                                   v
                  interaction points  {range, vel, direction, E-field}
                                   v
  +--------------------------- DEVICE MODEL ---------------------------+
  |  1 antenna       Rx gain map + Rx polarisation weight              |
  |  2 front end     transmit power, losses,                           |
  |                  P_Noise = k_B T B 10^(NF/10)   x NoiseScaling LUTs|
  |  3 Range-Doppler-Processing   2-D DFT -> RD map per virtual receiver|
  |  4 Range-Doppler-Filtering    OS-CFAR  (+ 2-D Peak-Finder)          |
  |      |                                                             |
  |      +--> VRx mode: output complex amplitudes per VRx, stop here    |
  |      |                                                             |
  |  5 Angular-Processing         spatial DFT over a lambda/2 ULA       |
  |  6 Angular-Filtering          CA-CFAR (+ Peak-Finder)               |
  |  7 Peak-Interpolation         parabolic, Eq. 527-529                |
  +--------------------------------+-----------------------------------+
                                   v
                    point cloud: coordinates, power [dBm], velocity
```

The data cube layout is documented by the shipped GPU sample
(`Examples/.../GPUCodingInterface.h`):

```
index = index_range + index_velocity * n_r + index_vrx * n_r * n_v
```

---

## 4. The physics

### Ray patterns

* **Phi-Theta (static)** — equidistant angular grid; the actual ray count can be
  below `nRays` because of the equidistance constraint.
* **Fibonacci (static)** — spherical Fibonacci lattice; the irregularity
  "helps to reduce exaggerated interference effects that could occur for certain
  alignments between the sensor and planes in the 3D environment".
* **Dynamic** — the FoV is split into buckets of ≤ 1° side length
  (`120° × 20°` → 2400 buckets); a *scene scan* (one ray per bucket centre)
  and/or *instance detection* (bounding spheres, ≤ 25 m radius in Movie NX)
  marks buckets occupied; the rays are then spread evenly over the occupied
  buckets only. `nRays` is never exceeded.

### Propagation

* Each ray carries **polarisation and a complex amplitude**; the whole
  computation is complex-valued so that interference is exact.
* Reflection: **Fresnel** with the material's relative permittivity.
* Roughness: "a material dependent stochastic component on the normal of the
  surface is introduced" — driven by `Material.<i>.Radar.Scattering` [deg].
* Damping: free-space spreading plus the frequency/rain-rate/fog-visibility
  model shared with the object-list radar.
* **Multipath**: every interaction is tracked both outward (further tracing) and
  toward the receiver. Reception uses **both** the line-of-sight path (if
  available) and the **back-propagation** of the incoming ray. The latter is
  what produces ghost targets (a target mirrored in a guardrail appears at the
  mirrored angle) and NLOS detections.
* **Doppler** "incorporating all reflections relative to the transmitter";
  bodies moving relative to each other give **micro-Doppler**.
* Every IA point must be in the **far field** — the model is only valid with
  enough rays per solid angle.

### Polarisation (Eq. 521–524)

```
A_TE = sqrt(w_p)                     A_TM = sqrt(1 - w_p)
E_TM = (-sin(theta) cos(phi), -sin(theta) sin(phi), cos(theta))
E_TE = (-sin(phi),             cos(phi),            0)
```

`w_p` is `PolarizationTransmit` / `PolarizationReceive`.

> ⚠️ The manual prints Equation 524 with `theta` in place of `phi`. With `theta`
> the two vectors are not orthogonal, so we use the azimuthal unit vector
> `(-sin phi, cos phi, 0)`.

### Front end (Eq. 525)

```
P_Noise = k_B * T_Noise * B_Noise * 10^(NF/10)
```

modelled as white Gaussian noise, scaled by two look-up tables
`NoiseScalingRange(range)` and `NoiseScalingDopplerVel(velocity)` in dB, and
reduced by the logarithmic `SystemLosses`.

### Back end

* **Range-Doppler-Processing**: two DFTs over `n_Range`/`n_Vel` samples with
  `n_zero` zero paddings and a window (`Rectangular` or `Hann`), producing the
  RD map. "The Range-Doppler-Processing incorporates aliasing, leakage,
  limitations on resolution and accuracy errors corresponding to a real radar
  processing unit." Powers are normalised so the output "corresponds to the
  signal power at the end of the Front-End Hardware" — i.e. directly comparable
  with the radar equation.
* **OS-CFAR** (Eq. 526) on the RD map of a single virtual receiver:
  `n_lay` surrounding layers, `n_guard` guarded layers, order statistic
  `m = T_OS-CFAR · N_lay_tot`, and a scaling factor `S_OS-CFAR` [dB].
* **Angular-Processing**: spatial DFT over a uniform linear array of virtual
  receivers spaced `d = λ/2` ("to ensure an unambiguous field of view of 180°"),
  `n_ULA` samples + `n_zero_ULA` paddings, optionally in elevation as well.
* **Angular-Filtering**: cell-averaging CFAR (`n_lay`, `n_guard`, `S_CA-CFAR`).
* **Peak-Finder**: a detection is kept only if it is the local maximum inside
  `n_lay_peak` layers — applied after both CFAR stages.
* **Peak-Interpolation** (Eq. 527–529), parabolic fit through the peak bin
  `beta` and its neighbours `alpha`, `gamma`:

```
p    = 1/2 * (alpha - gamma) / (alpha - 2 beta + gamma)
y(p) = beta - 1/4 * p * (alpha - gamma)
P_peak = P_maxV * P_maxR / P_CFAR          (2-D case)
```

* **VRx mode**: no angular processing at all; the complex amplitudes of all
  virtual receivers are exposed for a custom angle estimator. Positions and
  complex weights come from the `VRxCalibration` table.

---

## 5. How detections are generated

A detection is a **cell of the data cube that survives filtering**, not an
object. One physical vehicle typically produces several detections (several
scattering centres), a guardrail produces a line of them, and multipath produces
detections at ranges/angles where nothing physically is.

Documented key effects:

* false negatives from low power and filtering;
* false positives from mirror targets, noise and aliasing;
* detection merging through finite resolution;
* micro-Doppler;
* multipath interference *increasing or decreasing* the signal strength;
* angular ambiguities determined by the virtual-receiver layout;
* high-resolution imaging resolving multiple scattering centres on large objects;
* detection of NLOS objects via multipath back-propagation.

---

## 6. Outputs

`tRadarRSI` (`Sensor_RadarRSI.h`):

* `nDetections`, `TimeFired` (timestamp of the ray-tracing job),
  `t_ext` / `rot_zyx_ext` (external sensor motion).
* Point-cloud mode — `tDetPoint[]`:
  `Coordinates[3]` (Cartesian x/y/z, or spherical r/azimuth/elevation),
  `Power` [dBm], `Velocity` [m/s] (relative radial).
* VRx mode — `tDetVRx[]`: `Range` [m], `Velocity` [m/s],
  `AmpVRx[nVRx]` complex amplitudes [mV].

Optionally the whole complex data cube is written to **HDF5**
(`Sensor.<n>.FileOutput`, `SensorCluster.<i>.HDF5OutPath`), with `config`,
`legend` and `timestamps` datasets per sensor group.

---

## 7. Important configurable parameters

Transceiver: `FoV`, `Range` (min/max ray propagation), `nRays`, `RayPattern`,
`Frequency`, `TransmitPower`, `PolarizationTransmit/Receive`, `SystemLosses`,
`Antenna.FName`.
Noise: `NoiseBandwidth` [MHz], `NoiseFigure`, `NoiseTemperature`,
`NoiseScaling`, `NoiseScalingRange`, `NoiseScalingDopplerVel`.
Discretisation: `Discretization.Range.Max/Samples/ZeroPaddings`,
`.DopplerVel.MinMax/Samples/ZeroPaddings`,
`.AzimuthAngle.Samples/ZeroPaddings`, `.ElevationAngle.*`,
`ProcessElevationAngle`, `WindowFunction`.
Filtering: `CfarFilter.RangeDoppler.{Layers,Guards,SNR,Threshold}`,
`CfarFilter.AzimuthAngle.{Layers,Guards,SNR}`, `CfarFilter.ElevationAngle.*`,
`PeakFinder.{RangeDoppler,AzimuthAngle,ElevationAngle}.Layers`.
Output: `OutputType` (Cartesian/Spherical/VRx), `nMaxDetections` (default 2000).
Dynamic pattern: `nHorizontalBuckets`, `nVerticalBuckets`,
`MaxNumberOfRaysPerBucket`, `OccupancyCheckMode`.

Shipped example (`Data/Vehicle/Examples/DemoCar_SensorRadarRSI`):
FoV 120°×20°, range 0.1…150 m, 87 230 rays (Fibonacci), 77 GHz, 30 dBm,
B = 200 MHz, NF = 5 dB, T = 300 K, RD map 512 × 256 over 0…200 m and ±40 m/s,
60 azimuth samples + 100 zero paddings, OS-CFAR (5, 2, 11 dB, 95 %),
CA-CFAR azimuth (8, 2, 12 dB), peak finder layers 2/2, Cartesian output.

---

## 8. Difference to Radar Model 1

| | Radar Sensor (Model 1) | Radar RSI (Model 2) |
|---|---|---|
| Purpose | ADAS **function** development | radar **component**/signal-processing development |
| Object model | oriented bounding box + RCS look-up table | full 3D triangle mesh + material |
| Physics | radar equation on the whole object | ray-traced polarimetric EM field |
| Output | tracked object list | point cloud / raw VRx amplitudes |
| Multipath, ghosts, micro-Doppler | not modelled (only statistical false positives) | modelled physically |
| Occlusion | bounding-box overlap fraction | ray tracing |
| Noise | thermal + road clutter, added to a single SNR | WGN in every cube cell + CFAR |
| Compute | CPU, real-time | GPU, ray tracing |
| Cycle handling | `CycleTime` + `Latency` | sensor cluster, external cycle control |

They are complementary, not alternatives.

---

## 9. What is documented vs. what we inferred

**Fully documented and reproduced**

The pipeline order and every stage's name and parameters; Equations 521–529;
the data-cube layout and interaction-point structure (from the shipped GPU
sample); the noise model; the CFAR definitions; the peak finder and the
parabolic interpolation; the ULA `d = λ/2` convention; the antenna/VRx file
format and its two different interpolation rules (linear for transmit, next
smaller sample point for receive); the ray-pattern definitions including all
dynamic-pattern rules; the output structures.

**Inferred**

| Item | What the manual says | What we implement |
|---|---|---|
| Scattering kernel | "an analytical solution to the Maxwell's equation", exact "other than the restriction to the far-field" | physical-optics ray tube: `E_s = j k (L² dΩ)/(2π) · Γ · Λ(q) · E_inc/L · e^{-jkr₂}/r₂` with the PO aperture factor `Λ(q) = sinc(q_t1 a/2)·sinc(q_t2 a/2)`. Verified to reproduce `σ_plate = 4πA²/λ²` (validation §1) |
| Roughness distribution | "material dependent stochastic component on the normal" | Gaussian tilt of the normal, σ = `Radar.Scattering` [deg] |
| Number of bounces | not stated | `max_bounces`, default 3 |
| RD-map synthesis | "two Discrete-Fourier-Transformations processing raw data" | analytic windowed-DTFT kernel per IA point ("spectral" mode, default) — mathematically the same as DFT-ing the beat signal, but far cheaper. A "binning" mode reproduces the shipped `Radar.cu` sample instead |
| Doppler beyond `DopplerVel.MinMax` | those are "the minimum and maximum unambiguous velocities" | folded back into the interval (`doppler_aliasing = True`); `Radar.cu` clamps instead |
| CFAR scaling direction | "a scaling factor is applied to the CUT" | `P_CUT > Z · 10^(S/10)` |
| Which VRx feeds the RD-CFAR | "of a single virtual receiver" | VRx index 0 |
| Splitting RD power over angular peaks | not stated | proportional to the angular spectrum |
| VRx amplitudes in millivolt | unit documented, reference impedance not | `V = A·sqrt(2·50 Ω)·1000` |
| Elevation angular processing | "the same steps as for the azimuthal angle" | 1-D azimuth only in the reference implementation; a 2-D virtual array is needed and is not part of the default parameterisation |
| Ray-pattern → FoV mapping for Fibonacci | "spherical fibonacci lattice" | equal-solid-angle in elevation, uniform in azimuth |

**Not accessible**

The GPU implementation itself (`IPGRsi.dll`, `Sensor_RadarRSI.o` in
`lib/libcarmaker.a`). The IPG EULA §4(ii) forbids reverse engineering, so **no
disassembly was performed**; everything above comes from the documentation, the
public headers, the shipped GPU Coding Interface sample source, and the ASCII
data files.
