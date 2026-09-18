# Radar Model 1 — CarMaker "Radar Sensor" (object-list HiFi radar)

CarMaker identifier: `Sensor.Param.<no>.Type = "Radar"`
C interface: `include/Vehicle/Sensor_Radar.h` (+ deprecated `Sensor_Radar_protected.h`)
Documentation: *Reference Manual → Sensors → Radar Sensor* (chapters 1060–1083 of the TOC)

---

## 1. What the model simulates

An **object-list radar**: for every relevant object in the scene the model decides
whether the radar would detect it and, if so, reports a tracked object entry
(range, angles, relative velocity, size, RCS, SNR, confidences).
It does **not** produce a point cloud and does not simulate a waveform.

> "The Radar Sensor generates a radar specific object list based on defined
> objects. It belongs to the class of HiFi Sensors and is designed especially for
> the needs of function development and testing."
> — *Sensors → Radar Sensor → Introduction*

The key design decision, stated explicitly in the manual, is that **detection is
a link-budget decision, not a geometric one**:

> "The detection of the Radar Sensor is based on a physical model. The detection
> is explicitly not based on geometrical analyses. The parameters concerning
> range are merely used to identify detection candidates and to keep the
> simulation efficient. The parameters concerning field of view are used to
> calculate the antenna gain map."
> — *Detection Based on Signal-to-Noise Ratio*, footnote to Equation 502

Intended use: ADAS/AD **function** development and testing (ACC, AEB, …), where
you need a realistic object list with realistic failure modes but not raw
signals. It runs on the CPU, in real time, at a configurable cycle time.

---

## 2. Inputs

| Input | Source in CarMaker | Notes |
|---|---|---|
| Sensor pose | `Sensor.<n>.pos`, `Sensor.<n>.rot` (z-y-x, deg), `Sensor.<n>.Mounting` | boresight = +x of the sensor frame |
| Sensor kinematics | vehicle state | needed for relative velocity |
| Object list | "Sensor Surrounding" shared dataset | traffic objects, geometry objects, trees (trunk), guide posts, guard-rail posts |
| Per-object bounding box | `Basics.Dimension = l w h` of the traffic-object template | the model's entire geometric object model |
| Per-object RCS map | `RCSMap.FName` in the traffic-object template | **an object without an RCS map is undetectable** |
| Antenna gain map | `Sensor.Param.<n>.Gain.FName` → `Data/Sensor/<file>` | 2-D map over (azimuth, elevation) |
| Environment | `Env.Temperature`, `Env.RainRate`, `Env.VisRangeInFog` | enter noise and damping |
| Sensor parameters | `Sensor.Param.<n>.*` | see §9 |

Object candidates are limited to the **observation area**: a wedge of
`FoV` h/v, between `Range_min` and `Range_max`.

---

## 3. Processing pipeline

```
 for every radar cycle (Sensor.<n>.CycleTime):

 (1) candidate selection      objects inside [Range_min, Range_max] and the FoV wedge
       |
 (2) per object geometry      bounding-box corners -> range, azimuth, elevation,
       |                      angular extent, incidence azimuth in the object frame
       |
 (3) antenna gain             G(az, el) from the gain map                      [Eq. 507/508]
       |
 (4) RCS                      RCS_lut(phi) (+ extended-object averaging)
       |                      x Swerling-1 fluctuation                          [Eq. 509]
       |                      x occlusion factor (1 - o_h' o_v')                [Eq. 510]
       |
 (5) damping                  L_atm = L_satm * L_rain * L_fog                   [Eq. 512]
       |
 (6) radar equation           S = P G^2 lambda^2 RCS / ((4pi)^3 r^4 L_A L_atm)  [Eq. 505]
       |
 (7) noise                    N = N_thermal + N_clutter(r)                      [Eq. 506]
       |
 (8) detection                SNR = S/N ; detected if SNR > SNR_min             [Eq. 502/503/504]
       |
 (9) merging                  objects closer than delta*R in ALL of
       |                      (range, velocity, azimuth) become one object      [Eq. 511/514]
       |
(10) false positives          mirror objects + Poisson clutter                   (optional)
       |
(11) measurement noise        x_hat = x + N(0, accuracy)                        [Eq. 513]
       |
(12) object management        ProbExist up/down, MeasStat, MaxNumObj, list sort
       |
(13) latency                  result becomes visible Latency[0..1] * CycleTime later
       v
     tOutQuants[] + tOutQuantsGlob
```

---

## 4. The equations (all quoted from the Reference Manual)

**Detection threshold** (Eq. 502, 503, 504)

```
detected  <=>  SNR > SNR_min
SNR       =    S / N
SNR_min   =    2 * ( erfc^-1(2 P_FA) - erfc^-1(2 P_Dmin) )^2
```

`P_FA = 10^(-ProbFalseAlarmIdx)` with `ProbFalseAlarmIdx ∈ [1,10]`,
`P_Dmin = ProbDetectMin`. With the defaults of `DemoCar_SensorRadar`
(`P_FA = 1e-6`, `P_Dmin = 0.5`) this gives `SNR_min = 22.595 = 13.54 dB`.

The output quantity `ProbDetect` is documented as "resulting from Equation 504",
i.e. the same relation solved for `P_D`:

```
P_D = 0.5 * erfc( erfc^-1(2 P_FA) - sqrt(SNR / 2) )
```

so that `P_D(SNR_min) = P_Dmin` exactly. *(The inversion itself is our reading of
the documented sentence; it is the only inversion consistent with Eq. 504.)*

**Signal strength** (Eq. 505) — the classical radar equation with `G_t = G_r`:

```
S = ( P * G^2 * lambda^2 * RCS ) / ( (4 pi)^3 * r^4 ) * 1 / ( L_A * L_atm )
```

**Thermal noise** (Eq. 506):

```
N_thermal = T_0 * F_R * k_B * B_n
```

**Antenna gain** (Eq. 507/508) — uniform rectangular aperture:

```
f(Theta, Phi) = sinc(pi nu_y) * sinc(pi nu_z)
nu_y = (a/lambda) sin(Theta) cos(Phi)
nu_z = (b/lambda) sin(Theta) sin(Phi)
```

`Theta` is the polar angle away from boresight and `Phi` the roll angle about it,
so in the azimuth/elevation frame of the manual's figure
`sin(Theta)cos(Phi) = cos(el)·sin(az)` and `sin(Theta)sin(Phi) = sin(el)`.
`a`, `b` follow from the 3 dB beam widths; the scan range steers the beam.

**RCS** (Eq. 509, 510, 511):

```
Swerling-1:   p(RCS) = (1/RCS_lut(phi)) * exp( -RCS / RCS_lut(phi) )
occlusion:    RCS     = (1 - o_h' o_v') * RCS_no
merging:      RCS     = (1/n) * sum_i RCS_i
```

> ⚠️ The printed Equation 509 in the manual omits the minus sign in the
> exponent. A Swerling-1 / Rayleigh-amplitude target requires it (otherwise the
> pdf is not normalisable), so we implement `exp(-RCS/RCS_lut)`.

**Damping** (Eq. 512): `L_atm = L_satm(f) * L_rain(f, Q_rain) * L_fog(f, d_vis)`
with the references *DESK EW Handbook*, *ITU-R P.838-3* and *Brooker*.
Coefficients are **not published** by IPG.

**Separability** (Eq. 514): `Delta_i = delta * R_i` for i ∈ {range, velocity, azimuth}.

**Measurement noise** (Eq. 513): `x_hat = x + Delta_x`, `Delta_x ~ N(0, accuracy)`,
"a value based on a confidence level of 1 sigma".

---

## 5. How targets are represented and detected

* Every object is an **oriented bounding box** plus an **RCS map**. The bounding-box
  **centre (`BBC`)** is the reference point for range/velocity/angles
  (*Object Distance*). The eight corners give the angular extent and the
  occlusion overlap (`Sensor_Radar_protected.h` documents `tPTraffic` with the
  eight vertices, the extreme left/right/top/bottom vertices, and
  `azimuth_min`/`azimuth_max` "needed for Range resolution").
* The **RCS map argument** is the azimuth of incidence *in the object's own body
  frame* (x along the object's longitudinal axis) — *Direction of Incidence*.
* Elevation dependence of RCS is **not** tabulated; it is represented by the
  Swerling-1 stochastic process (*Radar Cross Section*: "the elevation dependency
  is modeled by the stochastic process (Swerling type 1, with
  Rayleigh-Distribution)").
* **Extended objects**: because automotive resolution cells are comparable to
  object size, "the RCS map is evaluated over a range of angles depending on the
  distance, sensor resolution and orientation of the object". *(How the interval
  is constructed is not documented — see §10.)*
* Detection is then purely `SNR > SNR_min`. A "false negative" is simply an
  object whose fluctuating RCS puts it below the threshold in this cycle.

---

## 6. How measurements are computed

Ground truth relative to the sensor frame → noise added per Eq. 513 with
`AccuracyDistance`, `AccuracyAzimuth`, `AccuracySpeed` as 1σ.

Reported per object (`tOutQuants`):
`Dist, DistX, DistY, DistZ, theta, Vrel, VrelX, VrelY, ArelX, RelCourseAngle`,
plus `Length/Width` (+ classes 0…7), `RCS [dBm²]`, `SignalStrength [dBW]`,
`SNR [dB]`, `ProbDetect`, `ProbExist` (0…7), `ProbObst`, `DynProp`, `MeasStat`,
`GroupId`, `LeftInGroup`, `ObjId`.

`ProbObst = 1 − g / 1 m`, where `g` is the gap between the object bounding box
and the ego driving path (ego width + 0.3 m clearance, current trajectory
curvature radius) — *Object Quantities*.

`ProbExist` is an integer 0…7 mapped to 0 %, 25 %, 50 %, 75 %, 90 %, 99 %,
99.9 %, 100 %. It is initialised from the object's RCS-map `ProbExist` field,
incremented while the object is detected and decremented otherwise; the object
stays on the list until it reaches 0 — this is the model's **tracking memory**.

---

## 7. Modelled error sources

| Effect | Documented? | How |
|---|---|---|
| False negatives | yes | SNR below threshold (incl. Swerling fading) |
| Measurement noise | yes | Gaussian, 1σ = accuracy parameters |
| Latency + cycle time | yes | output age between `Δt_latency` and `Δt_latency + T_radar` |
| Object merging | yes | separability × resolution in range **and** velocity **and** azimuth |
| Occlusion | yes | bounding-box overlap × occluder `OcclusionFactor`, multiple occluders |
| Atmospheric / rain / fog damping | yes (references only) | `L_atm` |
| Road clutter | qualitatively | range-dependent noise added to thermal noise |
| False positives — mirror objects | qualitatively | reflection via guardrail/tunnel wall; higher probability for trucks; depends on sensor–object and object–wall distance; ID = `−ObjId`; objects off road/RoadMargin excluded |
| False positives — clutter | yes | count ~ Poisson(`ClutterObjMean`), distance and dimensions uniform; ID = `−2` |
| Multipath (as signal interference) | **no** | not part of this model — use Radar RSI |
| Micro-Doppler, ghosting geometry | **no** | not part of this model |

---

## 8. Outputs

`tRadarSensor` (`Sensor_Radar.h`):

* `ObjList[]` of `tOutQuants` — one entry per detected/held object (§6).
* `GlobalInf` of `tOutQuantsGlob` — `RolCount`, `nObj`, `RelvTgt`,
  `nLanesL/R`, `DistToLeft/RightBorder`, and the DVA-writable sensor parameters
  (`TransmitPower`, `NoiseBandWidth`, `NoiseFigure`, `Separability`, the three
  resolutions, the three accuracies, `Range_min`, `Range_max`).

---

## 9. Important configurable parameters

Observation area: `FoV`, `Range_min`, `Range_max`, `MaxNumObj` (≤ 2000).
Front end: `Frequency` [GHz], `TransmitPower` [dBm], `SystemLosses` [dB],
`NoiseBandWidth` [Hz], `NoiseFigure` [dB].
Threshold: `ProbDetectMin`, `ProbFalseAlarmIdx` (1…10 → 1e-1…1e-10).
Resolution/accuracy: `Resolution{Distance,Azimuth,Speed}`,
`Accuracy{Distance,Azimuth,Speed}`, `Separability`.
False positives: `FalsePosActive`, `ClutterObjMean` (default 5).
Antenna: `Gain.FName` → file with `ScanRange`, `BeamWidth`, `AntennaEff`,
`Azimuth`, `Elevation`, `Gain`.
Timing: `Sensor.<n>.CycleTime`, `CycleOffset`, `Latency` (× CycleTime).
Reproducibility: Test Run parameter `RandomSeed`.

Shipped example (`Data/Vehicle/Examples/DemoCar_SensorRadar`):
FoV 40°×30°, range 0.2…200 m, 77 GHz, 14 dBm, B = 25 kHz, NF = 4.8 dB,
`ProbDetectMin = 0.5`, `ProbFalseAlarmIdx = 6`, resolutions 1.8 m / 1.6° / 0.4 km/h,
accuracies 0.4 m / 0.1° / 0.1 km/h, `Separability = 1.5`, cycle 60 ms,
`Gain.FName = Radar_Default`.

---

## 10. What is documented vs. what we inferred

**Fully documented and reproduced exactly**

Equations 502–514, the `ProbFalseAlarmIdx` mapping, the antenna file and RCS
file formats, the object list quantities and their units, the length/width class
table, the `ProbExist` semantics, the object-ID conventions for false positives,
the latency/cycle-time behaviour, and the parameter list.

**Documented qualitatively — our implementation is an interpretation**

| Item | What the manual says | What we implement |
|---|---|---|
| Road clutter `N_clutter(r)` | uses "antenna gain G, transmitted power P, reflectivity of the street σ₀ as well as the sensor resolution", "range dependent", added to thermal noise per object | textbook beam-limited surface clutter: `A_c = r·Θ_az·ΔR/cos ψ`, `σ_c = σ₀·A_c`, then Eq. 505. `σ₀` is a config parameter (default −35 dB); **IPG does not publish it** |
| Extended-object RCS | "RCS map is evaluated over a range of angles depending on the distance, sensor resolution and orientation" | average the *linear* RCS over the incidence-azimuth interval spanned by the bounding box, widened to at least `ResolutionAzimuth` |
| `o_h'`, `o_v'` in Eq. 510 | primes are never defined; "evaluated conservatively"; multiple occluders and transparency considered | `o = coverage_fraction × occluder.OcclusionFactor`, combined across occluders by `max` |
| "RCS corrected for the total noise" | one sentence, no formula | `RCS_reported = RCS·(1 + 1/SNR)` (the receiver measures S+N) |
| `L_satm`, `L_rain`, `L_fog` | literature references only | published functional forms with configurable coefficients (77 GHz defaults); **absolute damping will not match IPG** |
| Mirror false positives | probability grows for trucks, falls with sensor–object and object–wall distance | `p = truck_boost · exp(−d_wall/5 m) · exp(−r/Range_max)`; wall objects must be tagged by the caller |
| `RelvTgt` | "dep. on Dist and ObstacleProbability" | nearest object with `ProbObst > 0` |
| Object-list ordering | not specified | sorted by range |
| Antenna obliquity factor | not mentioned | `(1 + cos Θ)/2` added to Eq. 507 — see the validation report; it reduces the mismatch against the shipped `Radar_Default` map from 1.17 dB to 0.12 dB mean |
| `DynProp` stationary vs. stopped | "based on velocity direction of target and sensor" | single-frame classification only; CarMaker evidently keeps history |
| `MeasStat = 2` | header says "object not measured", Reference Manual says "currently not used" | we use it for objects held by `ProbExist` |

**Not accessible at all**

The compiled model itself (`Sensor_Radar_Calc.o` inside `lib/libcarmaker.a`).
The IPG EULA §4(ii) forbids reverse engineering, so **no disassembly was
performed**; everything above comes from the documentation, the public C
headers, and the shipped ASCII data files.
