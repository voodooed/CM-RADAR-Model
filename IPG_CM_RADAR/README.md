# IPG CarMaker radar models — analysis and independent reference implementation

Independent, CarMaker-free Python reference implementations of the **two** radar
sensor models in IPG CarMaker 15.1, together with the analysis they are based on
and the validation that was possible without running CarMaker.

Analysed against a **CarMaker 15.1** installation (`carmaker/win64-15.1`).
All paths in this package are relative to its own root; the only external
path you may need is your CarMaker installation directory — see *Running it*.

> **New to these models? Start with
> [`docs/IPG_radar_models.md`](docs/IPG_radar_models.md).**
> It builds up what a radar does from scratch, then explains both CarMaker models
> in plain terms — no CarMaker knowledge assumed. ~25 minutes. Read it before the
> Reference Manual, which is a lookup reference and assumes the background.

---

## The two radar models

CarMaker 15.1 contains exactly two radar models — confirmed by the headers in
`include/Vehicle/`, the Reference Manual table of contents, and the files in
`Data/Sensor` (see `docs/source_inventory.md`). They are **not** variants of one
another; they sit at different integration levels and are documented, configured
and consumed separately.

| | **Model 1 — Radar Sensor** | **Model 2 — Radar RSI** |
|---|---|---|
| CarMaker type | `Sensor.Param.<n>.Type = "Radar"` | `… = "RadarRSI"` |
| Class | HiFi sensor | Raw Signal Interface |
| Purpose | ADAS/AD **function** development & testing | radar **component** / signal-processing development |
| Object model | oriented bounding box + azimuth RCS look-up table | full 3D triangle mesh + per-face material |
| Physics | radar equation per object, SNR threshold | ray-traced polarimetric EM propagation |
| Output | tracked **object list** | **point cloud** or raw virtual-receiver amplitudes |
| Multipath / ghosts / micro-Doppler | not modelled (statistical false positives only) | modelled physically |
| Compute | CPU, real time | GPU ray tracing |
| C header | `Sensor_Radar.h` | `Sensor_RadarRSI.h` |
| Our implementation | `radar_object_list/` | `radar_rsi/` |

---

## Directory layout

```
IPG_CM_RADAR/
├── README.md                       this file
├── docs/
│   ├── IPG_radar_models.md         << START HERE: radar 101 + both models
│   └── source_inventory.md         every IPG file/section used, and the method's limits
│
├── common/                         shared, dependency-free helpers
│   ├── infofile.py                 reader for IPG "Info File" format
│   ├── frames.py                   frames, z-y-x rotations, bounding boxes, az/el
│   ├── paths.py                    locating a CarMaker install ($CARMAKER_DIR)
│   └── rng.py                      deterministic RNG (normal / exponential / Poisson)
│
├── radar_object_list/              MODEL 1 — CarMaker "Radar Sensor"
│   ├── explanation.md              how the model works (concise, with citations)
│   ├── implementation_logic.md     implementation specification
│   ├── config.py                   every Sensor.Param.<n>.* key, + inferred parameters
│   ├── maps.py                     antenna gain map + RCS map (load IPG files / generate)
│   ├── atmosphere.py               L_satm, L_rain, L_fog  (Eq. 512)
│   ├── targets.py                  target representation, geometry, occlusion
│   ├── detection.py                Eq. 502-506, erfcinv, clutter noise
│   ├── model.py                    the model: RadarSensorModel.step()
│   └── example.py                  runnable demo
│
├── radar_rsi/                      MODEL 2 — CarMaker "Radar RSI"
│   ├── explanation.md
│   ├── implementation_logic.md
│   ├── config.py                   every Sensor.Param.<n>.* key, + inferred parameters
│   ├── scene.py                    materials (MaterialLib), triangle scene, primitives
│   ├── raypattern.py               PhiTheta / Fibonacci / Dynamic ray patterns
│   ├── transceiver.py              Tx/Rx gain maps + VRx array (IPG file / generated)
│   ├── propagation.py              shooting-and-bouncing rays -> interaction points
│   ├── device.py                   front end + back end (DFT, OS-CFAR, CA-CFAR, peaks)
│   ├── model.py                    the model: RadarRSIModel.calculate()
│   └── example.py                  runnable demo
│
└── validation/
    ├── validation_report.md        what could be validated, what could not, results
    ├── run_all.sh                  run everything, write validation/results/
    ├── validate_ipg_files.py       read every shipped radar data file
    ├── validate_antenna_map.py     internal antenna model vs. IPG's Radar_Default
    ├── validate_detection_math.py  Eq. 502-506 and the erfcinv table
    ├── validate_object_list_model.py  documented Model 1 behaviour
    ├── validate_rsi_device.py      Eq. 525-529, CFAR, DFT, angular processing
    ├── validate_rsi_propagation.py Model 2 physics vs. closed-form results
    └── results/                    captured output of all of the above
```

Read in this order:
`docs/IPG_radar_models.md` (the primer) →
`radar_object_list/explanation.md` → `radar_rsi/explanation.md` →
the two `implementation_logic.md` → `validation/validation_report.md`.
`docs/source_inventory.md` is a lookup: it maps every claim back to the IPG file
or manual section it came from.

---

## Running it

Requirements: Python ≥ 3.10. `common/`, `radar_object_list/` and
`validation/validate_detection_math.py` need **no third-party packages**;
`radar_rsi/` needs **numpy** (`pip install numpy`).

A CarMaker installation is **optional**. Without one everything still runs,
using the internal antenna model and a synthetic material/RCS set. With one, the
models are driven by IPG's own data files. Point at it either way:

```bash
# from the root of this package
cd <path-to>/IPG_CM_RADAR

# optional: your CarMaker installation - the directory containing Data/, include/, doc/
export CARMAKER_DIR=<path-to>/carmaker/<arch>-<version>     # e.g. .../carmaker/win64-15.1

# Model 1 demo, using CarMaker's own antenna and RCS maps
python3 radar_object_list/example.py --ipg "$CARMAKER_DIR"

# Model 2 demo, using CarMaker's own transceiver config and MaterialLib
python3 radar_rsi/example.py --ipg "$CARMAKER_DIR" --rays 4000

# everything, with results written to validation/results/
./validation/run_all.sh                 # picks up $CARMAKER_DIR
./validation/run_all.sh <install-dir>   # or pass it explicitly
```

Every script accepts `--ipg <dir>`, falls back to `$CARMAKER_DIR`, and skips
itself with a clear message if neither is available. Run scripts from the
package root (`python3 radar_object_list/example.py`, not from inside the
subdirectory) so the local imports resolve.

### Minimal usage

```python
from common.frames import Pose
from common.paths import default_carmaker_dir
from radar_object_list.config import RadarSensorConfig
from radar_object_list.maps import AntennaGainMap, RcsMap
from radar_object_list.model import RadarSensorModel
from radar_object_list.targets import Target

IPG     = default_carmaker_dir()                    # $CARMAKER_DIR, or hard-code a path
cfg     = RadarSensorConfig()                       # = DemoCar_SensorRadar
antenna = AntennaGainMap.load(f"{IPG}/Data/Sensor/Radar_Default")
rcs_car = RcsMap.load(f"{IPG}/Data/Sensor/RCS_Car", "RCS_Car")
model   = RadarSensorModel(cfg, antenna=antenna)

pose    = Pose.from_pos_rot_deg(cfg.pos, cfg.rot_deg)   # sensor pose in the world
targets = [Target(obj_id=16000001, pos=(64.2, 0.0, 0.75), vel=(22.0, 0.0, 0.0),
                  length=4.8, width=1.8, height=1.5, rcs_map=rcs_car)]

# step() applies CycleTime and Latency, so the first cycle becomes visible one
# latency later; call it every simulation step.  Use calculate_now() to run a
# single cycle immediately (tests, single-scene evaluation).
for k in range(20):
    out = model.step(t=k * 0.01, sensor_pose_world=pose,
                     sensor_vel_world=(25.0, 0, 0), targets=targets,
                     ego_speed=25.0)
for o in out.objects:
    print(o.obj_id, o.dist, o.vrel, o.rcs, o.snr, o.prob_exist)
```

```python
from common.frames import Pose
from common.paths import default_carmaker_dir
from radar_rsi.config import RadarRSIConfig
from radar_rsi.model import RadarRSIModel
from radar_rsi.scene import MaterialLib, Scene, make_box, make_plane
from radar_rsi.transceiver import TransceiverConfig

IPG   = default_carmaker_dir()
cfg   = RadarRSIConfig(fov_deg=(40, 10), n_rays=4000, range_samples=128,
                       doppler_samples=64, azimuth_samples=16)
tc    = TransceiverConfig.load(f"{IPG}/Data/Sensor/RadarRSI_Default")
mats  = MaterialLib.load(f"{IPG}/Data/Sensor/MaterialLib")
model = RadarRSIModel(cfg, transceiver=tc)

scene = Scene()
scene.add_body(make_plane(mats.get("asphalt")))
car = make_box(mats.get("metal"), 4.5, 1.8, 1.4, subdiv=3)
car.pose = Pose((45.0, 0.0, 0.7)); car.velocity = (-5.0, 0.0, 0.0)
scene.add_body(car)

out = model.calculate(0.0, scene, Pose.from_pos_rot_deg(cfg.pos, cfg.rot_deg))
for p in out.det_points:
    print(p.coordinates, p.power_dbm, p.velocity)
```

---

## Integration with a custom simulator

This package has no dependency on CarMaker at runtime and no dependency on any
particular host simulator. Both models are plain computational kernels: they
consume a description of the scene at a point in time and return sensor output.
This section specifies the interface contract a host simulator must satisfy.

### 1. Runtime requirements

| Component | Requirement |
|---|---|
| `common/`, `radar_object_list/` (Model 1) | Python ≥ 3.10, **standard library only** |
| `radar_rsi/` (Model 2) | Python ≥ 3.10 + **numpy** |
| Everything else | none — no compiler, no CarMaker libraries, no GPU, no build step |

The absence of third-party dependencies in Model 1 is verified: the model
executes correctly with `numpy` made unimportable.

A CarMaker installation is **not** required at runtime. It is used only as an
optional source of sensor characterisation data (§4).

### 2. Common integration contract

Both models are instantiated once per physical sensor and driven once per
sensor cycle. The host simulator is responsible for supplying, in the units and
frames below, the state of the sensor and of the scene.

#### 2.1 Coordinate frames

All frames are **right-handed** and follow the CarMaker convention:

```
   z (up)
   |
   |____ y (left)
  /
 x (forward)
```

| Quantity | Frame | Definition |
|---|---|---|
| World frame | — | any fixed right-handed frame with **z up** |
| Sensor frame | world | origin at the transceiver; **boresight is +x**, y left, z up |
| Object body frame | world | origin at the **bounding-box centre**; +x along the object's longitudinal axis |
| Azimuth `φ` | sensor | `atan2(d_y, d_x)` — zero on boresight, positive to the left |
| Elevation `θ` | sensor | `atan2(d_z, hypot(d_x, d_y))` — zero in the x-y plane, positive up |

Rotations supplied as Euler angles are applied in **z-y-x order**:
`R = Rz(rz) · Ry(ry) · Rx(rx)`. `Pose.from_pos_rot_deg(pos, (rx, ry, rz))`
constructs a pose from a position and such a triple in **degrees**.

> If the host simulator uses a different convention (for example z-forward /
> y-up, or a left-handed frame), the adapter — not the model — is responsible
> for the conversion. Supplying a mirrored frame produces mirrored azimuths,
> which is the single most common integration defect.

#### 2.2 Units

| Quantity | Unit |
|---|---|
| Position, dimensions, range | m |
| Velocity | m/s |
| Acceleration | m/s² |
| Simulation time `t` | s (monotonically non-decreasing) |
| Angles in `*Config` objects and data files | deg |
| Angles in all runtime interfaces and outputs | rad |
| RCS supplied to `RcsMap` | **m² (linear, not dB)** |
| Power outputs | dBm² (RCS), dBW (signal strength), dB (SNR), dBm (Model 2 detections) |

#### 2.3 Sensor pose and velocity

The models take the pose of the **sensor**, not of the carrying vehicle. The
adapter must compose the vehicle pose with the mounting pose, and must account
for the lever arm when the vehicle rotates:

```
sensor_pose  = vehicle_pose.compose(mount_pose)
r_lever      = vehicle_pose.dir_to_parent(mount_pose.t)
sensor_vel   = vehicle_velocity + cross(vehicle_angular_velocity, r_lever)
```

Using the vehicle's velocity directly for a sensor mounted away from the centre
of rotation introduces a radial-velocity error proportional to the yaw rate.

### 3. Model-specific contracts

#### 3.1 Model 1 — object-list radar

Entry point:

```python
RadarSensorModel.step(
    t,                    # float  [s]
    sensor_pose_world,    # Pose   world <- sensor
    sensor_vel_world,     # (vx, vy, vz)  [m/s], world frame
    targets,              # Sequence[Target]
    ego_speed=0.0,        # float  [m/s]   - used for DynProp and ProbObst
    ego_yaw_rate=0.0,     # float  [rad/s] - driving-path curvature for ProbObst
    ego_width=1.8,        # float  [m]     - driving-path width for ProbObst
    road_info=None,       # optional dict, see below
) -> RadarOutput
```

Each element of `targets` is a `Target` with the following required fields:

| Field | Type | Meaning |
|---|---|---|
| `obj_id` | `int` | stable, unique identifier; used for tracking continuity across cycles |
| `pos` | `(x, y, z)` | **bounding-box centre**, world frame [m] |
| `vel` | `(vx, vy, vz)` | world frame [m/s] |
| `acc` | `(ax, ay, az)` | world frame [m/s²]; only `ArelX` is reported |
| `yaw`, `pitch`, `roll` | `float` | orientation [rad], applied z-y-x |
| `length`, `width`, `height` | `float` | bounding-box extents [m] |
| `rcs_map` | `RcsMap` | **mandatory** — a target with `rcs_map=None` is undetectable by design |
| `detect_mask` | `bool` | `False` suppresses the target entirely (CarMaker's `DetectMask`) |
| `name` | `str` | optional; values `"guardrail"`, `"wall"`, `"tunnel"` additionally mark an object as a mirror surface for false-positive generation |

`obj_id` **must be stable across cycles**. The model maintains per-object
tracking state (`ProbExist`, `MeasStat`) keyed on it; regenerating identifiers
every frame disables the tracking behaviour and produces only `MeasStat = 1`.

`road_info` is an optional dictionary carrying road-derived quantities the model
does not compute itself:
`{"nLanesL", "nLanesR", "DistToLeftBorder", "DistToRightBorder"}`.
Omitting it leaves those output fields at zero and affects nothing else.

`RadarSensorModel(cfg, antenna=None, road_z=0.0)` — `road_z` is the road surface
height in the world frame, used for the grazing angle of the road-clutter term.

#### 3.2 Model 2 — raw-signal radar

Entry point:

```python
RadarRSIModel.calculate(
    t,                    # float  [s]
    scene,                # Scene
    sensor_pose,          # Pose   world <- sensor
    sensor_vel=(0, 0, 0), # (vx, vy, vz)  [m/s], world frame
) -> RadarRSIOutput
```

The host must supply **geometry**, not bounding boxes. A `Scene` is a collection
of `Body` objects; each `Body` carries triangles expressed in its **own body
frame**, plus the pose and motion of that frame:

```python
Body(name, pose, velocity, angular_velocity, triangles)
Triangle(v0, v1, v2, material)         # vertices in the body frame
RadarMaterial(name, permittivity, scattering_deg)
```

Requirements on the geometry, following CarMaker's documented *Object Model
Requirements*:

1. **Outward normals.** The normal is derived from the winding order as
   `normalize(cross(v1 − v0, v2 − v0))` and must point out of a closed surface.
2. **Sufficient tessellation.** Faces must be small enough to represent surface
   curvature and to form the corner reflectors that dominate a vehicle's
   signature. Geometric detail represented only by textures does not scatter.
3. **A material on every face.** `permittivity` drives the Fresnel coefficient
   (use ≥ 1e9 for metal); `scattering_deg` is the standard deviation of the
   stochastic surface-normal perturbation that models roughness.
4. **Per-body motion.** `velocity` and `angular_velocity` are in the world
   frame. Bodies moving relative to one another produce micro-Doppler; a single
   rigid body produces a single Doppler line.

`calculate()` invokes `scene.build()` internally, so an adapter may retain one
`Scene` instance and update `body.pose` / `body.velocity` in place each cycle.

There is currently **no mesh importer**. Geometry must be constructed through
`make_box()` / `make_plane()` or by emitting `Triangle` objects directly from
the host's mesh representation.

### 4. Required sensor characterisation data

Independently of the host simulator, three data sets characterise the sensor.
Each may be loaded from a CarMaker installation, authored by the integrator, or
generated internally.

| Data set | Used by | Source options |
|---|---|---|
| Antenna gain map, `G(φ, θ)` [dB] | Model 1 | `AntennaGainMap.load(<file>)` (CarMaker `Data/Sensor/Radar_Default`), or `AntennaGainMap.generate(beam_width_deg, scan_range_deg, antenna_eff)` |
| RCS map per object class, `RCS(φ)` [m²] | Model 1 | `RcsMap.load(<file>)` (CarMaker `RCS_Car`, `RCS_Truck`, `RCS_Pedestrian`, `RCS_Bicycle`, `RCS_GuardRailPost`), or `RcsMap(name, azim_rad, rcs_linear, prob_exist, occlusion_factor)` |
| Transceiver configuration (Tx/Rx gain maps + virtual-receiver array) | Model 2 | `TransceiverConfig.load(<file>)` (CarMaker `RadarRSI_Default`), or `TransceiverConfig.generate(...)` |
| Material table | Model 2 | `MaterialLib.load(<file>)` (CarMaker `MaterialLib`), or `MaterialLib.default()` |

**The RCS map is the only genuinely new artefact an integrator must provide for
Model 1.** It is a one-dimensional table of RCS against incidence azimuth in the
object's body frame, together with two scalars: the initial probability of
existence (0…7) and an occlusion factor (0 = transparent, 1 = opaque). The file
format is documented in `docs/source_inventory.md` and parsed by
`common/infofile.py`.

Both models are fully functional without any CarMaker data: the internal antenna
model and `MaterialLib.default()` provide physically reasonable substitutes, at
the cost of no longer reproducing IPG's measured characteristics.

### 5. Reference adapter — Model 1

```python
from common.frames import Pose, rot_zyx, vadd, vcross
from radar_object_list.config import RadarSensorConfig
from radar_object_list.maps import AntennaGainMap, RcsMap
from radar_object_list.model import RadarSensorModel
from radar_object_list.targets import Target


class ObjectListRadarAdapter:
    """Binds one RadarSensorModel instance to a host simulator."""

    def __init__(self, cfg: RadarSensorConfig, antenna: AntennaGainMap,
                 rcs_by_class: dict[str, RcsMap], road_z: float = 0.0):
        self.model = RadarSensorModel(cfg, antenna=antenna, road_z=road_z)
        self.rcs_by_class = rcs_by_class
        # sensor pose relative to the vehicle body frame
        self.mount = Pose.from_pos_rot_deg(cfg.pos, cfg.rot_deg)
        self.model.set_initial_prob_exist({})   # optionally per obj_id

    def tick(self, t, ego, actors):
        """ego, actors: host-simulator objects. Returns RadarOutput."""
        vehicle_pose = Pose(ego.position, rot_zyx(ego.roll, ego.pitch, ego.yaw))
        sensor_pose = vehicle_pose.compose(self.mount)

        # sensor-origin velocity including the lever arm
        r_lever = vehicle_pose.dir_to_parent(self.mount.t)
        sensor_vel = vadd(ego.velocity, vcross(ego.angular_velocity, r_lever))

        targets = [
            Target(obj_id=a.id,
                   pos=a.bbox_center, vel=a.velocity, acc=a.acceleration,
                   yaw=a.yaw, pitch=a.pitch, roll=a.roll,
                   length=a.length, width=a.width, height=a.height,
                   rcs_map=self.rcs_by_class.get(a.object_class),
                   detect_mask=a.detectable,
                   name=a.semantic_tag)
            for a in actors
        ]

        return self.model.step(
            t, sensor_pose, sensor_vel, targets,
            ego_speed=ego.speed, ego_yaw_rate=ego.yaw_rate,
            ego_width=ego.width)
```

Actors whose `object_class` has no RCS map yield `rcs_map=None` and are silently
undetectable — the documented CarMaker behaviour. Log this during bring-up.

### 6. Reference adapter — Model 2

```python
from common.frames import Pose, rot_zyx
from radar_rsi.model import RadarRSIModel
from radar_rsi.scene import Body, Scene, Triangle


class RawSignalRadarAdapter:
    def __init__(self, cfg, transceiver, material_lib, static_bodies):
        self.model = RadarRSIModel(cfg, transceiver=transceiver)
        self.mats = material_lib
        self.mount = Pose.from_pos_rot_deg(cfg.pos, cfg.rot_deg)
        self.scene = Scene()
        for b in static_bodies:            # road, buildings, barriers
            self.scene.add_body(b)
        self.dynamic: dict[int, Body] = {}

    def _body_for(self, actor):
        body = self.dynamic.get(actor.id)
        if body is None:
            body = Body(name=str(actor.id),
                        triangles=[Triangle(v0, v1, v2, self.mats.get(mat))
                                   for (v0, v1, v2, mat) in actor.mesh_faces()])
            self.dynamic[actor.id] = body
            self.scene.add_body(body)
        return body

    def tick(self, t, ego, actors):
        for a in actors:                    # update poses in place
            body = self._body_for(a)
            body.pose = Pose(a.position, rot_zyx(a.roll, a.pitch, a.yaw))
            body.velocity = a.velocity
            body.angular_velocity = a.angular_velocity

        vehicle_pose = Pose(ego.position, rot_zyx(ego.roll, ego.pitch, ego.yaw))
        sensor_pose = vehicle_pose.compose(self.mount)
        return self.model.calculate(t, self.scene, sensor_pose, ego.velocity)
```

Mesh extraction (`actor.mesh_faces()`) and the mapping from host material names
to `RadarMaterial` are host-specific and are the bulk of the integration effort.

### 7. Timing semantics

**Model 1** implements CarMaker's cycle-time and latency behaviour internally.
Call `step()` **every simulation step**, not every radar cycle: the model
decides when to recompute (`CycleTime`) and when the result becomes observable
(`Latency`). Output age therefore varies between `Latency` and
`Latency + CycleTime`, as in the reference. Use `calculate_now()` to bypass this
for unit tests or single-frame evaluation.

**Model 2** has **no internal cycle or latency handling**. `calculate()`
executes one complete sensor cycle on every call; the host is responsible for
invoking it at the intended rate and for applying any latency.

Both models are deterministic for a fixed `random_seed`, scenario and parameter
set. Instantiate one model object per physical sensor; instances hold
independent state and independent random streams.

### 8. Computational cost

Measured single-core, CPython 3.14, on the development machine. These are
indicative, not benchmarks.

| Model 1 (targets in the observation area) | per radar cycle |
|---|---|
| 1 target | 0.13 ms |
| 10 targets | 1.3 ms |
| 50 targets | 6.6 ms |
| 200 targets | 42 ms |

| Model 2 (110-triangle scene) | per sensor cycle |
|---|---|
| 2 000 rays | 0.75 s |
| 4 000 rays | 1.3 s |
| 8 000 rays | 2.6 s |

The occlusion and merging stages are pairwise, i.e. `O(n²)` in the number of
detected objects, but at realistic object counts the per-object link budget
dominates: measured scaling from 10 to 200 targets is close to linear
(exponent ≈ 1.2). Model 1 is therefore suitable for faster-than-real-time batch
execution; a 60 ms radar cycle with 50 targets consumes ~11 % of one core. Model 2 is **not real-time** in Python: cost scales
with `rays × triangles × bounces`, and the reference implementation uses
brute-force ray/triangle intersection. Production use requires the C++ port
described in `radar_rsi/implementation_logic.md` §10, with a BVH or hardware ray
tracing replacing `intersect_batch`.

### 9. Integration checklist

1. Confirm the host frame is right-handed with **z up**; write the conversion in
   the adapter if not.
2. Confirm the sensor boresight maps to **+x** of the sensor frame.
3. Verify a static target directly ahead returns `azimuth ≈ 0` and
   `DistY ≈ 0`. A sign error here mirrors the entire scene.
4. Verify a target approaching head-on returns a **negative** `Vrel`.
5. Verify `obj_id` is stable across cycles — check that `MeasStat` transitions
   `1 → 3` rather than remaining `1`.
6. Verify the lever-arm term: with a stationary vehicle yawing in place, a
   laterally offset sensor must report non-zero radial velocity for static
   objects.
7. Confirm every object class of interest has an RCS map assigned.
8. Compare the reported detection range for a known RCS against
   `RadarSensorModel.detection_range(rcs_m2)`.
9. Fix `random_seed` and confirm run-to-run reproducibility.

### 10. Constraints

* Road-derived outputs (`nLanesL/R`, `DistToLeft/RightBorder`) are not computed;
  supply them via `road_info` or ignore them.
* Model 1 represents objects as bounding boxes only and reports the
  bounding-box centre. Model 2 reports scattering surfaces; a detection from
  Model 2 lies on the nearest illuminated face, not at the object centre.
* Model 2 elevation angle processing is not implemented; the elevation
  coordinate of every detection is zero.
* Model 2 has no mesh importer and no HDF5 raw-cube export.
* Runtime changes of sensor parameters (CarMaker's Direct Variable Access) are
  not supported; configuration is read at construction.

---

## Validation summary

Full detail and raw numbers: `validation/validation_report.md` and
`validation/results/`.

* **All 9 radar-relevant CarMaker data files parse unmodified** — both models can
  be driven with IPG's own parameterisation.
* **Antenna model vs. IPG's shipped `Radar_Default` map**: null positions match to
  0.2°, main lobe to **0.012 dB mean**, whole map above −25 dB to **0.120 dB
  mean**. This confirms our reading of Equations 507/508; it also revealed an
  undocumented obliquity factor (see limitations).
* **Documented equations** (502–506, 509–514, 525–529) reproduced to machine
  precision, including the `ProbFalseAlarmIdx` → `erfc⁻¹` table and the
  self-consistency of Eq. 504 with its inversion `P_D(SNR_min) = P_Dmin`.
* **Documented behaviour of Model 1** — Swerling-1 statistics, measurement-noise
  σ, the separability merge boundary, Eq. 511 merged RCS, `ProbExist`
  hysteresis, latency, Poisson clutter with ID −2, mirror objects with ID
  −`ObjId`, the class table — all verified.
* **Model 2 physics** — reproduces the physical-optics plate RCS
  `σ = 4πA²/λ²` to **0.76 dB mean**, the 1/r⁴ law to 0.1 dB, Doppler exactly,
  and produces the documented road-multipath interference fringes (±2…4 dB).
* **Model 2 device model** — noise power exact, amplitude normalisation exact,
  the documented window trade-off (rect −16.9 dBc leakage / Hann −40.3 dBc),
  Doppler aliasing folding, monotone OS-CFAR false-alarm behaviour, exact
  parabolic peak interpolation, and azimuth estimation to ≤ 0.005°.

---

## Limitations and assumptions

**This is not a bit-exact clone of CarMaker, and cannot be.** The models ship
only as compiled code, the IPG EULA §4(ii) forbids reverse engineering (no
decompilation or disassembly was performed here, and nothing in this package is
derived from the binaries), and CarMaker was not executed to produce reference
logs. What is reproduced is the **documented algorithm**; where the
documentation stops, the gap is filled by an inference that is flagged in place.

Everything inferred is marked `[INFERRED]` in the code and listed in
§10 of `radar_object_list/explanation.md` and §9 of `radar_rsi/explanation.md`.
The ones that matter most:

1. **Road clutter noise** (Model 1) — documented qualitatively only. We use the
   textbook beam-limited surface-clutter model; `clutter_sigma0_db` (default
   −35 dB) is **not** an IPG value and dominates the detection range of small
   targets beyond ~100 m.
2. **Atmospheric / rain / fog damping** — IPG cites DESK, ITU-R P.838-3 and
   Brooker but publishes no coefficients. Ours are public 77 GHz literature
   values; absolute damping will differ.
3. **Antenna obliquity factor** `(1 + cos Θ)/2` (Model 1) — not documented;
   inferred from the residual against `Radar_Default`. Switchable. Loading the
   IPG map directly avoids the issue entirely.
4. **Extended-object RCS averaging**, the **RCS noise correction**, the
   **occlusion combination rule**, the **mirror-target probability law**, the
   object-list ordering and `RelvTgt` tie-breaking (Model 1) are interpretations
   of one-sentence descriptions.
5. **The PO scattering kernel** (Model 2) is our choice of "analytical solution
   to the Maxwell's equation"; it is validated against the plate RCS, but IPG's
   formulation is unknown. Bounce count and energy cut-off are free parameters.
6. **RD-map synthesis** (Model 2) uses the analytic windowed-DTFT kernel rather
   than generating and transforming beat-signal samples — equivalent for an
   ideal FMCW chain, but it does not model ADC quantisation or other front-end
   non-idealities IPG may include. A `"binning"` mode reproduces the shipped
   `Radar.cu` GPU sample instead.
7. **Elevation angular processing** (Model 2) is not implemented (needs a 2-D
   virtual array; off in the shipped example parameterisation).
8. **Random streams** differ from IPG's for the same seed; only the
   distributions match, so per-cycle noise realisations will never coincide.
9. Two apparent **typos in the CarMaker 15.1 Reference Manual** were found and
   are implemented in their physically consistent form: Equation 509 is printed
   without the minus sign in the exponential, and Equation 524 uses `theta`
   where `phi` belongs.
10. Road-dependent outputs (`nLanesL/R`, `DistToLeft/RightBorder`) need a road
    model and are taken as optional caller-supplied inputs.

## Design for a later C++/CARLA port

* Core algorithm, I/O, configuration and demo code are in separate modules.
  `RadarSensorModel._calculate` and `RadarRSIModel.calculate` are pure functions
  of (time, pose, scene/targets, model state).
* The algorithms use plain loops, plain floats and small fixed-size vectors;
  `dataclass` records map to `struct`, `dict` tracks to
  `std::unordered_map`, `common.rng.Rng` to `std::mt19937_64` plus the standard
  distributions.
* `radar_rsi/propagation.py` is written as *array of rays × loop over bounces* —
  the same shape as the GPU kernels CarMaker uses; only `intersect_batch` needs
  replacing by a BVH/Embree/OptiX query or CARLA's ray cast.
* The interaction-point and data-cube layouts deliberately match
  `tUserRadarInteractionPoint` and `tUserRadarCube` of CarMaker's GPU Coding
  Interface, so the device model can be swapped in either direction.
* Per-model porting notes: §16 of `radar_object_list/implementation_logic.md`
  and §10 of `radar_rsi/implementation_logic.md`.
