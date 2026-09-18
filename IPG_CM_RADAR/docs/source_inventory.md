# Source inventory — which IPG files this work is based on

Installation examined: **CarMaker 15.1**, `carmaker/win64-15.1`.
Paths below are relative to that installation directory. The sibling top-level
directories of the IPG tree (`bin`, `control`, `graph`, `instruments`,
`movienx`) contain no radar model material.

## Method and limits

* Everything below is **documentation, public C headers, shipped example source,
  or shipped ASCII data files**.
* The models themselves are compiled (`lib/libcarmaker.a` members
  `Sensor_Radar.o`, `Sensor_Radar_Calc.o`, `Sensor_Radar_Params.o`,
  `Sensor_RadarRSI.o`; plus `GUI/IPGRsi.dll`). The IPG EULA §4(ii)
  (`doc/IPG_EULA.txt`) forbids decompiling, disassembling and reverse
  engineering. **No decompilation or disassembly was performed, and no
  statement in this package is derived from the compiled binaries.** Where the
  documentation is silent, the gap is filled by an explicitly flagged inference.
* The Reference Manual is MadCap-generated HTML with the equations rendered as
  MathJax SVG. Equation text quoted in this package was recovered by decoding
  the glyph references in those SVGs.

## How the two radar models were identified

`include/Vehicle/` contains exactly two radar headers —
`Sensor_Radar.h` (+ `Sensor_Radar_protected.h`) and `Sensor_RadarRSI.h` — and
the Reference Manual table of contents contains exactly two radar chapters
(entries 1060–1083 "Radar Sensor" and 1137–1163 "Radar RSI"). `Data/Sensor`
correspondingly ships `Radar_Default` (antenna gain map for the object-list
radar) and `RadarRSI_Default` (transceiver config for the RSI radar).
No third radar model exists in this release.

---

## Documentation (Reference Manual)

Path prefix: `doc/ReferenceManual/Content/ReferenceManual/`

### Radar Sensor (Model 1)

| Section | File |
|---|---|
| Introduction | `Sensors/Sensor_Radar_Intro.htm` |
| Detection Based on Signal-to-Noise Ratio (Eq. 502) | `Sensors/Detection_Based_on_Signa.htm` |
| Signal-to-Noise Ratio (Eq. 503) | `Sensors/Signal_to_Noise_Ratio.htm` |
| Detection Threshold (Eq. 504) | `Sensors/Detection_Threshold.htm` |
| Signal Strength (Eq. 505) | `Sensors/Signal_Strength.htm` |
| Detection Noise (Eq. 506, road clutter) | `Sensors/Detection_Noise.htm` |
| Antenna Gain (Eq. 507, 508) | `Sensors/Antenna_Gain.htm` |
| Object Distance (BBC reference point) | `Sensors/Object_Distance.htm` |
| Radar Cross Section (Eq. 509, 510, 511) | `Sensors/Radar_Cross_Section.htm` |
| Damping (Eq. 512) | `Sensors/Damping_1.htm` |
| Measurement Noise (Eq. 513) | `Sensors/Measurement_Noise.htm` |
| Time Delay (cycle time, latency) | `Sensors/Time_Delay.htm` |
| Separability (Eq. 514) | `Sensors/Separability.htm` |
| False positives | `Sensors/False_positives.htm` |
| Sensor Quantities | `Sensors/Sensor_Quantities.htm` |
| Object Quantities (ProbExist, ProbObst, DynProp) | `Sensors/Object_Quantities.htm` |
| Sensor Parameters | `Sensors/Parameterization_of_Radar_Sensor.htm` |
| Antenna Gain Map File format | `Sensors/Antenna_Gain_Map_File.htm` |
| RCS map file format | `Sensors/RCS_map_files.htm` |
| Visualization | `Sensors/Visualization_4.htm` |
| UAQ — General | `User_Accessible_Quantities/Radar_Sensor_General.htm` |
| UAQ — Object List | `User_Accessible_Quantities/Object_List_2.htm` |

### Radar RSI (Model 2)

| Section | File |
|---|---|
| Introduction | `Sensors/Sensor_RadarRSI_Intro.htm` |
| Ray Pattern | `Sensors/Ray_Pattern_1.htm` |
| Static Ray Pattern | `Sensors/Static_Ray_Pattern_1.htm` |
| Dynamic Ray Pattern | `Sensors/Dynamic_Ray_Pattern_1.htm` |
| Wave Propagation Model | `Sensors/Wave_Propagation_Model.htm` |
| 3D Geometry Illumination | `Sensors/3D_Geometry_Illumination.htm` |
| Material-Dependent Reflection | `Sensors/Material_Dependent_Refle.htm` |
| Multipath Propagation | `Sensors/Multipath_Propagation.htm` |
| Doppler Shift | `Sensors/Doppler_Shift.htm` |
| Scattered Electromagnetic Fields | `Sensors/Scattered_Electromagneti.htm` |
| Device Model | `Sensors/Device_Model.htm` |
| Antenna Characteristics (Eq. 521–524) | `Sensors/Antenna_Characteristics.htm` |
| Front-End Hardware (Eq. 525) | `Sensors/Front_End_Hardware.htm` |
| Back-End Software (Eq. 526–529) | `Sensors/Back_End_Software.htm` |
| Object Model Requirements | `Sensors/Object_Model_Requirement_1.htm` |
| Output Quantities | `Sensors/Output_Quantities_4.htm` |
| Output as HDF5-File | `Sensors/Output_as_HDF5_File.htm` |
| Sensor Parameters | `Sensors/Sensor_Parameters_12.htm` |
| Elevation Angle Processing | `Sensors/Elevation_Angle_Processi.htm` |
| Antenna / VRx File format | `Sensors/Antenna___VRx_File.htm` |
| GPU Coding Interface | `Sensors/GPU_Coding_Interface_3.htm` |
| UAQ | `User_Accessible_Quantities/Radar_RSI.htm` |

### Shared / supporting

| Topic | File |
|---|---|
| Sensor Assembly (`pos`, `rot` z-y-x, `CycleTime`, `Latency`, `Mounting`) | `Sensors/Sensor_Assembly.htm` |
| Sensor Surrounding (candidate object set, +20 m tolerance, 0.2 m static update) | `Sensors/Sensor_Surrounding.htm` |
| Surrounding Parameters | `Sensors/Surrounding_Parameters.htm` |
| Cycle and Latency Parameterization | `Sensors/Cycle_and_Latency_Parame.htm` |
| Object ID numbering | `Object_ID.htm` |
| Traffic Object Parameters (`Basics.Dimension`, `RCSMap.FName`, `DetectMask`) | `Traffic/Traffic_Object_Parameter.htm` |

### User's Guide

`doc/UsersGuide/Content/UsersGuide/Vehicle_Model/Radar_Sensor.htm` —
GUI-level descriptions; confirms the detectable object kinds (traffic objects,
geometry objects, tree trunks within 10 m of the road border, guide posts,
guard-rail posts) and the meaning of every dialog parameter.

---

## C headers

| File | Used for |
|---|---|
| `include/Vehicle/Sensor_Radar.h` | `tOutQuants`, `tOutQuantsGlob`, `tRadarSensor`, API |
| `include/Vehicle/Sensor_Radar_protected.h` | `tPnt`, `ePntType` (bounding-box vertex order), `tPTraffic`, `tObsPnt`, `tRadObject` (`azimuth_min`/`max`, `OccPhi`, `OccTheta`, `AntennaGain`, `Damping`, `IdRunning`), `tRadarSensorProtected` (`Damping_km`, `erfcinvPFA`, `ScanRange`, `BeamWidth`, `ant_eff`). **Deprecated by IPG — support ends with CarMaker 16.0** |
| `include/Vehicle/Sensor_RadarRSI.h` | `tDetPoint`, `tCompOut`, `tDetVRx`, `tRadarRSI`, `tOutTypeRadarRSI`, API incl. `RadarRSI_GetnVRx`, `RadarRSI_TriggerCalc` |
| `include/Vehicle/Sensor_Object.h` | `tObjectSensorObj` — the ideal object sensor, for contrast |

## Shipped example source

| File | Used for |
|---|---|
| `Examples/GPUCodingInterface/src.GPUCodingInterface/Include/GPUCodingInterface.h` | `tUserRadarInteractionPoint` (`range`, `vel`, `direction`, `electricField`), `tUserRadarCube` **incl. the documented index formula**, `tUserRadarDataType` |
| `Examples/GPUCodingInterface/src.GPUCodingInterface/Samples/4_Radar/SignalGenerationAndProcessing/Radar.cu` | reference implementation of writing IA points into the radar cube and of a peak finder; fixes the range/velocity/angle ↔ bin mapping |
| `.../Samples/4_Radar/.../README.md` | describes the signal-generation / signal-processing split |

## Data files

| File | Contents |
|---|---|
| `Data/Sensor/Radar_Default` | `CarMaker-AntennaGainMap`, `BeamWidth 20 15`, `ScanRange 0 0`, `AntennaEff 1.0`, 181×181 gain matrix, peak 20.3655 dB |
| `Data/Sensor/RadarRSI_Default` | `CarMaker-RadarRSI_TranceiverConfig`, Tx/Rx 181×181 maps over 180°×180°, 16 VRx at ≈1.95 mm spacing (= λ/2 at 77 GHz) |
| `Data/Sensor/RCS_{Car,Truck,Pedestrian,Bicycle,GuardRailPost}` | `CarMaker-RCS Map`, azimuth LUT + `ProbExist` + `OcclusionFactor` |
| `Data/Sensor/MaterialLib` | `Material.<i>.Radar.Permittivity` / `.Radar.Scattering` for 22 materials |
| `Data/Vehicle/Examples/DemoCar_SensorRadar` | full Model 1 parameter set |
| `Data/Vehicle/Examples/DemoCar_SensorRadar_FalsePositives` | same with `FalsePosActive = 1` |
| `Data/Vehicle/Examples/DemoCar_SensorRadarRSI` | full Model 2 parameter set incl. the NoiseScaling LUTs |
| `Data/Vehicle/Examples/Demo_IPG_CompanyCar_MultiRSI` | multi-RSI sensor cluster example |
