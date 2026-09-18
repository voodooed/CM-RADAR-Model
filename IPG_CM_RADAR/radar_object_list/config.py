"""
Configuration of the CarMaker "Radar Sensor" (object-list HiFi radar).

Every field maps to a documented ``Sensor.Param.<no>.*`` key
(Reference Manual -> Sensors -> Radar Sensor -> Parameterization of Radar
Sensor -> Sensor Parameters) or to a documented ``Sensor.<no>.*`` assembly key
(Reference Manual -> Sensors -> Sensor Assembly).  The CarMaker key is given in
the comment next to each field.

Fields marked "[INFERRED]" have no documented counterpart; they parameterise a
part of the model that CarMaker documents only qualitatively.  Their defaults
are chosen to be physically reasonable, not to match IPG bit-for-bit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

# Physical constants
K_BOLTZMANN = 1.380649e-23          # J/K
SPEED_OF_LIGHT = 299_792_458.0      # m/s


@dataclass
class RadarSensorConfig:
    # -------------------------------------------------- assembly / mounting
    name: str = "RA00"                       # Sensor.<no>.name
    pos: Tuple[float, float, float] = (4.2, 0.0, 0.4)      # Sensor.<no>.pos  [m]
    rot_deg: Tuple[float, float, float] = (0.0, 0.0, 0.0)  # Sensor.<no>.rot  [deg], z-y-x
    cycle_time_ms: float = 60.0              # Sensor.<no>.CycleTime  [ms]
    cycle_offset_ms: float = 0.0             # Sensor.<no>.CycleOffset
    latency_factors: Tuple[float, float] = (1.0, 1.0)      # Sensor.<no>.Latency (x CycleTime)

    # ------------------------------------------------------ observation area
    fov_deg: Tuple[float, float] = (40.0, 30.0)   # Sensor.Param.<n>.FoV  h/v [deg]
    range_min: float = 0.2                        # Sensor.Param.<n>.Range_min [m]
    range_max: float = 200.0                      # Sensor.Param.<n>.Range_max [m]
    max_num_obj: int = 200                        # Sensor.Param.<n>.MaxNumObj (<= 2000)

    # ----------------------------------------------------------- radar front end
    frequency_ghz: float = 77.0              # Sensor.Param.<n>.Frequency  [GHz]
    transmit_power_dbm: float = 14.0         # Sensor.Param.<n>.TransmitPower [dBm]
    system_losses_db: float = 0.0            # Sensor.Param.<n>.SystemLosses  [dB]
    noise_bandwidth_hz: float = 25_000.0     # Sensor.Param.<n>.NoiseBandWidth [Hz]
    noise_figure_db: float = 4.8             # Sensor.Param.<n>.NoiseFigure    [dB]

    # ----------------------------------------------------- detection threshold
    prob_detect_min: float = 0.5             # Sensor.Param.<n>.ProbDetectMin
    prob_false_alarm_idx: int = 6            # Sensor.Param.<n>.ProbFalseAlarmIdx (1..10 -> 1e-1..1e-10)

    # -------------------------------------------- resolution / accuracy / merge
    accuracy_distance: float = 0.4           # Sensor.Param.<n>.AccuracyDistance  [m]      (1 sigma)
    accuracy_azimuth_deg: float = 0.1        # Sensor.Param.<n>.AccuracyAzimuth   [deg]    (1 sigma)
    accuracy_speed_kmh: float = 0.1          # Sensor.Param.<n>.AccuracySpeed     [km/h]   (1 sigma)
    resolution_distance: float = 1.8         # Sensor.Param.<n>.ResolutionDistance [m]
    resolution_azimuth_deg: float = 1.6      # Sensor.Param.<n>.ResolutionAzimuth  [deg]
    resolution_speed_kmh: float = 0.4        # Sensor.Param.<n>.ResolutionSpeed    [km/h]
    separability: float = 1.5                # Sensor.Param.<n>.Separability

    # ----------------------------------------------------------- false positives
    false_pos_active: bool = False           # Sensor.Param.<n>.FalsePosActive
    clutter_obj_mean: float = 5.0            # Sensor.Param.<n>.ClutterObjMean

    # ------------------------------------------------------------------ antenna
    gain_fname: Optional[str] = None         # Sensor.Param.<n>.Gain.FName (Data/Sensor/<name>)
    # Used only when gain_fname is None -> internal uniform-rectangular-aperture model
    beam_width_deg: Tuple[float, float] = (20.0, 15.0)   # AntennaGainMap: BeamWidth h/v
    scan_range_deg: Tuple[float, float] = (0.0, 0.0)     # AntennaGainMap: ScanRange h/v
    antenna_eff: float = 1.0                             # AntennaGainMap: AntennaEff

    # --------------------------------------------------------------- environment
    temperature_k: float = 293.15            # Env.Temperature (converted to K)
    rain_rate_mm_h: float = 0.0              # Env.RainRate            [mm/h]
    vis_range_fog_m: float = 1.0e6           # Env.VisRangeInFog       [m]

    # ================================================================ INFERRED
    # Road-clutter noise. Documented only qualitatively ("antenna gain G,
    # transmitted power P, reflectivity of the street sigma_0 as well as the
    # sensor resolution are used to compute a range dependent clutter noise").
    clutter_enabled: bool = True                 # [INFERRED] switch
    clutter_sigma0_db: float = -35.0             # [INFERRED] normalised road reflectivity sigma_0 [dB m^2/m^2]
    clutter_grazing_min_deg: float = 0.2         # [INFERRED] floor on the grazing angle

    # Clear-air / rain / fog specific attenuation.  CarMaker cites DESK EW
    # handbook, ITU-R P.838-3 and Brooker but publishes neither coefficients nor
    # closed forms; the values below reproduce the published functional forms at
    # 77 GHz and are user-overridable.
    atm_gamma_db_km: Optional[float] = None      # [INFERRED] clear-air one-way [dB/km]; None -> model
    rain_k: float = 0.9                          # [INFERRED] ITU form gamma = k * R^alpha  (77 GHz)
    rain_alpha: float = 0.74                     # [INFERRED]
    fog_k_db_km_per_gm3: float = 4.0             # [INFERRED] specific attenuation per g/m^3 at 77 GHz

    # Extended-object RCS averaging ("The RCS map is evaluated over a range of
    # angles depending on the distance, sensor resolution and orientation").
    extended_object_rcs: bool = True             # [INFERRED] switch
    # RCS reported back to the user is "corrected for the total noise".
    rcs_noise_correction: bool = True            # [INFERRED] switch

    # Probability-of-existence bookkeeping (documented behaviour, undocumented step size)
    prob_exist_max: int = 7
    prob_exist_step: int = 1                     # [INFERRED]

    # Obstacle probability (documented): safety gap and lateral clearance
    obstacle_safety_gap_m: float = 1.0
    obstacle_extra_clearance_m: float = 0.3

    # Reproducibility
    random_seed: int = 1

    # ------------------------------------------------------------------ derived
    @property
    def wavelength(self) -> float:
        return SPEED_OF_LIGHT / (self.frequency_ghz * 1e9)

    @property
    def prob_false_alarm(self) -> float:
        return 10.0 ** (-self.prob_false_alarm_idx)

    @property
    def transmit_power_w(self) -> float:
        return 10.0 ** ((self.transmit_power_dbm - 30.0) / 10.0)

    def validate(self) -> None:
        assert 0 < self.max_num_obj <= 2000, "MaxNumObj must be in [0, 2000]"
        assert 1 <= self.prob_false_alarm_idx <= 10, "ProbFalseAlarmIdx must be in [1, 10]"
        assert self.frequency_ghz > 0.0
        assert 0.0 < self.fov_deg[0] <= 180.0 and 0.0 < self.fov_deg[1] <= 180.0
        assert self.resolution_azimuth_deg >= 1e-5
        assert self.noise_bandwidth_hz > 0.0
        assert 0.0 < self.prob_detect_min < 1.0


@dataclass
class LengthWidthClassTable:
    """Documented class boundaries for ``LengthClass`` / ``WidthClass``
    (Reference Manual -> User Accessible Quantities -> Radar Sensor -> Object List):
        0 unknown, 1 <0.5 m, 2 <1 m, 3 <2 m, 4 <3 m, 5 <4 m, 6 <6 m, 7 exceeds."""

    bounds: Tuple[float, ...] = (0.5, 1.0, 2.0, 3.0, 4.0, 6.0)

    def classify(self, value: float) -> int:
        if value <= 0.0:
            return 0
        for i, b in enumerate(self.bounds):
            if value < b:
                return i + 1
        return 7


LENGTH_WIDTH_CLASSES = LengthWidthClassTable()


# ------------------------------------------------------------------ enums ----
class MeasStat:
    """``MeasStat`` output values (Reference Manual, Object List)."""
    NO_OBJECT = 0
    NEW_OBJECT = 1
    NOT_MEASURED = 2      # documented as "currently not used"
    MEASURED = 3


class DynProp:
    """``DynProp`` output values (Reference Manual, Object List)."""
    NOT_CLASSIFIED = 0
    STATIONARY = 1
    STOPPED = 2
    MOVING = 3
    ONCOMING = 4


# Object-ID conventions (Reference Manual -> Object ID, and -> False positives)
OBJ_ID_UNKNOWN = -1
OBJ_ID_CLUTTER = -2       # clutter false positives are always -2
OBJ_ID_EGO = 15_000_000
OBJ_ID_TRAFFIC_BASE = 16_000_000
