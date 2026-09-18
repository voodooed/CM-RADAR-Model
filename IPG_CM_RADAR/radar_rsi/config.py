"""
Configuration of the CarMaker **Radar RSI** (Raw Signal Interface radar).

Every field maps to a documented ``Sensor.Param.<no>.*`` key
(Reference Manual -> Sensors -> Radar RSI -> Parameterization of the Radar RSI
-> Sensor Parameters / Elevation Angle Processing) or to a ``Sensor.<no>.*``
assembly key.  Defaults are those of the shipped example vehicle
``Data/Vehicle/Examples/DemoCar_SensorRadarRSI``.

Fields marked "[INFERRED]" have no documented counterpart.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Tuple

K_BOLTZMANN = 1.380649e-23
SPEED_OF_LIGHT = 299_792_458.0


class RayPattern(str, Enum):
    PHI_THETA = "PhiTheta"
    FIBONACCI = "Fibonacci"
    DYNAMIC = "Dynamic"


class OutputType(str, Enum):
    CARTESIAN = "Cartesian"
    SPHERICAL = "Spherical"
    VRX = "VRx"


class WindowFunction(str, Enum):
    RECTANGULAR = "Rectangular"
    HANN = "Hann"


class OccupancyCheckMode(str, Enum):
    INSTANCES_AND_SCENE_SCAN = "InstancesAndSceneScan"
    ONLY_SCENE_SCAN = "OnlySceneScan"
    ONLY_INSTANCES = "OnlyInstances"


@dataclass
class DiscretizationAxis:
    """One axis of the Range-Doppler / angular data cube."""
    samples: int
    zero_paddings: int = 0
    minimum: float = 0.0
    maximum: float = 0.0

    @property
    def n_bins(self) -> int:
        return self.samples + self.zero_paddings


@dataclass
class CfarConfig:
    """OS-CFAR (Range-Doppler) or CA-CFAR (angular) parameters."""
    layers: int = 5           # CfarFilter.*.Layers   n_lay
    guards: int = 2           # CfarFilter.*.Guards   n_guard
    snr_db: float = 11.0      # CfarFilter.*.SNR      S_CFAR  [dB]
    threshold_pct: float = 95.0   # CfarFilter.RangeDoppler.Threshold  T_OS-CFAR [%]


@dataclass
class RadarRSIConfig:
    # -------------------------------------------------- assembly / mounting
    name: str = "RARS00"                                   # Sensor.<no>.name
    pos: Tuple[float, float, float] = (4.4, 0.0, 0.5)      # Sensor.<no>.pos [m]
    rot_deg: Tuple[float, float, float] = (0.0, 0.0, 0.0)  # Sensor.<no>.rot [deg] z-y-x
    cycle_time_ms: float = 60.0                            # sensor cluster cycle time

    # ------------------------------------------------------------ transceiver
    fov_deg: Tuple[float, float] = (120.0, 20.0)   # Sensor.Param.<n>.FoV h/v [deg]
    range_min_max: Tuple[float, float] = (0.1, 150.0)  # Sensor.Param.<n>.Range [m]
    n_rays: int = 20000                            # Sensor.Param.<n>.nRays
    ray_pattern: RayPattern = RayPattern.FIBONACCI  # Sensor.Param.<n>.RayPattern
    frequency_ghz: float = 77.0                    # Sensor.Param.<n>.Frequency [GHz]
    transmit_power_dbm: float = 30.0               # Sensor.Param.<n>.TransmitPower [dBm]
    polarization_transmit: float = 0.5             # Sensor.Param.<n>.PolarizationTransmit (w_p)
    polarization_receive: float = 0.5              # Sensor.Param.<n>.PolarizationReceive
    system_losses_db: float = 0.0                  # Sensor.Param.<n>.SystemLosses [dB]
    antenna_fname: Optional[str] = None            # Sensor.Param.<n>.Antenna.FName

    # ------------------------------------------------------------ noise model
    noise_bandwidth_mhz: float = 200.0             # Sensor.Param.<n>.NoiseBandwidth [MHz]
    noise_figure_db: float = 5.0                   # Sensor.Param.<n>.NoiseFigure [dB]
    noise_temperature_k: float = 300.0             # Sensor.Param.<n>.NoiseTemperature [K]
    noise_scaling: bool = True                     # Sensor.Param.<n>.NoiseScaling
    # Sensor.Param.<n>.NoiseScalingRange  (range [m], additional noise figure [dB])
    noise_scaling_range: List[Tuple[float, float]] = field(
        default_factory=lambda: [(0.0, 10.0), (10.0, 8.0), (20.0, 5.0),
                                 (40.0, 0.0), (80.0, -1.0)])
    # Sensor.Param.<n>.NoiseScalingDopplerVel  (velocity [m/s], additional NF [dB])
    noise_scaling_doppler: List[Tuple[float, float]] = field(
        default_factory=lambda: [(-20.0, 0.0), (0.0, 0.0), (20.0, 0.0)])

    # --------------------------------------------------------- discretization
    range_max: float = 200.0            # Discretization.Range.Max [m] (min fixed to 0)
    range_samples: int = 512            # Discretization.Range.Samples
    range_zero_paddings: int = 0        # Discretization.Range.ZeroPaddings
    doppler_min_max: Tuple[float, float] = (-40.0, 40.0)   # Discretization.DopplerVel.MinMax [m/s]
    doppler_samples: int = 256          # Discretization.DopplerVel.Samples
    doppler_zero_paddings: int = 0      # Discretization.DopplerVel.ZeroPaddings
    azimuth_samples: int = 60           # Discretization.AzimuthAngle.Samples  (= n_ULA)
    azimuth_zero_paddings: int = 100    # Discretization.AzimuthAngle.ZeroPaddings
    elevation_samples: int = 7          # Discretization.ElevationAngle.Samples
    elevation_zero_paddings: int = 10   # Discretization.ElevationAngle.ZeroPaddings
    process_elevation_angle: bool = False   # Sensor.Param.<n>.ProcessElevationAngle
    window_function: WindowFunction = WindowFunction.RECTANGULAR   # Sensor.Param.<n>.WindowFunction

    # ------------------------------------------------------------- filtering
    cfar_range_doppler: CfarConfig = field(
        default_factory=lambda: CfarConfig(layers=5, guards=2, snr_db=11.0, threshold_pct=95.0))
    cfar_azimuth: CfarConfig = field(
        default_factory=lambda: CfarConfig(layers=8, guards=2, snr_db=12.0))
    cfar_elevation: CfarConfig = field(
        default_factory=lambda: CfarConfig(layers=2, guards=1, snr_db=12.0))
    peak_finder_range_doppler_layers: int = 2   # PeakFinder.RangeDoppler.Layers
    peak_finder_azimuth_layers: int = 2         # PeakFinder.AzimuthAngle.Layers
    peak_finder_elevation_layers: int = 1       # PeakFinder.ElevationAngle.Layers
    peak_interpolation: bool = True             # documented behaviour (Eq. 527-529)

    # ---------------------------------------------------------------- output
    output_type: OutputType = OutputType.CARTESIAN   # Sensor.Param.<n>.OutputType
    n_max_detections: int = 2000                     # Sensor.Param.<n>.nMaxDetections

    # --------------------------------------------- dynamic ray pattern params
    n_horizontal_buckets: int = -1        # -1 -> default (FoV_h deg, 1 deg per bucket)
    n_vertical_buckets: int = -1          # -1 -> default (FoV_v deg, 1 deg per bucket)
    max_rays_per_bucket: int = -1         # -1 -> unlimited
    occupancy_check_mode: OccupancyCheckMode = OccupancyCheckMode.INSTANCES_AND_SCENE_SCAN

    # ---------------------------------------------------------- environment
    rain_rate_mm_h: float = 0.0
    vis_range_fog_m: float = 1.0e6

    # =============================================================== INFERRED
    max_bounces: int = 3            # [INFERRED] multipath depth of the ray tracer
    ray_energy_cutoff_db: float = -80.0   # [INFERRED] stop tracing weak rays
    # How each interaction point is written into the Range-Doppler map:
    #   "binning"  - nearest bin, as in Examples/GPUCodingInterface .../Radar.cu
    #   "spectral" - windowed DFT kernel spread over all bins (models leakage,
    #                aliasing and finite resolution as the manual describes)
    rd_synthesis: str = "spectral"        # [INFERRED]
    # Doppler outside [DopplerVel.MinMax] folds back into the interval, because
    # those are "the minimum and maximum unambiguous velocities" (Reference
    # Manual -> Back-End Software).  Set False to clamp instead, which is what
    # the shipped GPU-Coding-Interface sample Radar.cu does.
    doppler_aliasing: bool = True
    rd_kernel_halfwidth: int = 8          # [INFERRED] bins evaluated around a peak
    random_seed: int = 1

    # ------------------------------------------------------------- derived
    @property
    def wavelength(self) -> float:
        return SPEED_OF_LIGHT / (self.frequency_ghz * 1e9)

    @property
    def transmit_power_w(self) -> float:
        return 10.0 ** ((self.transmit_power_dbm - 30.0) / 10.0)

    @property
    def noise_bandwidth_hz(self) -> float:
        return self.noise_bandwidth_mhz * 1e6

    @property
    def n_range_bins(self) -> int:
        return self.range_samples + self.range_zero_paddings

    @property
    def n_doppler_bins(self) -> int:
        return self.doppler_samples + self.doppler_zero_paddings

    @property
    def n_azimuth_bins(self) -> int:
        return self.azimuth_samples + self.azimuth_zero_paddings

    @property
    def n_elevation_bins(self) -> int:
        return self.elevation_samples + self.elevation_zero_paddings

    @property
    def n_vrx(self) -> int:
        """Number of virtual receivers used by the internal angular processing.

        Reference Manual -> Back-End Software -> Angular-Processing:
        "The FMCW-model incorporates uniform-linear arrays (ULA) of virtual
         receivers ... spaced d = lambda/2 ... The number of virtual receivers
         n_ULA (samples) as well as the number of zero-paddings n_zero^ULA are
         used to determine the angular spectrum using a spacial DFT."
        """
        return self.azimuth_samples

    def validate(self) -> None:
        assert self.n_rays > 0
        assert self.range_max > 0.0
        assert self.range_samples > 1 and self.doppler_samples > 1
        assert self.azimuth_samples > 1
        assert self.doppler_min_max[1] > self.doppler_min_max[0]
        assert 0.0 <= self.polarization_transmit <= 1.0
        assert 0.0 <= self.polarization_receive <= 1.0
        assert self.n_max_detections > 0
