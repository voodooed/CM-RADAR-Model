"""
Propagation damping for the CarMaker Radar Sensor.

Reference Manual -> Sensors -> Radar Sensor -> Detection Based on Signal-to-Noise
Ratio -> Damping (Equation 512):

        L_atm = L_satm * L_rain * L_fog

with
  * ``L_satm(f)``            clear-air attenuation, cited to
                             "DESK, EW Class. Electronic Warfare and Radar
                             Systems Engineering Handbook"
  * ``L_rain(f, Q_rain)``    rain attenuation, cited to ITU-R P.838-3, using
                             ``Env.RainRate``
  * ``L_fog(f, d_vis)``      fog attenuation, cited to
                             "BROOKER, Graham. Introduction to sensors for
                             ranging and imaging", using ``Env.VisRangeInFog``

CarMaker publishes only the *references*, not the coefficients.  The functions
below therefore reproduce the published functional forms with coefficients that
are correct in order of magnitude at 77 GHz; every coefficient is a
configuration parameter (see ``RadarSensorConfig``).  This is **[INFERRED]** -
absolute damping values will not match CarMaker exactly.

The protected header ``Sensor_Radar_protected.h`` documents
``double Damping_km;  // sum of atm., rain and fog damping [dB/km]`` and, per
object, ``double Damping;  // (Absolute) sum of atm., rain and fog damping [dB]``.
We therefore model a single specific attenuation in dB/km and convert it to a
per-object loss over the **two-way** path ``2 * r``.
"""

from __future__ import annotations

import math


def clear_air_db_per_km(frequency_ghz: float) -> float:
    """One-way clear-air specific attenuation [dB/km].

    Piecewise log-log interpolation through commonly published sea-level values
    for a standard atmosphere (oxygen + water-vapour), 7.5 g/m^3, 20 degC.
    At 77 GHz this yields ~0.45 dB/km, the value usually quoted for automotive
    radar.  [INFERRED]
    """
    # (GHz, dB/km) support points of the standard-atmosphere attenuation curve
    pts = ((1.0, 0.007), (10.0, 0.01), (22.0, 0.19), (35.0, 0.09),
           (60.0, 15.0), (77.0, 0.45), (94.0, 0.40), (140.0, 1.8), (300.0, 30.0))
    f = max(frequency_ghz, pts[0][0])
    if f >= pts[-1][0]:
        return pts[-1][1]
    for i in range(len(pts) - 1):
        f0, a0 = pts[i]
        f1, a1 = pts[i + 1]
        if f0 <= f <= f1:
            w = (math.log(f) - math.log(f0)) / (math.log(f1) - math.log(f0))
            return math.exp(math.log(a0) * (1 - w) + math.log(a1) * w)
    return pts[-1][1]


def rain_db_per_km(rain_rate_mm_h: float, k: float, alpha: float) -> float:
    """One-way rain attenuation [dB/km] using the ITU-R P.838-3 power law

        gamma_R = k * R**alpha

    ``k`` and ``alpha`` are frequency (and polarisation) dependent; defaults in
    ``RadarSensorConfig`` (k=0.9, alpha=0.74) reproduce the widely quoted
    ~10 dB/km at 25 mm/h for 77 GHz.  [INFERRED coefficients]
    """
    if rain_rate_mm_h <= 0.0:
        return 0.0
    return k * (rain_rate_mm_h ** alpha)


def fog_liquid_water_content(vis_range_m: float) -> float:
    """Liquid water content [g/m^3] from visual range [m].

    Kruse's empirical relation, as used in Brooker's textbook:
        V [km] = 0.024 / M**0.65      ->     M = (0.024 / V)**(1/0.65)
    """
    if vis_range_m <= 0.0 or vis_range_m > 1e5:
        return 0.0
    v_km = vis_range_m / 1000.0
    return (0.024 / v_km) ** (1.0 / 0.65)


def fog_db_per_km(vis_range_m: float, k_db_km_per_gm3: float) -> float:
    """One-way fog attenuation [dB/km] = K_l(f, T) * M.  [INFERRED coefficient]"""
    m = fog_liquid_water_content(vis_range_m)
    if m <= 0.0:
        return 0.0
    return k_db_km_per_gm3 * m


def total_damping_db_per_km(frequency_ghz: float,
                            rain_rate_mm_h: float,
                            vis_range_fog_m: float,
                            rain_k: float, rain_alpha: float,
                            fog_k: float,
                            atm_override: float | None = None) -> float:
    """``Damping_km`` of ``Sensor_Radar_protected.h``: one-way sum in dB/km."""
    atm = clear_air_db_per_km(frequency_ghz) if atm_override is None else atm_override
    return atm + rain_db_per_km(rain_rate_mm_h, rain_k, rain_alpha) \
               + fog_db_per_km(vis_range_fog_m, fog_k)


def two_way_damping_db(damping_db_per_km: float, range_m: float) -> float:
    """``Damping`` of ``Sensor_Radar_protected.h``: absolute two-way loss [dB]."""
    return damping_db_per_km * 2.0 * range_m / 1000.0
