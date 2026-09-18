"""
Detection mathematics of the CarMaker Radar Sensor.

All formulas in this module are taken verbatim from the Reference Manual
(-> Sensors -> Radar Sensor -> Detection Based on Signal-to-Noise Ratio):

  Eq. 502   detection  <=>  SNR > SNR_min
  Eq. 503   SNR = S / N
  Eq. 504   SNR_min = 2 * ( erfc^-1(2 P_FA) - erfc^-1(2 P_Dmin) )^2
  Eq. 505   S = P G^2 lambda^2 RCS / ((4 pi)^3 r^4) * 1 / (L_A L_atm)
  Eq. 506   N_thermal = T_0 F_R k_B B_n

The only non-documented part is the road-clutter contribution to N (see
``clutter_noise``), which the manual describes qualitatively only.
"""

from __future__ import annotations

import math
from radar_object_list.config import K_BOLTZMANN


# ===========================================================================
#  inverse complementary error function
# ===========================================================================
def erfcinv(y: float) -> float:
    """Inverse of ``math.erfc`` for y in (0, 2).

    Rational approximation (Giles / Acklam style) refined by two Newton steps
    against ``math.erfc``; accurate to ~1e-14 over the range needed here.
    """
    if y <= 0.0:
        return float("inf")
    if y >= 2.0:
        return float("-inf")
    if y == 1.0:
        return 0.0
    # initial guess from the inverse normal CDF: erfc^-1(y) = -Phi^-1(y/2)/sqrt(2)
    x = -_norm_ppf(0.5 * y) / math.sqrt(2.0)
    for _ in range(3):
        err = math.erfc(x) - y
        # d/dx erfc(x) = -2/sqrt(pi) * exp(-x^2)
        d = -2.0 / math.sqrt(math.pi) * math.exp(-x * x)
        if d == 0.0:
            break
        x -= err / d
    return x


def _norm_ppf(p: float) -> float:
    """Acklam's inverse standard-normal CDF (|rel err| < 1.15e-9)."""
    a = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
    b = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00)
    plow, phigh = 0.02425, 1.0 - 0.02425
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    if p > phigh:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)


# Precomputed erfc^-1(2 * 10^-idx) for the ten documented values of
# ``ProbFalseAlarmIdx`` (1..10 -> P_FA = 1e-1 .. 1e-10).  Index 0 is unused.
ERFCINV_PFA = tuple([0.0] + [erfcinv(2.0 * 10.0 ** (-k)) for k in range(1, 11)])


# ===========================================================================
#  threshold and probability of detection
# ===========================================================================
def snr_min(prob_false_alarm: float, prob_detect_min: float) -> float:
    """Reference Manual Equation 504 (linear, not dB)."""
    return 2.0 * (erfcinv(2.0 * prob_false_alarm) - erfcinv(2.0 * prob_detect_min)) ** 2


def prob_detect(snr: float, prob_false_alarm: float) -> float:
    """Probability of detection, i.e. Equation 504 solved for ``P_D``.

    From  SNR = 2 (erfc^-1(2 P_FA) - erfc^-1(2 P_D))^2   and  P_D > P_FA:
          erfc^-1(2 P_D) = erfc^-1(2 P_FA) - sqrt(SNR / 2)
          P_D            = 0.5 * erfc( erfc^-1(2 P_FA) - sqrt(SNR / 2) )

    The Reference Manual lists ``ProbDetect`` as an output quantity "resulting
    from Equation 504", so this inversion is the documented relation.
    """
    if snr <= 0.0:
        return 0.0
    x = erfcinv(2.0 * prob_false_alarm) - math.sqrt(0.5 * snr)
    return 0.5 * math.erfc(x)


# ===========================================================================
#  radar equation and noise
# ===========================================================================
def signal_strength(tx_power_w: float, gain_db: float, wavelength: float,
                    rcs_m2: float, range_m: float,
                    system_losses_db: float, damping_db: float) -> float:
    """Reference Manual Equation 505 - received power [W].

    ``gain_db`` is the *one-way* antenna gain; the equation uses G^2 because
    "Gain of the antenna G_t = G_r" (transmit and receive share one map).
    ``damping_db`` is the absolute two-way atmospheric loss L_atm [dB],
    ``system_losses_db`` is L_A [dB].
    """
    if range_m <= 0.0 or rcs_m2 <= 0.0:
        return 0.0
    g = 10.0 ** (gain_db / 10.0)
    l_a = 10.0 ** (system_losses_db / 10.0)
    l_atm = 10.0 ** (damping_db / 10.0)
    num = tx_power_w * g * g * wavelength * wavelength * rcs_m2
    den = ((4.0 * math.pi) ** 3) * (range_m ** 4) * l_a * l_atm
    return num / den


def thermal_noise(temperature_k: float, noise_figure_db: float,
                  noise_bandwidth_hz: float) -> float:
    """Reference Manual Equation 506 - N_thermal [W]."""
    f_r = 10.0 ** (noise_figure_db / 10.0)
    return temperature_k * f_r * K_BOLTZMANN * noise_bandwidth_hz


def clutter_noise(range_m: float, tx_power_w: float, gain_db: float,
                  wavelength: float, sigma0_db: float,
                  azimuth_resolution_rad: float, range_resolution_m: float,
                  sensor_height_m: float, system_losses_db: float,
                  damping_db: float, grazing_min_rad: float) -> float:
    """Road-clutter noise power [W] at a given range.  **[INFERRED]**

    The Reference Manual (-> Detection Noise) only states:
      "noise emerging from road clutter N_clutter is modeled.  The antenna gain
       G, transmitted power P, reflectivity of the street sigma_0 as well as the
       sensor resolution are used to compute a range dependent clutter noise.
       This clutter noise is evaluated and added to the thermal noise for every
       object with respect to its distance to the sensor."

    We implement the textbook beam-limited surface-clutter model that uses
    exactly those four ingredients:

        A_c(r)     = r * Theta_az * dR / cos(psi)          illuminated ground patch
        sigma_c(r) = sigma_0 * A_c(r)                      clutter RCS
        N_clutter  = radar equation evaluated with sigma_c

    with the grazing angle ``psi = atan(h_sensor / r)`` (floored at
    ``grazing_min_rad``).  ``sigma_0`` is a configuration parameter
    (``clutter_sigma0_db``); CarMaker does not publish its value.
    """
    if range_m <= 0.0:
        return 0.0
    psi = max(math.atan2(max(sensor_height_m, 0.0), range_m), grazing_min_rad)
    area = range_m * azimuth_resolution_rad * range_resolution_m / max(math.cos(psi), 1e-3)
    sigma_c = (10.0 ** (sigma0_db / 10.0)) * area
    return signal_strength(tx_power_w, gain_db, wavelength, sigma_c, range_m,
                           system_losses_db, damping_db)


def to_db(x: float, floor_db: float = -400.0) -> float:
    return 10.0 * math.log10(x) if x > 0.0 else floor_db


def to_dbw(x: float) -> float:
    """``SignalStrength`` is documented with unit dBW."""
    return to_db(x)
