#!/usr/bin/env python3
"""
Round-trip test of the IPG file readers against every radar-relevant file
shipped with CarMaker (``<install>/Data/Sensor``).

This validates the *interface* side of the reimplementation: the models can be
driven with CarMaker's own parameterisation, unmodified.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common.infofile import InfoFile                              # noqa: E402
from common.paths import (default_carmaker_dir, missing_message,   # noqa: E402
                          short)
from radar_object_list.maps import AntennaGainMap, RcsMap         # noqa: E402
from radar_rsi.scene import MaterialLib                           # noqa: E402
from radar_rsi.transceiver import TransceiverConfig               # noqa: E402



def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ipg", default=default_carmaker_dir(),
                    help="CarMaker installation dir (default: $CARMAKER_DIR)")
    args = ap.parse_args()
    if not args.ipg:
        print(missing_message())
        return
    d = os.path.join(args.ipg, "Data", "Sensor")
    if not os.path.isdir(d):
        print(f"[skip] {d} not found")
        return

    print(f"scanning {short(d, 4)}\n")
    files = sorted(f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f)))
    n_ok = n_skip = n_err = 0
    for name in files:
        path = os.path.join(d, name)
        try:
            ident = InfoFile.load(path).str("FileIdent", "")
        except Exception as exc:                       # noqa: BLE001
            print(f"  {name:<24} unreadable: {exc}")
            n_err += 1
            continue

        try:
            if ident.startswith("CarMaker-AntennaGainMap"):
                m = AntennaGainMap.load(path)
                print(f"  {name:<24} AntennaGainMap   {len(m.az)}x{len(m.el)} samples, "
                      f"BeamWidth {m.beam_width}, peak {m.peak_db:.2f} dB")
            elif ident.startswith("CarMaker-RCS Map"):
                m = RcsMap.load(path, name)
                print(f"  {name:<24} RCS Map          {len(m.azim)} samples, "
                      f"ProbExist {m.prob_exist}, OcclusionFactor "
                      f"{m.occlusion_factor:.2f}, "
                      f"max {10 * math.log10(max(m.rcs_lin)):.1f} dBm^2")
            elif "RadarRSI_TranceiverConfig" in ident:
                t = TransceiverConfig.load(path)
                print(f"  {name:<24} RadarRSI config  Tx {t.tx.n_az}x{t.tx.n_el} "
                      f"(peak {t.tx.peak_db:.2f} dB), Rx {t.rx.n_az}x{t.rx.n_el} "
                      f"(peak {t.rx.peak_db:.2f} dB), "
                      f"{len(t.vrx) if t.vrx else 0} VRx")
            elif ident.startswith("MaterialLib"):
                m = MaterialLib.load(path)
                print(f"  {name:<24} MaterialLib      {len(m.materials)} materials "
                      f"(e.g. asphalt eps={m.get('asphalt').permittivity}, "
                      f"scatter={m.get('asphalt').scattering_deg} deg)")
            else:
                print(f"  {name:<24} (not radar related: FileIdent '{ident}')")
                n_skip += 1
                continue
            n_ok += 1
        except Exception as exc:                       # noqa: BLE001
            print(f"  {name:<24} FAILED: {exc}")
            n_err += 1

    print(f"\n  {n_ok} radar file(s) read successfully, {n_skip} skipped, "
          f"{n_err} error(s)")


if __name__ == "__main__":
    main()
