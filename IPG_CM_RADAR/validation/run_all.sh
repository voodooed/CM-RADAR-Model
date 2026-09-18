#!/usr/bin/env bash
# Run every validation script and store its output under validation/results/.
#
# Usage:  ./validation/run_all.sh [<CarMaker install dir>]
#
# The install dir may also come from $CARMAKER_DIR.  Without either, the scripts
# that need IPG data files skip themselves and the rest still run.
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IPG="${1:-${CARMAKER_DIR:-}}"
OUT="$ROOT/validation/results"
mkdir -p "$OUT"
cd "$ROOT"

# only pass --ipg when we actually have a path
if [ -n "$IPG" ]; then IPG_ARG=(--ipg "$IPG"); else IPG_ARG=(); fi
export CARMAKER_DIR="$IPG"

run () {
    local name="$1"; shift
    echo "=== $name ==="
    "$@" > "$OUT/$name.txt" 2>&1
    local rc=$?
    if [ $rc -eq 0 ]; then
        echo "    ok   -> validation/results/$name.txt"
    else
        echo "    FAIL (exit $rc) -> validation/results/$name.txt"
    fi
    grep -c "FAIL" "$OUT/$name.txt" | sed 's/^/    FAIL markers: /'
}

run ipg_files          python3 validation/validate_ipg_files.py "${IPG_ARG[@]}"
run antenna_map        python3 validation/validate_antenna_map.py "${IPG_ARG[@]}"
run detection_math     python3 validation/validate_detection_math.py
run object_list_model  python3 validation/validate_object_list_model.py
run rsi_device         python3 validation/validate_rsi_device.py
run rsi_propagation    python3 validation/validate_rsi_propagation.py
run example_object_list python3 radar_object_list/example.py "${IPG_ARG[@]}"
run example_rsi        python3 radar_rsi/example.py "${IPG_ARG[@]}" --rays 4000

echo
echo "All outputs in $OUT"
