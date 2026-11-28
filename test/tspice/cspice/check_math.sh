#!/usr/bin/env bash
# Check that the program BIN takes from the C math library (libm.so) only functions
# whose result is exact (IEEE 754 requires it of them), so that its results are the same
# with every version of the library. Exit status 1, with the names, if it takes another
# one, or if BIN is not linked dynamically (then the check cannot see what it takes).
#
#     test/tspice/cspice/check_math.sh BIN
set -euo pipefail
export LC_ALL=C

[ $# -eq 1 ] || { echo "usage: check_math.sh BIN" >&2; exit 1; }
bin=$1
[ -f "$bin" ] || { echo "$bin: no such file" >&2; exit 1; }
if ! deps=$(ldd "$bin" 2>&1) || grep -q 'statically linked' <<< "$deps"; then
    echo "$bin is not linked dynamically: the check cannot see the functions that it takes" >&2
    exit 1
fi
libm=$(awk '$1 ~ /^libm\.so/ {print $3}' <<< "$deps")
[ -n "$libm" ] || exit 0
exact='sqrt|fabs|floor|ceil|trunc|round|roundeven|fmod|fma|fmin|fmax|fe[a-z]+'   # fe...: the floating-point environment
taken=$(comm -12 \
    <(nm -D --undefined-only "$bin" | awk '{print $2}' | sed 's/@.*//' | sort -u) \
    <(nm -D --defined-only "$libm" | awk '{print $3}' | sed 's/@.*//' | sort -u) |
    grep -vxE "$exact" || true)
[ -z "$taken" ] || { echo "$bin takes $(echo $taken) from $libm" >&2; exit 1; }
