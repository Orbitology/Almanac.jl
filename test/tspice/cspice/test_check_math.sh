#!/usr/bin/env bash
# Tests of check_math.sh on two small programs: one that takes cos from the C library
# (refused), and one that takes only fmod (accepted).
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
fail=0

prog() {   # NAME FUNCTION: a program that calls FUNCTION on a value it cannot fold
    printf '#include <math.h>\n#include <stdio.h>\nvolatile double x = 0.5;\nint main(void) { printf("%%.17g\\n", %s); return 0; }\n' "$2" > "$tmp/$1.c"
    gcc -O2 -o "$tmp/$1" "$tmp/$1.c" -lm
}

prog takes_cos 'cos(x)'
if bash "$here/check_math.sh" "$tmp/takes_cos" 2> "$tmp/err"; then
    echo "FAIL: a program that takes cos passed"; fail=1
elif ! grep -qw cos "$tmp/err"; then
    echo "FAIL: the refusal does not name cos: $(cat "$tmp/err")"; fail=1
fi

prog takes_fmod 'fmod(x, 0.3)'
bash "$here/check_math.sh" "$tmp/takes_fmod" || { echo "FAIL: a program that takes only fmod was refused"; fail=1; }

[ "$fail" = 0 ] && echo "check_math: 2 tests passed"
exit "$fail"
