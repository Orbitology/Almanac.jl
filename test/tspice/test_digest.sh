#!/usr/bin/env bash
# Tests of digest_recording.py on small fake recordings: the same content in two
# directories gives the same digest, a changed value in a log gives another digest, a
# changed output.txt in a snapshot directory gives the same digest, and so does a changed
# output of a failed call.
#
#     bash test/tspice/test_digest.sh
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
fail=0

# fake DIR VALUE OUTPUT LEVEL BAD: a recording of the family F_X with three calls. The
# first loads a kernel (an `S` line, and a snapshot directory with the kernel and the
# output of the run); the second has the input VALUE; the third fails, with the outputs
# BAD and, in a binary file named by the hash of its values, a long array of BAD. LEVEL
# is the compression level of the log.
fake() {
    local d=$1/logs bin
    mkdir -p "$d/F_X.log.files/1" "$d/F_X.log.files/bin"
    printf 'DAF/SPK \001\002\003' > "$d/F_X.log.files/1/k.bsp"
    printf '%s\n' "$3" > "$d/F_X.log.files/1/output.txt"
    ln -s 1 "$d/F_X.log.files/2"
    bin=$d/F_X.log.files/bin/$(printf '%s' "$5" | sha256sum | cut -c1-16).d
    printf '%s' "$5" | gzip > "$bin.gz"
    {
        printf '@ FURNSH 1 F\nS %s\n<file c 1\n|k.bsp\n! F \n' "$d/F_X.log.files/1"
        printf '@ VPACK 1 F\n<x d 1\n %s\n! F \n' "$2"
        printf '@ GETFOV 1 F\n<instid i 1\n-1000\nB %s\n>bounds d -1\n>shape c 1\n|%s\n! T SPICE(FRAMEMISSING)\n' "$bin.gz" "$5"
    } | gzip "-$4" > "$d/F_X.log.gz"
    printf 'F_X\tpassed\nF_Y\tpassed\n' > "$d/status.tsv"
}
digest() { python3 -I "$here/digest_recording.py" "$1" > "$1.tsv"; }

# the same content: other directory, compression and order of the families in status.tsv
fake "$tmp/a" 1.50000000000000000E+000 "run 1" 1 RECTANGLE
fake "$tmp/b" 1.50000000000000000E+000 "run 1" 9 RECTANGLE
printf 'F_Y\tpassed\nF_X\tpassed\n' > "$tmp/b/logs/status.tsv"
digest "$tmp/a"
digest "$tmp/b"
cmp -s "$tmp/a.tsv" "$tmp/b.tsv" || { echo "FAIL: the same recording in two directories gives two digests"; diff "$tmp/a.tsv" "$tmp/b.tsv" || true; fail=1; }

fake "$tmp/c" 1.50000000000000022E+000 "run 1" 1 RECTANGLE
digest "$tmp/c"
cmp -s "$tmp/a.tsv" "$tmp/c.tsv" && { echo "FAIL: a changed value in a log gives the same digest"; fail=1; }

fake "$tmp/d" 1.50000000000000000E+000 "run 2" 1 RECTANGLE
digest "$tmp/d"
cmp -s "$tmp/a.tsv" "$tmp/d.tsv" || { echo "FAIL: a changed output.txt gives another digest"; diff "$tmp/a.tsv" "$tmp/d.tsv" || true; fail=1; }

fake "$tmp/e" 1.50000000000000000E+000 "run 1" 1 'POLYGON?'
digest "$tmp/e"
cmp -s "$tmp/a.tsv" "$tmp/e.tsv" || { echo "FAIL: a changed output of a failed call gives another digest"; diff "$tmp/a.tsv" "$tmp/e.tsv" || true; fail=1; }

[ "$fail" = 0 ] && echo "digest_recording: 4 tests passed"
exit "$fail"
