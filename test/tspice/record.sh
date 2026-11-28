#!/usr/bin/env bash
# Record the calls that NAIF's test suite (tspice) makes to SPICELIB, for replay.jl.
#
#     test/tspice/record.sh OUTDIR [FAMILY ...]
#
# Builds the recording program of tspice_c, NAIF's test suite of CSPICE N0067, with
# cspice/build_cspice.sh (gcc; cspice/instrument_c.py writes the wrappers of the routines
# of routines.txt), and runs each f2c-translated family in its own directory with the log
# OUTDIR/logs/F_XXX.log.gz. Nothing here is committed to Almanac.jl: NAIF's archive and
# its build, the logs and the copies of the kernels that the families write stay in
# OUTDIR.
#
# A full recording (no FAMILY) writes OUTDIR/logs/manifest.tsv: the source, the SHA-256
# of the archive, the hashes of the files that made the program (the build script, the
# wrappers, the logger and CORE-MATH), and the list of the families; run.jl checks a
# recording against it. Each family adds its status to OUTDIR/logs/status.tsv. The exit
# status is 1 if a family failed or timed out.
#
# Environment: JOBS (parallel families, default: the CPU threads, at most one per 2 GB of
# available memory), TIMEOUT (seconds per family, default 1800), TSPICE_ALL=1 (log every
# call; by default a routine without state, or a query between two changes of state, is
# logged once per distinct input, and a routine without state at most TSPICE_MAXPURE
# times per family, default 20000; see cspice/tspice_log.h).
set -euo pipefail

[ $# -ge 1 ] || { echo "usage: record.sh OUTDIR [FAMILY ...]" >&2; exit 1; }
here=$(cd "$(dirname "$0")" && pwd)
out=$(mkdir -p "$1" && cd "$1" && pwd)
shift
memjobs=$(awk '/MemAvailable/ {print int($2 / 2097152)}' /proc/meminfo 2>/dev/null || echo 0)
jobs=${JOBS:-$(n=$(nproc); [ "$memjobs" -ge 1 ] && [ "$memjobs" -lt "$n" ] && n=$memjobs; echo "$n")}
echo "recording with $jobs parallel families"
timeout=${TIMEOUT:-1800}

exe=$(JOBS=$jobs bash "$here/cspice/build_cspice.sh" "$out" | tail -1)
makers="cspice/build_cspice.sh cspice/check_math.sh cspice/instrument_c.py cspice/tables.py
cspice/tspice_log.c cspice/signatures.tsv routines.txt
$(cd "$here" && find cspice/crmath -type f | LC_ALL=C sort)"

mkdir -p "$out/logs" "$out/work"
if [ $# -gt 0 ]; then
    families=$(printf '%s\n' "$@" | tr 'a-z' 'A-Z' | sed 's/^\(F_\)\{0,1\}/F_/')
else
    families=$("$exe" --list)
    # a full recording replaces the manifest and the status of every family
    {
        printf 'source\tcspice\n'
        printf 'tspice_c.tar\t%s\n' "$(sha256sum "$out/tspice_c.tar" | cut -d' ' -f1)"
        for f in $makers; do
            printf '%s\t%s\n' "$f" "$(sha256sum "$here/$f" | cut -d' ' -f1)"
        done
        printf 'TSPICE_ALL\t%s\n' "${TSPICE_ALL:-0}"
        printf 'TSPICE_MAXPURE\t%s\n' "${TSPICE_MAXPURE:-20000}"
        printf '%s\n' $families | sed 's/^/family\t/'
    } > "$out/logs/manifest.tsv"
    : > "$out/logs/status.tsv"
fi

# The log is written into a FIFO that gzip compresses to F_XXX.log.gz.
run_family() {
    local fam=$1 wd="$out/work/$1" log="$out/logs/$1.log" ok=0 gz
    rm -rf "$wd" "$log" "$log.gz" "$log.files"
    mkdir -p "$wd"
    mkfifo "$log"
    gzip -1 < "$log" > "$log.gz" &
    gz=$!
    (cd "$wd" && TSPICE_LOG="$log" timeout "$timeout" "$exe" "$fam" > output.txt 2>&1) || ok=$?
    # a family that wrote nothing never opened the FIFO: let gzip see its end
    kill -0 "$gz" 2>/dev/null && timeout 5 sh -c ": > '$log'" 2>/dev/null
    wait "$gz"
    rm -f "$log"
    # one short line per family: `>>` of it is atomic between the parallel families
    if [ "$ok" = 0 ]; then
        echo "$fam passed"
        printf '%s\tpassed\n' "$fam" >> "$out/logs/status.tsv"
    else
        [ "$ok" = 124 ] && why=timeout || why="exit $ok"
        echo "$fam FAILED ($why, see $wd/output.txt)"
        printf '%s\tfailed\t%s\n' "$fam" "$why" >> "$out/logs/status.tsv"
        return 1
    fi
}
export -f run_family
export out exe timeout TSPICE_ALL TSPICE_MAXPURE
status=0
printf '%s\n' $families | xargs -P "$jobs" -I{} bash -c 'run_family {}' || status=1
[ "$status" = 0 ] || echo "some families failed: see $out/logs/status.tsv"
exit "$status"
