#!/usr/bin/env bash
# Build the recording program of tspice_c, NAIF's test suite of CSPICE N0067, for
# record.sh.
#
#     test/tspice/cspice/build_cspice.sh OUTDIR
#
# Downloads NAIF's archive (below) into OUTDIR/tspice_c.tar unless it is there, checks
# its SHA-256, extracts it into OUTDIR/tspice_c, writes the wrappers and the driver with
# instrument_c.py into OUTDIR/gen, compiles CSPICE, its test utilities and the test
# families with NAIF's flags (gcc 12 or later), redirects the calls of the test program
# to the wrappers (objcopy), and links OUTDIR/tspice_c_rec. The last line of the output
# is the path of the program. Nothing here is committed to Almanac.jl. Environment: JOBS
# (parallel compilations, default 1).
#
# The program records the same calls with every compiler and C library, apart from the
# outputs that CSPICE leaves partly written in a call that fails. Two flags are added to
# NAIF's: -ftrivial-auto-var-init=zero starts each automatic variable at zero, because
# the test programs log some values that they never set (scratch arrays, unused inputs,
# outputs of a failed call), and -ffp-contract=off forbids the contraction of a product
# and a sum into one fused operation. The transcendental functions (sin, cos, exp, ...)
# come from CORE-MATH (crmath/), whose results are correctly rounded, because those of
# the C library differ between its versions in the last bit; after the link,
# check_math.sh refuses a program that takes another inexact function from the C
# library.
#
# The script also builds the oracle of the replay, OUTDIR/oracle/lib/libcspice.so and
# libcsupport.so: the CSPICE of the archive with the same flags and CORE-MATH, which
# run.jl has SPICE.jl load in place of the library of CSPICE_jll, so that the class of a
# call does not depend on the C library of the machine that replays either. Beside them,
# OUTDIR/oracle/lib/libcrmath.so holds the functions of CORE-MATH alone, under their own
# names (cr_sin, cr_cos, ...), which the replay calls in place of Julia's sin, cos, ...
# (crmath.jl), so that Almanac computes with the functions of the recording.
set -euo pipefail

[ $# -eq 1 ] || { echo "usage: build_cspice.sh OUTDIR" >&2; exit 1; }
here=$(cd "$(dirname "$0")" && pwd)
out=$(mkdir -p "$1" && cd "$1" && pwd)
url=https://naif.jpl.nasa.gov/pub/naif/misc/tspice/N0067/PC_Linux_64bit/tspice_c.tar
sha=373f3213cabab2717d21fc0110da703a713184d6cf73454372d0ae5a18c4adde
jobs=${JOBS:-1}
tar="$out/tspice_c.tar"
src="$out/tspice_c"

if ! echo "$sha  $tar" | sha256sum -c --quiet 2>/dev/null; then
    echo "downloading $url" >&2
    curl -fsSL --retry 3 -o "$tar.part" "$url"
    mv "$tar.part" "$tar"
    echo "$sha  $tar" | sha256sum -c --quiet || { echo "the archive does not have the SHA-256 $sha" >&2; exit 1; }
fi
# a new extraction also discards the objects: the sources keep their archive times, so
# the objects of an earlier archive would look newer than them
if [ ! -f "$src/.extracted" ] || [ "$src/.extracted" -ot "$tar" ]; then
    rm -rf "$src" "$out/obj"
    tar -xf "$tar" -C "$out"
    touch "$src/.extracted"
fi

python3 -B "$here/instrument_c.py" "$src" "$here/../routines.txt" "$out/gen" 2> "$out/gen.log" ||
    { cat "$out/gen.log" >&2; exit 1; }
grep -v '^mark ' "$out/gen.log" >&2 || true
cp "$here/tspice_log.c" "$here/tspice_log.h" "$out/gen/"

# NAIF's flags (mkprodct.csh of each directory) and the two flags above; -w: NAIF's
# sources are not ours to fix
cflags="-c -ansi -m64 -O2 -ffp-contract=off -ftrivial-auto-var-init=zero -DNON_UNIX_STDIO"
# the objects of other flags or of another gcc are discarded
stamp="gcc $(gcc -dumpfullversion) $cflags"
if [ "$(cat "$out/obj/.flags" 2>/dev/null)" != "$stamp" ]; then
    rm -rf "$out/obj"
    mkdir -p "$out/obj"
    printf '%s\n' "$stamp" > "$out/obj/.flags"
fi
compile() {   # DIR NAME [FLAG]: the .c files of DIR into OUTDIR/obj/NAME, then OUTDIR/obj/NAME.a
    local dir=$1 name=$2 flag=${3:-} o="$out/obj/$2"
    mkdir -p "$o"
    find "$dir" -maxdepth 1 -name '*.c' -print0 |
        xargs -0 -P "$jobs" -I{} bash -c 'f={}; b=$(basename "$f" .c); [ "$0/$b.o" -nt "$f" ] ||
            gcc '"$cflags $flag"' -w -I"$1" -I"$2" -o "$0/$b.o" "$f"' "$o" "$dir" "$src/include"
    rm -f "$out/obj/$name.a"
    ar crs "$out/obj/$name.a" "$o"/*.o
}
# CORE-MATH (gnu11: it uses asm); libm.c defines sin, cos, ... as its functions
crmath() {   # NAME [FLAG]: the functions of crmath/ into OUTDIR/obj/NAME
    local o="$out/obj/$1" f d
    mkdir -p "$o"
    for f in "$here"/crmath/*/*.c "$here/crmath/libm.c"; do
        d=$(basename "$(dirname "$f")")
        gcc -c -std=gnu11 -m64 -O2 -ffp-contract=off -w ${2:-} -I"$(dirname "$f")" \
            -o "$o/${d}_$(basename "$f" .c).o" "$f"
    done
}
echo "compiling CSPICE and the test families with $jobs jobs" >&2
compile "$src/src/cspice" cspice
compile "$src/src/csupport" csupport
compile "$src/src/tutils_c" tutils_c
compile "$src/src/tspice_c" tspice_c
# the generated files, warnings on
mkdir -p "$out/obj/gen"
for f in driver tlog_wrap tspice_log; do
    inc="-I$out/gen -I$src/include"
    [ "$f" = driver ] && inc="$inc -I$src/src/tutils_c"
    gcc $cflags -Wall $inc -o "$out/obj/gen/$f.o" "$out/gen/$f.c"
done
# the calls of the test program (the families, NAIF's test and support utilities) go to
# the wrappers; CSPICE's calls between its own routines do not
for lib in tspice_c tutils_c csupport; do
    objcopy --redefine-syms="$out/gen/wrap.syms" "$out/obj/$lib.a" "$out/obj/$lib.wrapped.a"
done
crmath crmath
gcc -m64 -o "$out/tspice_c_rec" "$out"/obj/gen/{driver,tlog_wrap,tspice_log}.o \
    -Wl,--start-group "$out"/obj/{tspice_c,tutils_c,csupport}.wrapped.a "$out/obj/cspice.a" -Wl,--end-group \
    "$out"/obj/crmath/*.o -lm
bash "$here/check_math.sh" "$out/tspice_c_rec" || exit 1

# the oracle: the two libraries that CSPICE_jll provides, built from the archive, each
# with CORE-MATH inside, where -Bsymbolic binds the library to it and --exclude-libs
# keeps it from the process that loads the library; erf_, erfc_, derf_ and derfc_ (erf
# and erfc of the f2c library, which CSPICE does not call) are left out
echo "compiling the oracle" >&2
compile "$src/src/cspice" cspice_pic -fPIC
compile "$src/src/csupport" csupport_pic -fPIC
crmath crmath_pic -fPIC
rm -f "$out/obj/crmath_pic.a"
ar crs "$out/obj/crmath_pic.a" "$out"/obj/crmath_pic/*.o
objs=()
for o in "$out"/obj/cspice_pic/*.o; do
    case ${o##*/} in erf_.o | erfc_.o | derf_.o | derfc_.o) ;; *) objs+=("$o") ;; esac
done
mkdir -p "$out/oracle/lib"
gcc -shared -Wl,-soname,libcspice.so -Wl,-Bsymbolic -o "$out/oracle/lib/libcspice.so" "${objs[@]}" \
    -Wl,--whole-archive "$out/obj/crmath_pic.a" -Wl,--no-whole-archive -Wl,--exclude-libs,ALL -lm
gcc -shared -Wl,-soname,libcsupport.so -Wl,-Bsymbolic -Wl,-rpath,'$ORIGIN' -o "$out/oracle/lib/libcsupport.so" \
    "$out"/obj/csupport_pic/*.o "$out/obj/crmath_pic.a" -Wl,--exclude-libs,ALL -L"$out/oracle/lib" -lcspice -lm
# the functions of CORE-MATH for the replay, without libm.c (crmath_libm.o), so that the
# library exports cr_sin, cr_cos, ... and not sin, cos, ...: one cr_ function for each
# directory of crmath/, and nothing else but the symbols that the linker defines
objs=()
for o in "$out"/obj/crmath_pic/*.o; do
    [ "${o##*/}" = crmath_libm.o ] || objs+=("$o")
done
gcc -shared -Wl,-soname,libcrmath.so -Wl,-Bsymbolic -o "$out/oracle/lib/libcrmath.so" "${objs[@]}" -lm
want=$(for d in "$here"/crmath/*/; do echo "cr_$(basename "$d")"; done | LC_ALL=C sort)
have=$(nm -D --defined-only "$out/oracle/lib/libcrmath.so" | awk '{print $3}' | sed 's/@.*//' |
    grep -vxE '_init|_fini|_edata|_end|__bss_start' | LC_ALL=C sort -u)
[ "$have" = "$want" ] || {
    echo "libcrmath.so exports $(echo $have), not $(echo $want)" >&2
    exit 1
}
for lib in libcspice libcsupport libcrmath; do
    bash "$here/check_math.sh" "$out/oracle/lib/$lib.so" || exit 1
done
echo "$out/tspice_c_rec"
