# NAIF's test suite (tspice)

This directory replays against Almanac.jl the calls that tspice, NAIF's test suite, makes
to SPICELIB, and checks each difference against CSPICE. SPICELIB is here the Fortran
library translated to C by f2c (`spkezr_`, ...), which the test programs call; CSPICE is
its C interface (`spkezr_c`, ...), which SPICE.jl calls and some of whose routines
(`vperp_c`, `nplnpt_c`, ...) are written in C, so that they can differ from SPICELIB in
the last bits. Both come from the same archive of NAIF. The run is not part of the
default test run.

```sh
julia --project=test test/tspice/run.jl
```

## Recording

The recording runs the 400 families of tspice translated to C by f2c, from NAIF's test
suite of CSPICE N0067
([`tspice_c.tar`](https://naif.jpl.nasa.gov/pub/naif/misc/tspice/N0067/PC_Linux_64bit/tspice_c.tar)),
against the library of the same archive. It needs gcc 12 or later, binutils (`objcopy`,
`nm`), `cmp`, Python 3 and curl. `cspice/build_cspice.sh` downloads the archive (checked
by its SHA-256), writes the wrappers with `cspice/instrument_c.py` and builds the
program; only the calls of the test program (the families and NAIF's test and support
utilities) go through the wrappers. `cspice/tables.py` lists how each routine is
recorded, and `cspice/tspice_log.h` describes the log. `cspice/signatures.tsv` holds the
kind, direction and dimensions of the arguments of each recorded routine, which the C
prototypes do not give: after a change of `routines.txt`, write it again with
`cspice/signatures.py` from the Fortran sources of SPICELIB (`tspice.tar` in the same
directory of NAIF).

The first run builds the program, runs the families and records their calls into
`~/.cache/almanac-tspice` (about 2 GB, 15 to 40 minutes), with a manifest of the
families and of what made the recording. Later runs only replay the recording, in
parallel processes. A recording whose manifest names a source other than tspice_c is
refused: record again with `--record`.

To keep the recording small, a call of a routine that has no state is kept once for
each distinct set of inputs, and a query is kept once between two changes of the loaded
kernels or of the settings. `TSPICE_ALL=1` keeps every call.

Apart from the outputs that SPICELIB leaves partly written in a call that fails, a
recording depends only on the sources of the harness and of NAIF's archive, so that a
recording made on one machine is the same as one made on another. The transcendental
functions of the C library (`sin`, `cos`, `exp`, ...) differ between its versions in the
last bit, which changes the recorded values and, through the rule above, the calls that
are kept; the program therefore takes them from CORE-MATH, whose results are correctly
rounded, and `cspice/check_math.sh` refuses a program that takes from the C math library
a function other than those whose results IEEE 754 defines exactly (`sqrt`, `fma`,
`floor`, ...). `cspice/crmath/` holds the sources of
[CORE-MATH](https://gitlab.inria.fr/core-math/core-math) at commit
`b464d7d0280f4c5091b8d2fff80b6cb03ab1d37f` (MIT license, `cspice/crmath/LICENSE`) for
the functions that the library calls. Some values that the test programs log are never set
(scratch arrays, unused inputs, outputs of a failed call), and the program is compiled
with `-ftrivial-auto-var-init=zero`, which starts each automatic variable at zero
rather than with what the stack held. The values that SPICELIB writes into an output
before a call fails can still differ between two machines, and the digest below leaves
them out.

Two recordings `A` and `B` are the same when their digests are equal:

```sh
diff <(python3 -I test/tspice/digest_recording.py A) <(python3 -I test/tspice/digest_recording.py B)
```

The digest lists the SHA-256 of each file of `logs`: for a log, that of its calls as the
replay reads them, without the outputs of the calls that failed; for another compressed
file, that of its decompressed content; and in each case without the directory in which
the recording was made. It leaves out the files that record the run, not its calls: the
output and the test log of each family (`output.txt` and `cspice*.log` in the snapshot
directories), which hold its dates, and the files `.last.*` of the logger. The workflow
`tspice-reproducible` records the full suite twice in parallel, on ubuntu-latest and in
a Fedora container with another gcc and another version of the C library, and fails if
the two digests differ; it runs on a pull request that changes the recorder, and by
hand.

A difference from SPICELIB is checked against an oracle: CSPICE, which
`cspice/build_cspice.sh` builds from the same archive into `oracle/lib` of the recording
directory with the flags and the math of the recording program. SPICE.jl loads it in
place of the library of CSPICE_jll, so that the class of a call does not depend on the C
library of the machine that replays either. To load it, `run.jl` writes an override of
the artifact of CSPICE_jll into `oracle_depot` of the recording directory and adds that
depot at the end of the depot path of the replay processes, and `replay.jl` refuses to
classify with another library. The oracle is reset before each family, together with the
watchers of its kernel pool and the body names that `boddef` defined. CSPICE is asked
once for each call, and not at all where its answer would not be that of CSPICE to the
call of the test: for a routine that SPICE.jl replaces with Julia code (`vnorm`,
`vzerog`, ...), for an error that SPICELIB makes because of the size of an argument
(`room`, `maxn`, the length of an output string), which the replay does not give
CSPICE, and for CVPOOL, whose answer depends on calls that CSPICE does not see.

The replay gives Almanac.jl the same functions. Julia's own `sin`, `cos`, `exp`, ... differ
from the correctly rounded results in the last bit for some arguments, and an output of
Almanac.jl that depends on one of them then differs from SPICELIB's for a reason outside
Almanac.jl. Before it loads Almanac.jl, `replay.jl` therefore includes `crmath.jl`, which
replaces the `Float64` methods of `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `exp`,
`log`, `sinh`, `cosh`, `tanh`, `sincos`, `atan(y, x)` and `x^y` (for `Float64` `x` and
`y`) with calls of
`oracle/lib/libcrmath.so` (the functions of CORE-MATH alone, which
`cspice/build_cspice.sh` builds beside the oracle); each of these methods throws a
`DomainError` where Julia's method throws one. The outputs are then compared at
tolerance zero, GF windows included, so that a change of the last bit in the arithmetic
of Almanac.jl is a mismatch. Almanac.jl also calls functions of Julia that are not in
this list (`log10`, `cospi`, `x^n` for an integer `n`), and the replay compares their
results as it compares any other, so that a difference that one of them causes is a
mismatch too.

`--julia-math` replays with Julia's functions, as a user of Almanac.jl has them, and
shows how far those results are from SPICELIB's. It compares an output within 1e-14
absolute (in the unit of the output) plus 1e-12 relative to the larger magnitude, and GF
windows within twice the convergence tolerance of the search. Its mismatches are not
compared with the baseline, which is made with CORE-MATH, and the run exits 0 unless it
stops. `TSPICE_RTOL` and `TSPICE_ATOL` replace the relative and the absolute tolerance,
also without `--julia-math`; a run with either of them set does not compare its
mismatches with the baseline either.

The recording also copies each file that a test closes after writing it, so that the
replay compares the file that Almanac.jl wrote with SPICELIB's, by content (arrays and
summaries, DAS data, EK entries), at the close. `files.tsv` (in `results/worker*`) lists
each such file and whether it was compared, or why not. Past 256 MB of copies in a
family, only the copy of the last close of each file is kept: a test that appends to
one file and closes it after each segment has its file compared once, with every
segment, at the last close. An input array too long for the log goes to a compressed
binary file beside it (at most `TSPICE_MAXBIN` bytes, default 256 MB, per family).

The log also names the file of each handle of a file open for reading. A test gets some
handles from routines outside the recording (the DSK segment search of ZZDSKBSR), and
the replay finds them by their file, if Almanac.jl has it loaded.

Windows and sets (cells) are logged as their size, their cardinality and their elements.
The GF searches that take procedures (GFEVNT, GFFOVE, GFOCCE) are replayed with the step
that GFSSTP or the last search set and with Almanac.jl's `gfrefn`, which are the
procedures that every tspice call passes. Not replayed:

- GFUDB, GFUDS, UDDC and UDDF, whose functions are written in the test program;
- writers without a counterpart in Almanac.jl (the SPKSnn subset writers, DAFWDR and
  other record-level writers). A file that one of them wrote is read from
  the copy of SPICELIB's file, and its later writes are skipped;
- calls on a file that the test opened with routines outside the recording;
- DASCUD and the DAS calls of F_DASCUD after it: Almanac.jl keeps a DAS file in memory
  and writes its directories at the close, so it has no step that updates them;
- calls of the writers that Almanac.jl's arguments cannot express, such as a degree
  that disagrees with the length of the data (Almanac.jl takes the degree from the
  records).

`TSPICE_WHY=1` writes `unmapped.tsv` (in `results/worker*`): for each family and routine,
the calls skipped for want of a file, and why the replay has none.

## When a run stops or fails

A run without a list of families is the full check.

It **stops** (exit code 2), without a comparison with the baseline, if a family is
missing from the recording, failed to record or fails to replay to its end (a log that
cannot be read, for example), if the manifest names another source than
tspice_c, or if the recording has no oracle or no `oracle/lib/libcrmath.so`
(`--julia-math` needs the oracle but not `libcrmath.so`).

It **fails** (exit code 1) if:

- a mismatch is not in the baseline `baseline_ids.tsv.gz`;
- a routine has more calls `skipped`, `near` or `spread` than `coverage.tsv` records, or
  fewer calls compared.

Each replayed call has one of these classes, and each class but the first two is a
mismatch:

- `match`: Almanac.jl gives SPICELIB's outputs to the last bit (with `--julia-math`,
  within the tolerances) or throws SPICELIB's error kind;
- `cspice`: CSPICE differs from SPICELIB, and Almanac.jl gives CSPICE's outputs to the
  last bit (with `--julia-math`, within 1e-9 relative) or its error kind;
- `near`: the outputs are within 1e-9 relative of SPICELIB's;
- `spread`: the outputs differ by more, SPICELIB and CSPICE differ from each other, and
  Almanac.jl is no farther from CSPICE than twice SPICELIB's distance;
- `value`: the outputs differ by more;
- `error`: one side signals an error and the other does not, or the error kinds differ;
- `adapter`: the replay cannot make the call, a fault of `replay.jl` and not of
  Almanac.jl.

Each mismatch is identified by family, routine, class and a hash of its inputs, so a new
mismatch fails even when another one went away, and a call that moves from `match` to
`near` fails too, but not one that moves to `cspice`.

The mismatches and the `cspice` calls are written to `results/replay.tsv` and
`results/ids.tsv`.

`julia --project=test test/tspice/selftest.jl [DIR]` checks the gate itself, on F_LTIME
and F_POOL of a copy of the recording `DIR` (in a few minutes). It requires the run to
fail (exit 1) when it removes a mismatch from the baseline, moves each output of
Almanac.jl by one unit in the last place, or changes the kind of its errors, and to stop
(exit 2) when it removes the oracle or `oracle/lib/libcrmath.so`, changes the source in
the manifest, deletes a log, damages a log or marks a family as failed. It also requires
that F_TANGPT3 has the same mismatches after F_DWPOOL in one process as alone, which
shows that the reset of the oracle leaves no state behind.

## The baseline

`baseline_ids.tsv.gz` lists the mismatches that the gate accepts. `baseline.tsv` counts
the `value`, `error` and `adapter` calls of each routine that has some, and gives each
such routine a label (column `class`) and a note, which `--update-baseline` keeps; a
`near` or `spread` call has an identity in `baseline_ids.tsv.gz` and no row in
`baseline.tsv`. The labels are:

- `documented`: a difference that a docstring or the manual states;
- `cspice`: Almanac.jl follows CSPICE where SPICELIB differs, in calls that the oracle
  does not class as `cspice` (the error kind of a CSPICE wrapper for an empty name, for
  example);
- `gap`: a call that the replay cannot reproduce, such as a read of a file that the
  test wrote with routines outside the recording;
- `unexplained`.

## Options

- `--record` records again.
- `--julia-math` replays with Julia's functions (see above), without the comparison with
  the baseline and the coverage check.
- `--dir DIR` takes the recording, the oracle and the results from `DIR` (default
  `~/.cache/almanac-tspice`); `A` and `B` of the digest and `DIR` of `selftest.jl` are
  such directories.
- `--update-baseline` makes the full run the new baseline (`baseline.tsv`,
  `baseline_ids.tsv.gz` and `coverage.tsv`) and keeps the labels and notes.
- `--baseline-from DIR` does the same from the `results` directory of a full run made
  elsewhere (the `tspice` workflow uploads one), without a replay.
- A list of families (`F_SPK08 F_CK06`) restricts the run to them, without the coverage
  check.
