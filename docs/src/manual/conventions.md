# Conventions

Almanac.jl keeps the names, the units and the numerical results of CSPICE N0067, and
its argument order with the few exceptions below. It departs from CSPICE where the C
interface relies on hidden state, on output arguments or on arguments that only manage
C memory.

Every docstring links to the NAIF page of the routine it ports: the `_c` page of
CSPICE, or the SPICELIB page for a Fortran routine that CSPICE does not wrap. The NAIF
page is the full specification. The docstring gives the Julia signature, the units and
what differs.

## Names, arguments and results

- **Names.** The `_c` suffix is dropped: `spkezr_c` is [`spkezr`](@ref). A call that
  changes one of its arguments ends in `!`, as [`wninsd!`](@ref).
- **Results** are returned, never written to output arguments. Several results come
  back as a tuple: `state, lt = spkezr(...)`.
- **Arrays.** Vectors and matrices are returned as `SVector` and `SMatrix` of
  StaticArrays.jl. Inputs can be any `AbstractVector` or `AbstractMatrix`.
- **Strings** are returned as a `String`, so no call takes a string length.
- **Memory arguments.** Counts, sizes, `room` and workspaces, which serve only C
  memory, are absent: [`gdpool`](@ref) takes only the name of the variable, and
  [`gfdist`](@ref) has no `nintvls`.
- **Indices** count from 1, where CSPICE counts from 0: window intervals, EK rows.
- **Cells.** A CSPICE cell becomes a `Vector`, kept sorted and without duplicates when
  it is a set, or a [`Window`](@ref). Both grow as needed and have no declared size. The calls that fill one ([`spkobj!`](@ref), [`spkcov!`](@ref), [`ckobj!`](@ref),
  [`ckcov!`](@ref), [`pckfrm!`](@ref), [`pckcov!`](@ref), [`dskobj!`](@ref),
  [`dsksrf!`](@ref)) take it as the first argument; CSPICE takes it last. SPICE.jl
  also puts it first, except in `ckobj!` and `ckcov!`.

## Contexts in place of global state

Loaded kernels, the kernel pool and settings such as [`timdef`](@ref) belong to a
[`SpiceContext`](@ref), not to the process. A call that needs them takes a context as
an optional first argument, written `[ctx,]` in the reference, and uses
[`default_context`](@ref) without one. `[ctx | ks,]` means that the call also accepts a
`KernelSet`, the immutable snapshot of a context that `Almanac.kernelset(ctx)` returns.
[Contexts and threads](@ref) shows contexts in use.

## Errors and `found` flags

- A call that fails throws a [`SpiceError`](@ref), in place of setting the error flag
  of CSPICE. Its `kind` is the CSPICE short message as a lower-case symbol
  (`SPICE(INVALIDDEGREE)` is `:invaliddegree`), and its `msg` is the text. Conditions
  that CSPICE does not check have kinds of Almanac.jl, such as `:badarg`.
- Missing SPK data throw the subtype [`SpkInsufficientDataError`](@ref), whose fields
  name the bodies and the epoch, and say whether the missing body has data at other
  epochs.
- A lookup that has a `found` flag in CSPICE returns `nothing` when nothing is found:
  [`bodn2c`](@ref) for an unknown name, [`ckgp`](@ref) without pointing data. A call
  that must return a value throws instead, as [`gdpool`](@ref) for a missing variable.

## File handles and files open for write

- A file handle is an object, a [`DAFHandle`](@ref) or a [`DASHandle`](@ref), passed
  to every call. No call acts on a hidden "current" file.
- A file open for write can also be read, by [`dafopr`](@ref), [`dasopr`](@ref),
  [`ekopr`](@ref) or [`furnsh`](@ref). CSPICE refuses this (`RWCONFLICT`). It is safe
  because a DAF array is complete before its summary is written, and a DAS file changes
  on disk only at [`dascls`](@ref).
- A loaded kernel can be opened for write, for instance by [`spkopa`](@ref). A context
  sees only the segments written before it loaded the file, so it sees new segments
  only after it loads the file again. A second write open of one
  file still throws `SpiceError(:fileopenconflict)`.

## No state between calls

CSPICE keeps buffers and the results of earlier searches between calls, for speed.
Almanac.jl keeps none: a query is a function of the loaded kernels and of its
arguments. This has three consequences.

- **Threads.** Queries take no lock, and kernels can be loaded while other threads
  query; see [Contexts and threads](@ref).
- **Settings.** [`timdef`](@ref), [`tsetyr`](@ref), [`tparch`](@ref),
  [`gfstol`](@ref) and [`dskstl`](@ref) change a context, not the process. The value
  set by `gfstol` is the default of the `tol` keyword of every search.
- **Repeatable answers.** Some CSPICE results depend on an earlier call: `ckgp` right
  after `ckgpav` on overlapping segments, the dates of a text kernel loaded after
  certain `str2et` calls, or the four padding bytes of a binary PCK segment summary,
  which CSPICE fills from the last SPK or CK summary it handled. Almanac.jl always
  returns what CSPICE returns in a new process.

## Porting from SPICE.jl

SPICE.jl already returns results, ends mutating calls in `!`, returns `nothing` for a
failed lookup and counts from 1. A program written for SPICE.jl changes in these
places:

- Results are static arrays, which cannot be changed in place.
- The error type is `Almanac.SpiceError`, with its `kind` field, not
  `SPICE.SpiceError`.
- Windows are a [`Window`](@ref), not a cell.
- Calls take an optional context.
- The memory arguments that SPICE.jl keeps are gone: `start` and `room` of `gdpool`,
  `nintvls` of `gfdist`, `maxn` of `bodvrd`.
- A call that must return a value throws where SPICE.jl returns `nothing`, as `gdpool`
  for a missing variable.
- Where a wrapper of SPICE.jl 0.2.5 has a defect, Almanac.jl gives the result of
  CSPICE: [`ckw02`](@ref) and [`ckw03`](@ref) write every record (SPICE.jl passes one
  fewer); [`pckw02`](@ref) has no degree argument and takes the degree from the record
  length (SPICE.jl passes one three times too large); [`spkw02`](@ref) and
  [`spkw03`](@ref) count the records from the right dimension; [`dskv02`](@ref) and
  [`dskp02`](@ref) count `start` from 1, as the routines they port; [`termpt`](@ref)
  returns only the points found; [`et2utc`](@ref) returns the whole string, which
  SPICE.jl can cut when `prec` is small.
- [`et2utc`](@ref) and [`etcal`](@ref) return a [`FixedString`](@ref Almanac.FixedString),
  an `AbstractString` stored inline, so that they allocate nothing; `String(s)` gives a
  `String`.
