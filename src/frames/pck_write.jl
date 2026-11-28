# Binary PCK writer: `pckopn`, `pckw02`, `pck03b`/`pck03a`/`pck03e`, `pckw20` and
# `pckcls`. A binary PCK is a DAF of the type
# PCK whose summaries hold 2 doubles (the coverage) and 5 integers (frame class ID of the
# body-fixed frame, base frame, data type, first and last address). `pckw02` makes the
# input checks of SPICELIB `pckw02` in the same order and stores the records of
# `spkw02` (`_seg_cheb_records`), so the file is equal byte for byte to the one that
# CSPICE writes from the same input.

"""
    _pck_irfnum(frame) -> Int

ID of the built-in inertial frame `frame` (a name or an ID), or 0. A name is compared
as `irfnum_` compares it: letter case and blanks do not count.
"""
function _pck_irfnum(frame::SegFrame)
    frame isa Integer && return 1 <= frame <= N_INERTIAL_FRAMES ? Int(frame) : 0
    key = uppercase(replace(frame, ' ' => ""))
    for (i, def) in enumerate(INERTIAL_FRAME_DEFS)   # invariant: no frame before i has the name
        def[1] == key && return i
    end
    return 0
end

"""
    pckopn(path, ifname, ncomch=0) -> DAFHandle

Create the binary PCK file `path` and open it for writing, as
`pckopn_c`.

`ifname` is the internal file name and `ncomch` is the number of comment characters to
make room for. The file is on disk when the call returns, byte for byte the file that
`pckopn_c` creates. An existing file is never replaced. Add segments with
[`pckw02`](@ref) and close the handle with [`pckcls`](@ref). Load the file with
[`furnsh`](@ref) only after the close.

Throws `SpiceError(:emptystring)` if `path` or `ifname` is empty and `:blankfilename` if
`path` is blank, as CSPICE, and `:fileopenfail` if the file exists or cannot be created.

# Examples
Write a segment for the frame IAU_MOON (frame class ID 301) whose angles are constant,
then load the file:
```jldoctest
julia> path = joinpath(mktempdir(), "example.bpc");

julia> h = pckopn(path, "example PCK");

julia> cdata = [[0.5, 0.25, 1.0]];    # phi, delta, w (rad), degree 0

julia> pckw02(h, 301, "J2000", 0.0, 100.0, "moon", 100.0, cdata, 0.0)

julia> pckcls(h)

julia> ctx = SpiceContext();

julia> furnsh(ctx, path)

julia> pxform(ctx, "J2000", "IAU_MOON", 50.0) ≈ eul2m(1.0, 0.25, 0.5, 3, 1, 3)
true
```

See also [`dafac`](@ref) to add comments.

# References

- [NAIF documentation of `pckopn_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pckopn_c.html)
"""
pckopn(path::AbstractString, ifname::AbstractString, ncomch::Integer = 0) =
    _seg_opn("pckopn", "PCK", 5, path, ifname, ncomch)

"""
    pckops(ifname="", ncomch=0) -> DAFHandle

Open a new PCK file that lives in memory, for writing: [`pckopn`](@ref) without the
path. No file is made. The writers (such as [`pckw02`](@ref)) take the handle as they take
one of a file and write the same bytes. [`furnsh`](@ref) loads its segments as a kernel,
and [`save_kernel`](@ref) writes them to a file; [`pckcls`](@ref) discards them, so load or
save before the close. This call
has no CSPICE counterpart; see [`dafops`](@ref).

Throws the errors of [`dafops`](@ref).
"""
pckops(ifname::AbstractString = "", ncomch::Integer = 0) = _seg_ops("PCK", 5, ifname, ncomch)

"""
    pckw02(handle, clssid, frame, first, last, segid, intlen, cdata, btime) -> Nothing

Write a binary PCK segment of type 2 (Chebyshev polynomials for three Euler angles, on
intervals of equal length) to the file of `handle`, as
`pckw02_c`.

- `handle`: a [`DAFHandle`](@ref) from [`pckopn`](@ref) or [`pckops`](@ref).
- `clssid`: the frame class ID of the body-fixed frame, the `classid` of
  [`frinfo`](@ref), for example 301 for IAU_MOON.
- `frame`: the base frame, a built-in inertial frame given by name or by frame ID (1 to
  21, which [`bltfrm`](@ref) lists). Letter case and blanks in the name do not count.
- `first`, `last`: the coverage of the segment, seconds past J2000 TDB.
- `segid`: the segment identifier, at most 40 printable characters.
- `intlen`, `btime`: the length of every interval (s) and the start of the first interval
  (seconds past J2000 TDB).
- `cdata`: `n` records, one for each interval. All records have the same length, a
  multiple of 3: the coefficients (rad) of the angle `phi`, then of `delta`, then of
  `w`, a third of the record for each, so the polynomial degree is a third of the
  length minus 1. The matrix from the base frame to the body-fixed frame is
  `R3(w) * R1(delta) * R3(phi)`, which is `eul2m(w, delta, phi, 3, 1, 3)`.

`cdata` is a vector of vectors or a matrix with one record in each column, as in
[`spkw02`](@ref).

!!! warning "Difference from CSPICE"
    The call has no arguments for the record count `n` and the degree, which CSPICE and
    SPICE.jl take: both follow from the shape of `cdata`.

The segment is on disk when the call returns. A call that throws writes nothing.

Throws `SpiceError`:
- `:emptystring` if `frame` or `segid` is empty;
- `:dafnosuchhandle` if the handle is closed, `:dafinvalidaccess` if it is open for
  reading;
- `:numcoeffsnotpos` if `cdata` has no record (an empty vector of vectors or a matrix
  without columns);
- `:badarg` if the records differ in length, or if their length is not a positive
  multiple of 3;
- `:intlennotpos` unless `intlen > 0`;
- `:invalidrefframe` if `frame` is not a built-in inertial frame;
- `:baddescrtimes` if `first > last`, or unless [`first`, `last`] is inside [`btime`,
  `btime + n * intlen`] (no tolerance);
- `:idstringtoolong` or `:nonprintablechars` if `segid` is not valid.

See [`pckopn`](@ref) for an example.

# References

- [NAIF documentation of `pckw02_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pckw02_c.html)
"""
function pckw02(
        h::DAFHandle, clssid::Integer, frame::SegFrame, first::Real, last::Real,
        segid::AbstractString, intlen::Real, cdata, btime::Real
    )
    who = "pckw02"
    cd, n, ninrec, first, last, intlen, btime = _seg_cheb_args(
        who, "Euler angle", 3, h, frame, first, last, segid, intlen, cdata, btime
    )
    intlen <= 0.0 && throw(SpiceError(:intlennotpos, "$who: the interval length is not positive: $intlen"))
    refcod = _pck_irfnum(frame)
    refcod == 0 && throw(
        SpiceError(
            :invalidrefframe,
            "$who: the reference frame $frame is not supported"
        )
    )
    first > last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is greater than the segment end time $last"
        )
    )
    first < btime && throw(
        SpiceError(
            :baddescrtimes, "$who: the segment descriptor " *
                "start time $first is less than the beginning time of the segment data $btime"
        )
    )
    ltime = btime + n * intlen
    last > ltime && throw(
        SpiceError(
            :baddescrtimes, "$who: the segment descriptor end " *
                "time $last is greater than the end time of the segment data $ltime"
        )
    )
    _seg_chckid(who, segid, "PCK")
    _daf_add_array!(
        h, Float64[first, last], Int[clssid, refcod, 2, 0, 0], segid,
        _seg_cheb_records(cd, n, ninrec, intlen, btime)
    )
    return nothing
end

"""
    PCK03Writer

A binary PCK segment of type 3 that is open for the addition of records.
[`pck03b`](@ref) returns it, and it is the first argument of [`pck03a`](@ref) and
[`pck03e`](@ref).
"""
mutable struct PCK03Writer
    const handle::DAFHandle
    const descr::Vector{Float64}
    const segid::String
    const ncoef::Int
    const packets::Vector{Float64}      # (reference, packet) pairs, as `sgwfpk` writes them
    const refs::Vector{Float64}
    open::Bool
end

Base.show(io::IO, w::PCK03Writer) = _seg_show_writer(io, w)

"""
    pck03b(handle, segid, body, frame, first, last, chbdeg) -> PCK03Writer

Begin a binary PCK segment of type 3 (Chebyshev polynomials for RA, DEC and W and for
the angular velocity, on intervals of unequal length) in the file of `handle`, as
`pck03b`
of SPICELIB, and return its writer. CSPICE has no C wrapper for it.

`body` is the frame class ID of the body-fixed frame, as `clssid` of [`pckw02`](@ref),
and `frame` the base frame, a built-in inertial frame given by name or by ID. `first`
and `last` are the coverage (seconds past J2000 TDB) and `chbdeg` the degree of the
polynomials. Add records with [`pck03a`](@ref) and end the segment with
[`pck03e`](@ref).

Throws `SpiceError`:
- `:invalidargument` if `chbdeg < 0`;
- `:barycenteridcode` if `0 <= body <= 9`;
- `:invalidrefframe` if `frame` is not a built-in inertial frame;
- `:baddescrtimes` unless `first < last`;
- `:dafnosuchhandle` or `:dafinvalidaccess` if the file is not open for writing;
- `:segidtoolong` or `:nonprintingchars` (not `:nonprintablechars` as in
  [`ckw04b`](@ref)) if `segid` is not valid.

# Examples
The frame ITRF93 (class ID 3000) turning at 1 degree per second about the z axis of
J2000, from 0 to 100 s: one record of degree 1 with RA = -90° and DEC = 90°, so that the
rotation is about the z axis, W = t degrees, and the angular velocity (0, 0, 1) degrees/s:
```jldoctest
julia> path = joinpath(mktempdir(), "type3.bpc");

julia> h = pckopn(path, "type 3 example");

julia> w = pck03b(h, "spin", 3000, "J2000", 0.0, 100.0, 1)
PCK03Writer("spin", 0 records)

julia> record = [50.0, 50.0, -90.0, 0.0, 90.0, 0.0, 50.0, 50.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];

julia> pck03a(w, record, [0.0])

julia> pck03e(w)

julia> pckcls(h)

julia> ctx = SpiceContext();

julia> furnsh(ctx, path)

julia> pxform(ctx, "J2000", "ITRF93", 90.0) ≈ [0 1 0; -1 0 0; 0 0 1]
true
```

# References

- [NAIF documentation of `pck03b`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/pck03b.html)
"""
function pck03b(
        h::DAFHandle, segid::AbstractString, body::Integer, frame::SegFrame, first::Real,
        last::Real, chbdeg::Integer
    )
    who = "pck03b"
    chbdeg < 0 && throw(
        SpiceError(
            :invalidargument, "$who: the degree of the Chebyshev polynomial was " *
                "negative, $chbdeg; the degree of the polynomial must be greater than or equal to zero"
        )
    )
    0 <= body <= 9 && throw(
        SpiceError(
            :barycenteridcode, "$who: you have attempted to create a segment for a " *
                "barycenter, and the PCK system does not support this"
        )
    )
    refcod = _pck_irfnum(frame)
    refcod == 0 && throw(SpiceError(:invalidrefframe, "$who: the reference frame $frame is not supported"))
    first, last = Float64(first), Float64(last)
    first >= last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is at or after the segment stop time $last"
        )
    )
    _daf_check_write(h, who)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid, :nonprintingchars)
    descr = dafps(Float64[first, last], Int[body, refcod, 3, 0, 0])
    return PCK03Writer(h, descr, String(segid), Int(chbdeg) + 1, Float64[], Float64[], true)
end

"""
    pck03a(writer, ncsets, coeffs, epochs) -> Nothing
    pck03a(writer, coeffs, epochs) -> Nothing

Add `ncsets` records to the type 3 segment `writer`, as
`pck03a`
of SPICELIB. CSPICE has no C wrapper for it. Without `ncsets`, the number of records is
`length(epochs)`.

A record has `2 + 6 * (chbdeg + 1)` numbers: the midpoint (seconds past J2000 TDB) and
the radius (s) of its interval, then the coefficients of RA, DEC and W (degrees) and of
the three components of the angular velocity (degrees/s, in the body-fixed frame). The
matrix from the base frame to the body-fixed frame is `R3(W) * R1(90° - DEC) *
R3(90° + RA)`. `coeffs` holds the records one after the other: a vector, or a matrix
with one record in each column. `epochs[i]` is the start of the interval of record `i`,
in seconds past J2000 TDB.

Throws `SpiceError`:
- `:segmentnotopen` if `pck03e` ended the writer;
- `:invalidargument` unless `ncsets > 0`;
- `:dafnosuchhandle` if the handle is closed;
- `:badarg` if `coeffs` or `epochs` does not have the size that `ncsets` requires;
- `:unorderedrefs` unless `epochs` is strictly increasing and starts after the last epoch
  of the call before.

See [`pck03b`](@ref) for an example.

# References

- [NAIF documentation of `pck03a`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/pck03a.html)
"""
function pck03a(w::PCK03Writer, ncsets::Integer, coeffs, epochs)
    who = "pck03a"
    _seg_check_open(who, w)
    ncsets <= 0 && throw(
        SpiceError(
            :invalidargument, "$who: the number of coefficient sets and epochs to be " *
                "added to the PCK segment was not positive; its value was $ncsets"
        )
    )
    _daf_check_write(w.handle, who)
    _seg_add_cheb_packets!(who, w.packets, w.refs, w.ncoef, Int(ncsets), coeffs, epochs)
    return nothing
end
pck03a(w::PCK03Writer, coeffs, epochs) = pck03a(w, length(epochs), coeffs, epochs)

"""
    pck03e(writer) -> Nothing

End the type 3 segment `writer` and write it to its file, as
`pck03e`
of SPICELIB. CSPICE has no C wrapper for it.

[`pck03b`](@ref) and [`pck03a`](@ref) do not change the file, so a segment that is never
ended is not in the file. The segment is on disk when the call returns, with the data
that SPICELIB writes. The writer cannot be used again.

Throws `SpiceError(:segmentnotopen)` if the writer was ended before, and
`SpiceError(:dafnosuchhandle)` if the handle was closed after [`pck03b`](@ref); the
segment is then lost.

# References

- [NAIF documentation of `pck03e`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/pck03e.html)
"""
pck03e(w::PCK03Writer) = _seg_cheb_end!("pck03e", w)

"""
    pckw20([ctx,] handle, clssid, frame, first, last, segid, intlen, cdata, ascale,
           tscale, initjd, initfr) -> Nothing

Write a binary PCK segment of type 20 (Chebyshev polynomials for the rates of three
Euler angles, on intervals of equal length) to the file of `handle`, as
`pckw20` of SPICELIB.

Each record of `cdata` holds `polydg + 1` coefficients of the rate of `phi`, then the
value of `phi` at the midpoint of the interval, then the same for `delta` and for `w`
(the angles of [`pckw02`](@ref)). The angles are in units of `ascale` radians and the
rates in units of `ascale` radians per `tscale` s. `intlen` is the length of every
interval in days, and `initjd + initfr` is the TDB Julian date of the start of the first
interval, as in [`spkw20`](@ref). Different from [`pckw02`](@ref), whose base frame must
be a built-in inertial frame, `frame` is any frame that `ctx` knows, by name or ID, as
`pckw20_` of CSPICE looks it up with `namfrm`. See [`pckw02`](@ref) for the other
arguments.

The segment is on disk when the call returns. A call that throws writes nothing.

Throws `SpiceError`:
- `:dafnosuchhandle` if the handle is closed, `:dafinvalidaccess` if it is open for
  reading;
- `:invalidcount` if `cdata` has no record;
- `:badarg` if the records differ in length, or if their length is not 3 * (degree + 2);
- `:invaliddegree` if the degree is greater than 50;
- `:intlennotpos` unless `intlen > 0`;
- `:invalidrefframe` if `frame` is a name that `ctx` does not know;
- `:baddescrtimes` unless `first < last`;
- `:nonpositivescale` unless `ascale > 0` and `tscale > 0`;
- `:coveragegap` unless [`first`, `last`] is inside the intervals of the data, with the
  tolerance of [`spkw20`](@ref);
- `:idstringtoolong` or `:nonprintablechars` if `segid` is not valid.

# References

- [NAIF documentation of `pckw20`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/pckw20.html)
"""
function pckw20(
        ctx::SegContext, h::DAFHandle, clssid::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, intlen::Real, cdata, ascale::Real, tscale::Real,
        initjd::Real, initfr::Real
    )
    refcod, data = _seg_cheb20(
        "pckw20", "PCK", ctx, h, frame, first, last, segid, intlen, cdata, ascale, tscale,
        initjd, initfr
    )
    _daf_add_array!(h, Float64[first, last], Int[clssid, refcod, 20, 0, 0], segid, data)
    return nothing
end
pckw20(
    h::DAFHandle, clssid::Integer, frame::SegFrame, first::Real, last::Real,
    segid::AbstractString, intlen::Real, cdata, ascale::Real, tscale::Real, initjd::Real,
    initfr::Real
) = pckw20(
    default_context(), h, clssid, frame, first, last, segid, intlen,
    cdata, ascale, tscale, initjd, initfr
)

"""
    pckcls(handle) -> Nothing

Close the binary PCK file of `handle`, as
`pckcls_c`. The call does nothing if the handle
is closed already. The close writes no data: every write call puts its bytes on disk
before it returns. Load the file with [`furnsh`](@ref) after the close. For a handle in memory
(from [`pckops`](@ref)), the close discards the bytes: load it with [`furnsh`](@ref) or
save it with [`save_kernel`](@ref) before the close.

Throws `SpiceError(:nosegmentsfound)` if the file is open for writing and holds no
segment. The handle then stays open; close it with [`dafcls`](@ref).

# References

- [NAIF documentation of `pckcls_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pckcls_c.html)
"""
pckcls(h::DAFHandle) = _seg_close(h, "PCK")
