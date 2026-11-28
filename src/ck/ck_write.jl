# CK writers: `ckopn`, `ckcls` and the segment writers `ckw01`, `ckw02`, `ckw03`,
# `ckw04b`/`ckw04a`/`ckw04e`, `ckw05`, `ckw06`. Each writer makes the input checks of the
# SPICELIB routine in the same order and writes the descriptor, the pointing data, the
# directories and the trailer in the layout of CSPICE, so the data area is equal bit for
# bit to the one that CSPICE writes from the same input.
#
# Argument lists follow CSPICE without the counts and flags that follow from the input:
# the record count is the length of the input, and `avflag` of `ckw01` and `ckw03` is
# `avvs !== nothing`. Quaternions are SPICE quaternions (scalar first). A table of
# quaternions, angular velocities or packets is a vector of vectors, a matrix with one
# row of the table in each column, or a flat vector. Times are encoded spacecraft clock
# ticks. See `spk_write.jl` for the frame argument, the optional first argument `ctx`
# and the error kinds.

"""
    ckopn(path, ifname, ncomch=0) -> DAFHandle

Create the CK file `path` and open it for writing, as
`ckopn_c`.

`ifname` is the internal file name and `ncomch` is the number of comment characters to
make room for. The file is on disk when the call returns. An existing file is never
replaced. Add segments with the `ckwNN` writers and close the handle with
[`ckcls`](@ref). Load the file with [`furnsh`](@ref) only after the close.

Throws `SpiceError(:emptystring)` if `path` or `ifname` is empty and `:blankfilename` if
`path` is blank, as CSPICE, and `:fileopenfail` if the file exists or cannot be created.

# Examples
Write a type 3 segment with two pointing records and one interpolation interval:
```jldoctest
julia> path = joinpath(mktempdir(), "example.bc");

julia> h = ckopn(path, "example CK");

julia> quats = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]];

julia> ckw03(h, 0.0, 10.0, -77000, "J2000", "attitude", [0.0, 10.0], quats, nothing, [0.0])

julia> ckcls(h)

julia> getfat(path)
("DAF", "CK")
```

See also [`ckw01`](@ref), [`ckw02`](@ref), [`ckw03`](@ref), [`ckw05`](@ref),
[`ckw06`](@ref), [`dafac`](@ref) to add comments.

# References

- [NAIF documentation of `ckopn_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckopn_c.html)
"""
ckopn(path::AbstractString, ifname::AbstractString, ncomch::Integer = 0) =
    _seg_opn("ckopn", "CK", 6, path, ifname, ncomch)

"""
    ckops(ifname="", ncomch=0) -> DAFHandle

Open a new CK file that lives in memory, for writing: [`ckopn`](@ref) without the
path. No file is made. The writers (such as [`ckw03`](@ref)) take the handle as they take
one of a file and write the same bytes. [`furnsh`](@ref) loads its segments as a kernel,
and [`save_kernel`](@ref) writes them to a file; [`ckcls`](@ref) discards them, so load or
save before the close. CSPICE has no such call; see [`dafops`](@ref).

Throws the errors of [`dafops`](@ref).
"""
ckops(ifname::AbstractString = "", ncomch::Integer = 0) = _seg_ops("CK", 6, ifname, ncomch)

"""
    ckcls(handle) -> Nothing

Close the CK file of `handle`, as
`ckcls_c`.
The call does nothing if the handle is closed already. The close writes no data: every
write call puts its bytes on disk before it returns. For a handle in memory
(from [`ckops`](@ref)), the close discards the bytes: load it with [`furnsh`](@ref) or
save it with [`save_kernel`](@ref) before the close.

Throws `SpiceError(:nosegmentsfound)` if the file is open for writing and holds no
segment. The handle then stays open; close it with [`dafcls`](@ref).

# References

- [NAIF documentation of `ckcls_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckcls_c.html)
"""
ckcls(h::DAFHandle) = _seg_close(h, "CK")

function _ck_add!(
        h::DAFHandle, inst::Integer, refcod::Int, dtype::Int, avflag::Bool,
        begtim::Real, endtim::Real, segid::AbstractString, data::Vector{Float64}
    )
    _daf_add_array!(
        h, Float64[begtim, endtim], Int[inst, refcod, dtype, avflag ? 1 : 0, 0, 0],
        segid, data
    )
    return nothing
end

"""
    _ck_avvs(who, avvs, nrec) -> (flat::Vector{Float64}, avflag::Bool)

Angular velocities for `nrec` records: none if `avvs` is `nothing`, else `avvs` must
hold `nrec` vectors. Another count throws `SpiceError(:badarg)`.
"""
function _ck_avvs(who::AbstractString, avvs, nrec::Int)
    avvs === nothing && return Float64[], false
    av, n, _ = _seg_flat(who, "avvs", avvs, 3)
    n == nrec || throw(
        SpiceError(
            :badarg,
            "$who: `avvs` must have $nrec elements (or be `nothing`) but has $n"
        )
    )
    return av, true
end

"""
    _ck_nonzero_quat(who, flat, a, i, mini = 0) -> Nothing

Throw `SpiceError(:zeroquaternion)` if the quaternion `flat[a:a+3]`, number `i` (of
mini-segment `mini` if `mini > 0`), is zero.
"""
@inline function _ck_nonzero_quat(
        who::AbstractString, flat::Vector{Float64}, a::Int, i::Int, mini::Int = 0
    )
    (flat[a] == 0.0 && flat[a + 1] == 0.0 && flat[a + 2] == 0.0 && flat[a + 3] == 0.0) &&
        _ck_zero_quat(who, i, mini)
    return nothing
end

@noinline _ck_zero_quat(who, i, mini) = throw(
    SpiceError(
        :zeroquaternion, "$who: the quaternion at index $i" *
            (mini > 0 ? " within mini-segment $mini" : "") * " has magnitude zero"
    )
)

"The quaternions and the angular velocities (see `_ck_avvs`) of `nrec` records of `ckw01` and `ckw03`; another count of quaternions throws `SpiceError(:badarg)`."
function _ck_quats_avvs(who::AbstractString, quats, avvs, nrec::Int)
    q, nq, _ = _seg_flat(who, "quats", quats, 4)
    nq == nrec || throw(SpiceError(:badarg, "$who: `quats` must have $nrec elements but has $nq"))
    av, avflag = _ck_avvs(who, avvs, nrec)
    return q, av, avflag
end

"""
    _ck_ids_times(who, ctx, ref, segid, t) -> Int

The checks of `ckw01` and `ckw03` after the descriptor times, in their order: the frame
(whose ID it returns), the segment ID, and times `t` that are not negative
(`SpiceError(:invalidsclktime)`) and strictly increasing.
"""
function _ck_ids_times(who::AbstractString, ctx, ref, segid::AbstractString, t::Vector{Float64})
    refcod = _seg_refcode(who, ctx, ref)
    _seg_idchars(who, segid)
    _seg_idlen(who, segid)
    t[1] < 0.0 && throw(SpiceError(:invalidsclktime, "$who: the first SCLKDP time $(t[1]) is negative"))
    _seg_increasing(who, t, "SCLKDP")
    return refcod
end

"Throws `SpiceError(:zeroquaternion)` if one of the `n` quaternions that start every `step` doubles in `flat` is zero."
function _ck_nonzero_quats(who::AbstractString, flat::Vector{Float64}, n::Int, step::Int)
    for i in 1:n
        _ck_nonzero_quat(who, flat, (i - 1) * step + 1, i)
    end
    return nothing
end

"Quaternions and, if `avflag`, angular velocities of `nrec` records interleaved as CK types 1 and 3 store them."
function _ck_pointing(q::Vector{Float64}, av::Vector{Float64}, nrec::Int, avflag::Bool)
    avflag || return copy(q)
    data = Vector{Float64}(undef, 7nrec)
    for i in 1:nrec                      # invariant: records 1 .. i-1 are in data
        copyto!(data, 7(i - 1) + 1, q, 4(i - 1) + 1, 4)
        copyto!(data, 7(i - 1) + 5, av, 3(i - 1) + 1, 3)
    end
    return data
end

"""
    ckw01([ctx,] handle, begtim, endtim, inst, ref, segid, sclkdp, quats, avvs) -> Nothing

Write a CK segment of type 1 (discrete pointing) to the file of `handle`, as
`ckw01_c`.

# Arguments
The arguments `ctx`, `handle`, `begtim`, `endtim`, `inst`, `ref` and `segid` have the same
meaning in every `ckwNN` writer. [`ckw06`](@ref) names the coverage `first`, `last`.
- `ctx`: the context in which a frame name is resolved (default: the default context).
- `handle`: a [`DAFHandle`](@ref) from [`ckopn`](@ref) or [`ckops`](@ref).
- `begtim`, `endtim`: the coverage of the segment, encoded spacecraft clock ticks
  ([`sce2c`](@ref) converts a time to them).
- `inst`: the NAIF ID of the instrument or structure.
- `ref`: the reference frame, a name or a frame ID. A name is resolved in `ctx`;
  built-in frames need no kernel.
- `segid`: the segment identifier, at most 40 printable characters.
- `sclkdp`: the times of the pointing records, encoded spacecraft clock ticks.
- `quats`: one SPICE quaternion (scalar first) for each record. [`q2m`](@ref) of it is
  the matrix that rotates vectors from `ref` to the frame of `inst`, the `cmat` of
  [`ckgp`](@ref).
- `avvs`: for each record, the angular velocity of the frame of `inst` relative to
  `ref`, in `ref` (rad/s); or `nothing`.

The number of records is `length(sclkdp)`. A table such as `quats` is a vector of vectors
(the form of SPICE.jl), a matrix with one record in each column, or a flat vector.

The segment has angular velocity if `avvs` holds one vector for each record, and none if
`avvs` is `nothing`; this takes the place of the `avflag` argument of CSPICE. Another
number of vectors throws `SpiceError(:badarg)`.

Every `ckwNN` writer puts its segment on disk when the call returns, with the data that
CSPICE writes from the same input. A call that throws writes nothing.

Throws `SpiceError`. Every `ckwNN` writer throws:
- `:emptystring` if `ref` or `segid` is empty, except `ckw06`, which has no C function in
  CSPICE N0067;
- `:dafnosuchhandle` if the handle is closed, `:dafinvalidaccess` if it is open for
  reading;
- `:invalidrefframe` if the frame is not known;
- `:nonprintablechars` or `:segidtoolong` if `segid` is not valid.

`ckw01` also throws:
- `:invalidnumrec` if there is no record;
- `:invaliddescrtime` unless `begtim <= sclkdp[1]` and `endtim >= sclkdp[end]`;
- `:invalidsclktime` if `sclkdp[1] < 0`;
- `:timesoutoforder` unless `sclkdp` is strictly increasing;
- `:zeroquaternion` if a quaternion is zero.

See also [`ckw02`](@ref), [`ckw03`](@ref), [`ckw05`](@ref), [`ckw06`](@ref), and
[`ckopn`](@ref) for an example.

# References

- [NAIF documentation of `ckw01_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckw01_c.html)
"""
function ckw01(
        ctx::SegContext, h::DAFHandle, begtim::Real, endtim::Real, inst::Integer,
        ref::SegFrame, segid::AbstractString, sclkdp, quats, avvs
    )
    who = "ckw01"
    ref isa AbstractString && _chkfstr(who, "string `ref`", ref)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    t = _seg_vec(who, "sclkdp", sclkdp, -1)
    nrec = length(t)
    q, av, avflag = _ck_quats_avvs(who, quats, avvs, nrec)
    nrec <= 0 && throw(
        SpiceError(
            :invalidnumrec,
            "$who: $nrec is an invalid number of pointing instances for type 1"
        )
    )
    begtim > t[1] && throw(
        SpiceError(
            :invaliddescrtime, "$who: the first d.p. component " *
                "of the descriptor is invalid: DCD(1) = $begtim and SCLKDP(1) = $(t[1])"
        )
    )
    endtim < t[nrec] && throw(
        SpiceError(
            :invaliddescrtime, "$who: the second d.p. component " *
                "of the descriptor is invalid: DCD(2) = $endtim and SCLKDP(NREC) = $(t[nrec])"
        )
    )
    refcod = _ck_ids_times(who, ctx, ref, segid, t)
    _ck_nonzero_quats(who, q, nrec, 4)
    data = _ck_pointing(q, av, nrec, avflag)
    append!(data, t)
    for i in 1:((nrec - 1) ÷ 100)              # directory: midpoints between records 100i and 100i+1
        push!(data, (t[100i] + t[100i + 1]) / 2.0)
    end
    push!(data, Float64(nrec))
    _ck_add!(h, inst, refcod, 1, avflag, begtim, endtim, segid, data)
    return nothing
end
ckw01(
    h::DAFHandle, begtim::Real, endtim::Real, inst::Integer, ref::SegFrame,
    segid::AbstractString, sclkdp, quats, avvs
) =
    ckw01(default_context(), h, begtim, endtim, inst, ref, segid, sclkdp, quats, avvs)

"""
    ckw02([ctx,] handle, begtim, endtim, inst, ref, segid, start, stop, quats, avvs, rates)
        -> Nothing

Write a CK segment of type 2 (pointing with a constant angular velocity on intervals), as
`ckw02_c`.

Record `i` is valid from `start[i]` to `stop[i]` (encoded spacecraft clock ticks).
`quats[i]` is its quaternion, `avvs[i]` its angular velocity (rad/s) and `rates[i]` its
number of seconds per tick. The number of records is `length(start)`. See
[`ckw01`](@ref) for the other arguments and for the errors of every writer.

Throws `SpiceError`:
- `:invalidnumrec` if there is no record;
- `:invaliddescrtime` unless `begtim <= start[1]` and `endtim >= stop[end]`;
- `:invalidsclktime` if `start[1] < 0`;
- `:degenerateinterval` unless `stop[i] > start[i]`;
- `:timesoutoforder` unless `start` is strictly increasing;
- `:badstoptime` if `stop[i-1] > start[i]`;
- `:zeroquaternion` if a quaternion is zero.

# References

- [NAIF documentation of `ckw02_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckw02_c.html)
"""
function ckw02(
        ctx::SegContext, h::DAFHandle, begtim::Real, endtim::Real, inst::Integer,
        ref::SegFrame, segid::AbstractString, start, stop, quats, avvs, rates
    )
    who = "ckw02"
    ref isa AbstractString && _chkfstr(who, "string `ref`", ref)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    t0 = _seg_vec(who, "start", start, -1)
    nrec = length(t0)
    t1 = _seg_vec(who, "stop", stop, nrec)
    rt = _seg_vec(who, "rates", rates, nrec)
    q, nq, _ = _seg_flat(who, "quats", quats, 4)
    av, na, _ = _seg_flat(who, "avvs", avvs, 3)
    (nq == nrec && na == nrec) || throw(
        SpiceError(
            :badarg,
            "$who: `quats` and `avvs` must have $nrec elements but have $nq and $na"
        )
    )
    nrec <= 0 && throw(
        SpiceError(
            :invalidnumrec,
            "$who: $nrec is an invalid number of pointing instances for type 2"
        )
    )
    begtim > t0[1] && throw(
        SpiceError(
            :invaliddescrtime, "$who: the first d.p. component " *
                "of the descriptor is invalid: DCD(1) = $begtim and START(1) = $(t0[1])"
        )
    )
    endtim < t1[nrec] && throw(
        SpiceError(
            :invaliddescrtime, "$who: the second d.p. component " *
                "of the descriptor is invalid: DCD(2) = $endtim and STOP(NREC) = $(t1[nrec])"
        )
    )
    refcod = _seg_refcode(who, ctx, ref)
    _seg_idchars(who, segid)
    _seg_idlen(who, segid)
    t0[1] < 0.0 && throw(SpiceError(:invalidsclktime, "$who: the first START time $(t0[1]) is negative"))
    for i in 1:nrec                      # invariant: intervals 1 .. i-1 are valid and ordered
        t1[i] <= t0[i] && throw(
            SpiceError(
                :degenerateinterval, "$who: the STOP time is " *
                    "less than or equal to the START time for interval number $i"
            )
        )
        i == 1 && continue
        t0[i] <= t0[i - 1] && throw(
            SpiceError(
                :timesoutoforder, "$who: the START times are " *
                    "not strictly increasing: START($(i - 1)) = $(t0[i - 1]) and START($i) = $(t0[i])"
            )
        )
        t1[i - 1] > t0[i] && throw(
            SpiceError(
                :badstoptime, "$who: the STOP time for " *
                    "interval $(i - 1) is greater than the following START time"
            )
        )
    end
    _ck_nonzero_quats(who, q, nrec, 4)
    data = Vector{Float64}(undef, 8nrec)
    for i in 1:nrec                      # invariant: records 1 .. i-1 are in data
        copyto!(data, 8(i - 1) + 1, q, 4(i - 1) + 1, 4)
        copyto!(data, 8(i - 1) + 5, av, 3(i - 1) + 1, 3)
        data[8i] = rt[i]
    end
    append!(data, t0)
    append!(data, t1)
    for i in 1:((nrec - 1) ÷ 100)
        push!(data, (t1[100i] + t0[100i + 1]) / 2.0)
    end
    _ck_add!(h, inst, refcod, 2, true, begtim, endtim, segid, data)
    return nothing
end
ckw02(
    h::DAFHandle, begtim::Real, endtim::Real, inst::Integer, ref::SegFrame,
    segid::AbstractString, start, stop, quats, avvs, rates
) =
    ckw02(
    default_context(), h, begtim, endtim, inst, ref, segid, start, stop, quats, avvs,
    rates
)

"""
    ckw03([ctx,] handle, begtim, endtim, inst, ref, segid, sclkdp, quats, avvs, starts)
        -> Nothing

Write a CK segment of type 3 (pointing with linear interpolation inside intervals), as
`ckw03_c`.

`sclkdp` holds the times of the `nrec` records (encoded spacecraft clock ticks), `quats`
their quaternions and `starts` the start times of the interpolation intervals. Inside
an interval the pointing is interpolated between adjacent records; between the last
record of an interval and the start of the next there is none. `avvs` is as for
[`ckw01`](@ref), which also gives the other arguments and the errors of every writer.

Throws `SpiceError`:
- `:invalidnumrec` if there is no record, `:invalidnumint` if there is no interval;
- `:invaliddescrtime` unless `begtim <= sclkdp[1]` and `endtim >= sclkdp[end]`;
- `:invalidsclktime` if `sclkdp[1] < 0`;
- `:timesoutoforder` unless `sclkdp` and `starts` are strictly increasing;
- `:timesdontmatch` unless `starts[1] == sclkdp[1]`;
- `:invalidstarttime` unless every start time is equal to one of the `sclkdp`;
- `:zeroquaternion` if a quaternion is zero.

See [`ckopn`](@ref) for an example.

# References

- [NAIF documentation of `ckw03_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckw03_c.html)
"""
function ckw03(
        ctx::SegContext, h::DAFHandle, begtim::Real, endtim::Real, inst::Integer,
        ref::SegFrame, segid::AbstractString, sclkdp, quats, avvs, starts
    )
    who = "ckw03"
    ref isa AbstractString && _chkfstr(who, "string `ref`", ref)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    t = _seg_vec(who, "sclkdp", sclkdp, -1)
    nrec = length(t)
    st = _seg_vec(who, "starts", starts, -1)
    nints = length(st)
    q, av, avflag = _ck_quats_avvs(who, quats, avvs, nrec)
    nrec <= 0 && throw(
        SpiceError(
            :invalidnumrec,
            "$who: $nrec is an invalid number of pointing instances for type 3"
        )
    )
    nints <= 0 && throw(
        SpiceError(
            :invalidnumint,
            "$who: $nints is an invalid number of interpolation intervals for type 3"
        )
    )
    begtim > t[1] && throw(
        SpiceError(
            :invaliddescrtime, "$who: the segment begin time is " *
                "greater than the time of the first pointing instance: DCD(1) = $begtim and " *
                "SCLKDP(1) = $(t[1])"
        )
    )
    endtim < t[nrec] && throw(
        SpiceError(
            :invaliddescrtime, "$who: the segment end time is " *
                "less than the time of the last pointing instance: DCD(2) = $endtim and " *
                "SCLKDP($nrec) = $(t[nrec])"
        )
    )
    refcod = _ck_ids_times(who, ctx, ref, segid, t)
    t[1] != st[1] && throw(
        SpiceError(
            :timesdontmatch, "$who: the start time of the first " *
                "interval $(st[1]) and the time of the first pointing instance $(t[1]) are not the same"
        )
    )
    _seg_increasing(who, st, "STARTS")
    index = 0
    for i in 1:nints                     # invariant: starts 1 .. i-1 match records up to `index`
        match = false
        while !match && index < nrec
            index += 1
            match = st[i] == t[index]
        end
        match || throw(
            SpiceError(
                :invalidstarttime,
                "$who: interval start time number $i is invalid: STARTS($i) = $(st[i])"
            )
        )
    end
    _ck_nonzero_quats(who, q, nrec, 4)
    data = _ck_pointing(q, av, nrec, avflag)
    append!(data, t)
    _seg_directory!(data, t, (nrec - 1) ÷ 100)
    append!(data, st)
    _seg_directory!(data, st, (nints - 1) ÷ 100)
    push!(data, Float64(nints), Float64(nrec))
    _ck_add!(h, inst, refcod, 3, avflag, begtim, endtim, segid, data)
    return nothing
end
ckw03(
    h::DAFHandle, begtim::Real, endtim::Real, inst::Integer, ref::SegFrame,
    segid::AbstractString, sclkdp, quats, avvs, starts
) =
    ckw03(
    default_context(), h, begtim, endtim, inst, ref, segid, sclkdp, quats, avvs,
    starts
)

"""
    CKW04Writer

A CK segment of type 4 that is open for the addition of records. [`ckw04b`](@ref)
returns it, and it is the first argument of [`ckw04a`](@ref) and [`ckw04e`](@ref).

SPICELIB keeps this state in the table of the generic segment writer `sgseqw`; here the
caller owns it.
"""
mutable struct CKW04Writer
    const handle::DAFHandle
    const begtim::Float64
    const inst::Int
    const refcod::Int
    const avflag::Bool
    const segid::String
    const packets::Vector{Float64}      # (reference, size, packet) triples, as `sgwvpk` writes them
    const sizes::Vector{Int}            # the stored packet sizes
    const refs::Vector{Float64}
    open::Bool
end

Base.show(io::IO, w::CKW04Writer) = _seg_show_writer(io, w)

"""
    ckw04b([ctx,] handle, begtim, inst, ref, avflag, segid) -> CKW04Writer

Begin a CK segment of type 4 (Chebyshev polynomials for the quaternion and the angular
velocity, on intervals of unequal length) in the file of `handle`, as
`ckw04b`
of SPICELIB, and return its writer. CSPICE has no C wrapper for it.

`begtim` is the begin time of the segment (encoded spacecraft clock ticks), `inst` the
NAIF ID of the instrument, `ref` the reference frame (a name resolved in `ctx`, or an
ID) and `avflag` whether the segment holds angular velocity. Add records with
[`ckw04a`](@ref) and end the segment with [`ckw04e`](@ref), which gives the end time.

Throws `SpiceError`:
- `:invalidrefframe` if `ref` is not a known frame;
- `:nonprintablechars` or `:segidtoolong` if `segid` is not valid;
- `:dafnosuchhandle` or `:dafinvalidaccess` if the file is not open for writing.

# Examples
One record of degree 0 (constant pointing) on the interval 0 to 200 ticks:
```jldoctest
julia> path = joinpath(mktempdir(), "type4.bc");

julia> h = ckopn(path, "type 4 example");

julia> w = ckw04b(h, 0.0, -77000, "J2000", false, "attitude")
CKW04Writer("attitude", 0 records)

julia> record = [100.0, 100.0, 1, 1, 1, 1, 0, 0, 0, 1.0, 0.0, 0.0, 0.0];

julia> ckw04a(w, [length(record)], record, [0.0])

julia> ckw04e(w, 200.0)

julia> ckcls(h)

julia> ctx = SpiceContext();

julia> furnsh(ctx, path)

julia> cmat, clkout = ckgp(ctx, -77000, 50.0, 0.0, "J2000");

julia> cmat == [1 0 0; 0 1 0; 0 0 1], clkout
(true, 50.0)
```

# References

- [NAIF documentation of `ckw04b`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/ckw04b.html)
"""
function ckw04b(
        ctx::SegContext, h::DAFHandle, begtim::Real, inst::Integer, ref::SegFrame,
        avflag::Bool, segid::AbstractString
    )
    who = "ckw04b"
    refcod = _seg_refcode(who, ctx, ref)
    _seg_idchars(who, segid)
    _seg_idlen(who, segid)
    _daf_check_write(h, who)
    return CKW04Writer(
        h, Float64(begtim), Int(inst), refcod, avflag, String(segid), Float64[], Int[],
        Float64[], true
    )
end
ckw04b(
    h::DAFHandle, begtim::Real, inst::Integer, ref::SegFrame, avflag::Bool,
    segid::AbstractString
) = ckw04b(default_context(), h, begtim, inst, ref, avflag, segid)

"""
    ckw04a(writer, npkts, pktsiz, pktdat, sclkdp) -> Nothing
    ckw04a(writer, pktsiz, pktdat, sclkdp) -> Nothing

Add `npkts` records to the type 4 segment `writer`, as
`ckw04a`
of SPICELIB. CSPICE has no C wrapper for it. Without `npkts`, the number of records is
`length(sclkdp)`.

Record `i` has `pktsiz[i]` numbers: the midpoint and the radius of its interval (ticks),
the seven coefficient counts of the quaternion components q0 to q3 and of the three
angular velocity components, then the Chebyshev coefficients, those of q0 first (the
angular velocity in rad/s).
`pktdat` holds the records one after the other. `sclkdp[i]` is the start of the interval
of record `i` (ticks). The inputs are not changed (SPICELIB changes `pktsiz` and
`pktdat`). [`ckgr04`](@ref) returns a record in this layout.

Throws `SpiceError`:
- `:segmentnotopen` if `ckw04e` ended the writer;
- `:badarg` if `pktsiz`, `pktdat` or `sclkdp` does not have the size that `npkts` and
  `pktsiz` require;
- `:invalidargument` unless `0 < pktsiz[i] <= 142`;
- `:dafnosuchhandle` if the handle is closed;
- `:numpacketsnotpos` unless `npkts > 0`;
- `:unorderedrefs` unless `sclkdp` is strictly increasing and starts after the last
  start of the call before;
- `:nonpospacketsize` if `pktsiz[i] <= 6`: SPICELIB packs the seven counts into one
  number, so the stored record has `pktsiz[i] - 6` numbers.

See [`ckw04b`](@ref) for an example.

# References

- [NAIF documentation of `ckw04a`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/ckw04a.html)
"""
function ckw04a(w::CKW04Writer, npkts::Integer, pktsiz, pktdat, sclkdp)
    who = "ckw04a"
    _seg_check_open(who, w)
    n = max(Int(npkts), 0)
    sz = _seg_vec(who, "pktsiz", pktsiz, n)
    for (k, s) in enumerate(sz)
        isinteger(s) || throw(SpiceError(:badarg, "$who: `pktsiz[$k]` = $s is not an integer"))
    end
    dat = _seg_vec(who, "pktdat", pktdat, Int(sum(sz; init = 0.0)))
    ref = _seg_vec(who, "sclkdp", sclkdp, n)
    for k in 1:n
        sz[k] <= 0 && throw(
            SpiceError(
                :invalidargument, "$who: the number of coefficient sets and epochs " *
                    "in the $k data packet (record) was not positive; its value was $(Int(sz[k]))"
            )
        )
        sz[k] >= 143 && throw(
            SpiceError(
                :invalidargument, "$who: the total size of the $k data packet (record) " *
                    "is greater than the maximum allowed type 4 record size 142; its value was $(Int(sz[k]))"
            )
        )
    end
    _daf_check_write(w.handle, who)
    npkts <= 0 && throw(
        SpiceError(
            :numpacketsnotpos, "$who: the number of packets to store is not positive; " *
                "the value supplied was $npkts"
        )
    )
    _seg_ordered_refs(who, w.refs, ref)
    for k in 1:n
        sz[k] - 6 <= 0 && throw(
            SpiceError(
                :nonpospacketsize, "$who: the packet size for packet $k was not " *
                    "positive; it had a value of $(Int(sz[k]) - 6)"
            )
        )
    end
    pos = 0
    for k in 1:n                    # invariant: packets 1:k-1 are in `w.packets`
        size = Int(sz[k])
        # The seven counts become one double, base 128 (`zzck4i2d`); the stored packet
        # is the first `size - 6` words of (mid, rad, counts, coefficients).
        push!(w.packets, ref[k], size - 6, dat[pos + 1])
        size >= 8 && push!(w.packets, dat[pos + 2])
        if size >= 9
            packed = 0.0
            x = 1.0
            for kk in 1:7
                packed += _cint(dat[pos + 2 + kk]) * x
                x *= 128.0
            end
            push!(w.packets, packed)
            append!(w.packets, @view dat[(pos + 10):(pos + size)])
        end
        push!(w.sizes, size - 6)
        pos += size
    end
    append!(w.refs, ref)
    return nothing
end
ckw04a(w::CKW04Writer, pktsiz, pktdat, sclkdp) = ckw04a(w, length(sclkdp), pktsiz, pktdat, sclkdp)

"""
    ckw04e(writer, endtim) -> Nothing

End the type 4 segment `writer` with the end time `endtim` (ticks) and write it to its
file, as
`ckw04e`
of SPICELIB. CSPICE has no C wrapper for it.

[`ckw04b`](@ref) and [`ckw04a`](@ref) do not change the file, so a segment that is never
ended is not in the file. The segment is on disk when the call returns, with the data that
SPICELIB writes. The writer cannot be used again.

Throws `SpiceError(:segmentnotopen)` if the writer was ended before, and
`SpiceError(:dafnosuchhandle)` if the handle was closed after [`ckw04b`](@ref); the
segment is then lost.

# References

- [NAIF documentation of `ckw04e`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/ckw04e.html)
"""
function ckw04e(w::CKW04Writer, endtim::Real)
    who = "ckw04e"
    _seg_check_open(who, w)
    _daf_check_write(w.handle, who)
    npkt = length(w.refs)
    data = copy(w.packets)
    # Packet directory of `sgwes`: the offset of each packet after its (reference,
    # size) pair, and one past the last.
    addr = 3.0
    for s in w.sizes
        push!(data, addr)
        addr += s + 2
    end
    push!(data, addr)
    append!(data, w.refs)
    _seg_directory!(data, w.refs, (npkt - 1) ÷ 100)
    total = sum(w.sizes; init = 0)
    pdrbas = total + 2npkt
    refbas = pdrbas + npkt + 1
    # Meta data of `sgwes` for packets of variable size with explicit reference values:
    # CONBAS, NCON, RDRBAS, NRDR, RDRTYP, REFBAS, NREF, PDRBAS, NPDR, PDRTYP, PKTBAS,
    # NPKT, RSVBAS, NRSV, PKTSZ, PKTOFF, NMETA.
    append!(
        data, Float64[
            0, 0, refbas + npkt, (npkt - 1) ÷ 100, 3, refbas, npkt, pdrbas, npkt + 1, 1,
            0, npkt, 0, 0, maximum(w.sizes; init = 0), 2, 17,
        ]
    )
    _ck_add!(w.handle, w.inst, w.refcod, 4, w.avflag, w.begtim, endtim, w.segid, data)
    w.open = false
    return nothing
end

"""
    ckw05([ctx,] handle, subtyp, degree, begtim, endtim, inst, ref, avflag, segid, sclkdp,
          packts, rate, nints, starts) -> Nothing

Write a CK segment of type 5 (Hermite or Lagrange interpolation of attitude packets), as
`ckw05_c`.

`subtyp` selects the interpolation and the content of a packet. `q` is a SPICE
quaternion (scalar first), `dq` its derivative, `av` the angular velocity in `ref`
(rad/s) and `dav` its derivative (rad/s²):

| `subtyp` | Interpolation | Packet                 | Values |
|:---------|:--------------|:-----------------------|:-------|
| 0        | Hermite       | `q`, `dq`              | 8      |
| 1        | Lagrange      | `q`                    | 4      |
| 2        | Hermite       | `q`, `dq`, `av`, `dav` | 14     |
| 3        | Lagrange      | `q`, `av`              | 7      |

The other arguments:
- `degree`: the degree of the interpolating polynomial. The window size is the number of
  packets that one interpolation uses: `degree + 1` for Lagrange and `(degree + 1) ÷ 2`
  for Hermite.
- `avflag`: `true` if the segment has angular velocity.
- `sclkdp`: the times of the packets `packts`, encoded spacecraft clock ticks.
- `rate`: the nominal clock rate, seconds per tick.
- `nints`, `starts`: the first `nints` elements of `starts` are the start times of the
  interpolation intervals. No interpolation uses packets of two intervals.

See [`ckw01`](@ref) for the remaining arguments and for the errors of every writer.

Throws `SpiceError`:
- `:toofewpackets` if there is no packet, `:invalidnumints` unless `nints > 0`;
- `:invalidsclktime` if `sclkdp[1] < 0`;
- `:timesoutoforder` unless `sclkdp` and `starts` are strictly increasing;
- `:invalidstarttime` unless every start time is equal to one of the `sclkdp`;
- `:notsupported` if `subtyp` is not in `0:3`;
- `:zeroquaternion` if a quaternion is zero;
- `:invaliddegree` unless `1 <= degree <= 23` and the window size is even: the degree is
  3 modulo 4 for a Hermite subtype, and odd for a Lagrange subtype;
- `:baddescrtimes` if `endtim < begtim`;
- `:emptysegment` if no time of `sclkdp` is inside [`begtim`, `endtim`];
- `:invalidvalue` if `rate` is zero.

# References

- [NAIF documentation of `ckw05_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ckw05_c.html)
"""
function ckw05(
        ctx::SegContext, h::DAFHandle, subtyp::Integer, degree::Integer,
        begtim::Real, endtim::Real, inst::Integer, ref::SegFrame, avflag::Bool,
        segid::AbstractString, sclkdp, packts, rate::Real, nints::Integer, starts
    )
    who = "ckw05"
    ref isa AbstractString && _chkfstr(who, "string `ref`", ref)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    t = _seg_vec(who, "sclkdp", sclkdp, -1)
    n = length(t)
    allst = _seg_vec(who, "starts", starts, -1)
    nints <= length(allst) || throw(
        SpiceError(
            :badarg,
            "$who: `starts` has $(length(allst)) elements, fewer than nints = $nints"
        )
    )
    n < 1 && throw(
        SpiceError(
            :toofewpackets,
            "$who: at least 1 packet is required for CK type 5; number of packets supplied: $n"
        )
    )
    nints <= 0 && throw(
        SpiceError(
            :invalidnumints,
            "$who: $nints is an invalid number of interpolation intervals for type 5"
        )
    )
    st = allst[1:nints]
    refcod = _seg_refcode(who, ctx, ref)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    t[1] < 0.0 && throw(SpiceError(:invalidsclktime, "$who: the first SCLKDP time $(t[1]) is negative"))
    _seg_increasing(who, t, "SCLKDP")
    _seg_increasing(who, st, "STARTS")
    for i in 1:nints
        k = searchsortedfirst(t, st[i])
        (k <= n && t[k] == st[i]) || throw(
            SpiceError(
                :invalidstarttime,
                "$who: interval start time number $i is invalid: STARTS($i) = $(st[i])"
            )
        )
    end
    (0 <= subtyp <= 3) || throw(SpiceError(:notsupported, "$who: CK type 5 subtype <$subtyp> is not supported"))
    packsz = _CK5_PACKSZ[subtyp + 1]
    winsiz = isodd(subtyp) ? degree + 1 : (degree + 1) ÷ 2
    pk, np, _ = _seg_flat(who, "packts", packts, packsz)
    np == n || throw(SpiceError(:badarg, "$who: `packts` must have $n elements but has $np"))
    _ck_nonzero_quats(who, pk, n, packsz)
    (degree < 1 || degree > 23) && throw(
        SpiceError(
            :invaliddegree, "$who: the " *
                "interpolating polynomials have degree $degree; the valid degree range is [1, 23]"
        )
    )
    isodd(winsiz) && throw(
        SpiceError(
            :invaliddegree, "$who: the interpolating polynomials " *
                "have degree $degree; for CK type 5, the degree must be equivalent to 3 mod 4 " *
                "for Hermite interpolation and odd for Lagrange interpolation"
        )
    )
    endtim < begtim && throw(
        SpiceError(
            :baddescrtimes,
            "$who: descriptor bounds are non-increasing: $begtim:$endtim"
        )
    )
    i = searchsortedfirst(t, begtim) - 1     # number of time tags less than begtim
    i == n && throw(
        SpiceError(
            :emptysegment,
            "$who: all time tags are less than segment start time $begtim"
        )
    )
    t[i + 1] > endtim && throw(
        SpiceError(
            :emptysegment, "$who: no time tags lie between " *
                "the segment start time $begtim and segment end time $endtim"
        )
    )
    rate == 0 && throw(SpiceError(:invalidvalue, "$who: the SCLK rate RATE was zero"))
    data = append!(pk, t)
    _seg_directory!(data, t, (n - 1) ÷ 100)
    append!(data, st)
    _seg_directory!(data, st, (nints - 1) ÷ 100)
    push!(data, Float64(rate), Float64(subtyp), Float64(winsiz), Float64(nints), Float64(n))
    _ck_add!(h, inst, refcod, 5, avflag, begtim, endtim, segid, data)
    return nothing
end
ckw05(
    h::DAFHandle, subtyp::Integer, degree::Integer, begtim::Real, endtim::Real,
    inst::Integer, ref::SegFrame, avflag::Bool, segid::AbstractString, sclkdp, packts,
    rate::Real, nints::Integer, starts
) =
    ckw05(
    default_context(), h, subtyp, degree, begtim, endtim, inst, ref, avflag, segid,
    sclkdp, packts, rate, nints, starts
)

# The mini-segments of CK type 6 for `_seg_minis` (see `_SPK19_MINI`). As in SPICELIB,
# the last epoch of a mini-segment is compared with the start of its interval: a CK
# type 6 mini-segment may end before its interval, which leaves a coverage gap.
const _CK06_MINI = (
    pktszs = _CK5_PACKSZ, maxsub = 3, maxdeg = 23, endbound = false,
    kname = "CK type 6", count = "mini-segment", bounds = "mini-segment interval",
    ivl = "mini-segment interval", union = "mini-segment intervals", times = "sclkdp",
)

"""
    _ck06_quats(who, pk, i, pktbeg, n, pktsiz, subtyp) -> Nothing

The quaternion checks of `ckw06` on the `n` packets of mini-segment `i`, which start at
`pk[pktbeg]`: `:zeroquaternion` for a zero quaternion, and for a Hermite subtype
`:badquatsign` if a quaternion is farther from its predecessor than its negative is.
"""
function _ck06_quats(
        who::AbstractString, pk::Vector{Float64}, i::Int, pktbeg::Int, n::Int,
        pktsiz::Int, subtyp::Int
    )
    for j in 1:n
        a = pktbeg + pktsiz * (j - 1)
        _ck_nonzero_quat(who, pk, a, j, i)
        if j > 1 && iseven(subtyp)
            p = a - pktsiz
            prevq = SVector(pk[p], pk[p + 1], pk[p + 2], pk[p + 3])
            q = SVector(pk[a], pk[a + 1], pk[a + 2], pk[a + 3])
            _m_vdistg(prevq, -q) < _m_vdistg(prevq, q) && throw(
                SpiceError(
                    :badquatsign, "$who: the " *
                        "quaternion at index $j within mini-segment $i is farther than its " *
                        "negative from its predecessor at index $(j - 1)"
                )
            )
        end
    end
    return nothing
end

"""
    ckw06([ctx,] handle, inst, ref, avflag, first, last, segid, npkts, subtps, degres,
          packts, rates, sclkdp, ivlbds, sellst) -> Nothing

Write a CK segment of type 6 (mini-segments of Hermite or Lagrange interpolation), as
`ckw06`
of SPICELIB. CSPICE has no C wrapper for it.

A mini-segment is a set of packets with its own subtype, degree and clock rate, valid on
its own interval. The number of mini-segments is `length(npkts)`. The subtypes, their
packets and the window size are those of [`ckw05`](@ref).

- `first`, `last`: the coverage of the segment, encoded spacecraft clock ticks. They
  take the place of `begtim` and `endtim`.
- `avflag`: `true` if the segment has angular velocity.
- `npkts[i]`, `subtps[i]`, `degres[i]`, `rates[i]`: the number of packets, the subtype,
  the degree and the clock rate (seconds per tick) of mini-segment `i`.
- `packts`: the flat vector of all packets, in order.
- `sclkdp`: the flat vector of the times of all packets, encoded spacecraft clock ticks.
- `ivlbds`: the `length(npkts) + 1` bounds of the intervals of the mini-segments, ticks.
- `sellst`: if `true`, a time on a bound of two intervals uses the later interval; else
  it uses the earlier one.

See [`ckw01`](@ref) for the other arguments and for the errors of every writer.

Throws `SpiceError`:
- `:baddescrtimes` if `first > last`;
- `:invalidcount` if there is no mini-segment;
- `:boundsoutoforder` unless `ivlbds` is strictly increasing;
- `:coveragegap` unless [`first`, `last`] is inside [`ivlbds[1]`, `ivlbds[end]`].

For a mini-segment it throws:
- `:toofewpackets` if it has fewer than two packets;
- `:invalidsubtype` if its subtype is not in `0:3`;
- `:invalidsclkrate` unless its rate is positive;
- `:invaliddegree` if its degree is not in `1:23`;
- `:badwindowsize` if its window size is odd;
- `:timesoutoforder` unless its times are strictly increasing;
- `:boundsdisagree` if its first time is after the start of its interval, or its last
  time is before that start. The times can end before the end of the interval, which
  leaves a gap in the coverage;
- `:zeroquaternion` if a quaternion is zero;
- `:badquatsign`, for a Hermite subtype, if a quaternion is farther from its
  predecessor than its negative is.

# References

- [NAIF documentation of `ckw06`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/ckw06.html)
"""
function ckw06(
        ctx::SegContext, h::DAFHandle, inst::Integer, ref::SegFrame, avflag::Bool,
        first::Real, last::Real, segid::AbstractString, npkts::AbstractVector{<:Integer},
        subtps::AbstractVector{<:Integer}, degres::AbstractVector{<:Integer}, packts, rates,
        sclkdp, ivlbds, sellst::Bool
    )
    who = "ckw06"
    _daf_check_write(h, who)
    nmini = _seg_mini_count(who, npkts, subtps, degres)
    pk = _seg_vec(who, "packts", packts, -1)
    t = _seg_vec(who, "sclkdp", sclkdp, -1)
    rt = _seg_vec(who, "rates", rates, nmini)
    bd = _seg_vec(who, "ivlbds", ivlbds, nmini < 1 ? -1 : nmini + 1)
    pktcheck = (i, pktbeg, n, pktsiz, subtyp) -> _ck06_quats(who, pk, i, pktbeg, n, pktsiz, subtyp)
    refcod, data = _seg_minis(
        who, _CK06_MINI, ctx, ref, first, last, segid, npkts, subtps, degres, pk, t, rt,
        bd, sellst, pktcheck
    )
    _ck_add!(h, inst, refcod, 6, avflag, first, last, segid, data)
    return nothing
end
ckw06(
    h::DAFHandle, inst::Integer, ref::SegFrame, avflag::Bool, first::Real, last::Real,
    segid::AbstractString, npkts::AbstractVector{<:Integer},
    subtps::AbstractVector{<:Integer}, degres::AbstractVector{<:Integer}, packts, rates,
    sclkdp, ivlbds, sellst::Bool
) =
    ckw06(
    default_context(), h, inst, ref, avflag, first, last, segid, npkts, subtps,
    degres, packts, rates, sclkdp, ivlbds, sellst
)
