# SPK writers: `spkopn`, `spkopa`, `spkcls` and the segment writers `spkwNN`. Each
# writer makes the input checks of the SPICELIB routine in the same order, packs the
# descriptor as `spkpds` does and writes the data, the directory and the trailer of the
# segment in the layout of CSPICE, so the data area is equal bit for bit to the one
# that CSPICE writes from the same input.
#
# Argument lists follow SPICE.jl. A table of states or packets is a vector of vectors
# (the form of SPICE.jl), a matrix with one state in each column, or a flat vector.
# The reference frame is a frame name, which is resolved against the kernel set of
# the optional first argument `ctx` (a `SpiceContext` or a `KernelSet`; default: the
# default context; built-in frames need no kernel), or a frame ID. Each writer has one
# method that takes `ctx` and one without it that passes `default_context()`. Errors
# have the short message of CSPICE in lower case as kind, for example
# `SpiceError(:invaliddegree, …)` for SPICE(INVALIDDEGREE).

"Reference frame of a segment writer: a frame name or a frame ID."
const SegFrame = Union{AbstractString, Integer}
"Where a segment writer resolves a frame name."
const SegContext = Union{SpiceContext, KernelSet}

_seg_ks(ctx::SpiceContext) = kernelset(ctx)
_seg_ks(ks::KernelSet) = ks

"""
    _seg_refcode(who, ctx, frame) -> Int

Frame ID of `frame` (a name resolved in `ctx`, or an ID). Throws
`SpiceError(:invalidrefframe)` if the name is unknown or the ID is 0.
"""
function _seg_refcode(who::AbstractString, ctx::SegContext, frame::SegFrame)
    code = frame isa Integer ? Int(frame) : _namfrm(_seg_ks(ctx), frame)
    code == 0 && throw(
        SpiceError(
            :invalidrefframe,
            "$who: the reference frame $frame is not supported"
        )
    )
    return code
end

function _seg_idlen(who::AbstractString, segid::AbstractString)
    _daf_lastnb(segid) > 40 && throw(
        SpiceError(
            :segidtoolong,
            "$who: segment identifier contains more than 40 characters"
        )
    )
    return nothing
end

"""
    _seg_idchars(who, segid, kind=:nonprintablechars) -> Nothing

Throws `SpiceError(kind)` if a character of `segid` before its trailing blanks is not
printable ASCII (codes 32 to 126). The segment writers of SPICELIB throw
`:nonprintablechars`; the generic segment writer (`sgbwfs`, used by `spkw10`, `spk14b`
and `pck03b`) throws `:nonprintingchars`.
"""
function _seg_idchars(who::AbstractString, segid::AbstractString, kind::Symbol = :nonprintablechars)
    cu = codeunits(String(segid))
    for i in 1:_daf_lastnb(segid)
        (cu[i] < 0x20 || cu[i] > 0x7e) && throw(
            SpiceError(
                kind,
                "$who: the segment identifier contains a nonprinting character at location $i"
            )
        )
    end
    return nothing
end

"The checks of `chckid` for a segment identifier of the kernel type `what`: `:idstringtoolong`, then `:nonprintablechars`."
function _seg_chckid(who::AbstractString, segid::AbstractString, what::AbstractString = "SPK")
    _daf_lastnb(segid) > 40 && throw(
        SpiceError(
            :idstringtoolong,
            "$who: the $what segment identifier is invalid; its last non-blank character " *
                "is at position $(_daf_lastnb(segid)); the maximum allowed length is 40"
        )
    )
    _seg_idchars(who, segid)
    return nothing
end

"""
    _seg_flat(who, what, x, rowlen) -> (flat::Vector{Float64}, n::Int, rowlen::Int)

The table `x` of `n` rows of `rowlen` doubles as one vector, row after row (the memory
layout that the C routine takes). `x` is a matrix with one row of the table in each
column, a vector of vectors, or a flat vector. With `rowlen == 0` the row length is
taken from `x` (not possible for a flat vector). Throws `SpiceError(:badarg)` if the
shape does not fit.
"""
function _seg_flat(who::AbstractString, what::AbstractString, x, rowlen::Int)
    if x isa AbstractMatrix{<:Real}
        m = size(x, 1)
        (rowlen == 0 || m == rowlen) || throw(
            SpiceError(
                :badarg,
                "$who: `$what` must have $rowlen rows (one table row in each column), it has $m"
            )
        )
        return _vec_of(Float64, x), size(x, 2), m
    elseif x isa AbstractVector{<:Real}
        rowlen > 0 || throw(
            SpiceError(
                :badarg,
                "$who: `$what` must be a matrix or a vector of vectors"
            )
        )
        length(x) % rowlen == 0 || throw(
            SpiceError(
                :badarg,
                "$who: the length $(length(x)) of `$what` is not a multiple of $rowlen"
            )
        )
        return _vec_of(Float64, x), length(x) ÷ rowlen, rowlen
    elseif x isa AbstractVector
        n = length(x)
        m = rowlen == 0 ? (n == 0 ? 0 : length(x[1])) : rowlen
        flat = Vector{Float64}(undef, m * n)
        for (k, row) in enumerate(x)     # invariant: rows 1 .. k-1 are copied
            length(row) == m || throw(
                SpiceError(
                    :badarg,
                    "$who: element $k of `$what` has $(length(row)) values, expected $m"
                )
            )
            for j in 1:m
                flat[(k - 1) * m + j] = Float64(row[j])
            end
        end
        return flat, n, m
    end
    throw(SpiceError(:badarg, "$who: `$what` must be a matrix or a vector of vectors"))
end

"The collection `v` as `Vector{Float64}`; throws `SpiceError(:badarg)` unless it has `n` elements."
function _seg_fixed(who::AbstractString, what::AbstractString, v, n::Int)
    length(v) == n || throw(SpiceError(:badarg, "$who: `$what` must have $n elements but has $(length(v))"))
    return _vec_of(Float64, v)
end

"`_seg_fixed` for a vector of numbers (`n < 0`: any length); throws `SpiceError(:badarg)` for anything else."
function _seg_vec(who::AbstractString, what::AbstractString, x, n::Int)
    x isa AbstractVector{<:Real} || throw(SpiceError(:badarg, "$who: `$what` must be a vector of numbers"))
    return n < 0 ? _vec_of(Float64, x) : _seg_fixed(who, what, x, n)
end

"""
    _seg_increasing(who, t, what="") -> Nothing

Throws `SpiceError(:timesoutoforder)` unless `t` is strictly increasing. `what` names the
times in the message of a CK writer (`"SCLKDP"`, `"STARTS"`); empty, the message speaks
of epochs, as the SPK writers do.
"""
function _seg_increasing(who::AbstractString, t::Vector{Float64}, what::AbstractString = "")
    for i in 2:length(t)                 # invariant: t[1:i-1] is strictly increasing
        t[i] <= t[i - 1] && throw(SpiceError(:timesoutoforder, _seg_order_msg(who, what, t, i)))
    end
    return nothing
end

function _seg_order_msg(who::AbstractString, what::AbstractString, t::Vector{Float64}, i::Int)
    isempty(what) &&
        return "$who: epoch $(t[i]) having index $i is not greater than its predecessor $(t[i - 1])"
    return "$who: the $what times are not strictly increasing: " *
        "$what($(i - 1)) = $(t[i - 1]) and $what($i) = $(t[i])"
end

"""
    _seg_ordered_refs(who, refs, ep) -> Nothing

The checks of `sgwfpk` and `sgwvpk` on the reference values `ep` of the packets added
to a generic segment whose packets before have the references `refs`: throws
`SpiceError(:unorderedrefs)` unless `ep[1] > refs[end]` and `ep` is strictly increasing.
"""
function _seg_ordered_refs(who::AbstractString, refs::Vector{Float64}, ep::Vector{Float64})
    (!isempty(refs) && !isempty(ep) && refs[end] >= ep[1]) && throw(
        SpiceError(
            :unorderedrefs, "$who: reference values are out of order; the value " *
                "$(ep[1]) is not greater than the last reference value $(refs[end]) of the " *
                "packets added before"
        )
    )
    for i in 2:length(ep)                # invariant: ep[1:i-1] is strictly increasing
        ep[i - 1] >= ep[i] && throw(
            SpiceError(
                :unorderedrefs, "$who: reference values are out of order; the value " *
                    "$(ep[i - 1]) was found to be out of order for index $(i - 1)"
            )
        )
    end
    return nothing
end

"""
    spkpds([ctx,] body, center, frame, type, first, last) -> Vector{Float64}

Return the packed descriptor of an SPK segment (five doubles), as
`spkpds_c`.

`body` and `center` are NAIF IDs, `frame` is a frame name or a frame ID, `type` is the SPK
data type, and `first`, `last` are the coverage in seconds past J2000 TDB. A frame name is
resolved in `ctx` (default: the default context); built-in frames need no kernel. The two
address components of the descriptor are zero; [`dafena`](@ref) sets them.

Throws `SpiceError`:
- `:emptystring` if `frame` is an empty string;
- `:barycenterephem` if `body == 0`;
- `:bodyandcentersame` if `body == center`;
- `:invalidrefframe` if the frame is not known;
- `:baddescrtimes` unless `first < last`;
- `:unknownspktype` unless `1 <= type <= 1000`.

# References

- [NAIF documentation of `spkpds_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkpds_c.html)
"""
function spkpds(
        ctx::SegContext, body::Integer, center::Integer, frame::SegFrame,
        type::Integer, first::Real, last::Real
    )
    _chkfstr("spkpds", "string `frame`", frame)       # CHKFSTR of spkpds_c
    body == 0 && throw(
        SpiceError(
            :barycenterephem, "spkpds: you've attempted to create " *
                "a segment for the solar system barycenter; this is not supported"
        )
    )
    body == center && throw(
        SpiceError(
            :bodyandcentersame, "spkpds: you've attempted to " *
                "create a segment for a body relative to itself; the body ID code was $body"
        )
    )
    refcod = _seg_refcode("spkpds", ctx, frame)
    first >= last && throw(
        SpiceError(
            :baddescrtimes, "spkpds: the segment start time " *
                "$first is at or after the segment stop time $last"
        )
    )
    (type <= 0 || type > 1000) && throw(
        SpiceError(
            :unknownspktype,
            "spkpds: the type $type is not supported within the SPK system"
        )
    )
    return dafps(Float64[first, last], Int[body, center, refcod, type, 0, 0])
end
spkpds(
    body::Integer, center::Integer, frame::SegFrame, type::Integer, first::Real,
    last::Real
) = spkpds(default_context(), body, center, frame, type, first, last)

function _spk_add!(
        h::DAFHandle, body::Integer, center::Integer, refcod::Int, dtype::Int,
        first::Real, last::Real, segid::AbstractString, data::Vector{Float64}
    )
    _daf_add_array!(h, Float64[first, last], Int[body, center, refcod, dtype, 0, 0], segid, data)
    return nothing
end

function _spk_add_pds!(
        h::DAFHandle, body::Integer, center::Integer, frame::SegFrame,
        dtype::Int, first::Real, last::Real, segid::AbstractString, data::Vector{Float64},
        ctx::SegContext
    )
    _daf_add_descr!(h, spkpds(ctx, body, center, frame, dtype, first, last), segid, data)
    return nothing
end

"`data` followed by every 100th epoch from the first `ndir * 100` epochs (the directory of the epoch table)."
function _seg_directory!(data::Vector{Float64}, epochs::AbstractVector{Float64}, ndir::Int)
    for i in 1:ndir
        push!(data, epochs[100i])
    end
    return data
end

"""
    _seg_coverage_gap(who, first, last, btime, ltime, tol) -> Nothing

Throws `SpiceError(:coveragegap)` if the descriptor interval [`first`, `last`] passes
the data interval [`btime`, `ltime`] by more than `tol`.
"""
function _seg_coverage_gap(
        who::AbstractString, first::Float64, last::Float64,
        btime::Float64, ltime::Float64, tol::Float64
    )
    first < btime - tol && throw(
        SpiceError(
            :coveragegap, "$who: the segment descriptor " *
                "start time $first is too much less than the beginning time of the segment " *
                "data $btime; the difference is $(btime - first) s, the tolerance is $tol s"
        )
    )
    last > ltime + tol && throw(
        SpiceError(
            :coveragegap, "$who: the segment descriptor " *
                "end time $last is too much greater than the end time of the segment data " *
                "$ltime; the difference is $(last - ltime) s, the tolerance is $tol s"
        )
    )
    return nothing
end

"Throws `SpiceError(:invaliddegree)` unless `1 <= degree <= 27`, and `:toofewstates` unless `n > degree`: the Lagrange checks of `spkw08` and `spkw09`."
function _seg_lagrange_degree(who::AbstractString, degree::Integer, n::Int)
    (degree < 1 || degree > 27) && throw(
        SpiceError(
            :invaliddegree, "$who: the " *
                "interpolating polynomials have degree $degree; the valid degree range is [1, 27]"
        )
    )
    n <= degree && throw(
        SpiceError(
            :toofewstates, "$who: at least $(degree + 1) states " *
                "are required to define a polynomial of degree $degree; number of states supplied: $n"
        )
    )
    return nothing
end

"""
    _seg_hermite_degree(who, dtype, degree, n) -> Nothing

The Hermite checks of `spkw12` and `spkw13` (SPK type `dtype`): `:invaliddegree` unless
`degree` is odd and in `1:27`, `:toofewstates` unless `n >= (degree + 1) ÷ 2`.
"""
function _seg_hermite_degree(who::AbstractString, dtype::Int, degree::Integer, n::Int)
    winsiz = (degree + 1) ÷ 2
    (degree < 1 || degree > 27) && throw(
        SpiceError(
            :invaliddegree, "$who: the " *
                "interpolating polynomials have degree $degree; the valid degree range is [1, 27]"
        )
    )
    iseven(degree) && throw(
        SpiceError(
            :invaliddegree, "$who: the interpolating " *
                "polynomials have degree $degree; for SPK type $dtype, the degree must be odd"
        )
    )
    n < winsiz && throw(
        SpiceError(
            :toofewstates, "$who: at least $winsiz states are " *
                "required to define a Hermite polynomial of degree $degree; number of states supplied: $n"
        )
    )
    return nothing
end

"The last checks of `spkw08` and `spkw12`: `:invalidstepsize` unless `step > 0`, then the coverage of the `n` states from `begtim` (`_seg_coverage_gap`)."
function _seg_step_coverage(
        who::AbstractString, first::Float64, last::Float64, begtim::Float64,
        step::Float64, n::Int
    )
    step <= 0.0 && throw(SpiceError(:invalidstepsize, "$who: the step size must be > 0 but was $step"))
    tol = max(abs(first), abs(last)) * 1.0e-13
    ltime = begtim + (n - 1) * step
    _seg_coverage_gap(who, first, last, begtim, ltime, tol)
    return nothing
end

"""
    _seg_epoch_bounds(who, ep, first, last; start=true, kind=:baddescrtimes) -> Nothing

The epoch checks of the SPK writers of types 1, 9, 13, 18 and 21: `:timesoutoforder`
unless `ep` is strictly increasing; if `start`, `:baddescrtimes` if `ep[1] > first`;
then `SpiceError(kind)` if `ep[end] < last`.
"""
function _seg_epoch_bounds(
        who::AbstractString, ep::Vector{Float64}, first::Real, last::Real;
        start::Bool = true, kind::Symbol = :baddescrtimes
    )
    n = length(ep)
    _seg_increasing(who, ep)
    start && ep[1] > first && throw(
        SpiceError(
            :baddescrtimes,
            "$who: segment start time $first precedes first epoch $(ep[1])"
        )
    )
    if ep[n] < last
        gap = kind === :coveragegap ? "segment has coverage gap: " : ""
        throw(SpiceError(kind, "$who: $(gap)segment end time $last follows last epoch $(ep[n])"))
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Open and close
# ---------------------------------------------------------------------------------------

"""
    spkopn(path, ifname, ncomch=0) -> DAFHandle

Create the SPK file `path` and open it for writing, as
`spkopn_c`.

`ifname` is the internal file name and `ncomch` is the number of comment characters to
make room for. The file is on disk when the call returns, byte for byte the file that
`spkopn_c` creates. An existing file is never replaced. Add segments with the `spkwNN`
writers and close the handle with [`spkcls`](@ref). Load the file with [`furnsh`](@ref)
only after the close.

Throws `SpiceError(:emptystring)` if `path` or `ifname` is empty and `:blankfilename` if
`path` is blank, as CSPICE, and `:fileopenfail` if the file exists or cannot be created.

# Examples
Write a type 2 segment that has one record of degree 1, then load the file:
```jldoctest
julia> path = joinpath(mktempdir(), "example.bsp");

julia> h = spkopn(path, "example SPK");

julia> cdata = [[1.0, 0.5, 2.0, 0.0, 3.0, 0.0]];    # coefficients of x, y, z

julia> spkw02(h, 499, 4, "J2000", 0.0, 100.0, "mars", 100.0, cdata, 0.0)

julia> spkcls(h)

julia> ctx = SpiceContext();

julia> furnsh(ctx, path)

julia> pos, lt = spkpos(ctx, 499, 100.0, "J2000", "NONE", 4);

julia> pos == [1.5, 2.0, 3.0]
true
```

See also [`spkopa`](@ref), [`spkw02`](@ref), [`spk14b`](@ref), [`dafac`](@ref) to add
comments.

# References

- [NAIF documentation of `spkopn_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkopn_c.html)
"""
spkopn(path::AbstractString, ifname::AbstractString, ncomch::Integer = 0) =
    _seg_opn("spkopn", "SPK", 6, path, ifname, ncomch)

"""
    _seg_opn(who, type, ni, path, ifname, ncomch) -> DAFHandle

`spkopn`, `ckopn` and `pckopn`: the checks of the C wrapper on `path` and `ifname`, then
a new DAF of the type `type` with 2 doubles and `ni` integers in a summary and room for
`ncomch` comment characters.
"""
function _seg_opn(
        who::AbstractString, type::AbstractString, ni::Int, path::AbstractString,
        ifname::AbstractString, ncomch::Integer
    )
    _chkfstr(who, "string `path`", path)
    _chkfstr(who, "string `ifname`", ifname)
    return dafonw(path, type, 2, ni, ifname, _seg_ncomr(ncomch))
end

"`spkops`, `ckops` and `pckops`: a new DAF in memory, as `_seg_opn` makes one on disk."
_seg_ops(type::AbstractString, ni::Int, ifname::AbstractString, ncomch::Integer) =
    dafops(type, 2, ni, ifname, _seg_ncomr(ncomch))

"The number of comment records that `ncomch` characters need, as `spkopn` computes it."
_seg_ncomr(ncomch::Integer) = ncomch > 0 ? (Int(ncomch) - 1) ÷ _DAF_CRECLEN + 1 : 0

"""
    spkops(ifname="", ncomch=0) -> DAFHandle

Open a new SPK file that lives in memory, for writing: [`spkopn`](@ref) without the
path. No file is made. The writers (such as [`spkw02`](@ref)) take the handle as they take
one of a file and write the same bytes. [`furnsh`](@ref) loads its segments as a kernel,
and [`save_kernel`](@ref) writes them to a file; [`spkcls`](@ref) discards them, so load or
save before the close. This call has no CSPICE counterpart; see [`dafops`](@ref).

Throws the errors of [`dafops`](@ref).
"""
spkops(ifname::AbstractString = "", ncomch::Integer = 0) = _seg_ops("SPK", 6, ifname, ncomch)

"""
    spkopa(path) -> DAFHandle

Open the existing SPK file `path` for writing, to append segments, as
`spkopa_c`.

The file must have the native byte order, and no other open handle may have it open
(see [`dafopw`](@ref)). Unload the file from every context before the call, and load it
again after [`spkcls`](@ref).

Throws `SpiceError(:emptystring)` if `path` is empty, `:blankfilename` if it holds only
blanks, `:filenotfound` if the file does not exist, `:fileisnotspk` if it is not a DAF
file of the type SPK, `:unsupportedbff` if it has the other byte order, and
`:fileopenconflict` if another handle has the file open.

# References

- [NAIF documentation of `spkopa_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkopa_c.html)
"""
function spkopa(path::AbstractString)
    _chk_filename("spkopa", path)
    _isfile(path) || throw(
        SpiceError(
            :filenotfound,
            "spkopa: the file $path is not recognized as an existing file"
        )
    )
    arch, kind = getfat(path)
    (arch == "DAF" && kind == "SPK") || throw(
        SpiceError(
            :fileisnotspk, "spkopa: the file " *
                "$path was not an SPK file; its architecture and type are $arch and $kind"
        )
    )
    return dafopw(path)
end

"Close `h` as `spkcls`/`ckcls` do: a file open for write must hold a segment (`:nosegmentsfound`)."
function _seg_close(h::DAFHandle, what::AbstractString)
    h.isopen || return nothing
    if h.writable
        dafbfs(h)
        daffna(h) || throw(
            SpiceError(
                :nosegmentsfound, "no segments were found in the " *
                    "$what file $(h.path); there must be at least one segment in the file when " *
                    "it is closed"
            )
        )
    end
    dafcls(h)
    return nothing
end

"""
    spkcls(handle) -> Nothing

Close the SPK file of `handle`, as
`spkcls_c`.
The call does nothing if the handle is closed already.

The close writes no data: every write call puts its bytes on disk before it returns,
except [`spk14b`](@ref) and [`spk14a`](@ref), so a type 14 segment that
[`spk14e`](@ref) did not end is not in the file. Load the file with [`furnsh`](@ref)
after the close. For a handle in memory (from [`spkops`](@ref)), the close discards the
bytes: load it with [`furnsh`](@ref) or save it with [`save_kernel`](@ref) before.

Throws `SpiceError(:nosegmentsfound)` if the file is open for writing and holds no
segment. The handle then stays open; close it with [`dafcls`](@ref).

# References

- [NAIF documentation of `spkcls_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcls_c.html)
"""
spkcls(h::DAFHandle) = _seg_close(h, "SPK")

# ---------------------------------------------------------------------------------------
# Types 2, 3 and 20: Chebyshev polynomials on equal intervals
# ---------------------------------------------------------------------------------------

"Checks and data that `spkw02` and `spkw03` share; `per` is 3 (positions) or 6 (states)."
function _spk_cheby(
        who::AbstractString, dtype::Int, per::Int, h::DAFHandle, body, center,
        frame, first, last, segid, intlen, cdata, btime, ctx
    )
    cd, n, ninrec, first, last, intlen, btime = _seg_cheb_args(
        who, "coordinate", per, h, frame, first, last, segid, intlen, cdata, btime
    )
    polydg = ninrec ÷ per - 1
    (polydg < 0 || polydg > 27) && throw(
        SpiceError(
            :invaliddegree, "$who: the " *
                "interpolating polynomials have degree $polydg; the valid degree range is [0, 27]"
        )
    )
    intlen <= 0.0 && throw(SpiceError(:intlennotpos, "$who: the interval length is not positive: $intlen"))
    refcod = _seg_refcode(who, ctx, frame)
    first > last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is greater than the segment end time $last"
        )
    )
    tol = max(abs(first), abs(last)) * 1.0e-13
    ltime = btime + n * intlen
    _seg_coverage_gap(who, first, last, btime, ltime, tol)
    _seg_chckid(who, segid)
    _spk_add!(
        h, body, center, refcod, dtype, first, last, segid,
        _seg_cheb_records(cd, n, ninrec, intlen, btime)
    )
    return nothing
end

"""
    _seg_cheb_args(who, what, per, h, frame, first, last, segid, intlen, cdata, btime)
        -> (cd, n, ninrec, first, last, intlen, btime)

The first checks of `spkw02`, `spkw03` and `pckw02`, which hold `per` series of `what`
coefficients in a record: the empty strings, the handle, the record count
(`:numcoeffsnotpos`) and the record length. Returns the flat coefficients `cd`, the
record count `n`, the record length `ninrec` and the times as `Float64`.
"""
function _seg_cheb_args(
        who::AbstractString, what::AbstractString, per::Int, h::DAFHandle, frame,
        first, last, segid, intlen, cdata, btime
    )
    frame isa AbstractString && _chkfstr(who, "string `frame`", frame)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    cd, n, ninrec = _seg_flat(who, "cdata", cdata, 0)
    first, last, intlen, btime = Float64(first), Float64(last), Float64(intlen), Float64(btime)
    n <= 0 && throw(
        SpiceError(
            :numcoeffsnotpos,
            "$who: the number of sets of $what coefficients is not positive: $n"
        )
    )
    ninrec % per == 0 && ninrec > 0 || throw(
        SpiceError(
            :badarg, "$who: every record of " *
                "`cdata` must hold $per * (degree + 1) coefficients, it holds $ninrec"
        )
    )
    return cd, n, ninrec, first, last, intlen, btime
end

"""
    _seg_cheb_records(cd, n, ninrec, intlen, btime) -> Vector{Float64}

Data area of a Chebyshev segment on `n` intervals of length `intlen` from `btime`, as
SPK types 2 and 3 and binary PCK type 2 store it: each record is the midpoint and the
radius of its interval followed by its `ninrec` coefficients from `cd`, and the trailer
is `btime`, `intlen`, the record size `ninrec + 2` and `n`.
"""
function _seg_cheb_records(
        cd::Vector{Float64}, n::Int, ninrec::Int, intlen::Float64,
        btime::Float64
    )
    @assert length(cd) == n * ninrec
    data = Vector{Float64}(undef, n * (ninrec + 2) + 4)
    radius = intlen / 2.0
    for i in 1:n                         # invariant: records 1 .. i-1 are in data
        base = (i - 1) * (ninrec + 2)
        data[base + 1] = btime + radius + (i - 1) * intlen
        data[base + 2] = radius
        copyto!(data, base + 3, cd, (i - 1) * ninrec + 1, ninrec)
    end
    data[end - 3] = btime
    data[end - 2] = intlen
    data[end - 1] = Float64(ninrec + 2)
    data[end] = Float64(n)
    return data
end

"""
    spkw02([ctx,] handle, body, center, frame, first, last, segid, intlen, cdata, btime)
        -> Nothing

Write an SPK segment of type 2 (Chebyshev polynomials for the position, on intervals of
equal length) to the file of `handle`, as
`spkw02_c`.

# Arguments
`ctx`, `handle`, `body`, `center`, `frame`, `first`, `last` and `segid` have the same
meaning in every `spkwNN` writer.
- `ctx`: the context in which a frame name is resolved (default: the default context).
- `handle`: a [`DAFHandle`](@ref) from [`spkopn`](@ref), [`spkopa`](@ref) or [`spkops`](@ref).
- `body`, `center`: NAIF IDs of the body and of its center of motion.
- `frame`: the reference frame, a name or a frame ID. A name is resolved in `ctx`;
  built-in frames need no kernel.
- `first`, `last`: the coverage of the segment, seconds past J2000 TDB.
- `segid`: the segment identifier, at most 40 printable characters.
- `intlen`, `btime`: the length of every interval (s) and the start of the first interval
  (seconds past J2000 TDB).
- `cdata`: one record for each interval, each with `3 * (polydg + 1)` coefficients (km):
  those of x, then of y, then of z. The number of records `n` and the degree `polydg`
  follow from its shape.

A table such as `cdata` is a vector of vectors (the form of SPICE.jl) or a matrix with one
record in each column. A table of states or of packets, whose record length the data type
fixes, can also be a flat vector; `cdata` cannot. The record count that CSPICE takes as an
argument is not an argument here: it follows from the input.

Every `spkwNN` writer puts its segment on disk when the call returns, with the data that
CSPICE writes from the same input. A call that throws writes nothing.

Throws `SpiceError`. Every `spkwNN` writer throws:
- `:emptystring` if `frame` or `segid` is empty, except `spkw01`, `spkw19` and `spkw21`,
  which have no C function in CSPICE N0067;
- `:dafnosuchhandle` if the handle is closed, `:dafinvalidaccess` if it is open for
  reading;
- `:invalidrefframe` if the frame is not known.

The other checks and their error kinds follow CSPICE writer by writer, so they differ
between writers:
- `spkw02`, `spkw03`, `spkw05`, `spkw08`, `spkw09`, `spkw18`, `spkw19` and `spkw20`
  accept `body == 0` and `body == center`; the others throw `:barycenterephem` and
  `:bodyandcentersame`.
- `spkw02`, `spkw03`, `spkw05`, `spkw18`, `spkw19` and `spkw20` accept `first == last`.
- A `segid` that is too long is `:idstringtoolong` in `spkw02`, `spkw03` and `spkw20`,
  and `:segidtoolong` in the others.

`spkw02` also throws:
- `:numcoeffsnotpos` if `cdata` has no record;
- `:invaliddegree` unless `0 <= polydg <= 27`;
- `:intlennotpos` unless `intlen > 0`;
- `:baddescrtimes` if `first > last`;
- `:coveragegap` unless [`first`, `last`] is inside [`btime`, `btime + n * intlen`], with
  the tolerance `1e-13 * max(abs(first), abs(last))`;
- `:idstringtoolong` or `:nonprintablechars` if `segid` is not valid.

See [`spkopn`](@ref) for an example.

See also [`spkw01`](@ref), [`spkw03`](@ref), [`spkw05`](@ref), [`spkw08`](@ref),
[`spkw09`](@ref), [`spkw10`](@ref), [`spkw12`](@ref), [`spkw13`](@ref), [`spkw15`](@ref),
[`spkw17`](@ref), [`spkw18`](@ref), [`spkw19`](@ref), [`spkw20`](@ref), [`spkw21`](@ref),
[`spk14b`](@ref).

# References

- [NAIF documentation of `spkw02_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw02_c.html)
"""
spkw02(
    ctx::SegContext, h::DAFHandle, body::Integer, center::Integer, frame::SegFrame,
    first::Real, last::Real, segid::AbstractString, intlen::Real, cdata, btime::Real
) =
    _spk_cheby("spkw02", 2, 3, h, body, center, frame, first, last, segid, intlen, cdata, btime, ctx)
spkw02(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real, last::Real,
    segid::AbstractString, intlen::Real, cdata, btime::Real
) =
    spkw02(default_context(), h, body, center, frame, first, last, segid, intlen, cdata, btime)

"""
    spkw03([ctx,] handle, body, center, frame, first, last, segid, intlen, cdata, btime)
        -> Nothing

Write an SPK segment of type 3 (Chebyshev polynomials for the position and the velocity,
on intervals of equal length), as
`spkw03_c`.

`cdata` holds one record for each interval, each with `6 * (polydg + 1)` coefficients:
those of x, y, z (km) and then of dx, dy, dz (km/s). The other arguments and the errors
are those of [`spkw02`](@ref).

# References

- [NAIF documentation of `spkw03_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw03_c.html)
"""
spkw03(
    ctx::SegContext, h::DAFHandle, body::Integer, center::Integer, frame::SegFrame,
    first::Real, last::Real, segid::AbstractString, intlen::Real, cdata, btime::Real
) =
    _spk_cheby("spkw03", 3, 6, h, body, center, frame, first, last, segid, intlen, cdata, btime, ctx)
spkw03(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real, last::Real,
    segid::AbstractString, intlen::Real, cdata, btime::Real
) =
    spkw03(default_context(), h, body, center, frame, first, last, segid, intlen, cdata, btime)

"""
    spkw20([ctx,] handle, body, center, frame, first, last, segid, intlen, cdata, dscale,
           tscale, initjd, initfr) -> Nothing

Write an SPK segment of type 20 (Chebyshev polynomials for the velocity, with the position
at the midpoint of each interval), as
`spkw20_c`.

`cdata` holds one record for each interval, each with `3 * (polydg + 2)` values: the
`polydg + 1` coefficients of the x velocity, then the x position at the midpoint of the
interval, then the same for y and for z. The positions are in units of `dscale` km and the
velocities in units of `dscale` km per `tscale` s; for data in AU and AU/day, `dscale` is
the km of one AU and `tscale` is 86400. `intlen` is the length of every interval in days,
not in seconds as in [`spkw02`](@ref). `initjd + initfr` is the TDB Julian date of the
start of the first interval, so the data cover [`b`, `e`] in seconds past J2000 TDB, with
`b = 86400 * (initjd + initfr - 2451545)` and `e = b + 86400 * n * intlen`. See
[`spkw02`](@ref) for the other arguments and for the errors of every writer.

Throws `SpiceError`:
- `:invalidcount` if `cdata` has no record;
- `:invaliddegree` if `polydg > 50`;
- `:intlennotpos` unless `intlen > 0`;
- `:baddescrtimes` if `first > last`;
- `:nonpositivescale` unless `dscale > 0` and `tscale > 0`;
- `:coveragegap` unless [`first`, `last`] is inside [`b`, `e`], with the tolerance
  `1e-13 * max(abs(b), abs(e))`;
- `:idstringtoolong` or `:nonprintablechars` if `segid` is not valid.

# References

- [NAIF documentation of `spkw20_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw20_c.html)
"""
function spkw20(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, intlen::Real,
        cdata, dscale::Real, tscale::Real, initjd::Real, initfr::Real
    )
    frame isa AbstractString && _chkfstr("spkw20", "string `frame`", frame)       # CHKFSTR of spkw20_c
    _chkfstr("spkw20", "string `segid`", segid)
    refcod, data = _seg_cheb20(
        "spkw20", "SPK", ctx, h, frame, first, last, segid, intlen, cdata, dscale, tscale,
        initjd, initfr
    )
    _spk_add!(h, body, center, refcod, 20, first, last, segid, data)
    return nothing
end

"""
    _seg_cheb20(who, kind, ctx, h, frame, first, last, segid, intlen, cdata, scale, tscale,
                initjd, initfr) -> (refcod, data)

The checks and the data of `spkw20` (`kind == "SPK"`) and `pckw20` (`kind == "PCK"`),
which differ in the test of the coverage (`pckw20` needs `first < last`, `spkw20` only
`first <= last`) and in the name of `scale` (distance or angle).
"""
function _seg_cheb20(
        who::AbstractString, kind::AbstractString, ctx::SegContext, h::DAFHandle,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, intlen::Real,
        cdata, scale::Real, tscale::Real, initjd::Real, initfr::Real
    )
    _daf_check_write(h, who)
    cd, n, ninrec = _seg_flat(who, "cdata", cdata, 0)
    first, last, intlen = Float64(first), Float64(last), Float64(intlen)
    scale, tscale, initjd, initfr = Float64(scale), Float64(tscale), Float64(initjd), Float64(initfr)
    n <= 0 && throw(
        SpiceError(
            :invalidcount,
            "$who: the number of sets of coordinate coefficients is not positive: $n"
        )
    )
    (ninrec % 3 == 0 && ninrec >= 6) || throw(
        SpiceError(
            :badarg, "$who: every record of " *
                "`cdata` must hold 3 * (degree + 2) values, it holds $ninrec"
        )
    )
    polydg = ninrec ÷ 3 - 2
    polydg > 50 && throw(
        SpiceError(
            :invaliddegree, "$who: the interpolating polynomials " *
                "have degree $polydg; the valid degree range is [0, 50]"
        )
    )
    intlen <= 0.0 && throw(SpiceError(:intlennotpos, "$who: the interval length is not positive: $intlen"))
    refcod = _seg_refcode(who, ctx, frame)
    if kind == "PCK"
        first >= last && throw(
            SpiceError(
                :baddescrtimes,
                "$who: the segment start time $first is not less than the segment end time $last"
            )
        )
    else
        first > last && throw(
            SpiceError(
                :baddescrtimes,
                "$who: the segment start time $first is greater than the segment end time $last"
            )
        )
    end
    what = kind == "PCK" ? "angle" : "distance"
    scale <= 0.0 && throw(SpiceError(:nonpositivescale, "$who: the $what scale is not positive: $scale"))
    tscale <= 0.0 && throw(SpiceError(:nonpositivescale, "$who: the time scale is not positive: $tscale"))
    btime = 86400.0 * (initjd - 2451545.0 + initfr)
    ltime = btime + n * intlen * 86400.0
    tol = max(abs(btime), abs(ltime)) * 1.0e-13
    _seg_coverage_gap(who, first, last, btime, ltime, tol)
    _seg_chckid(who, segid, kind)
    data = append!(cd, (scale, tscale, initjd, initfr, intlen, Float64(ninrec), Float64(n)))
    return refcod, data
end
spkw20(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, intlen::Real, cdata, dscale::Real, tscale::Real,
    initjd::Real, initfr::Real
) =
    spkw20(
    default_context(), h, body, center, frame, first, last, segid, intlen, cdata,
    dscale, tscale, initjd, initfr
)

# ---------------------------------------------------------------------------------------
# Types 5, 8, 9, 12, 13, 18: discrete states
# ---------------------------------------------------------------------------------------

"""
    spkw05([ctx,] handle, body, center, frame, first, last, segid, gm, states, epochs)
        -> Nothing

Write an SPK segment of type 5 (two-body propagation between discrete states), as
`spkw05_c`.

`states` holds `n` states (x, y, z in km; dx, dy, dz in km/s), `epochs` holds their `n`
epochs (seconds past J2000 TDB) and `gm` is the gravitational parameter of the central
body (km³/s²). See [`spkw02`](@ref) for the other arguments and for the errors of every
writer.

Throws `SpiceError`:
- `:nonpositivemass` unless `gm > 0`;
- `:numstatesnotpos` if there is no state;
- `:unorderedtimes` unless the epochs are strictly increasing;
- `:baddescrtimes` if `first > last`;
- `:nonprintablechars` or `:segidtoolong` if `segid` is not valid.

# References

- [NAIF documentation of `spkw05_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw05_c.html)
"""
function spkw05(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, gm::Real, states,
        epochs
    )
    who = "spkw05"
    frame isa AbstractString && _chkfstr(who, "string `frame`", frame)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    st, n, _ = _seg_flat(who, "states", states, 6)
    ep = _seg_vec(who, "epochs", epochs, n)
    gm <= 0 && throw(SpiceError(:nonpositivemass, "$who: GM = $gm; non-positive gravitational parameter"))
    refcod = _seg_refcode(who, ctx, frame)
    n <= 0 && throw(
        SpiceError(
            :numstatesnotpos,
            "$who: the number of states and epochs is not positive: $n"
        )
    )
    for i in 2:n
        ep[i] <= ep[i - 1] && throw(SpiceError(:unorderedtimes, "$who: epoch $i is out of order"))
    end
    first > last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is greater than the segment end time $last"
        )
    )
    _seg_idchars(who, segid)
    _seg_idlen(who, segid)
    data = append!(st, ep)
    _seg_directory!(data, ep, n ÷ 100)
    push!(data, Float64(gm), Float64(n))
    _spk_add!(h, body, center, refcod, 5, first, last, segid, data)
    return nothing
end
spkw05(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, gm::Real, states, epochs
) =
    spkw05(
    default_context(), h, body, center, frame, first, last, segid, gm, states,
    epochs
)

"""
    spkw08([ctx,] handle, body, center, frame, first, last, segid, degree, states, epoch1,
           step) -> Nothing

Write an SPK segment of type 8 (Lagrange interpolation of states at equal time steps), as
`spkw08_c`.

`states` holds `n` states (km, km/s). State `k`, counted from 1, is at the epoch
`epoch1 + (k - 1) * step` (seconds past J2000 TDB; `step` in s). `degree` is the degree
of the interpolating polynomial. See [`spkw02`](@ref) for the other arguments and for the
errors of every writer.

Throws `SpiceError`:
- `:segidtoolong` or `:nonprintablechars` if `segid` is not valid;
- `:invaliddegree` unless `1 <= degree <= 27`;
- `:toofewstates` unless `n > degree`;
- `:baddescrtimes` unless `first < last`;
- `:invalidstepsize` unless `step > 0`;
- `:coveragegap` unless [`first`, `last`] is inside [`epoch1`, `epoch1 + (n - 1) * step`],
  with the tolerance `1e-13 * max(abs(first), abs(last))`.

# References

- [NAIF documentation of `spkw08_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw08_c.html)
"""
function spkw08(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, degree::Integer,
        states, epoch1::Real, step::Real
    )
    who = "spkw08"
    frame isa AbstractString && _chkfstr(who, "string `frame`", frame)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    st, n, _ = _seg_flat(who, "states", states, 6)
    first, last, begtim, step = Float64(first), Float64(last), Float64(epoch1), Float64(step)
    refcod = _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    _seg_lagrange_degree(who, degree, n)
    first >= last && throw(
        SpiceError(
            :baddescrtimes, "$who: the segment start time " *
                "$first is greater than or equal to the segment end time $last"
        )
    )
    _seg_step_coverage(who, first, last, begtim, step, n)
    data = push!(st, begtim, step, Float64(degree), Float64(n))
    _spk_add!(h, body, center, refcod, 8, first, last, segid, data)
    return nothing
end
spkw08(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, degree::Integer, states, epoch1::Real, step::Real
) =
    spkw08(
    default_context(), h, body, center, frame, first, last, segid, degree, states,
    epoch1, step
)

"""
    spkw09([ctx,] handle, body, center, frame, first, last, segid, degree, states, epochs)
        -> Nothing

Write an SPK segment of type 9 (Lagrange interpolation of states at unequal time steps),
as `spkw09_c`.

`states` holds `n` states (km, km/s) and `epochs` their `n` epochs (seconds past J2000
TDB). `degree` is the degree of the interpolating polynomial. See [`spkw02`](@ref) for
the other arguments and for the errors of every writer.

Throws `SpiceError`:
- `:baddescrtimes` unless `first < last`, `epochs[1] <= first` and `epochs[n] >= last`;
- `:nonprintablechars` or `:segidtoolong` if `segid` is not valid;
- `:invaliddegree` unless `1 <= degree <= 27`;
- `:toofewstates` unless `n > degree`;
- `:timesoutoforder` unless the epochs are strictly increasing.

# References

- [NAIF documentation of `spkw09_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw09_c.html)
"""
function spkw09(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, degree::Integer,
        states, epochs
    )
    who = "spkw09"
    frame isa AbstractString && _chkfstr(who, "string `frame`", frame)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    st, n, _ = _seg_flat(who, "states", states, 6)
    ep = _seg_vec(who, "epochs", epochs, n)
    refcod = _seg_refcode(who, ctx, frame)
    first >= last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is greater than the segment end time $last"
        )
    )
    _seg_idchars(who, segid)
    _seg_idlen(who, segid)
    _seg_lagrange_degree(who, degree, n)
    _seg_epoch_bounds(who, ep, first, last)
    data = append!(st, ep)
    _seg_directory!(data, ep, (n - 1) ÷ 100)
    push!(data, Float64(degree), Float64(n))
    _spk_add!(h, body, center, refcod, 9, first, last, segid, data)
    return nothing
end
spkw09(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, degree::Integer, states, epochs
) =
    spkw09(
    default_context(), h, body, center, frame, first, last, segid, degree, states,
    epochs
)

"""
    spkw12([ctx,] handle, body, center, frame, first, last, segid, degree, states, epoch1,
           step) -> Nothing

Write an SPK segment of type 12 (Hermite interpolation of states at equal time steps), as
`spkw12_c`.

`states` holds `n` states (km, km/s). State `k`, counted from 1, is at the epoch
`epoch1 + (k - 1) * step` (seconds past J2000 TDB; `step` in s). `degree` is the degree
of the interpolating polynomial. See [`spkw02`](@ref) for the other arguments and for the
errors of every writer.

Throws `SpiceError`:
- `:segidtoolong` or `:nonprintablechars` if `segid` is not valid;
- `:invaliddegree` unless `degree` is odd and in `1:27`;
- `:toofewstates` unless `n >= (degree + 1) ÷ 2`;
- `:baddescrtimes` unless `first < last`;
- `:invalidstepsize` unless `step > 0`;
- `:coveragegap` unless [`first`, `last`] is inside [`epoch1`, `epoch1 + (n - 1) * step`],
  with the tolerance `1e-13 * max(abs(first), abs(last))`;
- `:barycenterephem` if `body == 0`, `:bodyandcentersame` if `body == center`.

# References

- [NAIF documentation of `spkw12_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw12_c.html)
"""
function spkw12(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, degree::Integer,
        states, epoch1::Real, step::Real
    )
    who = "spkw12"
    frame isa AbstractString && _chkfstr(who, "string `frame`", frame)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    st, n, _ = _seg_flat(who, "states", states, 6)
    first, last, begtim, step = Float64(first), Float64(last), Float64(epoch1), Float64(step)
    winsiz = (degree + 1) ÷ 2
    _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    _seg_hermite_degree(who, 12, degree, n)
    first >= last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is greater than the segment end time $last"
        )
    )
    _seg_step_coverage(who, first, last, begtim, step, n)
    data = push!(st, begtim, step, Float64(winsiz - 1), Float64(n))
    _spk_add_pds!(h, body, center, frame, 12, first, last, segid, data, ctx)
    return nothing
end
spkw12(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, degree::Integer, states, epoch1::Real, step::Real
) =
    spkw12(
    default_context(), h, body, center, frame, first, last, segid, degree, states,
    epoch1, step
)

"""
    spkw13([ctx,] handle, body, center, frame, first, last, segid, degree, states, epochs)
        -> Nothing

Write an SPK segment of type 13 (Hermite interpolation of states at unequal time steps),
as `spkw13_c`.

`states` holds `n` states (x, y, z in km; dx, dy, dz in km/s) and `epochs` their `n`
epochs (seconds past J2000 TDB). `degree` is the degree of the interpolating polynomial.
See [`spkw02`](@ref) for the other arguments and for the errors of every writer.

Throws `SpiceError`:
- `:segidtoolong` or `:nonprintablechars` if `segid` is not valid;
- `:invaliddegree` unless `degree` is odd and in `1:27`;
- `:toofewstates` unless `n >= (degree + 1) ÷ 2`;
- `:baddescrtimes` unless `first < last`, `epochs[1] <= first` and `epochs[n] >= last`;
- `:timesoutoforder` unless the epochs are strictly increasing;
- `:barycenterephem` if `body == 0`, `:bodyandcentersame` if `body == center`.

# References

- [NAIF documentation of `spkw13_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw13_c.html)
"""
function spkw13(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, degree::Integer,
        states, epochs
    )
    who = "spkw13"
    frame isa AbstractString && _chkfstr(who, "string `frame`", frame)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    st, n, _ = _seg_flat(who, "states", states, 6)
    ep = _seg_vec(who, "epochs", epochs, n)
    winsiz = (degree + 1) ÷ 2
    _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    _seg_hermite_degree(who, 13, degree, n)
    first >= last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is greater than the segment end time $last"
        )
    )
    _seg_epoch_bounds(who, ep, first, last)
    data = append!(st, ep)
    _seg_directory!(data, ep, (n - 1) ÷ 100)
    push!(data, Float64(winsiz - 1), Float64(n))
    _spk_add_pds!(h, body, center, frame, 13, first, last, segid, data, ctx)
    return nothing
end
spkw13(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, degree::Integer, states, epochs
) =
    spkw13(
    default_context(), h, body, center, frame, first, last, segid, degree, states,
    epochs
)

"Subtype code of an SPK type 18 segment from the code or the SPICE.jl symbol `:S18TP0`, `:S18TP1`."
_spk18_subtype(s::Integer) = Int(s)
_spk18_subtype(s::Symbol) = s === :S18TP0 ? 0 : s === :S18TP1 ? 1 : -1

"""
    spkw18([ctx,] handle, subtyp, body, center, frame, first, last, segid, degree, packts,
           epochs) -> Nothing

Write an SPK segment of type 18 (Hermite or Lagrange interpolation of packets at unequal
time steps), as
`spkw18_c`.

`subtyp` selects the interpolation and the content of a packet (km, km/s, km/s²):

| `subtyp`        | Interpolation | Packet                                            |
|:----------------|:--------------|:--------------------------------------------------|
| `0` or `:S18TP0` | Hermite       | x, y, z, dx/dt, dy/dt, dz/dt, vx, vy, vz, dvx/dt, dvy/dt, dvz/dt |
| `1` or `:S18TP1` | Lagrange      | x, y, z, vx, vy, vz                               |

In subtype 0 the velocity `vx, vy, vz` need not be equal to the derivative of the
position. `packts` holds the `n` packets and `epochs` their `n` epochs (seconds past
J2000 TDB). `degree` is the degree of the interpolating polynomial. The window size is
the number of packets that one interpolation uses: `degree + 1` for Lagrange and
`(degree + 1) ÷ 2` for Hermite. See [`spkw02`](@ref) for the other arguments and for the
errors of every writer.

Throws `SpiceError`:
- `:invalidvalue` if `subtyp` is not known;
- `:segidtoolong` or `:nonprintablechars` if `segid` is not valid;
- `:invaliddegree` unless `1 <= degree <= 15` and the window size is even: the degree is
  3 modulo 4 for subtype 0, and odd for subtype 1;
- `:toofewstates` unless `n >= 2`;
- `:baddescrtimes` unless `first <= last`, `epochs[1] <= first` and `epochs[n] >= last`;
- `:timesoutoforder` unless the epochs are strictly increasing.

# References

- [NAIF documentation of `spkw18_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw18_c.html)
"""
function spkw18(
        ctx::SegContext, h::DAFHandle, subtyp::Union{Integer, Symbol}, body::Integer,
        center::Integer, frame::SegFrame, first::Real, last::Real, segid::AbstractString,
        degree::Integer, packts, epochs
    )
    who = "spkw18"
    frame isa AbstractString && _chkfstr(who, "string `frame`", frame)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    sub = _spk18_subtype(subtyp)
    (sub == 0 || sub == 1) || throw(
        SpiceError(
            :invalidvalue,
            "$who: unexpected SPK type 18 subtype requested: $subtyp"
        )
    )
    packsz = sub == 0 ? 12 : 6
    winsiz = sub == 0 ? (degree + 1) ÷ 2 : degree + 1
    pk, n, _ = _seg_flat(who, "packts", packts, packsz)
    ep = _seg_vec(who, "epochs", epochs, n)
    refcod = _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    (degree < 1 || degree > 15) && throw(
        SpiceError(
            :invaliddegree, "$who: the " *
                "interpolating polynomials have degree $degree; the valid degree range is [1, 15]"
        )
    )
    isodd(winsiz) && throw(
        SpiceError(
            :invaliddegree, "$who: the interpolating " *
                "polynomials have degree $degree; for SPK type 18, the degree must be " *
                "equivalent to 3 mod 4 for Hermite interpolation and odd for Lagrange interpolation"
        )
    )
    n < 2 && throw(
        SpiceError(
            :toofewstates,
            "$who: at least 2 packets are required for SPK type 18; number of packets supplied: $n"
        )
    )
    first > last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is greater than the segment end time $last"
        )
    )
    _seg_epoch_bounds(who, ep, first, last)
    data = append!(pk, ep)
    _seg_directory!(data, ep, (n - 1) ÷ 100)
    push!(data, Float64(sub), Float64(winsiz), Float64(n))
    _spk_add!(h, body, center, refcod, 18, first, last, segid, data)
    return nothing
end
spkw18(
    h::DAFHandle, subtyp::Union{Integer, Symbol}, body::Integer, center::Integer,
    frame::SegFrame, first::Real, last::Real, segid::AbstractString, degree::Integer,
    packts, epochs
) =
    spkw18(
    default_context(), h, subtyp, body, center, frame, first, last, segid, degree,
    packts, epochs
)

# ---------------------------------------------------------------------------------------
# Types 15 and 17: precessing conic elements
# ---------------------------------------------------------------------------------------

"""
    spkw15([ctx,] handle, body, center, frame, first, last, segid, epoch, tp, pa, p, ecc,
           j2flg, pv, gm, j2, radius) -> Nothing

Write an SPK segment of type 15 (precessing conic propagation), as
`spkw15_c`.

- `epoch`: the epoch of periapsis, seconds past J2000 TDB.
- `tp`: the trajectory pole, a vector along the orbital angular momentum.
- `pa`: a vector from the central body to the periapsis.
- `p`: the semi-latus rectum, km.
- `ecc`: the eccentricity.
- `j2flg`: which J2 corrections the propagation applies. 1: regression of the nodes
  only. 2: precession of the line of apsides only. 3: none. Any other value: both.
- `pv`: a vector along the north pole of the central body.
- `gm`, `j2`, `radius`: the gravitational parameter (km³/s²), the J2 coefficient and the
  equatorial radius (km) of the central body.

`tp`, `pa` and `pv` are in `frame`; only their directions are stored. See
[`spkw02`](@ref) for the other arguments and for the errors of every writer.

Throws `SpiceError`:
- `:badlatusrectum` unless `p > 0`, `:badeccentricity` if `ecc < 0`;
- `:nonpositivemass` unless `gm > 0`, `:badradius` if `radius < 0`;
- `:badvector` if `tp`, `pa` or `pv` is zero;
- `:badinitstate` unless the unit vectors of `tp` and `pa` are orthogonal within 1e-5;
- `:segidtoolong` or `:nonprintablechars` if `segid` is not valid;
- the error kinds of [`spkpds`](@ref).

# References

- [NAIF documentation of `spkw15_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw15_c.html)
"""
function spkw15(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, epoch::Real, tp,
        pa, p::Real, ecc::Real, j2flg::Real, pv, gm::Real, j2::Real, radius::Real
    )
    who = "spkw15"
    frame isa AbstractString && _chkfstr(who, "string `frame`", frame)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    tp = _seg_fixed(who, "tp", tp, 3)
    pa = _seg_fixed(who, "pa", pa, 3)
    pv = _seg_fixed(who, "pv", pv, 3)
    mypa = _m_hat(_m_v3(pa))
    mytp = _m_hat(_m_v3(tp))
    record = Float64[epoch, mytp..., mypa..., p, ecc, j2flg, _m_hat(_m_v3(pv))..., gm, j2, radius]
    p <= 0 && throw(SpiceError(:badlatusrectum, "$who: the semi-latus rectum must be positive; it was $p"))
    ecc < 0 && throw(SpiceError(:badeccentricity, "$who: the eccentricity must be non-negative; it was $ecc"))
    gm <= 0 && throw(SpiceError(:nonpositivemass, "$who: the mass of the central body must be positive; GM was $gm"))
    all(iszero, tp) && throw(SpiceError(:badvector, "$who: the trajectory pole vector had length zero"))
    all(iszero, pa) && throw(SpiceError(:badvector, "$who: the periapse vector had length zero"))
    all(iszero, pv) && throw(SpiceError(:badvector, "$who: the central pole vector had length zero"))
    radius < 0 && throw(SpiceError(:badradius, "$who: the central body radius was negative: $radius"))
    dot = mypa[1] * mytp[1] + mypa[2] * mytp[2] + mypa[3] * mytp[3]
    abs(dot) > 1.0e-5 && throw(
        SpiceError(
            :badinitstate,
            "$who: the periapsis and trajectory pole vectors are not orthogonal (dot product $dot)"
        )
    )
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    _spk_add_pds!(h, body, center, frame, 15, first, last, segid, record, ctx)
    return nothing
end
spkw15(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, epoch::Real, tp, pa, p::Real, ecc::Real, j2flg::Real,
    pv, gm::Real, j2::Real, radius::Real
) =
    spkw15(
    default_context(), h, body, center, frame, first, last, segid, epoch, tp, pa, p,
    ecc, j2flg, pv, gm, j2, radius
)

"""
    spkw17([ctx,] handle, body, center, frame, first, last, segid, epoch, eqel, rapol,
           decpol) -> Nothing

Write an SPK segment of type 17 (equinoctial elements), as
`spkw17_c`.

`eqel` holds the nine equinoctial elements at `epoch` (seconds past J2000 TDB). `rapol`
and `decpol` are the right ascension and the declination (radians) in `frame` of the pole
of the central body, which is the pole of the reference plane of the elements. See
[`spkw02`](@ref) for the other arguments and for the errors of every writer.

| Index | Element                                                       |
|:------|:--------------------------------------------------------------|
| 1     | semi-major axis, km                                           |
| 2     | `e * sin(argp + node)`                                        |
| 3     | `e * cos(argp + node)`                                        |
| 4     | mean longitude `m0 + argp + node` at `epoch`, rad             |
| 5     | `tan(inc / 2) * sin(node)`                                    |
| 6     | `tan(inc / 2) * cos(node)`                                    |
| 7     | rate of the longitude of periapsis, rad/s                     |
| 8     | rate of the mean longitude, rad/s                             |
| 9     | rate of the longitude of the ascending node, rad/s            |

Here `e` is the eccentricity, `inc` the inclination, `argp` the argument of periapsis,
`node` the longitude of the ascending node and `m0` the mean anomaly.

Throws `SpiceError`:
- `:badsemiaxis` unless the semi-major axis `eqel[1]` is positive;
- `:badeccentricity` if the eccentricity `sqrt(eqel[2]^2 + eqel[3]^2)` is greater than
  0.9;
- `:segidtoolong` or `:nonprintablechars` if `segid` is not valid;
- the error kinds of [`spkpds`](@ref).

# References

- [NAIF documentation of `spkw17_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkw17_c.html)
"""
function spkw17(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, epoch::Real, eqel,
        rapol::Real, decpol::Real
    )
    who = "spkw17"
    frame isa AbstractString && _chkfstr(who, "string `frame`", frame)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    el = _seg_fixed(who, "eqel", eqel, 9)
    record = Float64[epoch, el..., rapol, decpol]
    a = record[2]
    ecc = sqrt(record[3] * record[3] + record[4] * record[4])
    a <= 0.0 && throw(SpiceError(:badsemiaxis, "$who: the semimajor axis must be positive; it was $a"))
    ecc > 0.9 && throw(
        SpiceError(
            :badeccentricity,
            "$who: the eccentricity $ecc is greater than 0.9"
        )
    )
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    _spk_add_pds!(h, body, center, frame, 17, first, last, segid, record, ctx)
    return nothing
end
spkw17(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, epoch::Real, eqel, rapol::Real, decpol::Real
) =
    spkw17(
    default_context(), h, body, center, frame, first, last, segid, epoch, eqel,
    rapol, decpol
)

# ---------------------------------------------------------------------------------------
# Types 1 and 21: modified difference arrays
# ---------------------------------------------------------------------------------------

"""
    spkw01([ctx,] handle, body, center, frame, first, last, segid, dlines, epochs)
        -> Nothing

Write an SPK segment of type 1 (modified difference arrays), as
`spkw01`
of SPICELIB. CSPICE has no C wrapper for it.

`dlines` holds `n` difference lines, each the 71 doubles of one modified difference
array: a vector of vectors, a matrix with one line in each column, or a flat vector of
`71n` doubles. `epochs` holds the final epoch of each line (seconds past J2000 TDB).
See [`spkw02`](@ref) for the other arguments and for the errors of every writer.

Throws `SpiceError`:
- `:segidtoolong` or `:nonprintablechars` if `segid` is not valid;
- `:invalidcount` if there is no difference line;
- `:baddescrtimes` unless `first < last` and `epochs[n] >= last`;
- `:timesoutoforder` unless the epochs are strictly increasing;
- `:barycenterephem` if `body == 0`, `:bodyandcentersame` if `body == center`.

# References

- [NAIF documentation of `spkw01`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/spkw01.html)
"""
function spkw01(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, dlines, epochs
    )
    who = "spkw01"
    _daf_check_write(h, who)
    dl, n, _ = _seg_flat(who, "dlines", dlines, 71)
    ep = _seg_vec(who, "epochs", epochs, n)
    _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    n < 1 && throw(
        SpiceError(
            :invalidcount,
            "$who: the difference line count was $n; the count must be at least one"
        )
    )
    first >= last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is greater than the segment end time $last"
        )
    )
    _seg_epoch_bounds(who, ep, first, last; start = false)
    data = append!(dl, ep)
    _seg_directory!(data, ep, n ÷ 100)
    push!(data, Float64(n))
    _spk_add_pds!(h, body, center, frame, 1, first, last, segid, data, ctx)
    return nothing
end
spkw01(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, dlines, epochs
) =
    spkw01(default_context(), h, body, center, frame, first, last, segid, dlines, epochs)

"""
    spkw21([ctx,] handle, body, center, frame, first, last, segid, dlines, epochs)
        -> Nothing

Write an SPK segment of type 21 (extended modified difference arrays), as
`spkw21`
of SPICELIB. CSPICE has no C wrapper for it.

`dlines` holds `n` difference lines of equal size `dlsize`: a matrix with one line in each
column, or a vector of vectors. `n` and `dlsize` follow from its shape. `epochs` holds
the final epoch of each line (seconds past J2000 TDB). See [`spkw02`](@ref) for the other
arguments and for the errors of every writer.

Throws `SpiceError`:
- `:difflinetoolarge` if `dlsize > 111`, `:difflinetoosmall` if `dlsize < 71`;
- `:segidtoolong` or `:nonprintablechars` if `segid` is not valid;
- `:invalidcount` if there is no difference line;
- `:baddescrtimes` unless `first < last`;
- `:timesoutoforder` unless the epochs are strictly increasing;
- `:coveragegap` if `epochs[n] < last`;
- `:zerostep` if a step size inside a line is zero;
- `:barycenterephem` if `body == 0`, `:bodyandcentersame` if `body == center`.

# References

- [NAIF documentation of `spkw21`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/spkw21.html)
"""
function spkw21(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, dlines, epochs
    )
    who = "spkw21"
    _daf_check_write(h, who)
    dl, n, dlsize = _seg_flat(who, "dlines", dlines, 0)
    ep = _seg_vec(who, "epochs", epochs, n)
    dlsize > 111 && throw(
        SpiceError(
            :difflinetoolarge,
            "$who: the input difference line size is $dlsize, the maximum supported is 111"
        )
    )
    dlsize < 71 && throw(
        SpiceError(
            :difflinetoosmall,
            "$who: the input difference line size is $dlsize, the minimum supported is 71"
        )
    )
    _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    n < 1 && throw(
        SpiceError(
            :invalidcount,
            "$who: the difference line count was $n; the count must be at least one"
        )
    )
    first > last && throw(
        SpiceError(
            :baddescrtimes,
            "$who: the segment start time $first is greater than the segment end time $last"
        )
    )
    _seg_epoch_bounds(who, ep, first, last; start = false, kind = :coveragegap)
    maxdim = (dlsize - 11) ÷ 4
    kqmloc = 4 * maxdim + 8
    for i in 1:n                         # invariant: lines 1 .. i-1 have no zero step
        base = (i - 1) * dlsize
        kqmax1 = round(Int, dl[base + kqmloc], RoundNearestTiesAway)
        for j in 2:(kqmax1 - 1)
            dl[base + j] == 0.0 && throw(
                SpiceError(
                    :zerostep, "$who: step size was zero " *
                        "at step size vector index $(j - 1) within difference line $i"
                )
            )
        end
    end
    data = append!(dl, ep)
    _seg_directory!(data, ep, n ÷ 100)
    push!(data, Float64(maxdim), Float64(n))
    _spk_add_pds!(h, body, center, frame, 21, first, last, segid, data, ctx)
    return nothing
end
spkw21(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, dlines, epochs
) =
    spkw21(default_context(), h, body, center, frame, first, last, segid, dlines, epochs)

# ---------------------------------------------------------------------------------------
# Type 19: mini-segments (ESOC/DDID, piecewise)
# ---------------------------------------------------------------------------------------

const _SPK19_PKTSZS = (12, 6, 6)

# What differs between the mini-segment writers `spkw19` and `ckw06` (`_CK06_MINI` in
# ck_write.jl), apart from the clock rates and the quaternion checks of `ckw06`: the
# packet sizes, the largest subtype and degree, the words of the messages, and the
# bound that the last epoch of a mini-segment is compared with (`ckw06` of SPICELIB
# compares it with the start of the interval, not with its end).
const _SPK19_MINI = (
    pktszs = _SPK19_PKTSZS, maxsub = 2, maxdeg = 27, endbound = true,
    kname = "SPK type 19", count = "mini-segment/interval", bounds = "interval",
    ivl = "interpolation interval", union = "interpolation intervals", times = "epochs",
)

"The number of mini-segments, `length(npkts)`; throws `SpiceError(:badarg)` unless `subtps` and `degres` have that length."
function _seg_mini_count(who::AbstractString, npkts, subtps, degres)
    n = length(npkts)
    (length(subtps) == n && length(degres) == n) || throw(
        SpiceError(
            :badarg,
            "$who: `npkts`, `subtps` and `degres` must have the same length"
        )
    )
    return n
end

"""
    _seg_minis(who, par, ctx, frame, first, last, segid, npkts, subtps, degres, pk, ep,
               rt, bd, sellst, pktcheck) -> (refcod, data)

The checks and the data of a segment of mini-segments, as `spkw19` and `ckw06` of
SPICELIB make them after their argument checks: `par` is `_SPK19_MINI` or `_CK06_MINI`,
`pk`, `ep` and `bd` the packets, epochs and interval bounds, `rt` the clock rates of
`ckw06` (`nothing` for `spkw19`), and `pktcheck(i, pktbeg, n, pktsiz, subtyp)` the
checks of the packets of mini-segment `i` that follow its bound checks (`nothing`: none).
"""
function _seg_minis(
        who::AbstractString, par::NamedTuple, ctx::SegContext, frame::SegFrame,
        first::Real, last::Real, segid::AbstractString, npkts, subtps, degres,
        pk::Vector{Float64}, ep::Vector{Float64}, rt::Union{Nothing, Vector{Float64}},
        bd::Vector{Float64}, sellst::Bool, pktcheck
    )
    nintvl = length(npkts)
    last < first && throw(
        SpiceError(
            :baddescrtimes, "$who: segment start time is " *
                "$first; stop time is $last; bounds must be in nondecreasing order"
        )
    )
    refcod = _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    nintvl < 1 && throw(
        SpiceError(
            :invalidcount,
            "$who: $(par.count) count was $nintvl; this count must be positive"
        )
    )
    for i in 1:nintvl
        bd[i] >= bd[i + 1] && throw(
            SpiceError(
                :boundsoutoforder, "$who: $(par.bounds) bounds at " *
                    "indices $i and $(i + 1) are $(bd[i]) and $(bd[i + 1]); they must be strictly increasing"
            )
        )
    end
    (first < bd[1] || last > bd[nintvl + 1]) && throw(
        SpiceError(
            :coveragegap, "$who: " *
                "segment coverage [$first, $last] must be contained within the union of the " *
                "$(par.union) [$(bd[1]), $(bd[nintvl + 1])]"
        )
    )
    winsizes = zeros(Int, nintvl)
    pktend = 0
    eepix = 0
    for i in 1:nintvl                    # invariant: mini-segments 1 .. i-1 are valid
        n = Int(npkts[i])
        n < 2 && throw(
            SpiceError(
                :toofewpackets, "$who: at least 2 packets are " *
                    "required for $(par.kname); $n were supplied in mini-segment $i"
            )
        )
        subtyp = Int(subtps[i])
        (subtyp < 0 || subtyp > par.maxsub) && throw(
            SpiceError(
                :invalidsubtype,
                "$who: unexpected $(par.kname) subtype $subtyp found in mini-segment $i"
            )
        )
        pktsiz = par.pktszs[subtyp + 1]
        winsiz = isodd(subtyp) ? degres[i] + 1 : (degres[i] + 1) ÷ 2    # odd: Lagrange
        (rt !== nothing && rt[i] <= 0.0) && throw(
            SpiceError(
                :invalidsclkrate,
                "$who: SCLK rate in mini-segment $i was $(rt[i]); the rate must be positive"
            )
        )
        (degres[i] < 1 || degres[i] > par.maxdeg) && throw(
            SpiceError(
                :invaliddegree, "$who: the " *
                    "interpolating polynomials of mini-segment $i have degree $(degres[i]); the " *
                    "valid degree range is [1, $(par.maxdeg)]"
            )
        )
        isodd(winsiz) && throw(
            SpiceError(
                :badwindowsize, "$who: the interpolating " *
                    "polynomials of mini-segment $i have window size $winsiz and degree $(degres[i])"
            )
        )
        winsizes[i] = winsiz
        pktbeg = pktend + 1
        pktend = pktbeg - 1 + n * pktsiz
        bepix = eepix + 1
        eepix = bepix - 1 + n
        (length(ep) >= eepix && length(pk) >= pktend) || throw(
            SpiceError(
                :badarg,
                "$who: `$(par.times)` or `packts` is shorter than the packet counts require"
            )
        )
        for k in bepix:(eepix - 1)
            ep[k] >= ep[k + 1] && throw(
                SpiceError(
                    :timesoutoforder, "$who: in mini-segment " *
                        "$i, epoch $(ep[k]) having index $k is greater than or equal to its successor"
                )
            )
        end
        ep[bepix] > bd[i] && throw(
            SpiceError(
                :boundsdisagree, "$who: $(par.ivl) " *
                    "$i start time $(bd[i]) precedes mini-segment's first epoch $(ep[bepix])"
            )
        )
        endbd = par.endbound ? bd[i + 1] : bd[i]
        ep[eepix] < endbd && throw(
            SpiceError(
                :boundsdisagree, "$who: $(par.ivl) " *
                    "$i end time $endbd exceeds mini-segment's last epoch $(ep[eepix])"
            )
        )
        pktcheck === nothing || pktcheck(i, pktbeg, n, pktsiz, subtyp)
    end
    data = Float64[]
    pktend = 0
    eepix = 0
    pointers = Float64[]
    segend = 0
    for i in 1:nintvl                    # invariant: mini-segments 1 .. i-1 are in data
        pktsiz = par.pktszs[subtps[i] + 1]
        n = Int(npkts[i])
        append!(data, view(pk, (pktend + 1):(pktend + n * pktsiz)))
        pktend += n * pktsiz
        bepix = eepix + 1
        eepix = bepix - 1 + n
        append!(data, view(ep, bepix:eepix))
        _seg_directory!(data, view(ep, bepix:eepix), (n - 1) ÷ 100)
        rt === nothing || push!(data, rt[i])
        push!(data, Float64(subtps[i]), Float64(winsizes[i]), Float64(n))
        minisz = n * (pktsiz + 1) + (n - 1) ÷ 100 + (rt === nothing ? 3 : 4)
        push!(pointers, Float64(segend + 1))
        segend += minisz
    end
    append!(data, bd)
    _seg_directory!(data, bd, nintvl ÷ 100)
    append!(data, pointers)
    push!(data, Float64(segend + 1), sellst ? 1.0 : -1.0, Float64(nintvl))
    return refcod, data
end

"""
    spkw19([ctx,] handle, body, center, frame, first, last, segid, npkts, subtps, degres,
           packts, epochs, ivlbds, sellst) -> Nothing

Write an SPK segment of type 19 (mini-segments of Hermite or Lagrange interpolation), as
`spkw19`
of SPICELIB. CSPICE has no C wrapper for it.

A mini-segment is a set of packets with its own subtype and degree, valid on its own
interval. The number of mini-segments is `length(npkts)`.

- `npkts[i]`, `subtps[i]`, `degres[i]`: the number of packets, the subtype and the degree
  of mini-segment `i`.
- `packts`: the flat vector of all packets, in order.
- `epochs`: the flat vector of the epochs of all packets, seconds past J2000 TDB.
- `ivlbds`: the `length(npkts) + 1` bounds of the intervals of the mini-segments,
  seconds past J2000 TDB.
- `sellst`: if `true`, an epoch on a bound of two intervals uses the later interval;
  else it uses the earlier one.

Subtypes 0 and 1 have the packets of [`spkw18`](@ref). Subtype 2 is Hermite
interpolation of packets x, y, z, vx, vy, vz (km, km/s). The window size is defined in
[`spkw18`](@ref). See [`spkw02`](@ref) for the other arguments and for the errors of
every writer.

Throws `SpiceError`:
- `:baddescrtimes` if `first > last`;
- `:segidtoolong` or `:nonprintablechars` if `segid` is not valid;
- `:invalidcount` if there is no mini-segment;
- `:boundsoutoforder` unless `ivlbds` is strictly increasing;
- `:coveragegap` unless [`first`, `last`] is inside [`ivlbds[1]`, `ivlbds[end]`].

For a mini-segment it throws:
- `:toofewpackets` if it has fewer than two packets;
- `:invalidsubtype` if its subtype is not in `0:2`;
- `:invaliddegree` if its degree is not in `1:27`;
- `:badwindowsize` if its window size is odd;
- `:timesoutoforder` unless its epochs are strictly increasing;
- `:boundsdisagree` unless its epochs cover its interval.

# References

- [NAIF documentation of `spkw19`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/spkw19.html)
"""
function spkw19(
        ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString,
        npkts::AbstractVector{<:Integer}, subtps::AbstractVector{<:Integer},
        degres::AbstractVector{<:Integer}, packts, epochs, ivlbds, sellst::Bool
    )
    who = "spkw19"
    _daf_check_write(h, who)
    nintvl = _seg_mini_count(who, npkts, subtps, degres)
    pk = _seg_vec(who, "packts", packts, -1)
    ep = _seg_vec(who, "epochs", epochs, -1)
    bd = _seg_vec(who, "ivlbds", ivlbds, nintvl < 1 ? -1 : nintvl + 1)
    refcod, data = _seg_minis(
        who, _SPK19_MINI, ctx, frame, first, last, segid, npkts, subtps, degres, pk, ep,
        nothing, bd, sellst, nothing
    )
    _spk_add!(h, body, center, refcod, 19, first, last, segid, data)
    return nothing
end
spkw19(
    h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, npkts::AbstractVector{<:Integer},
    subtps::AbstractVector{<:Integer}, degres::AbstractVector{<:Integer}, packts, epochs,
    ivlbds, sellst::Bool
) =
    spkw19(
    default_context(), h, body, center, frame, first, last, segid, npkts, subtps,
    degres, packts, epochs, ivlbds, sellst
)

"""
    SPK14Writer

An SPK segment of type 14 that is open for the addition of records. [`spk14b`](@ref)
returns it, and it is the first argument of [`spk14a`](@ref) and [`spk14e`](@ref), which
writes the segment.
"""
mutable struct SPK14Writer
    const handle::DAFHandle
    const descr::Vector{Float64}
    const segid::String
    const ncoef::Int
    const packets::Vector{Float64}      # (reference, packet) pairs, as `sgwfpk` writes them
    const refs::Vector{Float64}
    open::Bool
end

Base.show(io::IO, w::SPK14Writer) = _seg_show_writer(io, w)

"The `show` of a segment writer `w` (`SPK14Writer`, `PCK03Writer`, `CKW04Writer`): its type, `segid`, record count and whether it was ended."
_seg_show_writer(io::IO, w) = print(
    io, nameof(typeof(w)), "(\"", w.segid, "\", ",
    length(w.refs), " records", w.open ? "" : ", ended", ")"
)

"Throws `SpiceError(:segmentnotopen)` if the segment writer `w` was ended."
_seg_check_open(who::AbstractString, w) =
    w.open || throw(SpiceError(:segmentnotopen, "$who: the segment '$(w.segid)' was ended"))

"""
    spk14b([ctx,] handle, segid, body, center, frame, first, last, chbdeg) -> SPK14Writer

Begin an SPK segment of type 14 (Chebyshev polynomials for the position and the velocity,
on intervals of unequal length), as
`spk14b_c`,
and return its writer.

`chbdeg` is the degree of the polynomials. The other arguments are those of
[`spkw02`](@ref), but `segid` comes before `body`, as in CSPICE. Add records with
[`spk14a`](@ref) and end the segment with [`spk14e`](@ref).

Throws `SpiceError`:
- `:emptystring` if `frame` or `segid` is empty;
- `:dafnosuchhandle` or `:dafinvalidaccess` if the file is not open for writing;
- `:invalidargument` if `chbdeg < 0`;
- `:segidtoolong` or `:nonprintingchars` (the kind of `sgbwfs`) if `segid` is not valid;
- the error kinds of [`spkpds`](@ref).

# Examples
Two records of degree 1, each on an interval of 100 s:
```jldoctest
julia> path = joinpath(mktempdir(), "type14.bsp");

julia> h = spkopn(path, "type 14 example");

julia> w = spk14b(h, "probe", -99, 399, "J2000", 0.0, 200.0, 1)
SPK14Writer("probe", 0 records)

julia> record(mid, x) = [mid, 50.0, x, 1.0, 0.0, 0.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0];

julia> spk14a(w, [record(50.0, 10.0); record(150.0, 11.0)], [0.0, 100.0])

julia> spk14e(w)

julia> w
SPK14Writer("probe", 2 records, ended)

julia> spkcls(h)

julia> ctx = SpiceContext();

julia> furnsh(ctx, path)

julia> state, lt = spkezr(ctx, -99, 150.0, "J2000", "NONE", 399);

julia> state == [11.0, 0.0, 0.0, 0.02, 0.0, 0.0]
true
```

# References

- [NAIF documentation of `spk14b_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spk14b_c.html)
"""
function spk14b(
        ctx::SegContext, h::DAFHandle, segid::AbstractString, body::Integer,
        center::Integer, frame::SegFrame, first::Real, last::Real, chbdeg::Integer
    )
    who = "spk14b"
    frame isa AbstractString && _chkfstr(who, "string `frame`", frame)
    _chkfstr(who, "string `segid`", segid)
    _daf_check_write(h, who)
    chbdeg < 0 && throw(
        SpiceError(
            :invalidargument, "$who: the degree of the " *
                "Chebyshev polynomial was negative, $chbdeg; the degree must be greater than or equal to zero"
        )
    )
    descr = spkpds(ctx, body, center, frame, 14, first, last)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid, :nonprintingchars)     # the kind of `sgbwfs`
    return SPK14Writer(h, descr, String(segid), Int(chbdeg) + 1, Float64[], Float64[], true)
end
spk14b(
    h::DAFHandle, segid::AbstractString, body::Integer, center::Integer, frame::SegFrame,
    first::Real, last::Real, chbdeg::Integer
) =
    spk14b(default_context(), h, segid, body, center, frame, first, last, chbdeg)

"""
    spk14a(writer, ncsets, coeffs, epochs) -> Nothing
    spk14a(writer, coeffs, epochs) -> Nothing

Add `ncsets` records to the type 14 segment `writer`, as
`spk14a_c`.
Without `ncsets`, the number of records is `length(epochs)`.

A record has `2 + 6 * (chbdeg + 1)` numbers: the midpoint (seconds past J2000 TDB) and
the radius (s) of its interval, then the coefficients of x, y, z (km) and of dx, dy, dz
(km/s). `coeffs` holds the records one after the other: a vector, or a matrix with one
record in each column. `epochs[i]` is the start of the interval of record `i`, in seconds
past J2000 TDB.

Throws `SpiceError`:
- `:segmentnotopen` if `spk14e` ended the writer, `:dafnosuchhandle` if the handle is
  closed;
- `:invalidargument` unless `ncsets > 0`;
- `:badarg` if `coeffs` or `epochs` does not have the size that `ncsets` requires;
- `:unorderedrefs` unless `epochs` is strictly increasing and starts after the last epoch
  of the call before.

See [`spk14b`](@ref) for an example.

# References

- [NAIF documentation of `spk14a_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spk14a_c.html)
"""
function spk14a(w::SPK14Writer, ncsets::Integer, coeffs, epochs)
    who = "spk14a"
    _seg_check_open(who, w)
    _daf_check_write(w.handle, who)
    ncsets <= 0 && throw(
        SpiceError(
            :invalidargument, "$who: the number of coefficient " *
                "sets and epochs to be added to the SPK segment was not positive; its value was $ncsets"
        )
    )
    _seg_add_cheb_packets!(who, w.packets, w.refs, w.ncoef, Int(ncsets), coeffs, epochs)
    return nothing
end
spk14a(w::SPK14Writer, coeffs, epochs) = spk14a(w, length(epochs), coeffs, epochs)

# The checks of `sgwfpk` on `n >= 1` records of `2 + 6 * ncoef` numbers, and their
# addition to `packets` and `refs`.
function _seg_add_cheb_packets!(
        who::String, packets::Vector{Float64}, refs::Vector{Float64}, ncoef::Int, n::Int,
        coeffs, epochs
    )
    pktsiz = 2 + 6 * ncoef
    coeffs isa AbstractArray{<:Real} ||
        throw(SpiceError(:badarg, "$who: `coeffs` must be an array of numbers"))
    cf = _seg_vec(who, "coeffs", vec(collect(Float64, coeffs)), n * pktsiz)
    ep = _seg_vec(who, "epochs", epochs, n)
    _seg_ordered_refs(who, refs, ep)
    for i in 1:n                    # invariant: records 1:i-1 are in `packets`, each after its epoch
        push!(packets, ep[i])
        append!(packets, @view cf[((i - 1) * pktsiz + 1):(i * pktsiz)])
    end
    append!(refs, ep)
    @assert length(packets) == length(refs) * (pktsiz + 1)
    return nothing
end

"""
    _seg_cheb_packets(ncoef, packets, refs) -> Vector{Float64}

The data of a generic segment that `sgwes` writes for SPK type 14 and binary PCK type
3: the constant `ncoef`, the (reference, packet) pairs `packets` of `2 + 6 * ncoef`
numbers each, the references `refs`, their directory and the meta data.
"""
function _seg_cheb_packets(ncoef::Int, packets::Vector{Float64}, refs::Vector{Float64})
    n = length(refs)
    pktsiz = 2 + 6 * ncoef
    data = Float64[ncoef]
    append!(data, packets)
    append!(data, refs)
    ndir = (n - 1) ÷ 100                 # 0 for n = 0 too (`÷` truncates)
    _seg_directory!(data, refs, ndir)
    pktbas = 1
    refbas = pktbas + n * (pktsiz + 1)
    # Meta data of `sgwes` for packets of fixed size with explicit reference values:
    # CONBAS, NCON, RDRBAS, NRDR, RDRTYP, REFBAS, NREF, PDRBAS, NPDR, PDRTYP, PKTBAS,
    # NPKT, RSVBAS, NRSV, PKTSZ, PKTOFF, NMETA.
    append!(
        data, Float64[
            0, 1, refbas + n, ndir, 3, refbas, n, 0, 0, 0, pktbas, n, 0, 0,
            pktsiz, 1, 17,
        ]
    )
    return data
end

"""
    spk14e(writer) -> Nothing

End the type 14 segment `writer` and write it to its file, as
`spk14e_c`.

[`spk14b`](@ref) and [`spk14a`](@ref) do not change the file, so a segment that is never
ended is not in the file. The segment is on disk when the call returns, with the data that
CSPICE writes. The writer cannot be used again.

Throws `SpiceError(:segmentnotopen)` if the writer was ended before, and
`SpiceError(:dafnosuchhandle)` if the handle was closed after [`spk14b`](@ref); the
segment is then lost.

# References

- [NAIF documentation of `spk14e_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spk14e_c.html)
"""
spk14e(w::SPK14Writer) = _seg_cheb_end!("spk14e", w)

"`spk14e` and `pck03e`: write the segment of the writer `w` (`SPK14Writer` or `PCK03Writer`) and end it."
function _seg_cheb_end!(who::AbstractString, w)
    _seg_check_open(who, w)
    _daf_check_write(w.handle, who)
    _daf_add_descr!(w.handle, w.descr, w.segid, _seg_cheb_packets(w.ncoef, w.packets, w.refs))
    w.open = false
    return nothing
end
