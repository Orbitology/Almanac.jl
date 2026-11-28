# Spacecraft clocks (SCLK), type 1, as `sc01.f` and its callers in CSPICE.
#
# An encoded SCLK value ("ticks") counts the ticks of the least significant clock field
# from the start of the first partition. A type 1 clock relates ticks to a parallel time
# system (TDB or TDT) by a piecewise linear function, given as coefficient records
# (ticks, parallel time, rate). The data of every clock are read from the kernel pool
# once, when the `KernelSet` snapshot is built; the conversions are pure functions of
# that snapshot.

"Delimiters that separate the fields of a clock string on input (`sctk01`)."
const _SCLK_DELIMS = ".:-, "
"Output delimiters by code 1 to 5 (`SCLK01_OUTPUT_DELIM`)."
const _SCLK_OUTPUT_DELIM = ('.', ':', '-', ',', ' ')
"Largest number of clock fields (`MXNFLD`)."
const _SCLK_MXNFLD = 10

# `_dnint` (the `d_nint` of f2c) is in `time.jl`.

"""
    _sclk_build(pool, sc) -> Union{SclkType1,Nothing}

Data of the type 1 clock `sc` from the pool variables with the suffix `_<-sc>`, or
`nothing` if a required variable is absent or malformed. Required: `SCLK01_COEFFICIENTS`
(a positive multiple of three values), `SCLK_PARTITION_START` and `SCLK_PARTITION_END`
(equal counts, at least one), `SCLK01_MODULI`, `SCLK01_OFFSETS`, `SCLK01_N_FIELDS` (1 to
10, not more than the moduli and the offsets) and `SCLK01_OUTPUT_DELIM`.
`SCLK01_TIME_SYSTEM` is optional and defaults to 1 (TDB), as in `zzscad01`.
"""
function _sclk_build(pool::AbstractDict, sc::Integer)
    suffix = "_$(-Int(sc))"
    num(name) = _pool_numbers(pool, name * suffix)
    coef = num("SCLK01_COEFFICIENTS")
    pstart = num("SCLK_PARTITION_START")
    pstop = num("SCLK_PARTITION_END")
    moduli = num("SCLK01_MODULI")
    offsets = num("SCLK01_OFFSETS")
    nf = num("SCLK01_N_FIELDS")
    del = num("SCLK01_OUTPUT_DELIM")
    any(isnothing, (coef, pstart, pstop, moduli, offsets, nf, del)) && return nothing
    (length(coef) >= 3 && length(coef) % 3 == 0) || return nothing
    (length(pstart) == length(pstop) && !isempty(pstart)) || return nothing
    (length(nf) == 1 && length(del) == 1) || return nothing
    # The integer items are read by `scli01_` (`gipool_`); a value outside its range
    # leaves the clock without data.
    nf1 = _gipool_first(nf, ""; strict = false)
    delcde = _gipool_first(del, ""; strict = false)
    (nf1 === nothing || delcde === nothing) && return nothing
    nfield = Int(nf1)
    (
        1 <= nfield <= _SCLK_MXNFLD && length(moduli) >= nfield &&
            length(offsets) >= nfield
    ) || return nothing
    ts = num("SCLK01_TIME_SYSTEM")
    timsys = (ts === nothing || isempty(ts)) ? Int32(1) : _gipool_first(ts, ""; strict = false)
    timsys === nothing && return nothing
    tikmsc = 1.0
    for i in nfield:-1:2              # invariant: tikmsc = product of moduli[i+1:nfield]
        tikmsc *= moduli[i]
    end
    mxtick = 0.0
    for i in eachindex(pstart)        # invariant: mxtick = rounded length of partitions 1:i-1
        mxtick = _dnint(pstop[i] - pstart[i] + mxtick)
    end
    return SclkType1(
        Int32(sc), coef, pstart, pstop, moduli, offsets, Int32(nfield),
        delcde, timsys, tikmsc, mxtick
    )
end

"""
    build_sclk(pool) -> Dict{Int32,SclkType1}

Type 1 clocks of the pool, by clock ID: one for each variable `SCLK_DATA_TYPE_<n>` whose
first value rounds to 1 (the clock ID is `-n`) whose `SCLK01_*` data are complete. A clock with
incomplete data is left out, so that loading a kernel never fails; a conversion with
that clock then throws `SpiceError`.
"""
function build_sclk(pool::AbstractDict)
    out = Dict{Int32, SclkType1}()
    prefix = "SCLK_DATA_TYPE_"
    for (name, value) in pool
        startswith(name, prefix) || continue
        n = tryparse(Int32, SubString(name, ncodeunits(prefix) + 1))
        n === nothing && continue
        # The type is the first value, rounded, as `gipool_(name, 1, 1, ...)` of `scty01_`.
        _gipool_first(value, name; strict = false) == 1 || continue
        clock = _sclk_build(pool, -n)
        clock === nothing || (out[-n] = clock)
    end
    return out
end

"""
    sctype([ctx,] sc) -> Int

Return the data type of the spacecraft clock `sc`, as
`sctype` of SPICELIB.

`sc` is the NAIF ID of the spacecraft, a negative number. The value is the one of
`SCLK_DATA_TYPE_<-sc>`. The conversion calls of the package support type 1.

Throws `SpiceError(:kernelvarnotfound)` if no SCLK kernel for `sc` is loaded.

# References

- [NAIF documentation of `sctype`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/sctype.html)
"""
function sctype(ks::KernelSet, sc::Integer)
    name = "SCLK_DATA_TYPE_$(-Int(sc))"
    v = _pool_numbers(ks.pool, name)
    (v === nothing || isempty(v)) && throw(
        SpiceError(
            :kernelvarnotfound,
            "kernel variable $name was not found in the kernel pool; no SCLK kernel is " *
                "loaded for the clock $sc"
        )
    )
    return Int(_gipool_first(v, name)::Int32)
end

"Throw the error that tells why the clock `sc` has no data in `ks.sclk`."
@noinline function _sclk_missing(ks::KernelSet, sc::Integer)
    t = sctype(ks, sc)
    t == 1 || throw(SpiceError(:notsupported, "clock type $t of the clock $sc is not supported"))
    suffix = "_$(-Int(sc))"
    ps = _pool_numbers(ks.pool, "SCLK_PARTITION_START" * suffix)
    pe = _pool_numbers(ks.pool, "SCLK_PARTITION_END" * suffix)
    # `zzscup01` signals NUMPARTSUNEQUAL for start and stop lists of different lengths
    (ps !== nothing && pe !== nothing && length(ps) != length(pe)) && throw(
        SpiceError(
            :numpartsunequal, "the clock $sc has $(length(ps)) partition start " *
                "times and $(length(pe)) partition end times"
        )
    )
    throw(
        SpiceError(
            :kernelvarnotfound, "the SCLK01_* or SCLK_PARTITION_* kernel " *
                "variables of the spacecraft clock $sc are absent or malformed"
        )
    )
end

"""
    _sclk(ks, sc) -> SclkType1

Data of the type 1 clock `sc`. Throws `SpiceError` with kind `:kernelvarnotfound` (no or
incomplete data), `:numpartsunequal` (partition start and end lists of different
lengths) or `:notsupported` (a clock type other than 1).
"""
@inline function _sclk(ks::KernelSet, sc::Integer)
    clock = get(ks.sclk, Int32(sc), nothing)
    clock === nothing && _sclk_missing(ks, sc)
    return clock
end

@noinline _sclk_range(sc, what) = throw(
    SpiceError(
        :valueoutofrange,
        "$what of the spacecraft clock $sc is outside the range of the SCLK kernel"
    )
)
@noinline _sclk_badrate(sc) = throw(
    SpiceError(
        :valueoutofrange,
        "invalid SCLK rate (not positive) for the spacecraft clock $sc"
    )
)
@noinline _sclk_badsys(sc, sys) = throw(
    SpiceError(
        :valueoutofrange,
        "invalid time system code $sys was found for the spacecraft clock $sc"
    )
)
@noinline _sclk_nolsk() = throw(
    SpiceError(
        :missingtimeinfo,
        "a spacecraft clock with parallel time system TDT needs a leap-second kernel"
    )
)

"""
    _sclk_record(coef, col, x) -> Int

Index of the coefficient record to use for the value `x` of column `col` (1 ticks,
2 parallel time): the last record whose value is `<= x`, and the first record if `x` is
before all, found with the bisection of `scte01`. Requires `length(coef) >= 3`.
"""
@inline function _sclk_record(coef::Vector{Float64}, col::Int, x::Float64)
    n = length(coef) ÷ 3
    @inbounds x < coef[3 * (n - 1) + col] || return n
    lower, upper = 1, n
    while lower < upper - 1     # invariant: record `lower` is <= x or lower == 1; record `upper` is > x
        middle = (lower + upper) ÷ 2
        if @inbounds x < coef[3 * (middle - 1) + col]
            upper = middle
        else
            lower = middle
        end
    end
    return lower
end

"`x` converted by `unitim_` from `insys` to `outsys` if the parallel time of `clock` is TDT (2); `x` as it is for TDB (1)."
@inline function _sclk_unitim(ks::KernelSet, clock::SclkType1, x::Real, insys, outsys)
    clock.timsys == 1 && return x
    clock.timsys == 2 || _sclk_badsys(clock.id, clock.timsys)
    l = ks.lsk
    l === nothing && _sclk_nolsk()
    return _unitim(l, x, insys, outsys)
end

"TDB seconds past J2000 of the parallel time `partim` of `clock` (`unitim_`), generic over the number type."
@inline _sclk_to_tdb(ks::KernelSet, clock::SclkType1, partim::Real) =
    _sclk_unitim(ks, clock, partim, _US_TDT, _US_TDB)

"Parallel time of `clock` for the TDB epoch `et` (`unitim_`, three fixed-point steps), generic over the number type."
@inline _sclk_from_tdb(ks::KernelSet, clock::SclkType1, et::Real) =
    _sclk_unitim(ks, clock, et, _US_TDB, _US_TDT)

"""
    _sct2e(ks, clock, sclkdp::T) -> T

ET of the encoded SCLK `sclkdp` (`scte01`). The function is piecewise linear, so a
`Taylor` or dual number passes through; the record is chosen with `primal`. Throws
`SpiceError(:valueoutofrange)` for ticks outside the clock range.
"""
@inline function _sct2e(ks::KernelSet, clock::SclkType1, sclkdp::T) where {T <: Real}
    coef = clock.coef
    p = primal(sclkdp)
    (p < coef[1] || p > clock.mxtick || isnan(p)) && _sclk_range(clock.id, "the encoded SCLK $p")
    k = 3 * (_sclk_record(coef, 1, p) - 1)
    @inbounds c1, c2, c3 = coef[k + 1], coef[k + 2], coef[k + 3]
    c3 <= 0.0 && _sclk_badrate(clock.id)
    rate = c3 / clock.tikmsc
    return _sclk_to_tdb(ks, clock, c2 + rate * (sclkdp - c1))
end

"""
    _sce2c(ks, clock, et::T) -> T

Continuous encoded SCLK of the ET `et` (`scec01`), piecewise linear and generic over the
number type as `_sct2e`. Throws `SpiceError(:valueoutofrange)` for an epoch outside the
clock range.
"""
@inline function _sce2c(ks::KernelSet, clock::SclkType1, et::T) where {T <: Real}
    coef = clock.coef
    partim = _sclk_from_tdb(ks, clock, et)
    p = primal(partim)
    (p < coef[2] || isnan(p)) && _sclk_range(clock.id, "the epoch $(primal(et))")
    k = 3 * (_sclk_record(coef, 2, p) - 1)
    @inbounds c1, c2, c3 = coef[k + 1], coef[k + 2], coef[k + 3]
    c3 <= 0.0 && _sclk_badrate(clock.id)
    rate = 1.0 / (c3 / clock.tikmsc)
    sclkdp = c1 + rate * (partim - c2)
    primal(sclkdp) > clock.mxtick && _sclk_range(clock.id, "the epoch $(primal(et))")
    return sclkdp
end

"""
    _sce2t(ks, clock, et) -> Float64

Discrete encoded SCLK of the ET `et` (`scet01`): the value of the first or last record
is extended to epochs before the first record, and the result is rounded to whole ticks.
Throws `SpiceError(:valueoutofrange)` if the rounded value is outside the clock range.
"""
@inline function _sce2t(ks::KernelSet, clock::SclkType1, et::Float64)
    coef = clock.coef
    partim = _sclk_from_tdb(ks, clock, et)
    lower = partim < coef[2] ? 1 : _sclk_record(coef, 2, partim)
    k = 3 * (lower - 1)
    @inbounds c1, c2, c3 = coef[k + 1], coef[k + 2], coef[k + 3]
    c3 <= 0.0 && _sclk_badrate(clock.id)
    rate = 1.0 / (c3 / clock.tikmsc)
    sclkdp = _dnint(c1 + rate * (partim - c2))
    (sclkdp < _dnint(coef[1]) || sclkdp > clock.mxtick || isnan(sclkdp)) &&
        _sclk_range(clock.id, "the epoch $et")
    return sclkdp
end

# --- numeric interface ------------------------------------------------------------------

@inline _sclk_num(x::Integer) = Float64(x)
@inline _sclk_num(x::Real) = x

"""
    sct2e([ctx,] sc, sclkdp) -> Real

Return the epoch, in seconds past J2000 TDB, of the encoded spacecraft clock value
`sclkdp` of the clock `sc`, as
`sct2e_c`.

`sc` is the NAIF ID of the spacecraft. `sclkdp` is in ticks, the counts of the last
field of the clock, from the start of the first partition; [`scencd`](@ref) makes it
from a clock string. A partition is a span on which the clock counts without a reset
or a jump. An SCLK kernel for `sc` must be loaded, and a leap-second kernel if the
parallel time system of the clock, the time scale to which the kernel relates the
ticks, is TDT (`SCLK01_TIME_SYSTEM_<-sc> = 2`). The result has the number type of
`sclkdp` (`Float64` for an integer), so a dual number keeps its derivative. The call
does not allocate.

Throws `SpiceError(:kernelvarnotfound)` if no SCLK kernel for `sc` is loaded or its
data are incomplete, `SpiceError(:numpartsunequal)` if its partition start and end
lists differ in length, and `SpiceError(:notsupported)` if the clock is not of type 1.
Throws `SpiceError(:valueoutofrange)` if `sclkdp` is outside the range of the
kernel and `SpiceError(:missingtimeinfo)` if a leap-second kernel is needed and not
loaded.

# Examples
```julia
furnsh("naif0012.tls")
furnsh("mission.tsc")
et = sct2e(-77, scencd(-77, "1/0123456789:123"))
et = scs2e(-77, "1/0123456789:123")      # the same in one call
```

See also [`sce2c`](@ref), [`scs2e`](@ref).

# References

- [NAIF documentation of `sct2e_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sct2e_c.html)
"""
sct2e(ks::KernelSet, sc::Integer, sclkdp::Real) =
    _sct2e(ks, _sclk(ks, sc), _sclk_num(sclkdp))

"""
    sce2c([ctx,] sc, et) -> Real

Return the continuous encoded spacecraft clock value of the clock `sc` at `et`, as
`sce2c_c`.

`et` is in seconds past J2000 TDB. The result is in ticks from the start of the first
partition and can have a fraction, which makes it continuous; [`sce2t`](@ref) gives
the discrete value, rounded to whole ticks. This is
the form that [`ckgp`](@ref) takes. The result has the number type of `et`, and the
call does not allocate.

Throws the errors of [`sct2e`](@ref), with `et` in the place of `sclkdp`.

# References

- [NAIF documentation of `sce2c_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sce2c_c.html)
"""
sce2c(ks::KernelSet, sc::Integer, et::Real) = _sce2c(ks, _sclk(ks, sc), _sclk_num(et))

"""
    sce2t([ctx,] sc, et) -> Float64

Return the discrete encoded spacecraft clock value of the clock `sc` at `et`, in
whole ticks, as
`sce2t_c`.

`et` is in seconds past J2000 TDB. A dual-number `et` is taken at its plain value: the
result is a step function and has no derivative.

Throws the errors of [`sce2c`](@ref).

# References

- [NAIF documentation of `sce2t_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sce2t_c.html)
"""
sce2t(ks::KernelSet, sc::Integer, et::Real) = _sce2t(ks, _sclk(ks, sc), primal(et))

"""
    scpart([ctx,] sc) -> (pstart::Vector{Float64}, pstop::Vector{Float64})

Return the start and stop times of the partitions of the clock `sc`, in ticks, as
`scpart_c`.

The two vectors are new and have one element for each partition; their length takes
the place of the `nparts` result of CSPICE. [`sct2e`](@ref) says what a partition is.

Also throws the SCLK kernel errors of [`sct2e`](@ref).

# References

- [NAIF documentation of `scpart_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/scpart_c.html)
"""
function scpart(ks::KernelSet, sc::Integer)
    clock = _sclk(ks, sc)
    return copy(clock.pstart), copy(clock.pstop)
end

# --- strings ----------------------------------------------------------------------------

"""
    _lparsm(list, delims, nmax) -> Vector{SubString{String}}

Items of `list` separated by any character of `delims`, at most `nmax`, as `lparsm_`:
leading blanks of an item are dropped, two delimiters in a row give a blank item, and a
list that ends with a delimiter has a last blank item. A blank list has one blank item.
A blank is a delimiter only if `delims` holds one; then a run of blanks counts as one
delimiter and joins a non-blank delimiter next to it. Without a blank in `delims` an
item keeps the blanks before the delimiter that ends it, as the item of `lparsm_` does.
"""
function _lparsm(list::AbstractString, delims::AbstractString, nmax::Int)
    s = String(rstrip(list, ' '))
    items = SubString{String}[]
    isempty(s) && return push!(items, SubString(s, 1, 0))
    c = codeunits(s)
    isdelim(k) = Char(c[k]) in delims
    eol = length(c)
    b = 1
    while b <= eol
        # invariant: `items` holds the items that end before position b
        while b <= eol && c[b] == UInt8(' ')
            b += 1
        end
        e = b
        while e <= eol && !isdelim(e)
            e += 1
        end
        stop = e - 1                     # last character of the item
        if e <= eol && c[e] == UInt8(' ')
            while e <= eol && c[e] == UInt8(' ')
                e += 1
            end
            e <= eol && !isdelim(e) && (e -= 1)
        end
        # `stop` is the last byte of the item; a character of several bytes ends there
        # when the text is not ASCII, and a string index is the first byte of a character.
        push!(items, SubString(s, b, stop < b ? stop : thisind(s, stop)))
        length(items) < nmax || return items
        b = e + 1
    end
    Char(c[eol]) in delims && length(items) < nmax && push!(items, SubString(s, 1, 0))
    return items
end

"Ticks of one count of each of the `nfield` fields of `clock` (`CMPTKS` of `sc01`)."
function _sclk_cmptks(clock::SclkType1)
    nfield = Int(clock.nfield)
    cmptks = ones(nfield)
    for i in (nfield - 1):-1:1          # invariant: cmptks[i+1:nfield] are final
        cmptks[i] = cmptks[i + 1] * clock.moduli[i + 1]
    end
    return cmptks
end

"""
    sctiks([ctx,] sc, clkstr) -> Float64

Return the tick count of the clock string `clkstr` of the clock `sc`, as
`sctiks_c`.

`clkstr` has no partition number. Its fields are given most significant first and are
separated by `.`, `:`, `-`, `,` or blanks. Each field counts from its offset, the
value of `SCLK01_OFFSETS` in the SCLK kernel (the value of the field when the count is
zero). A blank field counts as its offset. Only the first 30 characters of a field are
read, as in CSPICE. The result is a count inside one partition,
not an encoded clock value; see [`scencd`](@ref).

Throws `SpiceError(:invalidsclkstring)` if the string is blank, has too many fields,
has a field that is not a number, or has a field below its offset.
Also throws the SCLK kernel errors of [`sct2e`](@ref).

See also [`scfmt`](@ref).

# References

- [NAIF documentation of `sctiks_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sctiks_c.html)
"""
function sctiks(ks::KernelSet, sc::Integer, clkstr::AbstractString)
    _chkfstr("sctiks", "string `clkstr`", clkstr)       # CHKFSTR of sctiks_c
    ticks, msg = _sclk_ticks(_sclk(ks, sc), clkstr)
    isempty(msg) || throw(SpiceError(:invalidsclkstring, msg))
    return ticks
end

"""
    _sclk_ticks(clock, clkstr) -> (ticks, errmsg)

Tick count of the clock string `clkstr` (no partition number) of `clock`, and `""`; or
0 and the message of `scps01_` for a string that is blank, has too many fields, has a
field that is not a number or a field below its offset. A field is cut to 30 bytes, and
only the blank counts as white space, as in `scps01_`. [`sctiks`](@ref) and the time
values of EK queries share it.
"""
function _sclk_ticks(clock::SclkType1, clkstr::AbstractString)
    isempty(strip(clkstr, ' ')) &&
        return 0.0, "Non partition part of the input clock string is blank."
    cmptks = _sclk_cmptks(clock)
    items = _lparsm(clkstr, _SCLK_DELIMS, _SCLK_MXNFLD)
    n = length(items)
    n > clock.nfield && return 0.0, _repmc(
            _repmc(
                _repmc(
                    "Input clock string # has # fields; maximum for this spacecraft clock is #.",
                    clkstr
                ), n
            ), Int(clock.nfield)
        )
    ticks = 0.0
    for (i, field) in enumerate(items)      # invariant: ticks = sum over the fields 1:i-1
        item = _cut_bytes(field, 30)         # `scps01_` reads items of 30 characters
        value = clock.offsets[i]
        if !isempty(strip(item, ' '))
            x = isascii(item) ? _try_nparsd(strip(item, ' ')) : nothing
            x === nothing && return 0.0, _repmc(
                    _repmc("Could not parse SCLK component # from # as a number.", item),
                    clkstr
                )
            value = x
        end
        value -= clock.offsets[i]
        _dnint(value) < 0 && return 0.0, _repmc(
                _repmc(
                    _repmc("Component number #, # in the SCLK string  # is invalid.", i),
                    item
                ), clkstr
            )
        ticks += value * cmptks[i]
    end
    return ticks, ""
end

"""
    scfmt([ctx,] sc, ticks) -> String

Return the clock string of the tick count `ticks` of the clock `sc`, without a
partition number, as
`scfmt_c`.

`ticks` is rounded to a whole number and split into the fields of the clock. Each
field is written with leading zeros to the width of its largest value (a field of
1e14 or more with its first 14 significant digits only, as in CSPICE), with the
output delimiter of the clock between the fields: the character `.`, `:`, `-`, `,` or
blank that `SCLK01_OUTPUT_DELIM_<-sc>` selects with 1 to 5. This is the inverse of
[`sctiks`](@ref).

Throws `SpiceError(:valueoutofrange)` if `ticks` is negative.
Also throws the SCLK kernel errors of [`sct2e`](@ref).

# References

- [NAIF documentation of `scfmt_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/scfmt_c.html)
"""
function scfmt(ks::KernelSet, sc::Integer, ticks::Real)
    clock = _sclk(ks, sc)
    nfield = Int(clock.nfield)
    cmptks = _sclk_cmptks(clock)
    rem = _dnint(Float64(ticks))
    (rem < 0 || !isfinite(rem)) && throw(
        SpiceError(
            :valueoutofrange,
            "negative or non-finite value for SCLK ticks: $ticks"
        )
    )
    1 <= clock.delim <= 5 || throw(
        SpiceError(
            :valueoutofrange,
            "output delimiter code $(clock.delim) of the clock $sc is not in 1:5"
        )
    )
    delim = _SCLK_OUTPUT_DELIM[clock.delim]
    io = IOBuffer()
    for i in 1:nfield       # invariant: fields 1:i-1 are written, rem is what they leave
        value = if i < nfield
            v = trunc(rem / cmptks[i]) + clock.offsets[i]
            # `DMOD`: where rem / cmptks rounds up to a whole number, this is negative,
            # and its magnitude is printed, as by `scfm01`; the exact remainder is not
            rem = f2c_dmod(rem, cmptks[i])
            v
        else
            rem + clock.offsets[i]
        end
        width = trunc(Int, log10(clock.moduli[i] + clock.offsets[i] - 1.0 + 0.5)) + 1
        i > 1 && print(io, delim)
        print(io, lpad(_sclk_field_digits(value), width, '0'))
    end
    return String(take!(io))
end

"""
    _sclk_field_digits(value) -> String

The digits that `scfm01` writes for a field: the characters before the decimal point of
`dpstrf_(value, 30, 'F')` without its sign. `dpstrf` keeps at most 14 significant
digits, so a field of 1e14 or more is rounded to them and ends in zeros.
"""
function _sclk_field_digits(value::Float64)
    abs(value) < 1.0e14 && return string(trunc(Int, abs(value)))    # all digits are exact
    s = String(_inline(_dpstrf(_ByteBuf(), value, 30)))
    dot = something(findfirst(==('.'), s), ncodeunits(s) + 1)
    return s[2:(dot - 1)]
end

"""
    _sclk_partitions(clock, rounded = true) -> (pstart, pstop, ptotls)

Partition bounds of `clock` and the rounded running totals of the partition lengths.
`rounded` rounds the bounds to whole ticks first, as `scencd_` and `scdecd_` do;
`scpars_` takes them as they are.
"""
function _sclk_partitions(clock::SclkType1, rounded::Bool = true)
    pstart = rounded ? _dnint.(clock.pstart) : copy(clock.pstart)
    pstop = rounded ? _dnint.(clock.pstop) : copy(clock.pstop)
    ptotls = similar(pstart)
    acc = 0.0
    for i in eachindex(pstart)      # invariant: acc = rounded length of partitions 1:i-1
        acc = _dnint(acc + pstop[i] - pstart[i])
        ptotls[i] = acc
    end
    return pstart, pstop, ptotls
end

"""
    scencd([ctx,] sc, sclkch) -> Float64

Return the encoded spacecraft clock value of the clock string `sclkch` of the clock
`sc`, as
`scencd_c`.

`sclkch` has the form `"partition/fields"`, such as `"1/0123456789:123"`; the
partition number is optional. The result is in ticks from the start of the first
partition. Without a partition number, the first partition that holds the count is
taken.

Throws `SpiceError`:
- `:emptystring` if `sclkch` is empty;
- `:badpartnumber` if the partition number cannot be read or is not the number of a
  partition;
- `:notinpart` if the count is not inside the given partition;
- `:nopartition` if no partition holds the count;
- `:invalidsclkstring` if the fields cannot be read (see [`sctiks`](@ref)).

Also throws the SCLK kernel errors of [`sct2e`](@ref).

See also [`scdecd`](@ref), [`scs2e`](@ref).

# References

- [NAIF documentation of `scencd_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/scencd_c.html)
"""
function scencd(ks::KernelSet, sc::Integer, sclkch::AbstractString)
    _chkfstr("scencd", "string `sclkch`", sclkch)       # CHKFSTR of scencd_c
    sclkdp, kind, msg = _sclk_encode(_sclk(ks, sc), sc, sclkch, true)
    isempty(msg) || throw(SpiceError(kind, msg))
    return sclkdp
end

"""
    _sclk_encode(clock, sc, sclkch, scencd) -> (sclkdp, kind, errmsg)

Encoded clock value of the clock string `sclkch` (`"partition/fields"`, the partition
optional) of `clock`, the clock of the spacecraft `sc`, with `:none` and `""`; or 0, the error kind of `scencd_` and the
message of `scpars_`. The fields are read before the partition number, which is cut
toward zero (`nparsi_`). With `scencd` true the count and the partition bounds are
rounded to whole ticks, as in `scencd_`; with `scencd` false they are not, and a blank
string is an error of its own, as in `scpars_`, which the time values of EK queries use.
"""
function _sclk_encode(clock::SclkType1, sc::Integer, sclkch::AbstractString, scencd::Bool)
    pstart, pstop, ptotls = _sclk_partitions(clock, scencd)
    nparts = length(pstart)
    !scencd && isempty(strip(sclkch, ' ')) &&
        return 0.0, :invalidsclkstring, "Input spacecraft clock string is blank."
    pos = something(findfirst('/', sclkch), 0)
    badpart = "Unable to parse the partition number from SCLK string #."
    kind = :badpartnumber    # the kinds are named by `kind = ...` for the registry test
    pos == 1 && return 0.0, kind, _repmc(badpart, sclkch)
    fields = SubString(sclkch, pos == 0 ? 1 : nextind(sclkch, pos))
    ticks, msg = _sclk_ticks(clock, fields)
    isempty(msg) || return 0.0, :invalidsclkstring, msg
    scencd && (ticks = _dnint(ticks))
    part = 1
    if pos > 1
        part = _try_nparsi(SubString(sclkch, 1, prevind(sclkch, pos)))
        part === nothing && return 0.0, kind, _repmc(badpart, sclkch)
        (part <= 0 || part > nparts) && return 0.0, kind, _repmc(
                _repmc(
                    _repmc(
                        "Partition number # taken from SCLK string # is not in acceptable " *
                        "range 1 to #.", part
                    ), sclkch
                ), nparts
            )
        kind = :notinpart
        (ticks < pstart[part] || ticks > pstop[part]) && return 0.0, kind, _repmc(
                _repmc(
                    "SCLK count from # does not fall in the boundaries of partition number #.",
                    fields
                ), part
            )
    else
        while part <= nparts && (ticks < pstart[part] || ticks > pstop[part])
            part += 1                # invariant: partitions 1:part-1 do not hold `ticks`
        end
        kind = :nopartition
        part > nparts && return 0.0, kind, _repmc(
                _repmc(
                    "SCLK count # does not fall in the boundaries of any of the partitions " *
                    "for spacecraft #.", fields
                ), Int(sc)
            )
    end
    sclkdp = part > 1 ? ticks - pstart[part] + ptotls[part - 1] : ticks - pstart[part]
    return sclkdp, :none, ""
end

"""
    scdecd([ctx,] sc, sclkdp) -> String

Return the clock string `"partition/fields"` of the encoded spacecraft clock value
`sclkdp` of the clock `sc`, as
`scdecd_c`.

`sclkdp` is in ticks from the start of the first partition and is rounded to whole
ticks. This is the inverse of [`scencd`](@ref).

Throws `SpiceError(:valueoutofrange)` if the ticks are not inside a partition.
Also throws the SCLK kernel errors of [`sct2e`](@ref).

# References

- [NAIF documentation of `scdecd_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/scdecd_c.html)
"""
function scdecd(ks::KernelSet, sc::Integer, sclkdp::Real)
    clock = _sclk(ks, sc)
    pstart, _, ptotls = _sclk_partitions(clock)
    nparts = length(pstart)
    ticks = _dnint(Float64(sclkdp))
    part = ticks == ptotls[nparts] ? nparts : searchsortedlast(ptotls, ticks) + 1
    (ticks < 0 || part > nparts || isnan(ticks)) && throw(
        SpiceError(
            :valueoutofrange,
            "value for ticks, $ticks, is not inside any partition of the spacecraft clock $sc"
        )
    )
    ticks = part == 1 ? ticks + pstart[1] : ticks + pstart[part] - ptotls[part - 1]
    return string(part, '/', scfmt(ks, sc, ticks))
end

"""
    sce2s([ctx,] sc, et) -> String

Return the clock string `"partition/fields"` of the clock `sc` at `et` (seconds past
J2000 TDB), as
`sce2s_c`.

The result is [`scdecd`](@ref) of [`sce2t`](@ref), and the errors are theirs.

See also [`scs2e`](@ref).

# References

- [NAIF documentation of `sce2s_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sce2s_c.html)
"""
sce2s(ks::KernelSet, sc::Integer, et::Real) = scdecd(ks, sc, sce2t(ks, sc, et))

"""
    scs2e([ctx,] sc, sclkch) -> Float64

Return the epoch, in seconds past J2000 TDB, of the clock string `sclkch` of the
clock `sc`, as
`scs2e_c`.

The result is [`sct2e`](@ref) of [`scencd`](@ref), and the errors are theirs.

See also [`sce2s`](@ref).

# References

- [NAIF documentation of `scs2e_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/scs2e_c.html)
"""
function scs2e(ks::KernelSet, sc::Integer, sclkch::AbstractString)
    _chkfstr("scs2e", "string `sclkch`", sclkch)       # CHKFSTR of scs2e_c
    return sct2e(ks, sc, scencd(ks, sc, sclkch))
end

# --- forms on a context and on the default context ---------------------------------------

for f in (:sct2e, :sce2c, :sce2t, :sce2s, :scdecd, :scfmt)
    @eval $f(ctx::SpiceContext, sc::Integer, x::Real) = $f(kernelset(ctx), sc, x)
    @eval $f(sc::Integer, x::Real) = $f(default_context(), sc, x)
end
for f in (:scs2e, :scencd, :sctiks)
    @eval $f(ctx::SpiceContext, sc::Integer, s::AbstractString) = $f(kernelset(ctx), sc, s)
    @eval $f(sc::Integer, s::AbstractString) = $f(default_context(), sc, s)
end
for f in (:scpart, :sctype)
    @eval $f(ctx::SpiceContext, sc::Integer) = $f(kernelset(ctx), sc)
    @eval $f(sc::Integer) = $f(default_context(), sc)
end
