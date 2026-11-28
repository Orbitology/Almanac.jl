# Kernel-pool, name, frame-name and time interface on a context.

"Return the value of pool variable `name`; throws `SpiceError(:kernelvarnotfound)` if it is absent."
function pool_value(ks::KernelSet, name::AbstractString)
    v = get(ks.pool, strip(name), nothing)          # a `SubString` key: no copy of the name
    v === nothing && throw(SpiceError(:kernelvarnotfound, "kernel-pool variable $name is not loaded"))
    return v
end

"""
    _pool_get(ks, name, who, V, what) -> V

The stored values of the pool variable `name`, not a copy, after the checks of the C
wrapper `who`: `SpiceError(:emptystring)`, `:kernelvarnotfound`, and `:typemismatch`
unless the values are a `V` (`what` names that type in the message).
"""
function _pool_get(
        ks::KernelSet, name::AbstractString, who::AbstractString, ::Type{V},
        what::AbstractString
    ) where {V}
    _chkfstr(who, "variable name", name)          # CHKFSTR of the C wrapper
    v = pool_value(ks, name)
    v isa V || throw(SpiceError(:typemismatch, "kernel-pool variable $name is not $what"))
    return v::V
end

"""
    gdpool([ctx,] name) -> Vector{Float64}

Return the numeric values of the kernel-pool variable `name`, as
`gdpool_c`.

The result is a new vector with every value of the variable. Blanks around `name` are
ignored; letter case is significant.

!!! warning "Difference from CSPICE"
    The call has no `start` and `room` arguments and no `found` flag: an absent
    variable throws `SpiceError(:kernelvarnotfound)`, where CSPICE sets `found` to false
    and SPICE.jl returns `nothing`. Use [`expool`](@ref) or [`dtpool`](@ref) to test
    first.

Throws `SpiceError(:emptystring)` if `name` is the empty string,
`SpiceError(:kernelvarnotfound)` if the variable is not in the pool and
`SpiceError(:typemismatch)` if it holds character values.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pdpool(ctx, "BODY399_RADII", [6378.1366, 6378.1366, 6356.7519]);

julia> gdpool(ctx, "BODY399_RADII")
3-element Vector{Float64}:
 6378.1366
 6378.1366
 6356.7519
```

See also [`gipool`](@ref), [`gcpool`](@ref), [`bodvrd`](@ref), [`pdpool`](@ref).

# References

- [NAIF documentation of `gdpool_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gdpool_c.html)
"""
function gdpool(ks::KernelSet, name::AbstractString)
    return copy(_pool_get(ks, name, "gdpool", Vector{Float64}, "numeric"))
end

"""
    gipool([ctx,] name) -> Vector{Int}

Return the values of the numeric kernel-pool variable `name` rounded to the nearest
integers, a tie away from zero, as
`gipool_c`.

The result holds every value.

!!! warning "Difference from CSPICE"
    There is no `found` flag: an absent variable throws, as in [`gdpool`](@ref).

Throws `SpiceError`:
- `:emptystring` if `name` is the empty string;
- `:kernelvarnotfound` if the variable is not in the pool;
- `:typemismatch` if it holds character values;
- `:intoutofrange` if a value is outside the 32-bit integers or NaN.

# References

- [NAIF documentation of `gipool_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gipool_c.html)
"""
function gipool(ks::KernelSet, name::AbstractString)
    v = _pool_get(ks, name, "gipool", Vector{Float64}, "numeric")
    for (i, x) in enumerate(v)
        # `gipool` of CSPICE: DBLE(INTMIN()) <= x <= DBLE(INTMAX()), then NINT
        _gipool_value(x) === nothing && _gipool_range_error(name, i, x)
    end
    return Int[_nint(x) for x in v]
end

"""
    gcpool([ctx,] name) -> Vector{String}

Return the character values of the kernel-pool variable `name`, as
`gcpool_c`.

The result is a new vector with every value.

!!! warning "Difference from CSPICE"
    There is no `found` flag: an absent variable throws, as in [`gdpool`](@ref).

Throws `SpiceError(:emptystring)` if `name` is the empty string,
`SpiceError(:kernelvarnotfound)` if the variable is not in the pool and
`SpiceError(:typemismatch)` if it holds numeric values.

# References

- [NAIF documentation of `gcpool_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gcpool_c.html)
"""
function gcpool(ks::KernelSet, name::AbstractString)
    return copy(_pool_get(ks, name, "gcpool", Vector{String}, "character"))
end

"""
    _pool_key(prefix, mid, s1, s2="", upper=false) -> String

Return the pool name `prefix * mid * s1 * s2`, built in one buffer: `mid` is an integer
(its decimal digits) or a string (in upper case if `upper`).
"""
function _pool_key(
        prefix::String, mid::Union{Integer, AbstractString}, s1::AbstractString,
        s2::AbstractString = "", upper::Bool = false
    )
    if mid isa Integer
        n = abs(Int(mid))
        nmid = ndigits(n) + (mid < 0)
    else
        nmid = ncodeunits(mid)
    end
    np, n1, n2 = ncodeunits(prefix), ncodeunits(s1), ncodeunits(s2)
    buf = Base.StringVector(np + nmid + n1 + n2)
    copyto!(buf, 1, codeunits(prefix), 1, np)
    if mid isa Integer
        mid < 0 && (buf[np + 1] = UInt8('-'))
        # Loop invariant: the `k` last digits of `n` are written, the last at `np + nmid`.
        for k in 0:(ndigits(n) - 1)
            n, d = divrem(n, 10)
            buf[np + nmid - k] = UInt8('0') + d
        end
    else
        for k in 1:nmid
            c = codeunit(mid, k)
            buf[np + k] = upper ? _upper(c) : c
        end
    end
    copyto!(buf, np + nmid + 1, codeunits(s1), 1, n1)
    copyto!(buf, np + nmid + n1 + 1, codeunits(s2), 1, n2)
    return String(buf)
end

_body_key(code::Integer, item::AbstractString) = _pool_key("BODY", code, "_", strip(item))

# `bodvrd_c` signals SPICE(NOTRANSLATION) for a name without a code, where the state
# calls signal SPICE(IDCODENOTFOUND).
_bodvrd_code(::KernelSet, body::Integer) = Int32(body)
function _bodvrd_code(ks::KernelSet, body::AbstractString)
    c = _body_code_lookup(ks.name2code, body)
    c === nothing && throw(SpiceError(:notranslation, "bodvrd: the body name '$body' has no NAIF ID code"))
    return c
end

"""
    bodvrd([ctx,] body, item) -> Vector{Float64}

Return the values of the kernel-pool variable `BODY<id>_<item>` for `body`, as
`bodvrd_c`.

`body` is a body name, a string that holds an integer (`"399"`) or an integer NAIF ID
code. `item` is case sensitive, as in CSPICE; blanks around it are ignored. The values have the units of the kernel: km for
`"RADII"`, km^3/s^2 for `"GM"` in the NAIF kernels. There is no `maxn` argument.

Throws `SpiceError`:
- `:emptystring` if `body` or `item` is the empty string;
- `:notranslation` if the body name is not known;
- `:kernelvarnotfound` if the variable is not in the pool;
- `:typemismatch` if it holds character values.

# Examples
```julia
furnsh("pck00011.tpc")
a, b, c = bodvrd("EARTH", "RADII")    # km
```

See also [`bodvcd`](@ref), [`gdpool`](@ref).

# References

- [NAIF documentation of `bodvrd_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html)
"""
function bodvrd(ks::KernelSet, body::NameOrID, item::AbstractString)
    body isa AbstractString && _chkfstr("bodvrd", "body name", body)   # CHKFSTR of bodvrd_c,
    _chkfstr("bodvrd", "item", item)                                     # body then item
    return gdpool(ks, _body_key(_bodvrd_code(ks, body), item))
end

const _RADII_KEYS = String[string("BODY", i, "_RADII") for i in 1:999]

@noinline function _bad_radii(body, v, kind::Symbol = :invalidcount)
    v === nothing && throw(
        SpiceError(
            :kernelvarnotfound,
            "the variable BODY$(body)_RADII could not be found in the kernel pool"
        )
    )
    (v isa Vector{Float64} || v isa NTuple{3, Float64}) || throw(
        SpiceError(
            :typemismatch,
            "the variable BODY$(body)_RADII is not numeric"
        )
    )
    length(v) == 3 || throw(
        SpiceError(
            length(v) > 3 ? :arraytoosmall : kind,
            "the variable BODY$(body)_RADII has $(length(v)) values; 3 are required"
        )
    )
    throw(SpiceError(:badaxislength, "a radius of the body $body is not positive: $(collect(v))"))
end

"""
    _radii3(ks, body, kind) -> NTuple{3,Float64}

`zzbodvcd` for `RADII`: three values, with no test of their signs. Throws `SpiceError`
(`:kernelvarnotfound`, `:typemismatch`, `:arraytoosmall`, and `kind` for fewer than 3
values).
"""
@inline function _radii3(ks::KernelSet, body::Integer, kind::Symbol)
    v = get(
        ks.pool, 1 <= body <= 999 ? (@inbounds _RADII_KEYS[body]) :
            string("BODY", body, "_RADII"), nothing
    )
    v isa Vector{Float64} || _bad_radii(body, v, kind)
    length(v) == 3 || _bad_radii(body, v, kind)
    return @inbounds (v[1], v[2], v[3])
end

"""
    _positive_radii(ks, body) -> NTuple{3,Float64}

`zzgftreb`: the three radii of `body` from the kernel pool. Throws `SpiceError`
(`:kernelvarnotfound`, `:typemismatch`, `:invalidcount`, `:badaxislength`) unless
they are three positive values. No allocation for an ID in 1:999.
"""
@inline function _positive_radii(ks::KernelSet, body::Integer)
    a, b, c = _radii3(ks, body, :invalidcount)
    (a > 0.0 && b > 0.0 && c > 0.0) || _bad_radii(body, (a, b, c))
    return a, b, c
end

"""
    bods2c([ctx,] name) -> Union{Int,Nothing}

Return the NAIF ID code of the body `name`, or `nothing` if the name is not known, as
`bods2c_c`.

`name` is a body name or a string that holds an integer, such as `"399"`. Letter case
and blanks around the name are ignored. Names from the kernel pool and from
[`boddef`](@ref) have priority over the built-in names; only the first 36 characters
of a name count, as in CSPICE. `nothing` takes the place of the `found` flag of CSPICE.

Throws `SpiceError(:emptystring)` if `name` is the empty string.

# Examples
```jldoctest
julia> bods2c("Earth"), bods2c("-77"), bods2c("NO SUCH BODY")
(399, -77, nothing)
```

See also [`bodn2c`](@ref), [`bodc2n`](@ref).

# References

- [NAIF documentation of `bods2c_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bods2c_c.html)
"""
function bods2c(ks::KernelSet, name::AbstractString)
    _chkfstr("bods2c", "body name", name)         # CHKFSTR of bods2c_c
    c = _body_code_lookup(ks.name2code, name)
    return c === nothing ? nothing : Int(c)
end

@noinline _no_body_id(what, name) = throw(
    SpiceError(
        :idcodenotfound,
        "the $what, '$name', is not a recognized name for an ephemeris object"
    )
)

"The body code of `body` (a code, or a name or an integer string for `bods2c`); `what` names it in the error."
@inline _body_id(::KernelSet, id::Integer, ::AbstractString) = Int32(id)
function _body_id(ks::KernelSet, name::AbstractString, what::AbstractString)
    # The lookup of `bods2c` without the CHKFSTR of `bods2c_c`: an empty name that the C
    # function of the call does not check is a blank name for SPICELIB, which no body has.
    code = _body_code_lookup(ks.name2code, name)
    code === nothing && _no_body_id(what, name)
    return Int32(code)
end

"""
    bodn2c([ctx,] name) -> Union{Int,Nothing}

Return the NAIF ID code of the body name `name`, or `nothing` if the name is not
known, as
`bodn2c_c`.

Different from [`bods2c`](@ref), a string that holds an integer is not a name.
`nothing` takes the place of the `found` flag of CSPICE.

Throws `SpiceError(:emptystring)` if `name` is the empty string.

# Examples
```jldoctest
julia> bodn2c("MOON"), bodn2c("301")
(301, nothing)
```

# References

- [NAIF documentation of `bodn2c_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodn2c_c.html)
"""
function bodn2c(ks::KernelSet, name::AbstractString)
    _chkfstr("bodn2c", "body name", name)         # CHKFSTR of bodn2c_c
    c = _body_name_code(ks.name2code, name)
    return c === nothing ? nothing : Int(c)
end

"""
    bodc2n([ctx,] code) -> Union{String,Nothing}

Return the name of the body with the NAIF ID `code`, or `nothing` if the code has no
name, as
`bodc2n_c`.

A name from the kernel pool or from [`boddef`](@ref) has priority over the built-in
name.

# Examples
```jldoctest
julia> bodc2n(399), bodc2n(123456789)
("EARTH", nothing)
```

See also [`bodn2c`](@ref).

# References

- [NAIF documentation of `bodc2n_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodc2n_c.html)
"""
bodc2n(ks::KernelSet, code::Integer) = body_name(ks, code)

"""
    bodc2s([ctx,] code) -> String

Return the name of the body with the NAIF ID `code`, or the decimal text of `code` if
the code has no name, as
`bodc2s_c`.

The name is the one of [`bodc2n`](@ref).

# Examples
```jldoctest
julia> bodc2s(399), bodc2s(123456789)
("EARTH", "123456789")
```

See also [`bods2c`](@ref).

# References

- [NAIF documentation of `bodc2s_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodc2s_c.html)
"""
function bodc2s(ks::KernelSet, code::Integer)
    name = body_name(ks, code)
    return name === nothing ? string(code) : name
end

"""
    scn2id([ctx,] clknam) -> Union{Int,Nothing}

Return the NAIF ID code of the spacecraft clock `clknam`, or `nothing` if the name is
not known, as
`scn2id`
of SPICELIB. CSPICE has no C wrapper for it.

The name is usually the body name of the spacecraft with the word `SCLK` after it, as
in `"MGS SCLK"`. Only the first 32 characters count, as in SPICELIB; letter case is
ignored, the last `SCLK` among them is blanked (it need not be a separate word), and the
rest is looked up with [`bodn2c`](@ref), so a bare body name gives the same code.

# Examples
```jldoctest
julia> scn2id("CASSINI SCLK"), scn2id("cassini"), scn2id("NO SUCH CLOCK")
(-82, -82, nothing)
```

See also [`scid2n`](@ref).

# References

- [NAIF documentation of `scn2id`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/scn2id.html)
"""
function scn2id(ks::KernelSet, clknam::AbstractString)
    # `ucase` into a 32-character buffer, then the last "SCLK" becomes blanks.
    s = codeunits(String(clknam))
    b = UInt8[_upper(c) for c in view(s, 1:min(length(s), 32))]
    for i in (length(b) - 3):-1:1
        if b[i] == UInt8('S') && b[i + 1] == UInt8('C') && b[i + 2] == UInt8('L') &&
                b[i + 3] == UInt8('K')
            b[i:(i + 3)] .= UInt8(' ')
            break
        end
    end
    name = String(b)
    isempty(strip(name)) && return nothing
    return bodn2c(ks, name)
end

"""
    scid2n([ctx,] clkid) -> Union{String,Nothing}

Return the name of the spacecraft clock with the NAIF ID `clkid`, or `nothing` if the
code has no body name, as
`scid2n`
of SPICELIB. CSPICE has no C wrapper for it.

The name is the body name of [`bodc2n`](@ref) with `" SCLK"` after it.

# Examples
```jldoctest
julia> scid2n(-82), scid2n(123456789)
("CASSINI SCLK", nothing)
```

See also [`scn2id`](@ref).

# References

- [NAIF documentation of `scid2n`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/scid2n.html)
"""
function scid2n(ks::KernelSet, clkid::Integer)
    name = body_name(ks, clkid)
    return name === nothing ? nothing : string(rstrip(name), " SCLK")
end

"""
    namfrm([ctx,] name) -> Int

Return the ID of the frame `name`, or 0 if the frame is not known, as
`namfrm_c`.

Built-in frames and frames defined in the kernel pool are known. Throws
`SpiceError(:emptystring)` if `name` is empty.

# Examples
```jldoctest
julia> namfrm("J2000"), namfrm("eclipj2000"), namfrm("NO_SUCH_FRAME")
(1, 17, 0)
```

See also [`frmnam`](@ref), [`frinfo`](@ref).

# References

- [NAIF documentation of `namfrm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/namfrm_c.html)
"""
function namfrm(ks::KernelSet, name::AbstractString)
    _chkfstr("namfrm", "frame name", name)
    return _namfrm(ks, name)
end

"Return the frame ID of `name`, or 0, as `namfrm_` of SPICELIB (without the empty-string check of `namfrm_c`)."
function _namfrm(ks::KernelSet, name::AbstractString)
    # The lookups of `frame_id`, without its error for an unknown name.
    id = name_lookup(BUILTIN_FRAMENAME2ID, name)
    id === nothing && (id = name_lookup(ks.framename2id, name))
    return id === nothing ? 0 : Int(id)
end

"The frame ID of `ref`: an ID as it is, a name through `_namfrm` (0 if it is not known)."
@inline _frame_code(::KernelSet, ref::Integer) = Int32(ref)
_frame_code(ks::KernelSet, ref::AbstractString) = Int32(_namfrm(ks, ref))

"""
    frmnam([ctx,] id) -> String

Return the name of the frame with ID `id`, or `""` if the frame is not known, as
`frmnam_c`.

# Examples
```jldoctest
julia> frmnam(17), frmnam(-123456)
("ECLIPJ2000", "")
```

See also [`namfrm`](@ref).

# References

- [NAIF documentation of `frmnam_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/frmnam_c.html)
"""
function frmnam(ks::KernelSet, id::Integer)
    name = frame_name(ks, id)
    return name === nothing ? "" : name
end

# Context and default-context forms of the functions above.
for f in (:gdpool, :gipool, :gcpool, :bods2c, :bodn2c, :scn2id, :namfrm)
    @eval begin
        $f(ctx::SpiceContext, name::AbstractString) = $f(kernelset(ctx), name)
        $f(name::AbstractString) = $f(default_context(), name)
    end
end
for f in (:bodc2n, :bodc2s, :scid2n, :frmnam)
    @eval begin
        $f(ctx::SpiceContext, code::Integer) = $f(kernelset(ctx), code)
        $f(code::Integer) = $f(default_context(), code)
    end
end
bodvrd(ctx::SpiceContext, body::NameOrID, item::AbstractString) =
    bodvrd(kernelset(ctx), body, item)
bodvrd(body::NameOrID, item::AbstractString) = bodvrd(default_context(), body, item)
