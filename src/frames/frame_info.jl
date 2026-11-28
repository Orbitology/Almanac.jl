# Frame information (`frinfo`, `cidfrm`, `cnmfrm`, `ccifrm` of `framex.c`), the lists of
# frames (`bltfrm`, `kplfrm`), fixed-offset frame data (`tkfram`), body orientation
# (`tipbod`, `tisbod`), the contents of a binary PCK file (`pckfrm`, `pckcov`) and body
# constants (`bodfnd`, `bodvcd`). Every call that reads kernel data has a form on a
# `KernelSet`, on a `SpiceContext` and on the default context.

"""
    frinfo([ctx,] frcode) -> Union{Nothing,Tuple{Int,Int,Int}}

Return `(center, class, classid)` of the frame with ID `frcode`, as
`frinfo_c`, or
`nothing` if the frame is not known.

`center` is the NAIF ID of the body at the frame origin. `class` is 1 for inertial, 2
for PCK, 3 for CK, 4 for fixed offset, 5 for dynamic and 6 for switch frames;
`classid` is the ID of the frame inside its class. A built-in frame has priority over
a frame of the kernel pool. An ID in 13000:13999 without a definition is a PCK frame
of the Earth with the class ID `frcode - 10000`, as in CSPICE. `nothing` takes the
place of the `found` flag.

# Examples
```jldoctest
julia> frinfo(namfrm("IAU_EARTH")), frinfo(-999)
((399, 2, 399), nothing)
```

See also [`namfrm`](@ref), [`ccifrm`](@ref).

# References

- [NAIF documentation of `frinfo_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/frinfo_c.html)
"""
function frinfo(ks::KernelSet, frcode::Integer)
    f = _frame_def_or_nothing(ks, frcode)
    return f === nothing ? nothing : (Int(f.center), Int(f.class), Int(f.classid))
end

"""
    _frame_of_object(ks, value) -> Union{Nothing,Tuple{Int,String}}

Return the frame ID and frame name for the value of an `OBJECT_<body>_FRAME` variable,
as `cidfrm_` and `cnmfrm_` read it: a frame name (built-in, or with a `FRAME_<name>`
variable) or a frame ID (built-in, or with a `FRAME_<id>_NAME` variable). `nothing` if
the value does not identify a frame. A name is taken as `gcpool_` returns it into 32
characters and compared exactly (`zzhscchk_`): `iau_mars` or ` IAU_MARS` is not the
built-in frame `IAU_MARS`.
"""
function _frame_of_object(ks::KernelSet, value::PoolValue)
    isempty(value) && return nothing
    if value isa Vector{String}
        raw = value[1]
        ncodeunits(raw) > 32 && (raw = String(codeunits(raw)[1:32]))
        name = String(rstrip(raw, ' '))
        # The built-in names are in the form of `normalize_name`, so the exact match
        # is a match of the normalized name by a name that is already in that form.
        if normalize_name(name) == name
            code = get(BUILTIN_FRAMENAME2ID, name, nothing)
            code === nothing || return Int(code), name
        end
        # `gipool_(FRAME_<name>, 1, 8, ...)`, then the test of one value
        key = "FRAME_" * name
        fcode = _gipool_first(get(ks.pool, key, nothing), key; room = 8, single = true)
        return fcode === nothing ? nothing : (Int(fcode), name)
    end
    code = Int(_gipool_first(value, "OBJECT_<body>_FRAME")::Int32)
    f = get(BUILTIN_FRAMES, Int32(code), nothing)
    f === nothing || return code, f.name
    v = get(ks.pool, "FRAME_$(code)_NAME", nothing)
    (v isa Vector{String} && length(v) == 1) || return nothing
    return code, String(strip(v[1]))
end

"""
The built-in frames of `zzfdat_` as `(ID, name, center)`, in the order of the centers
and, for equal centers, of the table (the order vector that `orderi_` gives).
"""
const _FRAME_CENTER_ORDER = let all = Tuple{Int, String, Int}[]
    for (i, d) in enumerate(INERTIAL_FRAME_DEFS)
        push!(all, (i, String(d[1]), 0))
    end
    for f in BUILTIN_NONINERTIAL_FRAMES
        push!(all, (Int(f.id), f.name, Int(f.center)))
    end
    sort!(all; by = x -> x[3], alg = MergeSort)
end

"""
    cidfrm([ctx,] cent) -> Union{Nothing,Tuple{Int,String}}

Return `(frcode, frname)` of the frame that is associated with the body with ID
`cent`, as
`cidfrm_c`, or `nothing` if there is none.

The kernel-pool variable `OBJECT_<cent>_FRAME` (or `OBJECT_<name>_FRAME`) has
priority. Without it the answer is `IAU_EARTH` for the Earth, and otherwise a built-in
frame whose center is `cent`. For a center with several built-in frames the answer is
the frame that CSPICE gives. `nothing` takes the place of the `found` flag.

# Examples
```jldoctest
julia> cidfrm(499)
(10014, "IAU_MARS")
```

See also [`cnmfrm`](@ref).

# References

- [NAIF documentation of `cidfrm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/cidfrm_c.html)
"""
function cidfrm(ks::KernelSet, cent::Integer)
    v = get(ks.pool, _pool_key("OBJECT_", cent, "_FRAME"), nothing)
    if v === nothing
        name = body_name(ks, cent)
        name === nothing || (v = get(ks.pool, _pool_key("OBJECT_", name, "_FRAME", "", true), nothing))
    end
    v === nothing || return _frame_of_object(ks, v)
    cent == 399 && return 10013, "IAU_EARTH"
    # The binary search of `bschoi_` over the built-in frames in the order of their
    # centers. For a center with several frames (all inertial frames have the center 0)
    # the answer is the frame that this search meets first, as in CSPICE.
    left, right = 1, length(_FRAME_CENTER_ORDER)
    while left <= right     # invariant: a frame with the center `cent` is in left:right, if any
        i = (left + right) ÷ 2
        id, name, center = _FRAME_CENTER_ORDER[i]
        cent == center && return id, name
        cent < center ? (right = i - 1) : (left = i + 1)
    end
    return nothing
end

"""
    cnmfrm([ctx,] cname) -> Union{Nothing,Tuple{Int,String}}

Return `(frcode, frname)` of the frame that is associated with the body called
`cname`, as
`cnmfrm_c`, or `nothing` if there is none.

The kernel-pool variable `OBJECT_<cname>_FRAME` (or `OBJECT_<id>_FRAME`) has priority.
Without it the answer is the built-in frame `IAU_<CNAME>`, if it exists. `nothing`
takes the place of the `found` flag. Throws `SpiceError(:emptystring)` if `cname` is
empty.

# Examples
```jldoctest
julia> cnmfrm("Mars"), cnmfrm("NO SUCH BODY")
((10014, "IAU_MARS"), nothing)
```

See also [`cidfrm`](@ref).

# References

- [NAIF documentation of `cnmfrm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/cnmfrm_c.html)
"""
function cnmfrm(ks::KernelSet, cname::AbstractString)
    _chkfstr("cnmfrm", "string `cname`", cname)
    name = String(strip(cname))
    v = get(ks.pool, uppercase("OBJECT_$(name)_FRAME"), nothing)
    if v === nothing
        code = _body_name_code(ks.name2code, name)       # `bodn2c_`, without the check of `bodn2c_c`
        code === nothing || (v = get(ks.pool, "OBJECT_$(code)_FRAME", nothing))
    end
    v === nothing || return _frame_of_object(ks, v)
    frname = uppercase("IAU_" * name)
    code = get(BUILTIN_FRAMENAME2ID, normalize_name(frname), nothing)
    return code === nothing ? nothing : (Int(code), frname)
end

"""
    ccifrm([ctx,] frclss, clssid) -> Union{Nothing,Tuple{Int,String,Int}}

Return `(frcode, frname, center)` of the frame with the class `frclss` and the class
ID `clssid`, as
`ccifrm_c`, or `nothing` if there is none.

The built-in frames are searched first. If several frames of the kernel pool match,
the one with the smallest ID is returned; CSPICE returns the first one in the order
of its hash table. `nothing` takes the place of the `found` flag.

Throws `SpiceError(:invalidframedef)` if a `FRAME_*_CLASS` variable that the search
meets is not a number, or a matching definition lacks its `_CLASS_ID`, `_NAME` or
`FRAME_<name>` variable, and `SpiceError(:kernelvarnotfound)` if it lacks its
`FRAME_<code>_CENTER` variable.

# Examples
```jldoctest
julia> ccifrm(2, 499)
(10014, "IAU_MARS", 499)
```

See also [`frinfo`](@ref).

# References

- [NAIF documentation of `ccifrm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ccifrm_c.html)
"""
function ccifrm(ks::KernelSet, frclss::Integer, clssid::Integer)
    if frclss == 1 && 1 <= clssid <= N_INERTIAL_FRAMES
        return Int(clssid), INERTIAL_FRAME_DEFS[clssid][1], 0
    end
    for f in BUILTIN_NONINERTIAL_FRAMES
        f.class == frclss && f.classid == clssid && return Int(f.id), f.name, Int(f.center)
    end
    _ccifrm_check(ks, frclss, clssid)
    best = nothing
    for f in values(ks.frames)      # invariant: `best` is the match with the smallest ID so far
        (f.class == frclss && f.classid == clssid && !haskey(BUILTIN_FRAMES, f.id)) || continue
        (best === nothing || f.id < best.id) && (best = f)
    end
    return best === nothing ? nothing : (Int(best.id), best.name, Int(best.center))
end

"""
    _ccifrm_check(ks, frclss, clssid)

The checks of the pool search of `ccifrm_`, on the `FRAME_*_CLASS` variables in sorted
order up to the first that matches: a class that is not a number, or a matching class
without `_CLASS_ID`, a matching ID without `_NAME`, or a name without `FRAME_<name>`,
throws `SpiceError(:invalidframedef)`; a definition without `FRAME_<code>_CENTER`
throws `SpiceError(:kernelvarnotfound)`. CSPICE takes the order of its hash table, so it
can return a match before it meets a malformed definition that this check meets first.
"""
function _ccifrm_check(ks::KernelSet, frclss::Integer, clssid::Integer)
    bad(key, what) = throw(
        SpiceError(
            :invalidframedef, "invalid frame specification " *
                "found in kernel pool: frame class keyword is $key but $what"
        )
    )
    for key in sort!([k for k in keys(ks.pool) if occursin(r"^FRAME_.*_CLASS$", k)])
        v = ks.pool[key]
        (v isa Vector{Float64} && !isempty(v)) ||
            bad(key, "integer class was not associated with this keyword")
        _gipool_first(v, key) == frclss || continue
        id = get(ks.pool, key * "_ID", nothing)
        (id isa Vector{Float64} && !isempty(id)) ||
            bad(key, "associated integer class ID assignment was not found")
        _gipool_first(id, key * "_ID") == clssid || continue
        name = get(ks.pool, key[1:(end - 6)] * "_NAME", nothing)
        (name isa Vector{String} && !isempty(name)) ||
            bad(key, "associated frame name assignment was not found")
        code = get(ks.pool, "FRAME_" * strip(name[1]), nothing)
        (code isa Vector{Float64} && !isempty(code)) || throw(
            SpiceError(
                :invalidframedef,
                "invalid frame specification found in kernel pool: frame name is " *
                    "$(strip(name[1])) but associated frame ID assignment was not found"
            )
        )
        center = "FRAME_$(_gipool_first(code, "FRAME_" * strip(name[1])))_CENTER"
        haskey(ks.pool, center) || throw(
            SpiceError(
                :kernelvarnotfound,
                "the center of frame $(strip(name[1])), the kernel variable $center, is not in " *
                    "the kernel pool"
            )
        )
        return nothing
    end
    return nothing
end

"Throws `SpiceError(:badframeclass)` unless `frmcls` is -1 or a class 1 to 6, as `bltfrm_` and `kplfrm_`."
function _check_frame_class(who::AbstractString, frmcls::Integer)
    (frmcls == -1 || 1 <= frmcls <= 6) || throw(
        SpiceError(
            :badframeclass,
            "$who: the frame class $frmcls is not supported; it must be -1 (all) or 1 to 6"
        )
    )
    return nothing
end

"The built-in frames of `zzfdat_` as `(ID, class)`, in the order of the IDs."
const _BUILTIN_FRAME_CLASSES = sort!(
    vcat(
        [(i, 1) for i in eachindex(INERTIAL_FRAME_DEFS)],
        [(Int(f.id), Int(f.class)) for f in BUILTIN_NONINERTIAL_FRAMES]
    )
)

"""
    bltfrm(frmcls) -> Vector{Int}

Return the IDs of the built-in frames of the class `frmcls`, in increasing order, as
`bltfrm_c` of CSPICE.

`frmcls` is a frame class of [`frinfo`](@ref) (1 to 6), or -1 for every class. The
frames that the kernel pool defines are not included; [`kplfrm`](@ref) gives those.

Throws `SpiceError(:badframeclass)` if `frmcls` is not -1 or a class 1 to 6.

# Examples
```jldoctest
julia> length(bltfrm(1)), bltfrm(1)[1:3], bltfrm(3)
(21, [1, 2, 3], Int64[])
```

# References

- [NAIF documentation of `bltfrm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bltfrm_c.html)
"""
function bltfrm(frmcls::Integer)
    _check_frame_class("bltfrm", frmcls)
    return Int[id for (id, class) in _BUILTIN_FRAME_CLASSES if frmcls == -1 || class == frmcls]
end

"""
    kplfrm([ctx,] frmcls) -> Vector{Int}

Return the IDs of the frames of the class `frmcls` that the kernel pool defines, in
increasing order and without duplicates, as
`kplfrm_c` of CSPICE.

`frmcls` is a frame class of [`frinfo`](@ref) (1 to 6), or -1 for every class. As in
CSPICE, a frame is listed if the pool has

- `FRAME_<name> = <id>`, for the value `<name>` of a `FRAME_*_NAME` variable;
- `FRAME_<id>_NAME`;
- `FRAME_<id>_CLASS` or `FRAME_<name>_CLASS`.

The other variables of a frame definition are not checked, so an ID can be listed that
[`frinfo`](@ref) does not know. The built-in frames are not included unless the pool
defines them; [`bltfrm`](@ref) gives those.

Throws `SpiceError(:badframeclass)` if `frmcls` is not -1 or a class 1 to 6.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> lmpool(ctx, ["FRAME_TOPO = 1500001", "FRAME_1500001_NAME = 'TOPO'",
                    "FRAME_1500001_CLASS = 4", "FRAME_1500001_CLASS_ID = 1500001",
                    "FRAME_1500001_CENTER = 399"]);

julia> kplfrm(ctx, 4), kplfrm(ctx, 3)
([1500001], Int64[])
```

# References

- [NAIF documentation of `kplfrm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/kplfrm_c.html)
"""
function kplfrm(ks::KernelSet, frmcls::Integer)
    _check_frame_class("kplfrm", frmcls)
    pool = ks.pool
    # `gipool_` gives the first value of a numeric variable, rounded; names are 32
    # characters long in CSPICE, so a longer name is cut.
    getint(key) = (k = first(key, 32); _gipool_first(get(pool, k, nothing), k))
    ids = Int[]
    for (key, val) in pool
        _matchi(key, "FRAME_*_NAME") || continue
        (val isa Vector{String} && !isempty(val)) || continue
        name = strip(first(val[1], 32), ' ')
        isempty(name) && (name = " ")
        idcode = getint("FRAME_" * name)
        idcode === nothing && continue
        v = get(pool, "FRAME_$(idcode)_NAME", nothing)
        (v isa Vector{String} && !isempty(v)) || continue
        class = getint("FRAME_$(idcode)_CLASS")
        class === nothing && (class = getint("FRAME_$(name)_CLASS"))
        class === nothing && continue
        (frmcls == -1 || class == frmcls) && push!(ids, idcode)
    end
    return unique!(sort!(ids))
end

"""
    tkfram([ctx,] frcode) -> Tuple{SMatrix{3,3,Float64,9},Int}

Return `(rot, relative)` of the fixed-offset (class 4) frame with the class ID
`frcode`, as
`tkfram_c`. `frcode` is the value of `FRAME_<id>_CLASS_ID`.

`rot` rotates a vector from the frame to its relative frame, the frame with ID
`relative`; this is the meaning of `TKFRAME_<frcode>_MATRIX`. The definition comes
from the `TKFRAME_<frcode>_*` variables of the kernel pool, or from the
`TKFRAME_<name>_*` variables, where `<name>` is the name of the frame with the ID
`frcode`, as in CSPICE. A TK frame whose class ID is not its frame ID therefore needs
another frame with that ID, in CSPICE as here.

Throws `SpiceError`:
- `:zeroframeid` if `frcode == 0`;
- `:incompleteframe` if no frame has the ID `frcode`;
- `:competingframespec` if a variable is present in both forms;
- `:variablenotfound`, `:badvariablesize` or `:badvariabletype` if a variable that the
  definition needs is missing, has the wrong number of values or the wrong type;
- `:badframespec` if the relative frame is not known, `:badframespec2` if it is the
  frame itself;
- `:unknownframespec` if `TKFRAME_<frcode>_SPEC` is not `MATRIX`, `ANGLES` or
  `QUATERNION`;
- `:unitsnotrec` or `:incompatibleunits` if `TKFRAME_<frcode>_UNITS` is not an angle
  unit, and `:badaxisnumbers` if an axis is not 1, 2 or 3.

[`pxform`](@ref) and the other calls that pass through such a frame throw the same.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> lmpool(ctx, ["FRAME_TOPO = -1000", "FRAME_-1000_NAME = 'TOPO'",
                    "FRAME_-1000_CLASS = 4", "FRAME_-1000_CLASS_ID = -1000",
                    "FRAME_-1000_CENTER = 399", "TKFRAME_-1000_RELATIVE = 'J2000'",
                    "TKFRAME_-1000_SPEC = 'ANGLES'", "TKFRAME_-1000_ANGLES = (0, 0, 90)",
                    "TKFRAME_-1000_AXES = (1, 2, 3)", "TKFRAME_-1000_UNITS = 'DEGREES'"])

julia> rot, relative = tkfram(ctx, -1000);

julia> relative, round.(rot; digits=12) == [0 1 0; -1 0 0; 0 0 1]
(1, true)

julia> pdpool(ctx, "TKFRAME_-1000_AXES", [1.0, 2.0]);

julia> try tkfram(ctx, -1000) catch err; err.kind end
:badvariablesize
```

# References

- [NAIF documentation of `tkfram_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/tkfram_c.html)
"""
function tkfram(ks::KernelSet, frcode::Integer)
    r = frcode == 0 ? nothing : ccifrm(ks, 4, frcode)
    if r !== nothing
        tk = get(ks.tk, Int32(r[1]), nothing)
        tk === nothing || return tk.rot, Int(tk.relative)
    end
    tk = _tk_parse(ks.pool, ks.frames, frcode)    # throws the error of `tkfram_`
    return tk.rot, Int(tk.relative)
end

"ID of the built-in inertial frame `ref` (name or ID); throws `SpiceError(:irfnotrec)` otherwise, as `irfrot_`."
function _irf_id(ref::AbstractString)
    id = get(BUILTIN_FRAMENAME2ID, normalize_name(ref), Int32(0))
    1 <= id <= N_INERTIAL_FRAMES || throw(
        SpiceError(
            :irfnotrec,
            "the frame '$ref' is not a recognized inertial frame"
        )
    )
    return id
end
function _irf_id(ref::Integer)
    1 <= ref <= N_INERTIAL_FRAMES || throw(
        SpiceError(
            :irfnotrec,
            "the frame with ID $ref is not a recognized inertial frame"
        )
    )
    return Int32(ref)
end

"""
    _tis_jet(ks, ref, body, et) -> (R, dR)

Return the matrix from the inertial frame `ref` to the body-fixed frame of `body`
(a body ID, or the frame class ID of a binary PCK frame) at `et` and its time
derivative: `_body_rotation` (binary PCK first, text constants otherwise) times the
constant rotation from `ref` to J2000.
"""
function _tis_jet(ks::KernelSet, ref::NameOrID, body::Integer, et::Real)
    req = _irf_id(ref)
    t = check_epoch(et)
    # `tisbod_`: the matrix of the native frame of the PCK data, times `irfrot_` from
    # `req` to that frame if the two differ
    m = t isa Float64 ? _body_state(ks, Int32(body), t, req) :
        _body_rotation(ks, Int32(body), taylor_var(t, Val(1)), req)
    return map(x -> x.v, m), map(x -> x.d[1], m)
end

"""
    tipbod([ctx,] ref, body, et) -> SMatrix{3,3}

Return the matrix that rotates a position vector from the inertial frame `ref` to the
body-fixed frame of `body` at `et`, as
`tipbod_c`.

`ref` is the name or the ID of a built-in inertial frame. `body` is a NAIF body ID, or
the frame class ID of a binary PCK frame, the ID that its segments carry (31008 for
`MOON_PA_DE440`). `et` is in seconds past J2000 TDB. A binary PCK has priority over the
constants of a text PCK.

Throws `SpiceError`:
- `:irfnotrec` if `ref` is not a built-in inertial frame;
- `:badarg` if `et` is not finite;
- `:framedatanotfound` if no PCK data are loaded for the body;
- `:kernelvarnotfound` if the text constants of the body have `BODY<id>_PM` without
  `BODY<id>_POLE_RA` or `BODY<id>_POLE_DEC`;
- `:competingepochspec` or `:competingframespec` if the pool holds both forms of the
  epoch (`BODY<sys>_CONSTANTS_JED_EPOCH`, `BODY<sys>_CONSTS_JED_EPOCH`) or the frame
  variable of the body system;
- `:degreeoutofrange` if the system has phase angles and `BODY<sys>_MAX_PHASE_DEGREE`
  is outside 1:3; `:insufficientangles` if the body has more terms than angles.

See also [`tisbod`](@ref), [`pxform`](@ref).

# References

- [NAIF documentation of `tipbod_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/tipbod_c.html)
"""
tipbod(ks::KernelSet, ref::NameOrID, body::Integer, et::Real) = _tis_jet(ks, ref, body, et)[1]

"""
    tisbod([ctx,] ref, body, et) -> SMatrix{6,6}

Return the matrix that transforms a state from the inertial frame `ref` to the
body-fixed frame of `body` at `et`, as
`tisbod_c`.

The matrix is `[R 0; dR R]`, with `R` the matrix of [`tipbod`](@ref) and `dR` its
time derivative in 1/s. The arguments and the errors are those of `tipbod`.

!!! warning "Difference from CSPICE"
    For a body whose text PCK has cubic nutation and precession angles
    (`BODY<sys>_MAX_PHASE_DEGREE = 3`), `dR` differs from that of `tisbod_c`: CSPICE
    N0067 divides the rate of the cubic term by the square of the seconds per century,
    so its `dR` misses that term. This call gives the exact derivative.

See also [`sxform`](@ref).

# References

- [NAIF documentation of `tisbod_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/tisbod_c.html)
"""
function tisbod(ks::KernelSet, ref::NameOrID, body::Integer, et::Real)
    r, dr = _tis_jet(ks, ref, body, et)
    return vcat(hcat(r, zero(r)), hcat(dr, r))
end

"""
    pckfrm!(ids, pck) -> ids

Insert the frame class IDs of the frames that have data in the binary PCK file `pck`
into the vector `ids`, as
`pckfrm_c` does with a SPICE cell.

`ids` is used as a set: on return it is sorted and without duplicates. The file is
read from disk and need not be loaded.

!!! warning "Difference from CSPICE"
    The vector comes first, as in SPICE.jl; CSPICE has the cell last.

Throws the errors of [`pckfrm`](@ref).

# References

- [NAIF documentation of `pckfrm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pckfrm_c.html)
"""
function pckfrm!(ids::AbstractVector{<:Integer}, pck::AbstractString)
    _require_binary_kernel("pckfrm", pck, "PCK")
    sort!(unique!(ids))
    _each_summary(pck) do dc, ic, h
        k = searchsortedfirst(ids, ic[1])
        (k <= length(ids) && ids[k] == ic[1]) || insert!(ids, k, ic[1])
    end
    return ids
end

"""
    pckfrm(pck) -> Vector{Int}

Return the frame class IDs of the frames that have data in the binary PCK file `pck`,
sorted and without duplicates, as
`pckfrm_c`.

The file is read from disk and need not be loaded. An open [`DAFHandle`](@ref)
(in memory, or of a file) can take the place of the path: the call then reads the
segments that the handle holds now, without an array that `dafbna` began and `dafena`
did not end, and throws `SpiceError(:dafnosuchhandle)` if the handle is closed.

Throws `SpiceError(:filenotfound)` if the file does not exist, `SpiceError(:invalidformat)`
if it is a transfer file (a binary kernel in text form), `SpiceError(:invalidarchtype)` if
it is not a DAF, and `SpiceError(:invalidfiletype)` if it is a DAF of another kernel type.

See also [`pckcov`](@ref).

# References

- [NAIF documentation of `pckfrm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pckfrm_c.html)
"""
pckfrm(pck::AbstractString) = pckfrm!(Int[], pck)

"""
    pckcov!(cover, pck, idcode) -> cover

Insert the coverage of the frame with the class ID `idcode` in the binary PCK file
`pck` into the window `cover`, as
`pckcov_c`.

The result is the union of the old content of `cover` with the new intervals. See
[`pckcov`](@ref) for the intervals and the errors.

!!! warning "Difference from CSPICE"
    The window comes first, as in SPICE.jl; CSPICE has the cell last.

# References

- [NAIF documentation of `pckcov_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pckcov_c.html)
"""
function pckcov!(cover::Window, pck::AbstractString, idcode::Integer)
    _require_binary_kernel("pckcov", pck, "PCK")
    _each_summary(pck) do dc, ic, h
        ic[1] == idcode && wninsd!(cover, dc[1], dc[2])
    end
    return cover
end

"""
    pckcov(pck, idcode) -> Window

Return the coverage [`Window`](@ref) of the frame with the class ID `idcode` in the
binary PCK file `pck`, as
`pckcov_c`.

The window is the union of the intervals of the segments of that frame, in seconds
past J2000 TDB. It is empty if the frame has no data. The file is read from disk and
need not be loaded. An open [`DAFHandle`](@ref)
(in memory, or of a file) can take the place of the path: the call then reads the
segments that the handle holds now, without an array that `dafbna` began and `dafena`
did not end, and throws `SpiceError(:dafnosuchhandle)` if the handle is closed.

Throws the errors of [`pckfrm`](@ref).

See also [`pckcov!`](@ref), [`spkcov`](@ref).

# References

- [NAIF documentation of `pckcov_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pckcov_c.html)
"""
pckcov(pck::AbstractString, idcode::Integer) = pckcov!(Window(), pck, idcode)

"""
    bodfnd([ctx,] body, item) -> Bool

Return `true` if the kernel pool holds the variable `BODY<body>_<item>`, as
`bodfnd_c`.

`body` is a NAIF ID code. `item` is case sensitive; blanks around it are ignored. The
variable can be numeric or character.

See also [`bodvcd`](@ref).

# References

- [NAIF documentation of `bodfnd_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodfnd_c.html)
"""
bodfnd(ks::KernelSet, body::Integer, item::AbstractString) =
    haskey(ks.pool, _body_key(body, item))

"""
    bodvcd([ctx,] bodyid, item) -> Vector{Float64}

Return the values of the kernel-pool variable `BODY<bodyid>_<item>`, as
`bodvcd_c`.

`bodyid` is a NAIF ID code; [`bodvrd`](@ref) also takes a body name. `item` is case
sensitive. The result is a new vector with every value, in the units of the kernel;
there is no `maxn` argument.

Throws `SpiceError(:kernelvarnotfound)` if the variable is not in the pool and
`SpiceError(:typemismatch)` if it holds character values.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pdpool(ctx, "BODY399_RADII", [6378.1366, 6378.1366, 6356.7519]);

julia> bodvcd(ctx, 399, "RADII")
3-element Vector{Float64}:
 6378.1366
 6378.1366
 6356.7519
```

# References

- [NAIF documentation of `bodvcd_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvcd_c.html)
"""
bodvcd(ks::KernelSet, bodyid::Integer, item::AbstractString) =
    gdpool(ks, _body_key(bodyid, item))

# Context and default-context forms.
for f in (:frinfo, :cidfrm, :tkfram)
    @eval begin
        $f(ctx::SpiceContext, code::Integer) = $f(kernelset(ctx), code)
        $f(code::Integer) = $f(default_context(), code)
    end
end
cnmfrm(ctx::SpiceContext, cname::AbstractString) = cnmfrm(kernelset(ctx), cname)
cnmfrm(cname::AbstractString) = cnmfrm(default_context(), cname)
ccifrm(ctx::SpiceContext, frclss::Integer, clssid::Integer) = ccifrm(kernelset(ctx), frclss, clssid)
ccifrm(frclss::Integer, clssid::Integer) = ccifrm(default_context(), frclss, clssid)
for f in (:tipbod, :tisbod)
    @eval begin
        $f(ctx::SpiceContext, ref::NameOrID, body::Integer, et::Real) =
            $f(kernelset(ctx), ref, body, et)
        $f(ref::NameOrID, body::Integer, et::Real) = $f(default_context(), ref, body, et)
    end
end
bodfnd(ctx::SpiceContext, body::Integer, item::AbstractString) = bodfnd(kernelset(ctx), body, item)
bodfnd(body::Integer, item::AbstractString) = bodfnd(default_context(), body, item)
bodvcd(ctx::SpiceContext, bodyid::Integer, item::AbstractString) =
    bodvcd(kernelset(ctx), bodyid, item)
bodvcd(bodyid::Integer, item::AbstractString) = bodvcd(default_context(), bodyid, item)
kplfrm(ctx::SpiceContext, frmcls::Integer) = kplfrm(kernelset(ctx), frmcls)
kplfrm(frmcls::Integer) = kplfrm(default_context(), frmcls)
