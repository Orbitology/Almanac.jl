# Dynamic frames (frame class 5): the parameterized families of `zzdynfrm.c` and
# `zzdynrot.c`. The definitions are parsed once when the snapshot is built
# (`_dyn_parse`); the evaluation reads the parsed `DynFrame` only, is a pure function of
# the snapshot and does not allocate.
#
# Two modes, the ones of the CK frames (`ck.jl`):
#
# - `_CKExact()`: the rotation that `pxform` returns, as `zzdynrot_` computes it, in the
#   number type of the epoch. With a `Taylor` epoch the result holds the exact time
#   derivatives of that rotation.
# - `_CKState()`: the state transformation of `zzdynfrm_` at a `Float64` epoch, as a
#   matrix of `Taylor{1,Float64}` values whose derivative part is the derivative block of
#   `sxform`. It follows CSPICE where CSPICE does not differentiate the rotation: the
#   velocities of `spkezr`, the difference quotient of the unit velocity, the nutation
#   rates of `zzwahr_`, a constant vector with zero derivative, `ROTATION_STATE` and
#   `FREEZE_EPOCH`.
#
# Recursion: a dynamic frame of the snapshot `ks` is evaluated with the states and the
# frames of `ks.inner`, which holds the same data one level down. The innermost
# snapshot has `inner === nothing`, and a dynamic frame that is needed there throws
# `SpiceError(:recursiontoodeep)`. A definition can therefore refer to other dynamic
# frames, and those cannot, which is the limit of `zzfrmgt1_`.

"Family codes of `DynFrame.family`."
const _DYN_MEAN_EQ = Int8(1)
const _DYN_TRUE_EQ = Int8(2)
const _DYN_MEAN_ECL = Int8(3)
const _DYN_TWOVEC = Int8(4)
const _DYN_EULER = Int8(5)
const _DYN_PRODUCT = Int8(6)

"Number of snapshots below the one that the user holds: levels 1 and 2 of CSPICE."
const _DYN_LEVELS = 2

const _DYN_NOVECTOR = DynVector(0, 1, false, 0x00, 0, 0, 0, 0, 0, zero(SVector{3, Float64}))

# ---------------------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------------------

"""
    _dyn_lookup(pool, id, name, item, required) -> Union{Nothing,Tuple{String,PoolValue}}

Return the name and value of the frame definition variable `FRAME_<id>_<item>`, or of
`FRAME_<name>_<item>` if the first is absent, with the rules of `zzdynvac_`: a name
longer than 32 characters is not searched for. Throws
`SpiceError(:varnametoolong)` if both names are too long and
`SpiceError(:kernelvarnotfound)` if `required` and no variable exists; returns `nothing`
if not `required` and no variable exists.
"""
function _dyn_lookup(
        pool::AbstractDict, id::Integer, name::AbstractString,
        item::AbstractString, required::Bool
    )
    k1 = "FRAME_$(id)_$(item)"
    k2 = "FRAME_$(name)_$(item)"
    kv = _frame_key(pool, k1, k2)
    kv === nothing || return kv
    (length(k1) > _TK_MAXNAME && length(k2) > _TK_MAXNAME) && throw(
        SpiceError(
            :varnametoolong,
            "kernel variables $k1 and $k2 are longer than $_TK_MAXNAME characters"
        )
    )
    required || return nothing
    throw(
        SpiceError(
            :kernelvarnotfound, "at least one of the kernel variables $k1 or " *
                "$k2 is needed to define the frame $name, and neither is in the kernel pool"
        )
    )
end

@noinline _dyn_badtype(key, name, want) = throw(
    SpiceError(
        :badvariabletype,
        "the kernel variable $key of the definition of frame $name must be $want"
    )
)
@noinline _dyn_badsize(key, name, maxn, n) = throw(
    SpiceError(
        :badvariablesize,
        "the kernel variable $key of the definition of frame $name has $n values; at most " *
            "$maxn are allowed"
    )
)

"""
    _dyn_values(pool, id, name, item, T, maxn, required = true) -> Union{Vector{T}, Nothing}

Return the values of the definition item `item`, of element type `String` or `Float64`
and at most `maxn` of them, as `zzdynvac_`/`zzdynvad_` (required) and
`zzdynoac_`/`zzdynoad_` (optional) do. An optional item that is absent gives `nothing`.
"""
function _dyn_values(pool, id, name, item, ::Type{T}, maxn::Int, required::Bool = true) where {T}
    r = _dyn_lookup(pool, id, name, item, required)
    r === nothing && return nothing
    key, v = r
    v isa Vector{T} || _dyn_badtype(key, name, T === String ? "character" : "numeric")
    length(v) > maxn && _dyn_badsize(key, name, maxn, length(v))
    return v
end

_dyn_char(pool, id, name, item) = _dyn_values(pool, id, name, item, String, 1)[1]

"The single value of type `T` of an optional definition item, or `nothing`, as `zzdynoac_` and `zzdynoad_`."
function _dyn_optvalue(pool, id, name, item, ::Type{T}) where {T}
    v = _dyn_values(pool, id, name, item, T, 1, false)
    return v === nothing ? nothing : v[1]
end

"""
    _dyn_frameid(ks, id, name, item) -> Int32

Return the frame ID of a required definition item that holds a frame name or a frame
ID, as `zzdynfid_`. Throws `SpiceError(:notranslation)` for a name that is not a frame
name and not an integer.
"""
function _dyn_frameid(ks::KernelSet, id, name, item)
    key, v = _dyn_lookup(ks.pool, id, name, item, true)
    length(v) > 1 && _dyn_badsize(key, name, 1, length(v))
    v isa Vector{Float64} && return _gipool_first(v, key)::Int32
    code = _namfrm(ks, v[1])
    code != 0 && return Int32(code)
    n = tryparse(Int32, strip(v[1]))
    n === nothing && throw(
        SpiceError(
            :notranslation, "the kernel variable $key of the " *
                "definition of frame $name has the value '$(v[1])', which is not a frame name"
        )
    )
    return n
end

"""
    _dyn_bodyid(ks, id, name, item) -> Int32

Return the body ID of a required definition item that holds a body name or a body ID,
as `zzdynbid_`. Throws `SpiceError(:notranslation)` for an unknown name.
"""
function _dyn_bodyid(ks::KernelSet, id, name, item)
    key, v = _dyn_lookup(ks.pool, id, name, item, true)
    length(v) > 1 && _dyn_badsize(key, name, 1, length(v))
    v isa Vector{Float64} && return _gipool_first(v, key)::Int32
    code = bods2c(ks, v[1])
    code === nothing && throw(
        SpiceError(
            :notranslation, "the body name '$(v[1])' of the " *
                "kernel variable $key could not be translated to an ID code"
        )
    )
    return Int32(code)
end

"""
    _dyn_degrees(units) -> Float64

Return the number of degrees in one angle unit `units`, the factor of `convrt_`.
An angle `x` is converted to radians as `x * f / (180 / pi)`, the two operations of
`convrt_`, so that the result has the bits of CSPICE. Throws `SpiceError(:unitsnotrec)`
for a name that is not an angle unit, as `convrt_`.
"""
function _dyn_degrees(units::AbstractString)
    u = uppercase(strip(units))
    u == "RADIANS" && return 180.0 / pi
    u == "DEGREES" && return 1.0
    u == "ARCMINUTES" && return 0.016666666666666666
    u == "ARCSECONDS" && return 2.7777777777777778e-4
    u == "HOURANGLE" && return 15.0
    u == "MINUTEANGLE" && return 0.25
    u == "SECONDANGLE" && return 0.0041666666666666666
    throw(SpiceError(:unitsnotrec, "the angle unit '$units' of a dynamic frame is not recognized"))
end

"""
    _DynCorrError(msg)

The error of a position or near-point vector whose aberration correction is not one of
the nine values of `spkezr`. CSPICE checks that correction where it computes the vector,
so the kind depends on the computation: SPICE(SPKINVALIDOPTION) from `zzspkpa0_` in a
rotation (`zzdynrot_`), SPICE(INVALIDOPTION) from `zzvalcor_` in a state transformation
(`zzdynfrm_`).
"""
struct _DynCorrError <: Exception
    msg::String
end

"""
    _dyn_abcorr(s, kind) -> Aberration

Return the aberration correction `s` of an observer-target vector of `kind` 1
(position), 2 (velocity) or 3 (near point). `s` is cut to its first five bytes, the
length of the string in `zzdynrot_` and `zzdynfrm_`. A value other than the nine of
`spkezr` throws `SpiceError(:invalidoption)` for a velocity, which CSPICE checks with
`zzvalcor_` in both routines, and `_DynCorrError` for the other two kinds.
"""
function _dyn_abcorr(s::AbstractString, kind::Integer)
    try
        return Aberration(_cut_bytes(s, 5))
    catch e
        (e isa SpiceError && e.kind === :invalidoption && kind != 2) || rethrow()
        throw(_DynCorrError(e.msg))
    end
end

"""
    _dyn_constcorr(s) -> (Aberration, relativistic::Bool)

Return the aberration correction `s` of a constant vector, as `zzprscor_` parses it: the
nine values of `spkezr`, `S`, `XS`, and `RL`, `XRL`, which give `LT`, `XLT` with
`relativistic` true. `s` is cut to its first five bytes, as in `zzdynrot_`. Throws
`SpiceError(:invalidoption)` for another value and for a value with both light time and
stellar aberration (`LT+S`, `RL+S`, ...), as `zzdynrot_` and `zzdynfrm_` do.
"""
function _dyn_constcorr(s::AbstractString)
    u = _squeeze(_cut_bytes(s, 5))
    rel = u in ("RL", "RL+S", "XRL", "XRL+S")
    flags = u == "S" ? _AB_S : u == "XS" ? (_AB_S | _AB_X) :
        _abcorr_flags(_abcode(rel ? replace(u, "RL" => "LT") : u))
    flags == 0xff && throw(
        SpiceError(:invalidoption, "aberration correction '$s' of a constant vector is not recognized")
    )
    c = _aberration(flags)
    c.uselt && c.stellar && throw(
        SpiceError(
            :invalidoption, "aberration correction '$s' of a " *
                "constant vector: light time and stellar aberration are mutually exclusive"
        )
    )
    return c, rel
end

"Class and center of frame `id`; throws `SpiceError(:framedatanotfound)` if the frame is unknown."
function _dyn_frameinfo(ks::KernelSet, id::Int32, name)
    r = frinfo(ks, id)
    r === nothing && throw(
        SpiceError(
            :framedatanotfound, "the definition of frame " *
                "$name refers to the frame with ID $id, which has no frame specification"
        )
    )
    return Int32(r[2]), Int32(r[1])
end

"""
    _dyn_parse_vector(ks, fd, prefix) -> DynVector

Return the vector definition `FRAME_<id>_<prefix>*` (`prefix` is `PRI_` or `SEC_`)
of the two-vector frame `fd`, with the checks of `zzdynrot_` and `zzdynfrm_` in their
order: the correction of an observer-target vector is checked where CSPICE first uses
it, after the other items of the vector. Throws the `SpiceError` of the first check
that fails, or `_DynCorrError`.
"""
function _dyn_parse_vector(ks::KernelSet, fd::FrameDef, prefix::AbstractString)
    pool, id, name = ks.pool, fd.id, fd.name
    ax = _squeeze(_dyn_char(pool, id, name, prefix * "AXIS"))
    negate = startswith(ax, '-')
    ax = replace(ax, '-' => "", '+' => "")
    axis = ax == "X" ? 1 : ax == "Y" ? 2 : ax == "Z" ? 3 : 0
    axis == 0 && throw(
        SpiceError(
            :invalidaxis, "the definition of frame $name has the " *
                "axis '$ax'; the valid values are X, -X, Y, -Y, Z, -Z"
        )
    )
    vdef = _squeeze(_dyn_char(pool, id, name, prefix * "VECTOR_DEF"))
    zero3 = zero(SVector{3, Float64})
    if vdef == "OBSERVER_TARGET_POSITION"
        targ = _dyn_bodyid(ks, id, name, prefix * "TARGET")
        obs = _dyn_bodyid(ks, id, name, prefix * "OBSERVER")
        corr = _dyn_abcorr(_dyn_char(pool, id, name, prefix * "ABCORR"), 1)
        return DynVector(1, axis, negate, corr.flags, targ, obs, 1, 1, 0, zero3)
    elseif vdef == "OBSERVER_TARGET_VELOCITY"
        targ = _dyn_bodyid(ks, id, name, prefix * "TARGET")
        obs = _dyn_bodyid(ks, id, name, prefix * "OBSERVER")
        cs = _dyn_char(pool, id, name, prefix * "ABCORR")
        frid = _dyn_frameid(ks, id, name, prefix * "FRAME")
        corr = _dyn_abcorr(cs, 2)
        class, center = _dyn_frameinfo(ks, frid, name)
        return DynVector(2, axis, negate, corr.flags, targ, obs, frid, class, center, zero3)
    elseif vdef == "TARGET_NEAR_POINT"
        targ = _dyn_bodyid(ks, id, name, prefix * "TARGET")
        obs = _dyn_bodyid(ks, id, name, prefix * "OBSERVER")
        cs = _dyn_char(pool, id, name, prefix * "ABCORR")
        fix = cidfrm(ks, targ)
        fix === nothing && throw(
            SpiceError(
                :framedatanotfound, "the definition of " *
                    "frame $name needs the body-fixed frame of the target body $targ, and no " *
                    "frame is associated with this body"
            )
        )
        radii = bodvcd(ks, targ, "RADII")
        length(radii) == 3 || throw(
            SpiceError(
                :invalidcount,
                "only $(length(radii)) axes were found for ID $targ; three axes expected"
            )
        )
        all(>(0.0), radii) || throw(
            SpiceError(
                :badaxislength,
                "an axis of body $targ is negative or zero: BODY$(targ)_RADII = $radii"
            )
        )
        corr = _dyn_abcorr(cs, 3)
        class, center = _dyn_frameinfo(ks, Int32(fix[1]), name)
        return DynVector(
            3, axis, negate, corr.flags, targ, obs, fix[1], class, center,
            SVector{3, Float64}(radii)
        )
    elseif vdef == "CONSTANT"
        frid = _dyn_frameid(ks, id, name, prefix * "FRAME")
        cs = _dyn_optvalue(pool, id, name, prefix * "ABCORR", String)
        corr, rel = _dyn_constcorr(cs === nothing ? "NONE" : cs)
        class, center, obs = Int32(1), Int32(0), Int32(0)
        if corr.uselt
            class, center = _dyn_frameinfo(ks, frid, name)
            if class != 1
                obs = _dyn_bodyid(ks, id, name, prefix * "OBSERVER")
                # The light time to the center of the frame is that of `spkpos`, which
                # rejects a relativistic correction; in an inertial frame it is not used.
                rel && throw(
                    SpiceError(
                        :spkinvalidoption, "aberration correction '$cs' of a constant " *
                            "vector in frame $frid, which is not inertial: relativistic " *
                            "corrections are not supported"
                    )
                )
            end
        elseif corr.stellar
            obs = _dyn_bodyid(ks, id, name, prefix * "OBSERVER")
        end
        spec = _squeeze(_dyn_char(pool, id, name, prefix * "SPEC"))
        if spec == "RECTANGULAR"
            v = _dyn_values(pool, id, name, prefix * "VECTOR", Float64, 3)
            length(v) == 3 || _dyn_badsize("FRAME_$(id)_$(prefix)VECTOR", name, 3, length(v))
            dir = SVector{3, Float64}(v)
        elseif spec == "LATITUDINAL" || spec == "RA/DEC"
            f = _dyn_degrees(_dyn_char(pool, id, name, prefix * "UNITS"))
            lonitem, latitem = spec == "LATITUDINAL" ? ("LONGITUDE", "LATITUDE") : ("RA", "DEC")
            lon = _dyn_values(pool, id, name, prefix * lonitem, Float64, 1)[1] * f / (180.0 / pi)
            lat = _dyn_values(pool, id, name, prefix * latitem, Float64, 1)[1] * f / (180.0 / pi)
            slon, clon = sincos(lon)
            slat, clat = sincos(lat)
            dir = SVector{3, Float64}(clon * clat, slon * clat, slat)
        else
            throw(
                SpiceError(
                    :notsupported, "the definition of frame $name has the " *
                        "constant vector specification '$spec', which is not supported"
                )
            )
        end
        return DynVector(4, axis, negate, corr.flags, 0, obs, frid, class, center, dir)
    end
    throw(
        SpiceError(
            :notsupported, "the definition of the two-vector frame $name has " *
                "the vector definition '$vdef', which is not supported"
        )
    )
end

"""
    _dyn_parse(ks, fd) -> DynFrame

Return the definition of the class 5 frame `fd` from the pool of `ks`, with the
checks of `zzdynfrm_` in its order. The function does not throw: a definition that is
not valid gives a `DynFrame` with `errkind !== :none`, and an evaluation of the frame
throws that error, as CSPICE reports it at evaluation time.
"""
function _dyn_parse(ks::KernelSet, fd::FrameDef)
    try
        return _dyn_parse_checked(ks, fd)
    catch err
        if err isa _DynCorrError
            kind, state = :spkinvalidoption, :invalidoption
        else
            err isa SpiceError || rethrow()
            kind = state = err.kind
        end
        return DynFrame(
            fd.id, fd.name, fd.center, 0, 0, false, 0.0, false, false,
            _DYN_NOVECTOR, _DYN_NOVECTOR, 0.0, (Int32(0), Int32(0), Int32(0)), 0.0, 1.0,
            Float64[], Float64[], Float64[], Int32[], Int32[], kind, state, err.msg
        )
    end
end

"Like `_dyn_parse`, but throw the error of the first check that fails."
function _dyn_parse_checked(ks::KernelSet, fd::FrameDef)
    pool, id, name = ks.pool, fd.id, fd.name
    base = _dyn_frameid(ks, id, name, "RELATIVE")
    style = _dyn_char(pool, id, name, "DEF_STYLE")
    _eqstr(style, "PARAMETERIZED") || throw(
        SpiceError(
            :notsupported,
            "dynamic frame style '$style' (in the definition of frame $name) is not supported"
        )
    )
    family = _squeeze(_dyn_char(pool, id, name, "FAMILY"))
    fam = family == "MEAN_EQUATOR_AND_EQUINOX_OF_DATE" ? _DYN_MEAN_EQ :
        family == "TRUE_EQUATOR_AND_EQUINOX_OF_DATE" ? _DYN_TRUE_EQ :
        family == "MEAN_ECLIPTIC_AND_EQUINOX_OF_DATE" ? _DYN_MEAN_ECL :
        family == "TWO-VECTOR" ? _DYN_TWOVEC : family == "EULER" ? _DYN_EULER :
        family == "PRODUCT" ? _DYN_PRODUCT : Int8(0)
    ofdate = 1 <= fam <= 3
    freeze = _dyn_optvalue(pool, id, name, "FREEZE_EPOCH", Float64)
    frozen = freeze !== nothing
    rotsta = _dyn_optvalue(pool, id, name, "ROTATION_STATE", String)
    frozen && rotsta !== nothing && throw(
        SpiceError(
            :framedeferror, "the definition " *
                "of frame $name has both FREEZE_EPOCH and ROTATION_STATE; at most one is allowed"
        )
    )
    ofdate && !frozen && rotsta === nothing && throw(
        SpiceError(
            :framedeferror,
            "the definition of frame $name of the family $family has neither FREEZE_EPOCH " *
                "nor ROTATION_STATE; exactly one is required"
        )
    )
    inertial = false
    if rotsta !== nothing
        inertial = _eqstr(rotsta, "INERTIAL")
        inertial || _eqstr(rotsta, "ROTATING") || throw(
            SpiceError(
                :notsupported,
                "the definition of frame $name has ROTATION_STATE = '$rotsta'; the valid " *
                    "values are ROTATING and INERTIAL"
            )
        )
    end
    v1 = v2 = _DYN_NOVECTOR
    minsep = 0.0
    axes = (Int32(0), Int32(0), Int32(0))
    epoch, scale = 0.0, 1.0
    coef = (Float64[], Float64[], Float64[])
    from, to = Int32[], Int32[]
    if ofdate
        for (item, model, on) in (
                ("PREC_MODEL", "EARTH_IAU_1976", true),
                ("NUT_MODEL", "EARTH_IAU_1980", fam == _DYN_TRUE_EQ),
                ("OBLIQ_MODEL", "EARTH_IAU_1980", fam == _DYN_MEAN_ECL),
            )
            on || continue
            value = _dyn_char(pool, id, name, item)
            _eqstr(value, model) || throw(
                SpiceError(
                    :notsupported, "the definition " *
                        "of frame $name has $item = '$value', which is not recognized"
                )
            )
            fd.center == 399 || throw(
                SpiceError(
                    :invalidselection, "the definition of " *
                        "frame $name has the center $(fd.center) and $item = '$model'; this " *
                        "model is not applicable to that body"
                )
            )
        end
    elseif fam == _DYN_TWOVEC
        v1 = _dyn_parse_vector(ks, fd, "PRI_")
        v2 = _dyn_parse_vector(ks, fd, "SEC_")
        sep = _dyn_optvalue(pool, id, name, "ANGLE_SEP_TOL", Float64)
        minsep = sep === nothing ? 1.0e-3 : sep
        v1.axis == v2.axis && throw(
            SpiceError(
                :undefinedframe, "the definition of " *
                    "frame $name gives the same axis to both vectors"
            )
        )
    elseif fam == _DYN_EULER
        epoch = _dyn_values(pool, id, name, "EPOCH", Float64, 1)[1]
        key, a = _dyn_lookup(pool, id, name, "AXES", true)
        a isa Vector{Float64} || _dyn_badtype(key, name, "integer")
        length(a) == 3 || _dyn_badsize(key, name, 3, length(a))
        _gipool_first(a, key; room = 3)      # `zzdynvai_`: `gipool_` after the size check
        axes = map(x -> _gipool_value(x)::Int32, (a[1], a[2], a[3]))
        coef = ntuple(i -> copy(_dyn_values(pool, id, name, "ANGLE_$(i)_COEFFS", Float64, 20)), 3)
        scale = _dyn_degrees(_dyn_char(pool, id, name, "UNITS"))
        all(x -> 1 <= x <= 3, axes) || throw(
            SpiceError(
                :badaxisnumbers,
                "the Euler axes of frame $name are $axes; each must be 1, 2 or 3"
            )
        )
    elseif fam == _DYN_PRODUCT
        fnames = _dyn_values(pool, id, name, "FROM_FRAMES", String, 10)
        tnames = _dyn_values(pool, id, name, "TO_FRAMES", String, 10)
        length(fnames) == length(tnames) || throw(
            SpiceError(
                :badframecount,
                "the product frame $name has $(length(fnames)) \"from\" frames and " *
                    "$(length(tnames)) \"to\" frames; these counts must match"
            )
        )
        from = Int32[_namfrm(ks, s) for s in fnames]
        to = Int32[_namfrm(ks, s) for s in tnames]
        k = findfirst(==(0), vcat(from, to))
        k === nothing || throw(
            SpiceError(
                :unknownframe, "the product frame $name " *
                    "refers to the frame '$(vcat(fnames, tnames)[k])', which is not known"
            )
        )
    else
        throw(
            SpiceError(
                :notsupported, "dynamic frame family '$family' (in the " *
                    "definition of frame $name) is not supported"
            )
        )
    end
    return DynFrame(
        id, name, fd.center, base, fam, frozen, frozen ? freeze : 0.0, inertial,
        1 <= base <= N_INERTIAL_FRAMES, v1, v2, minsep, axes, epoch, scale, coef[1],
        coef[2], coef[3], from, to, :none, :none, ""
    )
end

"""
    KernelSet(loaded, daf, spk, pck, pool, name2code, code2name, frames, framename2id, tk,
        bodyrot, lsk, ck, sclk, ckmeta)

Return the snapshot with these fields, no dynamic or switch frame definitions, no
inner snapshot, no pool edits and no `boddef` pairs. `_dyn_levels` completes it.
"""
KernelSet(
    loaded, daf, spk, pck, pool, name2code, code2name, frames, framename2id, tk,
    bodyrot, lsk, ck, sclk, ckmeta
) = KernelSet(
    loaded, daf, spk, pck, pool, name2code,
    code2name, frames, framename2id, tk, bodyrot, lsk, ck, sclk, ckmeta,
    Dict{Int32, DynFrame}(), Dict{Int32, SwitchFrame}(), nothing, false, Tuple{String, Int32}[],
    DSKFile[], Dict{Int32, Vector{DSKSegment}}(), EKFile[], EKTable[]
)

"""
    _dyn_levels(base, edited, bodydefs) -> KernelSet

Return the snapshot with the data of `base`, the parsed definitions of its class 5
and class 6 frames, and the flags `edited` and `bodydefs`. If a dynamic frame exists,
the result has a chain of `_DYN_LEVELS` inner snapshots (`inner.ks`) with the same data,
the innermost with `inner === nothing`; otherwise `inner === nothing`. No data are
copied.
"""
function _dyn_levels(base::KernelSet, edited::Bool, bodydefs::Vector{Tuple{String, Int32}})
    dyn = Dict{Int32, DynFrame}()
    sw = Dict{Int32, SwitchFrame}()
    for fd in values(base.frames)
        haskey(BUILTIN_FRAMES, fd.id) && continue
        fd.class == 5 && (dyn[fd.id] = _dyn_parse(base, fd))
        fd.class == 6 && (sw[fd.id] = _sw_parse(base, fd))
    end
    ks = nothing
    # invariant: `ks` is the snapshot of the level below the one that is built next
    for _ in 0:(isempty(dyn) ? 0 : _DYN_LEVELS)
        ks = KernelSet(
            base.loaded, base.daf, base.spk, base.pck, base.pool,
            base.name2code, base.code2name, base.frames, base.framename2id, base.tk,
            base.bodyrot, base.lsk, base.ck, base.sclk, base.ckmeta, dyn, sw,
            ks === nothing ? nothing : DynLevel{KernelSet}(ks), edited, bodydefs,
            base.dskfiles, base.dsk, base.ekfiles, base.ektables
        )
    end
    return ks::KernelSet
end

# ---------------------------------------------------------------------------------------
# Small generic helpers
# ---------------------------------------------------------------------------------------

"Vector of `Taylor{1,Float64}` with the values `r` and the derivatives `v`."
@inline _dyn_pack(r::Vec3{Float64}, v::Vec3{Float64}) = map((a, b) -> Taylor((a, b)), r, v)

"`x` with a zero derivative part (state mode) or `x` itself (plain numbers)."
@inline _dyn_still(x::Taylor{1, Float64}) = Taylor((x.v, 0.0))
@inline _dyn_still(x::Float64) = x

@inline _dyn_lift(::Type{S}, x::Float64) where {S <: Real} = convert(S, x)
@inline _dyn_lift(::Type{S}, x::Taylor{1, Float64}) where {S <: Real} = convert(S, _dyn_still(x))

"""
    _dyn_xf(mode, ks, from, to, et) -> Mat3

Return the matrix from frame `from` to frame `to` at `et`: `frame_rotation` in the
exact mode (`zzrefch0_`), `_ck_state_rotation` in the state mode (`zzfrmch0_`).
"""
@inline _dyn_xf(::_CKExact, ks::KernelSet, from::Int32, to::Int32, et::T) where {T <: Real} =
    frame_rotation(ks, from, to, et)::Mat3{T}
@inline _dyn_xf(::_CKState, ks::KernelSet, from::Int32, to::Int32, et::Float64) =
    _ck_state_rotation(ks, from, to, et)::Mat3{Taylor{1, Float64}}

# ---------------------------------------------------------------------------------------
# Earth precession, nutation and obliquity models
# ---------------------------------------------------------------------------------------

"Radians per arcsecond, `rpd_() / 3600`."
const _DYN_ARCSEC = (pi / 180) / 3600.0

"""
    _dyn_prec(et) -> Mat3

Return the 1976 IAU precession matrix from J2000 to the mean equator and equinox
of `et` (TDB seconds), `[-z]_3 [theta]_2 [-zeta]_3` of `zzeprc76_`, in the number type
of `et`.
"""
@inline function _dyn_prec(et::S) where {S <: Real}
    zeta, z, theta = _dyn_prec_angles(et)
    return _m_eul2m_unchecked(-z, theta, -zeta, 3, 2, 3)
end

"""
    _dyn_prec_angles(et) -> (zeta, z, theta)

Return the precession angles of `zzeprc76_` (radians) at `et` (TDB seconds), with the
arithmetic of `zzeprc76_`, in the number type of `et`.
"""
@inline function _dyn_prec_angles(et::S) where {S <: Real}
    t = et / _SEC_PER_CENTURY
    zeta = t * (t * (t * 0.017998 + 0.30188) + 2306.2181) * _DYN_ARCSEC
    z = t * (t * (t * 0.018203 + 1.09468) + 2306.2181) * _DYN_ARCSEC
    theta = t * (t * (t * -0.041833 - 0.42665) + 2004.3109) * _DYN_ARCSEC
    return zeta, z, theta
end

"""
    _dyn_prec_rates(et) -> (dzeta, dz, dtheta)

Return the analytic rates (radians per second) of the angles of `_dyn_prec_angles`,
with the arithmetic of `zzeprc76_`.
"""
@inline function _dyn_prec_rates(et::S) where {S <: Real}
    t = et / _SEC_PER_CENTURY
    ts = 1.0 / _SEC_PER_CENTURY
    dzeta = ts * (t * (t * 3 * 0.017998 + 0.60375999999999996) + 2306.2181) * _DYN_ARCSEC
    dz = ts * (t * (t * 3 * 0.018203 + 2.1893600000000002) + 2306.2181) * _DYN_ARCSEC
    dtheta = ts * (t * (t * 3 * -0.041833 - 0.85329999999999995) + 2004.3109) * _DYN_ARCSEC
    return dzeta, dz, dtheta
end

"Return the mean obliquity of the ecliptic of date in radians (1980 IAU), as `zzmobliq_`."
@inline function _dyn_mobliq(et::S) where {S <: Real}
    t = et / 31557600.0 / 100.0
    return _DYN_ARCSEC * (t * (t * (t * 0.001813 - 5.9e-4) - 46.815) + 84381.448)
end

"""
    _DYN_WAHR

The 106 terms of the 1980 IAU nutation series as `zzwahr_` stores them: the five
multipliers of the fundamental arguments (l, l', F, D, Omega), then the longitude
coefficient and its rate, then the obliquity coefficient and its rate, in units of
1e-4 arcsecond and 1e-4 arcsecond per 1000 Julian years.
"""
const _DYN_WAHR = (
    (0, 0, 0, 0, 1, -171996, -1742, 92025, 89),
    (0, 0, 0, 0, 2, 2062, 2, -895, 5),
    (-2, 0, 2, 0, 1, 46, 0, -24, 0),
    (2, 0, -2, 0, 0, 11, 0, 0, 0),
    (-2, 0, 2, 0, 2, -3, 0, 1, 0),
    (1, -1, 0, -1, 0, -3, 0, 0, 0),
    (0, -2, 2, -2, 1, -2, 0, 1, 0),
    (2, 0, -2, 0, 1, 1, 0, 0, 0),
    (0, 0, 2, -2, 2, -13187, -16, 5736, -31),
    (0, 1, 0, 0, 0, 1426, -34, 54, -1),
    (0, 1, 2, -2, 2, -517, 12, 224, -6),
    (0, -1, 2, -2, 2, 217, -5, -95, 3),
    (0, 0, 2, -2, 1, 129, 1, -70, 0),
    (2, 0, 0, -2, 0, 48, 0, 1, 0),
    (0, 0, 2, -2, 0, -22, 0, 0, 0),
    (0, 2, 0, 0, 0, 17, -1, 0, 0),
    (0, 1, 0, 0, 1, -15, 0, 9, 0),
    (0, 2, 2, -2, 2, -16, 1, 7, 0),
    (0, -1, 0, 0, 1, -12, 0, 6, 0),
    (-2, 0, 0, 2, 1, -6, 0, 3, 0),
    (0, -1, 2, -2, 1, -5, 0, 3, 0),
    (2, 0, 0, -2, 1, 4, 0, -2, 0),
    (0, 1, 2, -2, 1, 4, 0, -2, 0),
    (1, 0, 0, -1, 0, -4, 0, 0, 0),
    (2, 1, 0, -2, 0, 1, 0, 0, 0),
    (0, 0, -2, 2, 1, 1, 0, 0, 0),
    (0, 1, -2, 2, 0, -1, 0, 0, 0),
    (0, 1, 0, 0, 2, 1, 0, 0, 0),
    (-1, 0, 0, 1, 1, 1, 0, 0, 0),
    (0, 1, 2, -2, 0, -1, 0, 0, 0),
    (0, 0, 2, 0, 2, -2274, -2, 977, -5),
    (1, 0, 0, 0, 0, 712, 1, -7, 0),
    (0, 0, 2, 0, 1, -386, -4, 200, 0),
    (1, 0, 2, 0, 2, -301, 0, 129, -1),
    (1, 0, 0, -2, 0, -158, 0, -1, 0),
    (-1, 0, 2, 0, 2, 123, 0, -53, 0),
    (0, 0, 0, 2, 0, 63, 0, -2, 0),
    (1, 0, 0, 0, 1, 63, 1, -33, 0),
    (-1, 0, 0, 0, 1, -58, -1, 32, 0),
    (-1, 0, 2, 2, 2, -59, 0, 26, 0),
    (1, 0, 2, 0, 1, -51, 0, 27, 0),
    (0, 0, 2, 2, 2, -38, 0, 16, 0),
    (2, 0, 0, 0, 0, 29, 0, -1, 0),
    (1, 0, 2, -2, 2, 29, 0, -12, 0),
    (2, 0, 2, 0, 2, -31, 0, 13, 0),
    (0, 0, 2, 0, 0, 26, 0, -1, 0),
    (-1, 0, 2, 0, 1, 21, 0, -10, 0),
    (-1, 0, 0, 2, 1, 16, 0, -8, 0),
    (1, 0, 0, -2, 1, -13, 0, 7, 0),
    (-1, 0, 2, 2, 1, -10, 0, 5, 0),
    (1, 1, 0, -2, 0, -7, 0, 0, 0),
    (0, 1, 2, 0, 2, 7, 0, -3, 0),
    (0, -1, 2, 0, 2, -7, 0, 3, 0),
    (1, 0, 2, 2, 2, -8, 0, 3, 0),
    (1, 0, 0, 2, 0, 6, 0, 0, 0),
    (2, 0, 2, -2, 2, 6, 0, -3, 0),
    (0, 0, 0, 2, 1, -6, 0, 3, 0),
    (0, 0, 2, 2, 1, -7, 0, 3, 0),
    (1, 0, 2, -2, 1, 6, 0, -3, 0),
    (0, 0, 0, -2, 1, -5, 0, 3, 0),
    (1, -1, 0, 0, 0, 5, 0, 0, 0),
    (2, 0, 2, 0, 1, -5, 0, 3, 0),
    (0, 1, 0, -2, 0, -4, 0, 0, 0),
    (1, 0, -2, 0, 0, 4, 0, 0, 0),
    (0, 0, 0, 1, 0, -4, 0, 0, 0),
    (1, 1, 0, 0, 0, -3, 0, 0, 0),
    (1, 0, 2, 0, 0, 3, 0, 0, 0),
    (1, -1, 2, 0, 2, -3, 0, 1, 0),
    (-1, -1, 2, 2, 2, -3, 0, 1, 0),
    (-2, 0, 0, 0, 1, -2, 0, 1, 0),
    (3, 0, 2, 0, 2, -3, 0, 1, 0),
    (0, -1, 2, 2, 2, -3, 0, 1, 0),
    (1, 1, 2, 0, 2, 2, 0, -1, 0),
    (-1, 0, 2, -2, 1, -2, 0, 1, 0),
    (2, 0, 0, 0, 1, 2, 0, -1, 0),
    (1, 0, 0, 0, 2, -2, 0, 1, 0),
    (3, 0, 0, 0, 0, 2, 0, 0, 0),
    (0, 0, 2, 1, 2, 2, 0, -1, 0),
    (-1, 0, 0, 0, 2, 1, 0, -1, 0),
    (1, 0, 0, -4, 0, -1, 0, 0, 0),
    (-2, 0, 2, 2, 2, 1, 0, -1, 0),
    (-1, 0, 2, 4, 2, -2, 0, 1, 0),
    (2, 0, 0, -4, 0, -1, 0, 0, 0),
    (1, 1, 2, -2, 2, 1, 0, -1, 0),
    (1, 0, 2, 2, 1, -1, 0, 1, 0),
    (-2, 0, 2, 4, 2, -1, 0, 1, 0),
    (-1, 0, 4, 0, 2, 1, 0, 0, 0),
    (1, -1, 0, -2, 0, 1, 0, 0, 0),
    (2, 0, 2, -2, 1, 1, 0, -1, 0),
    (2, 0, 2, 2, 2, -1, 0, 0, 0),
    (1, 0, 0, 2, 1, -1, 0, 0, 0),
    (0, 0, 4, -2, 2, 1, 0, 0, 0),
    (3, 0, 2, -2, 2, 1, 0, 0, 0),
    (1, 0, 2, -2, 0, -1, 0, 0, 0),
    (0, 1, 2, 0, 1, 1, 0, 0, 0),
    (-1, -1, 0, 2, 1, 1, 0, 0, 0),
    (0, 0, -2, 0, 1, -1, 0, 0, 0),
    (0, 0, 2, -1, 2, -1, 0, 0, 0),
    (0, 1, 0, 2, 0, -1, 0, 0, 0),
    (1, 0, -2, -2, 0, -1, 0, 0, 0),
    (0, -1, 2, 0, 1, -1, 0, 0, 0),
    (1, 1, 0, -2, 1, -1, 0, 0, 0),
    (1, 0, -2, 2, 0, -1, 0, 0, 0),
    (2, 0, 0, 2, 0, 1, 0, 0, 0),
    (0, 0, 2, 4, 2, -1, 0, 0, 0),
    (0, 1, 0, 1, 0, 1, 0, 0, 0),
)

"Fundamental arguments of `zzwahr_`: constant (degrees), rate (degrees per day), and the two higher coefficients, for l, l', F, D and Omega."
const _DYN_WAHR_ARGS = (
    (134.96298138888886, 13.064992947243136, 6.5192872572139397e-4, 3.6484365631332527e-7),
    (357.52772333333331, 0.98560028309377146, -1.201414483363923e-5, -6.8408185558748495e-8),
    (93.271910277777778, 13.229350240603848, -2.760338267929901e-4, 6.2707503428852773e-8),
    (297.85036305555559, 12.190749116495549, -1.4348262053484912e-4, 1.0831296046801845e-7),
    (125.04452222222223, -0.052953764841432813, 1.5522608272925558e-4, 4.5605457039165659e-8),
)

"Degrees per radian as `zzwahr_` computes it."
const _DYN_RADIAN = 180.0 / pi
"Divisor that converts the units of the nutation series to radians."
const _DYN_WAHR_FACTOR = _DYN_RADIAN * 3600.0 * 1.0e4

"""
    _dyn_wahr_args(et) -> (dd, t, ang)

Return the time arguments of `zzwahr_` at `et` (TDB seconds): `dd`, the days past
J2000 over 1e4, `t`, the Julian millennia, and `ang`, the five fundamental arguments
l, l', F, D, Omega in radians.
"""
@inline function _dyn_wahr_args(et::Real)
    dj = et / 86400.0
    dd = dj / 1.0e4
    dd2 = dd * dd
    t = dj / 365250.0
    ang = map(_DYN_WAHR_ARGS) do c
        f2c_dmod(c[1] + dj * c[2] + (c[3] + dd * c[4]) * dd2, 360.0) / _DYN_RADIAN
    end
    return dd, t, ang
end

"Coefficient `c0 + c1 * t` of a nutation term, in the type of `t`, as `zzwahr_` forms it (no product when `c1` is 0)."
@inline _dyn_wahr_coef(c0::Integer, c1::Integer, t::Real) = c1 != 0 ? c0 + c1 * t : zero(t) + c0

"""
    _dyn_wahr(et) -> (dpsi, deps)

Return the nutation in longitude and in obliquity (radians) of the 1980 IAU theory at
`et` (TDB seconds), with the operations of `zzwahr_`, in the number type of `et`. With
a `Taylor` epoch the derivatives are those of the series, including the slow change
of the coefficients that the rates of `zzwahr_` leave out.
"""
function _dyn_wahr(et::S) where {S <: Real}
    dd, t, ang = _dyn_wahr_args(et)
    dpsi = zero(S)
    deps = zero(S)
    for i in 1:106      # invariant: dpsi, deps hold the sums of the terms 1:i-1
        m = @inbounds _DYN_WAHR[i]
        arg = zero(S)
        for j in 1:5
            if m[j] != 0
                arg += m[j] * ang[j]
                arg = f2c_dmod(arg, _TWOPI)
            end
        end
        cl = _dyn_wahr_coef(m[6], m[7], t)
        ce = _dyn_wahr_coef(m[8], m[9], t)
        s, c = sincos(arg)
        dpsi += cl * s / _DYN_WAHR_FACTOR
        deps += ce * c / _DYN_WAHR_FACTOR
    end
    return dpsi, deps
end

"""
    _dyn_wahr_rates(et) -> (dpsi, deps, ddpsi, ddeps)

Return the nutation angles of `_dyn_wahr` and the rates that `zzwahr_` returns
(radians per second): the rates of the arguments, without the change of the
coefficients. Generic over the number type of `et`; with a `Float64` epoch the
arithmetic is that of `zzwahr_`.
"""
function _dyn_wahr_rates(et::S) where {S <: Real}
    dd, t, ang = _dyn_wahr_args(et)
    dddj = dd / 1.0e4
    rate = map(_DYN_WAHR_ARGS) do c
        f2c_dmod(c[2] + dddj * (c[3] * 2.0 + dd * 3.0 * c[4]), 360.0) / _DYN_RADIAN
    end
    dpsi = deps = ddpsi = ddeps = zero(S)
    for i in 1:106      # invariant: the four sums hold the terms 1:i-1
        m = @inbounds _DYN_WAHR[i]
        arg = argrt = zero(S)
        for j in 1:5
            if m[j] != 0
                arg += m[j] * ang[j]
                argrt += m[j] * rate[j]
                arg = f2c_dmod(arg, _TWOPI)
            end
        end
        cl = _dyn_wahr_coef(m[6], m[7], t)
        ce = _dyn_wahr_coef(m[8], m[9], t)
        s, c = sincos(arg)
        dpsi += cl * s / _DYN_WAHR_FACTOR
        deps += ce * c / _DYN_WAHR_FACTOR
        ddpsi += cl * c * argrt / _DYN_WAHR_FACTOR
        ddeps -= ce * s * argrt / _DYN_WAHR_FACTOR
    end
    return dpsi, deps, ddpsi / 86400.0, ddeps / 86400.0
end

"""
    _dyn_wahr_state(et::Float64) -> (dpsi, deps)

Return the nutation angles of `_dyn_wahr` as `Taylor{1,Float64}` values whose
derivative parts are the rates of `_dyn_wahr_rates`.
"""
@inline function _dyn_wahr_state(et::Float64)
    dpsi, deps, ddpsi, ddeps = _dyn_wahr_rates(et)
    return Taylor((dpsi, ddpsi)), Taylor((deps, ddeps))
end

"Nutation angles in the number type of the mode: `_dyn_wahr` or `_dyn_wahr_state`."
@inline _dyn_nutation(::_CKExact, et::Real) = _dyn_wahr(et)
@inline _dyn_nutation(::_CKState, et::Float64) = _dyn_wahr_state(et)

"""
    _dyn_ofdate(mode, ks, def, et) -> Mat3

Return the matrix from the of-date frame `def` (families 1 to 3) to its base frame
at `et`: the transpose of the precession matrix, of nutation times precession, or of
the obliquity rotation times precession, then the rotation from J2000 to the base.
"""
function _dyn_ofdate(mode::_CKExact, ks::KernelSet, def::DynFrame, et::T) where {T <: Real}
    m = _dyn_prec(et)
    if def.family == _DYN_TRUE_EQ
        dpsi, deps = _dyn_nutation(mode, et)
        mob = _dyn_mobliq(et)
        m = _m_mxm(_m_eul2m_unchecked(-mob - deps, -dpsi, mob, 1, 3, 1), m)
    elseif def.family == _DYN_MEAN_ECL
        m = _m_mxm(_m_rotate(_dyn_mobliq(et), 1), m)
    end
    x = typeof(m)(transpose(m))
    def.base == _J2000_FRAME && return x
    return _m_mxm(_dyn_xf(mode, ks, _J2000_FRAME, def.base, et), x)
end

# State mode, as `zzdynfrm_`: the state transformations of `zzeprc76_`, `zzenut80_` and
# of the obliquity rotation come from the angles and their analytic rates by `eul2xf_`,
# their products are `mxmg_` (`_sx_mxm`) and the inverse is `invstm_` (the transpose of
# each block). The derivative of the `Taylor` evaluation is the same function, rounded
# otherwise.
function _dyn_ofdate(mode::_CKState, ks::KernelSet, def::DynFrame, et::Float64)
    m = _dyn_prec_state(et)
    if def.family == _DYN_TRUE_EQ
        m = _sx_mxm(_dyn_enut80_state(et), m)
    elseif def.family == _DYN_MEAN_ECL
        mob, dmob = _dyn_mobliq_state(et)
        m = _sx_mxm(_dyn_eul2xf(0.0, 0.0, mob, 0.0, 0.0, dmob, 1, 3, 1), m)
    end
    x = Mat3{Taylor{1, Float64}}(transpose(m))
    def.base == _J2000_FRAME && return x
    return _sx_mxm(_dyn_xf(mode, ks, _J2000_FRAME, def.base, et), x)
end

"""
    _dyn_eul2xf(a3, a2, a1, d3, d2, d1, x3, x2, x1) -> Mat3{Taylor{1,Float64}}

Return the state transformation of `eul2xf_` for the angles `a3`, `a2`, `a1` about
the axes `x3`, `x2`, `x1` and their rates `d3`, `d2`, `d1`, as a matrix of
`Taylor{1,Float64}` (value: the rotation block; derivative: the derivative block).
"""
@inline function _dyn_eul2xf(a3, a2, a1, d3, d2, d1, x3::Integer, x2::Integer, x1::Integer)
    r, d = _m_blocks(eul2xf(SVector{6, Float64}(a3, a2, a1, d3, d2, d1), x3, x2, x1))
    return map((a, b) -> Taylor((a, b)), r, d)::Mat3{Taylor{1, Float64}}
end

"""
    _dyn_prec_state(et::Float64) -> Mat3{Taylor{1,Float64}}

Return the precession matrix of `_dyn_prec` with its derivative, with the arithmetic
of `zzeprc76_`: the angles and their analytic rates go to `eul2xf_`.
"""
function _dyn_prec_state(et::Float64)
    zeta, z, theta = _dyn_prec_angles(et)
    dzeta, dz, dtheta = _dyn_prec_rates(et)
    return _dyn_eul2xf(-z, theta, -zeta, -dz, dtheta, -dzeta, 3, 2, 3)
end

"Return the mean obliquity of `_dyn_mobliq` and its rate (radians per second), with the arithmetic of `zzmobliq_`, in the number type of `et`."
@inline function _dyn_mobliq_state(et::Real)
    t = et / 31557600.0 / 100.0
    mob = _dyn_mobliq(et)
    dmob = _DYN_ARCSEC * (t * (t * 3 * 0.001813 - 0.0011800000000000001) - 46.815) *
        (1.0 / (31557600.0 * 100.0))
    return mob, dmob
end

"""
    _dyn_enut80_state(et::Float64) -> Mat3{Taylor{1,Float64}}

Return the nutation matrix of the true of-date frame with its derivative, with the
arithmetic of `zzenut80_`: the angles of `zzwahr_` (`_dyn_wahr_state`) and
`zzmobliq_` with their rates go to `eul2xf_`.
"""
function _dyn_enut80_state(et::Float64)
    dpsi, deps = _dyn_wahr_state(et)
    mob, dmob = _dyn_mobliq_state(et)
    return _dyn_eul2xf(
        -mob - deps.v, -dpsi.v, mob, -dmob - deps.d[1], -dpsi.d[1], dmob, 1, 3, 1
    )
end

# ---------------------------------------------------------------------------------------
# States and positions of the vector definitions
# ---------------------------------------------------------------------------------------

"""
    _dyn_stelab(p, v) -> Vec3

Return the position `p` with the stellar aberration for the observer velocity `v`,
`cos(phi) p + |p| v_perp / c`, the rotation of `stelab_` in closed form (no inverse
sine), for any number type. Requires `p != 0` and `|v| < c`.
"""
@inline function _dyn_stelab(p, v)
    pm = _norm3(p)
    u = p / pm
    vp = v - _m_dot(v, u) * u
    c = sqrt(1 - _m_dot(vp, vp) / CLIGHT^2)
    return c * p + (pm / CLIGHT) * vp
end

"""
    _dyn_apo(ks, targ, et, ref, corr, obs) -> (r, lt)

Return the position of `targ` relative to `obs` in frame `ref` with the correction
`corr`, and the light time, as `zzspkzp0_` (`spkpos`). For a `Float64` epoch the values
are those of `spkpos`. For another number type the same algorithm runs on that type
(one light-time iteration for `LT`, three for `CN`, stellar aberration in closed form),
so the derivatives are those of the position that `spkpos` returns. No allocation.
"""
function _dyn_apo(
        ks::KernelSet, targ::Int32, et::T, ref::Int32, corr,
        obs::Int32
    ) where {T <: Real}
    # The body relative to itself: zero, and no square root of a zero series.
    targ == obs && return zero(Vec3{T}), zero(T)
    if corr.geometric
        r = spk_jet(ks, targ, et, ref, obs, Val(0))[1]::Vec3{T}
        return r, _m_norm(r) / CLIGHT
    end
    T === Float64 && return abcorr_position(ks, targ, et, ref, corr, obs)::Tuple{Vec3{T}, T}
    inertial = _is_inertial(ks, ref)
    iref = inertial ? ref : _J2000_FRAME
    nos = _aberration(corr.flags & ~_AB_S)
    robs, vobs = _ssb(ks, obs, et, iref)::NTuple{2, Vec3{T}}
    r, lt = _apo(ks, targ, et, iref, nos, robs, vobs)::Tuple{Vec3{T}, T}
    if corr.stellar && primal(_m_dot(r, r)) > 0
        r = _dyn_stelab(r, corr.transmit ? -vobs : +vobs)
    end
    inertial && return r, lt
    center = frame_def(ks, ref).center
    ltc = center == obs ? zero(lt) : center == targ ? lt :
        _apo(ks, center, et, iref, nos, robs, vobs)[2]::T
    return _m_mxv(_dyn_xf(_CKExact(), ks, _J2000_FRAME, ref, et + _ltsign(corr) * ltc), r), lt
end

"""
    _dyn_spkez(ks, targ, et, ref, corr, obs) -> (r, v, lt)

Return the position, velocity and light time of `targ` relative to `obs` in frame `ref`
with the correction `corr`, as `zzspkez0_` (`spkezr`), in the number type of `et`. The
velocity is the one of CSPICE (see `spkder` for its difference from the derivative of
the position). No allocation.
"""
function _dyn_spkez(
        ks::KernelSet, targ::Int32, et::T, ref::Int32, corr,
        obs::Int32
    ) where {T <: Real}
    if corr.geometric
        r, v = spk_jet(ks, targ, et, ref, obs, Val(1))::NTuple{2, Vec3{T}}
        return r, v, _m_norm(r) / CLIGHT
    end
    s, lt = abcorr_state(ks, targ, et, ref, corr, obs)::Tuple{SVector{6, T}, T}
    return Vec3{T}(s[1], s[2], s[3]), Vec3{T}(s[4], s[5], s[6]), lt
end

"""
    _dyn_nearpt(p, radii) -> Vec3

Return the point of the ellipsoid with the semi-axes `radii` (all positive) that is
nearest to the point `p`, in the number type of `p`: the near point of `nearpt`
(`_m_nearpt`), so the plain value is the one of CSPICE and each derivative order comes
from the stationarity condition. Postcondition: the result is on the ellipsoid to
rounding.
"""
@inline _dyn_nearpt(p::Vec3{S}, radii::Vec3{Float64}) where {S <: Real} =
    Vec3{S}(_m_nearpt(p, radii[1], radii[2], radii[3])[1])

"""
    _dyn_dnearp(s, radii, vd) -> Vec3{Taylor{1,Float64}}

Return the state of the near point on the ellipsoid with the semi-axes `radii` for
the state `s`, as `dnearp_` gives it to `zzdynfrm_` (`dnearp`: the velocity from
`zzdnpt_`). Throws `SpiceError(:degeneratecase)` where `dnearp_` finds no velocity of
the near point, as `zzdynfrm_` does.
"""
function _dyn_dnearp(s::Vec3{Taylor{1, Float64}}, radii::Vec3{Float64}, vd::DynVector)
    st = SVector(s[1].v, s[2].v, s[3].v, s[1].d[1], s[2].d[1], s[3].d[1])
    r = dnearp(st, radii[1], radii[2], radii[3])
    r === nothing && throw(
        SpiceError(
            :degeneratecase,
            "a defining vector is the near point on body $(vd.target) as seen from body " *
                "$(vd.observer); the state of this near point was not found (see dnearp)"
        )
    )
    n = r[1]
    return Vec3{Taylor{1, Float64}}(Taylor((n[1], n[4])), Taylor((n[2], n[5])), Taylor((n[3], n[6])))
end

"Epoch `et` moved by the light time `lt` in the direction of the correction `corr`, as `zzcorepc_`."
@inline _dyn_corepc(corr, et, lt) = corr.uselt ? et + _ltsign(corr) * lt : et

"The epoch `et` corrected by `corr` for the light time from the observer to the center of the frame of `vd` if `uselt` and the frame is not inertial; `et` otherwise."
@inline function _dyn_frame_epoch(ks::KernelSet, vd::DynVector, corr, et::Real, uselt::Bool)
    (uselt && vd.frclass != 1) || return et
    lt = _dyn_apo(ks, vd.frcenter, et, _J2000_FRAME, corr, vd.observer)[2]
    return _dyn_corepc(corr, et, lt)
end

"""
    _dyn_vector(mode, ks, vd, et) -> Vec3

Return the defining vector `vd` of a two-vector frame in J2000 at `et`, before the
sign of the axis is applied.

Exact mode (`zzdynrot_`), in the number type of `et`: the position of `spkpos`; the
unit vector of the velocity of `spkezr`, rotated from its frame at the epoch corrected
for the light time to the center of that frame; the vector from the observer to the
near point on the target ellipsoid, rotated from the body-fixed frame at the epoch
corrected for the light time; the constant vector rotated from its frame, with the
optional stellar aberration.

State mode (`zzdynfrm_`), `Taylor{1,Float64}` values: the state of `spkezr`; the unit
velocity and the central difference of the unit velocity over `±max(1, et * 2^-27)` s;
the state of the near point for the state of `spkezr`; the constant vector with zero
derivative before the rotation. A stellar aberration correction of a constant vector
changes the vector and leaves the derivative as it is, as in CSPICE.
"""
function _dyn_vector(mode::_CKExact, ks::KernelSet, vd::DynVector, et::T) where {T <: Real}
    corr = _aberration(vd.corr)
    if vd.kind == 1
        return _dyn_apo(ks, vd.target, et, _J2000_FRAME, corr, vd.observer)[1]::Vec3{T}
    elseif vd.kind == 2
        _, v, _ = _dyn_spkez(ks, vd.target, et, vd.frame, corr, vd.observer)
        u = _m_hat(v)
        vd.frame == _J2000_FRAME && return u::Vec3{T}
        vet = _dyn_frame_epoch(ks, vd, corr, et, corr.uselt)
        return _m_mxv(_dyn_xf(_CKExact(), ks, vd.frame, _J2000_FRAME, vet), u)::Vec3{T}
    elseif vd.kind == 3
        p, lt = _dyn_apo(ks, vd.target, et, vd.frame, corr, vd.observer)
        pobs = -p
        d = _dyn_nearpt(pobs, vd.vec) - pobs
        fet = _dyn_corepc(corr, et, lt)
        return _m_mxv(_dyn_xf(_CKExact(), ks, vd.frame, _J2000_FRAME, fet), d)::Vec3{T}
    end
    # A constant vector.
    fet = _dyn_frame_epoch(ks, vd, corr, et, corr.uselt && !corr.stellar)
    c = Vec3{T}(vd.vec)
    vd.frame == _J2000_FRAME || (c = _m_mxv(_dyn_xf(_CKExact(), ks, vd.frame, _J2000_FRAME, fet), c))
    if corr.stellar
        vobs = spk_jet(ks, vd.observer, et, _J2000_FRAME, Int32(0), Val(1))[2]
        vs = corr.transmit ? -vobs : +vobs
        # A plain epoch takes `stelab_` itself, as `zzdynrot_` does; another number type
        # the closed form, which has the same value up to rounding.
        c = T === Float64 ? _stelab(c, vs) : _dyn_stelab(c, vs)
    end
    return c::Vec3{T}
end

function _dyn_vector(mode::_CKState, ks::KernelSet, vd::DynVector, et::Float64)
    corr = _aberration(vd.corr)
    if vd.kind == 1
        r, v, _ = _dyn_spkez(ks, vd.target, et, _J2000_FRAME, corr, vd.observer)
        return _dyn_pack(r, v)
    elseif vd.kind == 2
        delta = max(1.0, et * 7.4505805969238281e-9)
        um = _m_hat(_dyn_spkez(ks, vd.target, et - delta, vd.frame, corr, vd.observer)[2])
        u0 = _m_hat(_dyn_spkez(ks, vd.target, et, vd.frame, corr, vd.observer)[2])
        up = _m_hat(_dyn_spkez(ks, vd.target, et + delta, vd.frame, corr, vd.observer)[2])
        s = _dyn_pack(u0, (0.5 / delta) * up + (-0.5 / delta) * um)
        vd.frame == _J2000_FRAME && return s
        vet = _dyn_frame_epoch(ks, vd, corr, et, corr.uselt)
        return _sx_mxv(_dyn_xf(_CKState(), ks, vd.frame, _J2000_FRAME, vet), s)
    elseif vd.kind == 3
        r, v, lt = _dyn_spkez(ks, vd.target, et, vd.frame, corr, vd.observer)
        sobs = _dyn_pack(-r, -v)
        d = _dyn_dnearp(sobs, vd.vec, vd) - sobs
        fet = _dyn_corepc(corr, et, lt)
        return _sx_mxv(_dyn_xf(_CKState(), ks, vd.frame, _J2000_FRAME, fet), d)
    end
    fet = _dyn_frame_epoch(ks, vd, corr, et, corr.uselt && !corr.stellar)
    c = _dyn_pack(vd.vec, zero(Vec3{Float64}))
    vd.frame == _J2000_FRAME || (c = _sx_mxv(_dyn_xf(_CKState(), ks, vd.frame, _J2000_FRAME, fet), c))
    if corr.stellar
        vobs = spk_jet(ks, vd.observer, et, _J2000_FRAME, Int32(0), Val(1))[2]
        p = _stelab(map(x -> x.v, c), corr.transmit ? -vobs : +vobs)
        c = map((a, b) -> Taylor((a, b.d[1])), p, c)
    end
    return c
end

@noinline function _dyn_throw_sep(def::DynFrame, sep, et)
    throw(
        SpiceError(
            :degeneratecase, "angular separation of the vectors that define " *
                "the two-vector frame $(def.name) is $sep radians at ET $(primal(et)); the " *
                "minimum allowed difference from 0 or pi is $(def.minsep) radians " *
                "(keyword ANGLE_SEP_TOL)"
        )
    )
end

"""
    _dyn_twovec(mode, ks, def, et) -> Mat3

Return the matrix from the two-vector frame `def` to its base frame at `et`. The
columns are the frame axes in J2000 as `zztwovxf_` and `twovec_` build them: the axis
of the primary vector along that vector, the axis of the secondary vector in the
half-plane of the two vectors. Throws `SpiceError(:degeneratecase)` if the sine of the
separation of the vectors is below the sine of `def.minsep`.
"""
function _dyn_twovec(mode, ks::KernelSet, def::DynFrame, et::T) where {T <: Real}
    a = _dyn_vector(mode, ks, def.v1, et)
    b = _dyn_vector(mode, ks, def.v2, et)
    def.v1.negate && (a = -a)
    def.v2.negate && (b = -b)
    sep = _m_vsep(map(primal, a), map(primal, b))
    sin(sep) < sin(def.minsep) && _dyn_throw_sep(def, sep, et)
    x = _dyn_twovec_axes(a, b, Int(def.v1.axis), Int(def.v2.axis))
    def.base == _J2000_FRAME && return x
    return _ck_mxm(mode, _dyn_xf(mode, ks, _J2000_FRAME, def.base, et), x)
end

"""
    _dyn_twovec_axes(a, b, i1, ip) -> Mat3

Return the matrix whose columns are the axes, in the frame of `a` and `b`, of the
frame with axis `i1` along `a` and axis `ip` in the half-plane of `a` and `b`, with the
arithmetic of CSPICE: for vectors with a derivative (`Taylor{1,Float64}`) that of
`zztwovxf_` (`dvhat_`, `ducrss_`), which `zzdynfrm_` calls, for any other number type
that of `twovec_` (`vhat_`, `ucrss_`), which `zzdynrot_` calls. Requires `i1 != ip`,
both in 1:3.
"""
@inline function _dyn_twovec_axes(a::Vec3{S}, b::Vec3{S}, i1::Int, ip::Int) where {S <: Real}
    i2 = i1 % 3 + 1
    c1 = _m_hat(a)
    if ip == i2
        c3 = _m_ucrss(a, b)
        c2 = _m_ucrss(c3, a)
    else
        c2 = _m_ucrss(b, a)
        c3 = _m_ucrss(a, c2)
    end
    return _dyn_axes_matrix(c1, c2, c3, i1)
end

@inline function _dyn_twovec_axes(
        a::Vec3{Taylor{1, Float64}}, b::Vec3{Taylor{1, Float64}}, i1::Int, ip::Int
    )
    pack(t) = SVector{6, Float64}(t[1].v, t[2].v, t[3].v, t[1].d[1], t[2].d[1], t[3].d[1])
    unpack(s) = Vec3{Taylor{1, Float64}}(
        Taylor((s[1], s[4])), Taylor((s[2], s[5])), Taylor((s[3], s[6]))
    )
    sa, sb = pack(a), pack(b)
    i2 = i1 % 3 + 1
    s1 = _m_dvhat(sa, false)
    if ip == i2
        s3 = _m_ducrss(sa, sb)
        s2 = _m_ducrss(s3, sa)
    else
        s2 = _m_ducrss(sb, sa)
        s3 = _m_ducrss(sa, s2)
    end
    return _dyn_axes_matrix(unpack(s1), unpack(s2), unpack(s3), i1)
end

"The matrix with column `i1` `c1`, the next column (cyclically) `c2` and the last `c3`."
@inline _dyn_axes_matrix(c1, c2, c3, i1::Int) =
    i1 == 1 ? hcat(c1, c2, c3) : i1 == 2 ? hcat(c3, c1, c2) : hcat(c2, c3, c1)

"""
    _dyn_euler(mode, def, et) -> Mat3

Return the matrix from the Euler frame `def` to its base frame at `et`,
`[a1]_x1 [a2]_x2 [a3]_x3` with the angles `a_i` polynomials in `et - def.epoch`.
"""
function _dyn_euler(mode, def::DynFrame, et::T) where {T <: Real}
    delta = _ck_epoch(mode, et) - def.epoch
    a1 = _horner(def.coef1, delta) * def.scale / _DYN_RADIAN
    a2 = _horner(def.coef2, delta) * def.scale / _DYN_RADIAN
    a3 = _horner(def.coef3, delta) * def.scale / _DYN_RADIAN
    return _dyn_euler_matrix(a1, a2, a3, def.axes)
end

@inline _dyn_euler_matrix(a1::S, a2::S, a3::S, axes::NTuple{3, Int32}) where {S <: Real} =
    _m_eul2m_unchecked(a1, a2, a3, axes[1], axes[2], axes[3])

# State mode: the angles and their rates (the Horner scheme of a `Taylor` epoch is the
# arithmetic of `polyds_`) go to `eul2xf_`, as in `zzdynfrm_`; the derivative of the
# `Taylor` product of the three rotations rounds otherwise.
@inline function _dyn_euler_matrix(
        a1::Taylor{1, Float64}, a2::Taylor{1, Float64}, a3::Taylor{1, Float64},
        axes::NTuple{3, Int32}
    )
    return _dyn_eul2xf(a1.v, a2.v, a3.v, a1.d[1], a2.d[1], a3.d[1], axes[1], axes[2], axes[3])
end

"""
    _dyn_product(mode, ks, def, et) -> Mat3

Return the matrix from the product frame `def` to its base frame at `et`: the
product, last factor on the left, of the matrices from `def.to[i]` to `def.from[i]`.
"""
function _dyn_product(mode, ks::KernelSet, def::DynFrame, et::T) where {T <: Real}
    x = _dyn_xf(mode, ks, (@inbounds def.to[1]), (@inbounds def.from[1]), et)
    for i in 2:length(def.to)       # invariant: x is the product of the factors 1:i-1
        x = _ck_mxm(mode, _dyn_xf(mode, ks, (@inbounds def.to[i]), (@inbounds def.from[i]), et), x)
    end
    return x
end

"""
    _dyn_inert(mode, ks, def, et, x) -> Mat3

Return `x` with the derivative part that `ROTATION_STATE = 'INERTIAL'` defines, as
`zzdynfrm_`: zero relative to J2000, so zero if the base is a built-in inertial frame
and otherwise the derivative of the rotation from J2000 to the base, applied to the
fixed orientation. The exact mode has no derivative part to replace and returns `x`.
"""
@inline _dyn_inert(::_CKExact, ks::KernelSet, def::DynFrame, et, x) = x
function _dyn_inert(
        ::_CKState, ks::KernelSet, def::DynFrame, et::Float64,
        x::Mat3{Taylor{1, Float64}}
    )
    def.inertial || return x
    def.base_inertial && return map(_dyn_still, x)
    # `zzdynfrm_`: the product with the transformation to J2000, its derivative block
    # zeroed, and the product with the inverse transformation (`mxmg_`, `invstm_`).
    b = _dyn_xf(_CKState(), ks, def.base, _J2000_FRAME, et)
    out = _sx_mtxm(b, map(_dyn_still, _sx_mxm(b, x)))
    return map((a, o) -> Taylor((a.v, o.d[1])), x, out)
end

"""
    _dyn_eval(mode, ks, def, et) -> Mat3

Return the matrix from the dynamic frame `def` to its base frame at `et`, evaluated
with the states and frames of `ks` (the inner snapshot of the one that holds `def`).
Requires `def.errkind === :none`. No allocation.
"""
function _dyn_eval(mode, ks::KernelSet, def::DynFrame, et::T) where {T <: Real}
    S = _ck_numtype(mode, T)
    fam = def.family
    x = fam <= _DYN_MEAN_ECL ? _dyn_ofdate(mode, ks, def, et)::Mat3{S} :
        fam == _DYN_TWOVEC ? _dyn_twovec(mode, ks, def, et)::Mat3{S} :
        fam == _DYN_EULER ? _dyn_euler(mode, def, et)::Mat3{S} :
        _dyn_product(mode, ks, def, et)::Mat3{S}
    return _dyn_inert(mode, ks, def, et, x)::Mat3{S}
end

# ---------------------------------------------------------------------------------------
# The recursion barrier
# ---------------------------------------------------------------------------------------
#
# The evaluation of a dynamic frame calls the state and frame functions of the package,
# and these call the evaluation of a dynamic frame again when a definition names one.
# CSPICE limits this recursion to two levels at run time, but the compiler sees a cycle
# without a limit in which the number types can nest (a series of series), and type
# inference then gives up on some calls of the cycle, also for kernels without any
# dynamic frame. The call from `_dyn_frame_step` to `_dyn_eval` therefore goes through
# a typed C function pointer: the argument and result types are exact, nothing is boxed,
# and the compiler does not follow the call. There is one pointer for each mode and
# number type, made on first use.

"Holder of the table of pointers; the field is atomic, so that a reader without the lock sees a complete table."
mutable struct _DynPointerTable
    @atomic table::IdDict{Any, Ptr{Cvoid}}
end

"Pointers of `_dyn_call`, by `Tuple{mode type, number type}`. The table is replaced, never mutated."
const _DYN_POINTERS = _DynPointerTable(IdDict{Any, Ptr{Cvoid}}())
const _DYN_POINTER_LOCK = ReentrantLock()

"C entry points of `_dyn_call`: `_dyn_eval` with the snapshot of the level `lev`."
_dyn_entry_exact(lev::DynLevel{KernelSet}, def::DynFrame, et::T) where {T <: Real} =
    _dyn_eval(_CKExact(), lev.ks, def, et)
_dyn_entry_state(lev::DynLevel{KernelSet}, def::DynFrame, et::Float64) =
    _dyn_eval(_CKState(), lev.ks, def, et)

"Make the C function pointer of `_dyn_eval` for the mode and number type of `F`."
@generated function _dyn_cfunction(::Type{Tuple{_CKExact, T}}) where {T}
    return :(
        @cfunction(
            _dyn_entry_exact, $(Mat3{T}),
            (Ref{DynLevel{KernelSet}}, Ref{DynFrame}, $T)
        )
    )
end
_dyn_cfunction(::Type{Tuple{_CKState, Float64}}) = @cfunction(
    _dyn_entry_state,
    Mat3{Taylor{1, Float64}}, (Ref{DynLevel{KernelSet}}, Ref{DynFrame}, Float64)
)

"""
    _dyn_pointer_slow(F) -> Ptr{Cvoid}

Make the pointer for `F`, publish a new table that has it, and return it. The
table is replaced under a lock, so a reader without the lock sees the old table or the
new one, both complete. The compiler does not follow the call that makes the pointer.
"""
@noinline function _dyn_pointer_slow(@nospecialize(F))
    return @lock _DYN_POINTER_LOCK begin
        p = get(@atomic(:acquire, _DYN_POINTERS.table), F, C_NULL)::Ptr{Cvoid}
        if p == C_NULL
            p = Base.invokelatest(_dyn_cfunction, F)::Ptr{Cvoid}
            table = copy(@atomic(:acquire, _DYN_POINTERS.table))
            table[F] = p
            @atomic :release _DYN_POINTERS.table = table
        end
        p
    end
end

"Pointer of `_dyn_call` for `F = Tuple{mode type, number type}`. No allocation after the first use."
@inline function _dyn_pointer(::Type{F}) where {F}
    p = get(@atomic(:acquire, _DYN_POINTERS.table), F, C_NULL)::Ptr{Cvoid}
    return p == C_NULL ? _dyn_pointer_slow(F) : p
end

"Forget the pointers of an earlier process (they are not valid in a loaded package image)."
function __init__()
    @atomic :release _DYN_POINTERS.table = IdDict{Any, Ptr{Cvoid}}()
    return nothing
end

"""
    _dyn_call(mode, lev, def, et) -> Mat3

Return `_dyn_eval(mode, lev.ks, def, et)` through the typed function pointer of the mode and of the number type of `et` (see "The recursion barrier"). A number type
that is not plain bits is evaluated with a direct call. No allocation.
"""
@generated function _dyn_call(
        ::_CKExact, lev::DynLevel{KernelSet}, def::DynFrame,
        et::T
    ) where {T <: Real}
    isbitstype(T) || return :(_dyn_eval(_CKExact(), lev.ks, def, et))
    return :(
        ccall(
            _dyn_pointer(Tuple{_CKExact, $T}), $(Mat3{T}),
            (Ref{DynLevel{KernelSet}}, Ref{DynFrame}, $T), lev, def, et
        )
    )
end
@inline _dyn_call(::_CKState, lev::DynLevel{KernelSet}, def::DynFrame, et::Float64) =
    ccall(
    _dyn_pointer(Tuple{_CKState, Float64}), Mat3{Taylor{1, Float64}},
    (Ref{DynLevel{KernelSet}}, Ref{DynFrame}, Float64), lev, def, et
)

@noinline _dyn_throw_missing(id) = throw(
    SpiceError(
        :framedatanotfound,
        "no definition of the dynamic frame $id is loaded"
    )
)
@noinline _dyn_throw_depth(def::DynFrame) = throw(
    SpiceError(
        :recursiontoodeep,
        "the dynamic frame $(def.name) is needed for the evaluation of a dynamic frame " *
            "that is itself needed for the evaluation of a dynamic frame; dynamic frames may " *
            "not be used at recursion level 1"
    )
)

# `::Type{S}` with a type parameter: Julia 1.10 does not specialize a closure that
# captures a type, and the call then allocates.
"Return the matrix `x` with each element in the number type `S` (`_dyn_lift`)."
@inline _dyn_liftmat(::Type{S}, x::Mat3) where {S} = Mat3{S}(ntuple(i -> _dyn_lift(S, x[i]), Val(9)))

"The kind of the error of the definition `def` that is not valid: that of `zzdynrot_` for a rotation, of `zzdynfrm_` for a state transformation."
_dyn_errkind(::_CKExact, def::DynFrame) = def.errkind
_dyn_errkind(::_CKState, def::DynFrame) = def.errkind_state

"""
    _dyn_frame_step(mode, ks, fd, et) -> (base::Int32, rot)

One link of a frame chain for the class 5 frame `fd`: its base frame (`RELATIVE`) and
the matrix from the frame to the base at `et`, in the number type of the mode. A
frozen frame is evaluated at its `FREEZE_EPOCH` and is constant. Throws the error of
an invalid definition, and `SpiceError(:recursiontoodeep)` at the recursion limit.
No allocation. Not inlined, so that the frame functions that hold the call stay small
for the frames of the other classes.
"""
@noinline function _dyn_frame_step(mode, ks::KernelSet, fd::FrameDef, et::T) where {T <: Real}
    def = get(ks.dyn, fd.id, nothing)
    def === nothing && _dyn_throw_missing(fd.id)
    def.errkind === :none || _throw_spice(_dyn_errkind(mode, def), def.errmsg)
    inner = ks.inner
    inner === nothing && _dyn_throw_depth(def)
    S = _ck_numtype(mode, T)
    if def.frozen
        x = _dyn_call(mode, inner, def, def.t0)
        return def.base, _dyn_liftmat(S, x)
    end
    return def.base, _dyn_call(mode, inner, def, et)::Mat3{S}
end

"""
    _dyn_chain_step(mode, ks, fd, et) -> (parent::Int32, rot)

One link of a frame chain for a frame of class 5 (`_dyn_frame_step`) or class 6
(`_sw_frame_step`), the only classes that `_chain_step` passes.
"""
@noinline function _dyn_chain_step(mode, ks::KernelSet, fd::FrameDef, et::T) where {T <: Real}
    fd.class == 5 && return _dyn_frame_step(mode, ks, fd, et)
    return _sw_frame_step(mode, ks, fd, et)
end
