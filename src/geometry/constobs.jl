# States for an observer or a target with a constant position or a constant velocity
# relative to its center of motion: `spkcpo`, `spkcpt`, `spkcvo`, `spkcvt` of CSPICE
# N0067, with their private routines `zzcvxsta`, `zzspkfzo`, `zzspkfzt`, `zzspkfao`,
# `zzspkfat`, `zzspkfap`, `zzspkflt` and `zzcorsxf`.
#
# The functions are generic over the epoch number type; all discrete decisions and the
# iteration counts use `primal`. No function allocates.

"""
    _GeoFixed

An object with a constant velocity relative to the body `ctr`: the saved values of
`zzcvssta`. `pos` and `vel` are in the frame `ref` at the epoch `epoch`.
"""
struct _GeoFixed
    pos::Vec3{Float64}
    vel::Vec3{Float64}
    epoch::Float64
    ctr::Int32
    ref::Int32
end

"""
    _geo_cvx(ks, f, et, ref) -> (r, v)

`zzcvxsta`: the state of the object `f` relative to its center at `et`, in the frame
`ref`: the position extrapolated with the constant velocity, transformed with the
state transformation from `f.ref` to `ref` at `et`.
"""
@inline function _geo_cvx(ks::KernelSet, f::_GeoFixed, et, ref::Int32)
    x = f.pos + (et - f.epoch) * f.vel
    v = f.vel - zero(x)
    f.ref == ref && return x, v
    return rotate_jet(frame_jet(ks, f.ref, ref, et, Val(1)), (x, v))
end

@inline function _geo_cvx_ssb(ks::KernelSet, f::_GeoFixed, et, iref::Int32)
    r, v = _geo_cvx(ks, f, et, iref)
    rc, vc = _ssb(ks, f.ctr, et, iref)
    return r + rc, v + vc
end

"""
    _geo_flt(ks, f, et, iref, corr, robs, vobs) -> (; r, v, lt, dlt)

`zzspkflt`: the state of the fixed target `f` relative to an observer with the
barycentric state `(robs, vobs)` at `et`, in the inertial frame `iref`, corrected for
light time as `corr` says, the light time and its rate. The iteration is that of
`spkltc`, with the barycentric state of the target from `_geo_cvx_ssb`.
"""
function _geo_flt(ks::KernelSet, f::_GeoFixed, et, iref::Int32, corr::Aberration, robs, vobs)
    rt, vt = _geo_cvx_ssb(ks, f, et, iref)
    r = rt - robs
    v = vt - vobs
    dist = _m_norm(r)
    lt = dist / CLIGHT
    primal(lt) == 0 && return (; r, v, lt, dlt = zero(lt))
    corr.uselt || return (; r, v, lt, dlt = _m_dot(r, v) / (dist * CLIGHT))
    s = _ltsign(corr)
    numitr = corr.converged ? 5 : 1
    lterr = 1.0
    i = 0
    # invariant: `lt` is the light time after `i` iterations
    while i < numitr && lterr > 1.0e-17
        epoch = et + s * lt
        rt, vt = _geo_cvx_ssb(ks, f, epoch, iref)
        r = rt - robs
        v = vt - vobs
        prvlt = lt
        lt = _m_norm(r) / CLIGHT
        lterr = abs(primal(lt) - primal(prvlt)) / max(1.0, abs(primal(epoch)))
        i += 1
    end
    a = 1 / (CLIGHT * _m_norm(r))
    b = _m_dot(r, v)
    c = _m_dot(r, vt)
    primal(s * c * a) > 1 - 1.0e-10 && _throw_range_rate()
    dlt = a * b / (1 - s * c * a)
    v = (1 + s * dlt) * vt - vobs
    return (; r, v, lt, dlt)
end

"""
    _geo_fa(ks, f::_GeoFixed, et, iref, corr, obs) -> (; r, v, lt, dlt)

`zzspkfat` with `zzspkfap`: the state of the fixed target `f` relative to the body
`obs` in the inertial frame `iref`, corrected for light time and stellar aberration.
The acceleration of the observer is the central difference of its velocity over ±1 s.
"""
@inline function _geo_fa(
        ks::KernelSet, f::_GeoFixed, et, iref::Int32, corr::Aberration,
        obs::Int32
    )
    robs, vobs, accobs = _obs_state(ks, obs, et, iref, corr)
    q = _geo_flt(ks, f, et, iref, corr, robs, vobs)
    corr.stellar || return q
    scorr, dscorr = _zzstelab(corr.transmit, accobs, vobs, q.r, q.v)
    return (; r = q.r + scorr, v = q.v + dscorr, q.lt, q.dlt)
end

"""
    _geo_fa(ks, targ, et, iref, corr, f::_GeoFixed) -> (; r, v, lt, dlt)

`zzspkfao`: the state of the body `targ` relative to the fixed observer `f` in the
inertial frame `iref`, corrected as `corr` says (`spkaps` with the barycentric state
of the observer; its acceleration is the central difference of the velocity over ±1 s).
"""
@inline function _geo_fa(
        ks::KernelSet, targ::Int32, et, iref::Int32, corr::Aberration,
        f::_GeoFixed
    )
    robs, vobs = _geo_cvx_ssb(ks, f, et, iref)
    accobs = zero(vobs)
    if corr.stellar
        vm = _geo_cvx_ssb(ks, f, et - 1.0, iref)[2]
        vp = _geo_cvx_ssb(ks, f, et + 1.0, iref)[2]
        accobs = 0.5 * vp - 0.5 * vm
    end
    q = _aps(ks, targ, et, iref, corr, robs, vobs, accobs)
    return (; q.r, q.v, q.lt, q.dlt)
end

"""
`zzcorsxf` applied to a state: `[R 0; (1 + s dlt) dR R] * (r, v)`, the derivative block
scaled element by element and the product formed as `mxvg_` forms it (`rotate_jet`).
"""
@inline function _geo_corsxf(R, dR, scale, r, v)
    pr, pv = rotate_jet((R, map(x -> x * scale, dR)), (r, v))
    return vcat(pr, pv)
end

@inline function _geo_frame_class(ks::KernelSet, ref::Int32)
    center, class, found = _geo_frinfo(ks, ref)
    found || throw(SpiceError(:unknownframe, "the frame with ID $ref is not recognized"))
    return center, class
end

# The geometric state, and the light time and its rate between the observer and the
# center `center` of a frame, for `_geo_fz`: one method for a fixed observer
# (`zzspkfzo`) and one for a fixed target (`zzspkfzt`).
@inline function _geo_fz_geometric(ks::KernelSet, targ::Int32, et, ref::Int32, f::_GeoFixed)
    ro, vo = _geo_cvx(ks, f, et, ref)
    rt, vt = spk_jet(ks, targ, et, ref, f.ctr, Val(1))
    r = rt - ro
    return vcat(r, vt - vo), _m_norm(r) / CLIGHT
end
@inline function _geo_fz_geometric(ks::KernelSet, f::_GeoFixed, et, ref::Int32, obs::Int32)
    rt, vt = _geo_cvx(ks, f, et, ref)
    rc, vc = spk_jet(ks, f.ctr, et, ref, obs, Val(1))
    r = rt + rc
    return vcat(r, vt + vc), _m_norm(r) / CLIGHT
end

@inline function _geo_fz_ltcent(
        ks::KernelSet, targ::Int32, et, corr::Aberration, f::_GeoFixed,
        center::Int32, q
    )
    center == targ && return q.lt, q.dlt
    robs, vobs = _geo_cvx_ssb(ks, f, et, Int32(1))
    qc = _ltc(ks, center, et, Int32(1), corr, robs, vobs)
    return qc.lt, qc.dlt
end
@inline function _geo_fz_ltcent(
        ks::KernelSet, ::_GeoFixed, et, corr::Aberration, obs::Int32,
        center::Int32, q
    )
    center == obs && return zero(q.lt), zero(q.dlt)
    robs, vobs = _ssb(ks, obs, et, Int32(1))
    qc = _ltc(ks, center, et, Int32(1), corr, robs, vobs)
    return qc.lt, qc.dlt
end

"""
    _geo_fz(ks, targ, et, ref, corr, obs) -> (state, lt)

`zzspkfzo` (the observer `obs` is a `_GeoFixed`) or `zzspkfzt` (the target `targ` is
one): the state of `targ` relative to `obs` in the frame `ref`. A frame that is not
inertial is evaluated at the epoch corrected for the light time between the observer
and the center of the frame.
"""
function _geo_fz(ks::KernelSet, targ, et, ref::Int32, corr::Aberration, obs)
    corr.geometric && return _geo_fz_geometric(ks, targ, et, ref, obs)
    center, class = _geo_frame_class(ks, ref)
    if class == Int32(1)
        q = _geo_fa(ks, targ, et, ref, corr, obs)
        return vcat(q.r, q.v), q.lt
    end
    q = _geo_fa(ks, targ, et, Int32(1), corr, obs)
    ltcent, dltctr = _geo_fz_ltcent(ks, targ, et, corr, obs, center, q)
    s = _ltsign(corr)
    R, dR = frame_jet(ks, Int32(1), ref, et + s * ltcent, Val(1))
    return _geo_corsxf(R, dR, 1 + s * dltctr, q.r, q.v), q.lt
end

const _GEO_AT_OBSERVER = 1
const _GEO_AT_TARGET = 2
const _GEO_AT_CENTER = 3

@noinline _geo_bad_refloc(refloc) = throw(
    SpiceError(
        :notsupported,
        "the output frame evaluation locus \"$refloc\" is not recognized; the values are " *
            "OBSERVER, TARGET and CENTER"
    )
)

"The locus code of `refloc` (`esrchc`: case and blanks are ignored); `SpiceError(:notsupported)` otherwise."
@inline function _geo_refloc(refloc::AbstractString)
    refloc == "OBSERVER" && return _GEO_AT_OBSERVER
    refloc == "TARGET" && return _GEO_AT_TARGET
    refloc == "CENTER" && return _GEO_AT_CENTER
    return _geo_refloc_slow(refloc)
end
@noinline function _geo_refloc_slow(refloc::AbstractString)
    _eqstr(refloc, "OBSERVER") && return _GEO_AT_OBSERVER
    _eqstr(refloc, "TARGET") && return _GEO_AT_TARGET
    _eqstr(refloc, "CENTER") && return _GEO_AT_CENTER
    return _geo_bad_refloc(refloc)
end

@noinline _geo_unknown_frame(ref) = throw(
    SpiceError(
        :unknownframe,
        "the reference frame $ref is not recognized"
    )
)

@inline function _geo_known_frame(ks::KernelSet, ref::NameOrID)
    id = _frame_code(ks, ref)
    _geo_frinfo(ks, id)[3] || _geo_unknown_frame(ref)
    return id
end

"""
    _geo_cv(ks, targ, et, orfcde, evltyp, corr, obs) -> (state, lt)

The computation of `spkcvo` (the observer `obs` is a `_GeoFixed`) or `spkcvt` (the
target `targ` is one) after the checks: the state of `targ` relative to `obs` in the
frame `orfcde`, evaluated at the locus `evltyp`.
"""
function _geo_cv(
        ks::KernelSet, targ, et, orfcde::Int32, evltyp::Int, corr::Aberration, obs
    )
    (!corr.uselt || evltyp == _GEO_AT_CENTER) && return _geo_fz(ks, targ, et, orfcde, corr, obs)
    if evltyp == _GEO_AT_OBSERVER
        st, lt = _geo_fz(ks, targ, et, Int32(1), corr, obs)
        orfcde == Int32(1) && return st, lt
        pr, pv = rotate_jet(frame_jet(ks, Int32(1), orfcde, et, Val(1)), _m_pv(st))
        return vcat(pr, pv), lt
    end
    q = _geo_fa(ks, targ, et, Int32(1), corr, obs)
    orfcde == Int32(1) && return vcat(q.r, q.v), q.lt
    s = _ltsign(corr)
    R, dR = frame_jet(ks, Int32(1), orfcde, et + s * q.lt, Val(1))
    return _geo_corsxf(R, dR, 1 + s * q.dlt, q.r, q.v), q.lt
end

"""
    _geo_cv_args(ks, abcorr, ctr, ctrwhat, body, bodywhat, outref, refloc)
        -> (corr, ctrcde, bodcde, orfcde, evltyp)

The first checks of `spkcvo`, `spkcpo`, `spkcvt` and `spkcpt`, in their order: the
correction, the center `ctr` of the constant state, the other body `body`, the output
frame `outref` and the locus `refloc`. The checks of the state and of its frame follow
in the caller.
"""
@inline function _geo_cv_args(
        ks::KernelSet, abcorr::AbcorrArg, ctr::NameOrID, ctrwhat::String,
        body::NameOrID, bodywhat::String, outref::NameOrID, refloc::AbstractString
    )
    corr = Aberration(abcorr)
    ctrcde = _body_id(ks, ctr, ctrwhat)
    bodcde = _body_id(ks, body, bodywhat)
    orfcde = _frame_code(ks, outref)
    orfcde == Int32(0) && _geo_unknown_frame(outref)
    evltyp = _geo_refloc(refloc)
    return corr, ctrcde, bodcde, orfcde, evltyp
end

@inline function _geo_state(s::AbstractVector)
    length(s) == 6 || throw(SpiceError(:badarg, "a state vector must have 6 elements"))
    @inbounds return Vec3{Float64}(s[begin], s[begin + 1], s[begin + 2]),
        Vec3{Float64}(s[begin + 3], s[begin + 4], s[begin + 5])
end

# --- public interface -------------------------------------------------------------------

"""
    spkcvo([ctx,] target, et, outref, refloc, abcorr, obssta, obsepc, obsctr, obsref)
        -> (state, lt)

Return the state of `target` relative to an observer that moves at a constant velocity
relative to the body `obsctr`, as
`spkcvo_c` of CSPICE. `state` is an `SVector{6}` (km,
km/s) in the frame `outref`; `lt` is the one-way light time in seconds.

`obssta` is the state of the observer relative to the center of `obsctr` (km, km/s, 6
elements) in the frame `obsref` at the epoch `obsepc`. The position at another epoch is
extrapolated with the constant velocity. `et` and `obsepc` are in seconds past J2000
TDB.

`refloc` tells at which epoch a rotating `outref` is evaluated when `abcorr` uses the
light time:

| `refloc` | Epoch of `outref` |
|---|---|
| `"OBSERVER"` | `et` |
| `"TARGET"` | `et` corrected for the light time to the target |
| `"CENTER"` | `et` corrected for the light time to the center of `outref` |

`abcorr` is a correction string of [`spkezr`](@ref) or an [`Aberration`](@ref).
Bodies and frames are names or integer IDs. `et` can be a `ForwardDiff.Dual`. The call
does not allocate when the bodies and the frames are integer IDs.

Throws `SpiceError`:
- `:invalidoption` if `abcorr` is not a valid correction;
- `:idcodenotfound` if a body name is not known;
- `:unknownframe` if `outref` or `obsref` is not known;
- `:notsupported` if `refloc` is not one of the three values;
- `:badarg` if `obssta` does not have 6 elements;
- the kinds of [`spkezr`](@ref) for SPK or frame data that are not loaded.

See also [`spkcpo`](@ref), [`spkcvt`](@ref), [`azlcpo`](@ref).

# References

- [NAIF documentation of `spkcvo_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcvo_c.html)
"""
function spkcvo(
        ks::KernelSet, target::NameOrID, et::Real, outref::NameOrID,
        refloc::AbstractString, abcorr::AbcorrArg, obssta::AbstractVector, obsepc::Real,
        obsctr::NameOrID, obsref::NameOrID
    )
    _chkfstr("spkcvo", "string `target`", target)       # CHKFSTR of spkcvo_c
    _chkfstr("spkcvo", "string `outref`", outref)
    _chkfstr("spkcvo", "string `refloc`", refloc)
    _chkfstr("spkcvo", "string `abcorr`", abcorr)
    _chkfstr("spkcvo", "string `obsctr`", obsctr)
    _chkfstr("spkcvo", "string `obsref`", obsref)
    corr, ctrcde, trgcde, orfcde, evltyp =
        _geo_cv_args(ks, abcorr, obsctr, "observer center", target, "target", outref, refloc)
    pos, vel = _geo_state(obssta)
    f = _GeoFixed(pos, vel, Float64(obsepc), ctrcde, _geo_known_frame(ks, obsref))
    return _geo_cv(ks, trgcde, check_epoch(et), orfcde, evltyp, corr, f)
end

"""
    spkcpo([ctx,] target, et, outref, refloc, abcorr, obspos, obsctr, obsref) -> (state, lt)

Return the state of `target` relative to an observer at a constant position on the
body `obsctr`, as
`spkcpo_c` of CSPICE. `state` is an `SVector{6}` (km, km/s) in the
frame `outref`; `lt` is the one-way light time in seconds.

`obspos` is the observer position relative to the center of `obsctr` (km, 3 elements)
in the frame `obsref`, for example a ground station in the body-fixed frame of its
planet. This is [`spkcvo`](@ref) with a zero velocity; the other arguments, the notes
and the errors are the same.

See also [`spkcpt`](@ref), [`azlcpo`](@ref).

# References

- [NAIF documentation of `spkcpo_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcpo_c.html)
"""
function spkcpo(
        ks::KernelSet, target::NameOrID, et::Real, outref::NameOrID,
        refloc::AbstractString, abcorr::AbcorrArg, obspos::AbstractVector, obsctr::NameOrID,
        obsref::NameOrID
    )
    _chkfstr("spkcpo", "string `target`", target)       # CHKFSTR of spkcpo_c
    _chkfstr("spkcpo", "string `outref`", outref)
    _chkfstr("spkcpo", "string `refloc`", refloc)
    _chkfstr("spkcpo", "string `abcorr`", abcorr)
    _chkfstr("spkcpo", "string `obsctr`", obsctr)
    _chkfstr("spkcpo", "string `obsref`", obsref)
    corr, ctrcde, trgcde, orfcde, evltyp =
        _geo_cv_args(ks, abcorr, obsctr, "observer center", target, "target", outref, refloc)
    f = _GeoFixed(_dsk_vec3(obspos), zero(Vec3{Float64}), 0.0, ctrcde, _geo_known_frame(ks, obsref))
    return _geo_cv(ks, trgcde, check_epoch(et), orfcde, evltyp, corr, f)
end

"""
    spkcvt([ctx,] trgsta, trgepc, trgctr, trgref, et, outref, refloc, abcorr, obsrvr)
        -> (state, lt)

Return the state, relative to the observer `obsrvr`, of a target that moves at a
constant velocity relative to the body `trgctr`, as
`spkcvt_c` of CSPICE. `state` is an
`SVector{6}` (km, km/s) in the frame `outref`; `lt` is the one-way light time in seconds.

`trgsta` is the state of the target relative to the center of `trgctr` (km, km/s, 6
elements) in the frame `trgref` at the epoch `trgepc`. The position at another epoch is
extrapolated with the constant velocity. `refloc`, `abcorr` and the notes on `et` and
allocation are those of [`spkcvo`](@ref).

Throws the errors of [`spkcvo`](@ref), with `trgsta` and `trgref` in the place of
`obssta` and `obsref`.

See also [`spkcpt`](@ref), [`spkcvo`](@ref).

# References

- [NAIF documentation of `spkcvt_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcvt_c.html)
"""
function spkcvt(
        ks::KernelSet, trgsta::AbstractVector, trgepc::Real, trgctr::NameOrID,
        trgref::NameOrID, et::Real, outref::NameOrID, refloc::AbstractString,
        abcorr::AbcorrArg, obsrvr::NameOrID
    )
    _chkfstr("spkcvt", "string `trgctr`", trgctr)       # CHKFSTR of spkcvt_c
    _chkfstr("spkcvt", "string `trgref`", trgref)
    _chkfstr("spkcvt", "string `outref`", outref)
    _chkfstr("spkcvt", "string `refloc`", refloc)
    _chkfstr("spkcvt", "string `abcorr`", abcorr)
    _chkfstr("spkcvt", "string `obsrvr`", obsrvr)
    corr, ctrcde, obscde, orfcde, evltyp =
        _geo_cv_args(ks, abcorr, trgctr, "target center", obsrvr, "observer", outref, refloc)
    pos, vel = _geo_state(trgsta)
    f = _GeoFixed(pos, vel, Float64(trgepc), ctrcde, _geo_known_frame(ks, trgref))
    return _geo_cv(ks, f, check_epoch(et), orfcde, evltyp, corr, obscde)
end

"""
    spkcpt([ctx,] trgpos, trgctr, trgref, et, outref, refloc, abcorr, obsrvr) -> (state, lt)

Return the state, relative to the observer `obsrvr`, of a target at a constant position
on the body `trgctr`, as
`spkcpt_c` of CSPICE. `state` is an `SVector{6}` (km, km/s) in
the frame `outref`; `lt` is the one-way light time in seconds.

`trgpos` is the target position relative to the center of `trgctr` (km, 3 elements) in
the frame `trgref`, for example a surface point in the body-fixed frame. This is
[`spkcvt`](@ref) with a zero velocity; the other arguments, the notes and the errors
are the same.

See also [`spkcpo`](@ref), [`spkcvt`](@ref).

# References

- [NAIF documentation of `spkcpt_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcpt_c.html)
"""
function spkcpt(
        ks::KernelSet, trgpos::AbstractVector, trgctr::NameOrID, trgref::NameOrID,
        et::Real, outref::NameOrID, refloc::AbstractString, abcorr::AbcorrArg,
        obsrvr::NameOrID
    )
    _chkfstr("spkcpt", "string `trgctr`", trgctr)       # CHKFSTR of spkcpt_c
    _chkfstr("spkcpt", "string `trgref`", trgref)
    _chkfstr("spkcpt", "string `outref`", outref)
    _chkfstr("spkcpt", "string `refloc`", refloc)
    _chkfstr("spkcpt", "string `abcorr`", abcorr)
    _chkfstr("spkcpt", "string `obsrvr`", obsrvr)
    corr, ctrcde, obscde, orfcde, evltyp =
        _geo_cv_args(ks, abcorr, trgctr, "target center", obsrvr, "observer", outref, refloc)
    f = _GeoFixed(_dsk_vec3(trgpos), zero(Vec3{Float64}), 0.0, ctrcde, _geo_known_frame(ks, trgref))
    return _geo_cv(ks, f, check_epoch(et), orfcde, evltyp, corr, obscde)
end

# Context and default-context forms.
spkcvo(
    ctx::SpiceContext, target::NameOrID, et::Real, outref::NameOrID, refloc::AbstractString,
    abcorr::AbcorrArg, obssta::AbstractVector, obsepc::Real, obsctr::NameOrID, obsref::NameOrID
) =
    spkcvo(kernelset(ctx), target, et, outref, refloc, abcorr, obssta, obsepc, obsctr, obsref)
spkcvo(
    target::NameOrID, et::Real, outref::NameOrID, refloc::AbstractString, abcorr::AbcorrArg,
    obssta::AbstractVector, obsepc::Real, obsctr::NameOrID, obsref::NameOrID
) =
    spkcvo(default_context(), target, et, outref, refloc, abcorr, obssta, obsepc, obsctr, obsref)
spkcpo(
    ctx::SpiceContext, target::NameOrID, et::Real, outref::NameOrID, refloc::AbstractString,
    abcorr::AbcorrArg, obspos::AbstractVector, obsctr::NameOrID, obsref::NameOrID
) =
    spkcpo(kernelset(ctx), target, et, outref, refloc, abcorr, obspos, obsctr, obsref)
spkcpo(
    target::NameOrID, et::Real, outref::NameOrID, refloc::AbstractString, abcorr::AbcorrArg,
    obspos::AbstractVector, obsctr::NameOrID, obsref::NameOrID
) =
    spkcpo(default_context(), target, et, outref, refloc, abcorr, obspos, obsctr, obsref)
spkcvt(
    ctx::SpiceContext, trgsta::AbstractVector, trgepc::Real, trgctr::NameOrID,
    trgref::NameOrID, et::Real, outref::NameOrID, refloc::AbstractString, abcorr::AbcorrArg,
    obsrvr::NameOrID
) =
    spkcvt(kernelset(ctx), trgsta, trgepc, trgctr, trgref, et, outref, refloc, abcorr, obsrvr)
spkcvt(
    trgsta::AbstractVector, trgepc::Real, trgctr::NameOrID, trgref::NameOrID, et::Real,
    outref::NameOrID, refloc::AbstractString, abcorr::AbcorrArg, obsrvr::NameOrID
) =
    spkcvt(default_context(), trgsta, trgepc, trgctr, trgref, et, outref, refloc, abcorr, obsrvr)
spkcpt(
    ctx::SpiceContext, trgpos::AbstractVector, trgctr::NameOrID, trgref::NameOrID, et::Real,
    outref::NameOrID, refloc::AbstractString, abcorr::AbcorrArg, obsrvr::NameOrID
) =
    spkcpt(kernelset(ctx), trgpos, trgctr, trgref, et, outref, refloc, abcorr, obsrvr)
spkcpt(
    trgpos::AbstractVector, trgctr::NameOrID, trgref::NameOrID, et::Real, outref::NameOrID,
    refloc::AbstractString, abcorr::AbcorrArg, obsrvr::NameOrID
) =
    spkcpt(default_context(), trgpos, trgctr, trgref, et, outref, refloc, abcorr, obsrvr)
