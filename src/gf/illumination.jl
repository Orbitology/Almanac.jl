# Geometry finder, illumination angles at a surface point: the quantity `ILLUMINATION
# ANGLE` of `gfevnt` and the call `gfilum` (CSPICE N0067: `gfilum`, `zzgfilu` with
# `zzgfilin`, `zzgfilgq`, `zzgfildc`, and `zzilusta`).
#
# The value is the angle that `illumg` gives. The test "is decreasing" uses the rate of
# the angle as `zzilusta` computes it: the states of the observer and of the source
# relative to the surface point, and the state of the surface normal, all in J2000,
# with the velocity of the source and the rotation rate of the normal multiplied by the
# rate of the epoch at the point. This is the formula of CSPICE, not the exact
# derivative of the value, so that the brackets of the extrema are those of CSPICE.

"Names of the illumination angles, in the order of `zzgfilu`."
const _GF_ANGNMS = ("PHASE", "INCIDENCE", "EMISSION")

"""
    GFIllumination

An illumination angle at the surface point `spoint` of `targ` (body-fixed frame `ref`),
for the source `illum` and the observer `obs` (`zzgfilu`). `idx` is 1 for the phase
angle, 2 for the incidence angle and 3 for the emission angle; `normal` is the outward
normal of the reference ellipsoid at `spoint`.
"""
struct GFIllumination <: GFQuantity
    ks::KernelSet
    targ::Int32
    illum::Int32
    obs::Int32
    corr::Aberration
    ref::Int32
    spoint::Vec3{Float64}
    normal::Vec3{Float64}
    idx::Int
end

"""
    GFIllumination(ks, method, angtyp, target, illmn, fixref, abcorr, obsrvr, spoint)

The checks of `zzgfilin`, in its order: `:idcodenotfound`, `:bodiesnotdistinct` (the
source can be the observer), `:invalidoption` and `:notsupported` (transmission) for the
correction, `:kernelvarnotfound`, `:invalidcount` and `:badaxislength` for the radii,
`:invalidmethod`, `:unknownframe` and `:invalidframe` for `fixref`, and `:notsupported`
for `angtyp`. `spoint` must have 3 elements (`:badarg`); it is not checked to be on the
surface.
"""
function GFIllumination(
        ks::KernelSet, method::AbstractString, angtyp::AbstractString,
        target::NameOrID, illmn::NameOrID, fixref::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, spoint::AbstractVector
    )
    targ = _body_id(ks, target, "target object")
    obs = _body_id(ks, obsrvr, "observer")
    illum = _body_id(ks, illmn, "illumination source")
    _gf_distinct(targ, obs, target, obsrvr)
    targ == illum && throw(
        SpiceError(
            :bodiesnotdistinct,
            "The target and illumination source must be distinct objects, but are not: " *
                "TARGET = $target; ILLUM = $illmn."
        )
    )
    corr = Aberration(abcorr)
    corr.transmit && throw(
        SpiceError(
            :notsupported, "Aberration correction was " *
                "$(string(corr)); transmission corrections are not allowed by this routine."
        )
    )
    a, b, c = _positive_radii(ks, targ)
    sp = _dsk_vec3(spoint)
    normal = _m_surfnm(a, b, c, sp)
    _ljucase(method) == "ELLIPSOID" || throw(
        SpiceError(
            :invalidmethod, "The only supported computation method is ELLIPSOID; the input " *
                "method was $method."
        )
    )
    ref = _frame_code(ks, fixref isa AbstractString ? strip(fixref, ' ') : fixref)
    center, _, found = _geo_frinfo(ks, ref)
    found || throw(
        SpiceError(
            :unknownframe, "Reference frame $fixref is not " *
                "recognized by the SPICE frame subsystem. Possibly a required frame definition " *
                "kernel has not been loaded."
        )
    )
    center == targ || throw(
        SpiceError(
            :invalidframe, "Reference frame $fixref is not " *
                "centered at the target body $target. The ID code of the frame center is $center."
        )
    )
    idx = something(findfirst(n -> _eqstr(angtyp, n), _GF_ANGNMS), 0)
    idx == 0 && throw(
        SpiceError(
            :notsupported,
            "Illumination angle type $angtyp is not recognized."
        )
    )
    return GFIllumination(ks, targ, illum, obs, corr, ref, sp, normal, idx)
end

_gf_label(::GFIllumination) = "Illumination angle pass # of #"

"The angle of `illumg` (`zzgfilgq`)."
@inline function _gf_value(q::GFIllumination, et::Float64)
    trgepc, srfvec, tipos = _geo_illum_states(
        q.ks, q.corr, q.targ, q.obs, q.illum, q.ref,
        q.spoint, et
    )
    obspos = -srfvec
    q.idx == 1 && return Float64(vsep(obspos, tipos))
    q.idx == 2 && return Float64(vsep(q.normal, tipos))
    return Float64(vsep(q.normal, obspos))
end

"""
    _gf_ilurate(q, et) -> Float64

Rate of the angle of `q` at `et`, as `zzilusta`: `dvsep` of two of the three states in
J2000 (observer and source relative to the surface point, and the unit normal). With
light time the rotation rate of the body-fixed frame and the velocity of the source
are multiplied by `1 - dlt`, where `dlt` is the rate of the light time between the
observer and the surface point.
"""
function _gf_ilurate(q::GFIllumination, et::Float64)
    ks = q.ks
    corr = q.corr
    f = _GeoFixed(q.spoint, _GF_ZERO3, 0.0, q.targ, q.ref)
    starg, lt = _geo_cv(ks, f, et, _J2000_FRAME, _GEO_AT_TARGET, corr, q.obs)
    etsurf = et + _geo_sign(corr) * lt
    srcsta, _ = _geo_cv(ks, q.illum, etsurf, _J2000_FRAME, _GEO_AT_OBSERVER, corr, f)
    R, dR = frame_jet(ks, q.ref, _J2000_FRAME, etsurf, Val(1))
    tr, tv = _m_pv(starg)
    sr, sv = _m_pv(srcsta)
    if corr.uselt
        dlt = vdot(tv, vhat(tr)) / clight()
        dR = (1.0 - dlt) * dR
        sv = (1.0 - dlt) * sv
    end
    n = Vec3{Float64}(vhat(q.normal))
    nr, nv = rotate_jet((R, dR), (n, _GF_ZERO3))     # `mxvg_` of the state of the normal
    nrmsta = vcat(nr, nv)
    obssta = vcat(-tr, -tv)
    src = vcat(sr, sv)
    q.idx == 1 && return Float64(dvsep(obssta, src))
    q.idx == 2 && return Float64(dvsep(nrmsta, src))
    return Float64(dvsep(nrmsta, obssta))
end

"The angle decreases if its rate is negative (`zzgfildc`)."
@inline _gf_decreasing(q::GFIllumination, et::Float64) = _gf_ilurate(q, et) < 0.0

"""
    gfilum([ctx,] method, angtyp, target, illmn, fixref, abcorr, obsrvr, spoint, relate,
           refval, adjust, step, cnfine; tol=gfstol(ctx), bail=nothing, report=nothing,
           threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which an illumination angle
in radians at the surface point `spoint` of `target` satisfies the relation `relate`
with `refval`, as
`gfilum_c`. The angle is the one of [`illumg`](@ref).

# Arguments
- `method`: must be `"ELLIPSOID"`, in any case. The surface normal is the one of the
  reference ellipsoid.
- `angtyp`: `"PHASE"` (the angle at the point between the directions to the source
  `illmn` and to the observer `obsrvr`), `"INCIDENCE"` (between the normal and the
  direction to the source) or `"EMISSION"` (between the normal and the direction to the
  observer). Case and blanks are ignored.
- `target`, `illmn`, `obsrvr`: body names or IDs. The target must differ from the other
  two; the source can be the observer.
- `fixref`: body-fixed frame of `spoint`. It must be centered on `target`.
- `abcorr`: `"NONE"` or a reception correction.
- `spoint`: a vector of 3 elements in km. It is taken as a point of the surface without
  a check.

The other arguments, the keywords and the result are those of [`gfdist`](@ref).

Throws `SpiceError(:invalidstep)`, `SpiceError(:invalidtolerance)`,
`SpiceError(:notrecognized)`, `SpiceError(:valueoutofrange)` and
`SpiceError(:invalidvalue)` as [`gfdist`](@ref) does, and
- `SpiceError(:idcodenotfound)` if a body name is not known, and
  `SpiceError(:bodiesnotdistinct)` if the target is the observer or the source;
- `SpiceError(:invalidoption)` if `abcorr` is not a correction, and
  `SpiceError(:notsupported)` if it is a transmission correction;
- `SpiceError(:badarg)` if `spoint` does not have 3 elements;
- `SpiceError(:invalidmethod)` if `method` is not `"ELLIPSOID"`;
- `SpiceError(:unknownframe)` if `fixref` is not known, and `SpiceError(:invalidframe)`
  if it is not centered on `target`;
- `SpiceError(:notsupported)` if `angtyp` is not one of the three.

Also throws the kinds of [`srfrec`](@ref) if the radii of `target` are missing or not
valid. The errors of [`spkcpt`](@ref) and [`spkcpo`](@ref) come from the epochs that
the search samples.

See also [`ilumin`](@ref), [`illumg`](@ref), [`gfpa`](@ref).

# References

- [NAIF documentation of `gfilum_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gfilum_c.html)
"""
function gfilum(
        ks::KernelSet, method::AbstractString, angtyp::AbstractString,
        target::NameOrID, illmn::NameOrID, fixref::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, spoint::AbstractVector, relate::AbstractString, refval::Real,
        adjust::Real, step::Real, cnfine::Window; tol::Real = _GF_CNVTOL, bail = nothing,
        report = nothing, threads::Bool = false
    )
    _chkfstr("gfilum", "string `method`", method)       # CHKFSTR of gfilum_c
    _chkfstr("gfilum", "string `angtyp`", angtyp)
    _chkfstr("gfilum", "string `target`", target)
    _chkfstr("gfilum", "string `illmn`", illmn)
    _chkfstr("gfilum", "string `fixref`", fixref)
    _chkfstr("gfilum", "string `abcorr`", abcorr)
    _chkfstr("gfilum", "string `obsrvr`", obsrvr)
    _chkfstr("gfilum", "string `relate`", relate)
    run = _gf_run(step, tol; bail, report, threads)
    _gf_relcode(relate)
    q = GFIllumination(ks, method, angtyp, target, illmn, fixref, abcorr, obsrvr, spoint)
    return _gf_search(q, relate, refval, adjust, cnfine, run)
end
gfilum(
    ctx::SpiceContext, method::AbstractString, angtyp::AbstractString, target::NameOrID,
    illmn::NameOrID, fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
    spoint::AbstractVector, relate::AbstractString, refval::Real, adjust::Real, step::Real,
    cnfine::Window; tol::Real = gfstol(ctx), kw...
) =
    gfilum(
    kernelset(ctx), method, angtyp, target, illmn, fixref, abcorr, obsrvr, spoint,
    relate, refval, adjust, step, cnfine; tol, kw...
)
gfilum(
    method::AbstractString, angtyp::AbstractString, target::NameOrID, illmn::NameOrID,
    fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, spoint::AbstractVector,
    relate::AbstractString, refval::Real, adjust::Real, step::Real, cnfine::Window; kw...
) =
    gfilum(
    default_context(), method, angtyp, target, illmn, fixref, abcorr, obsrvr, spoint,
    relate, refval, adjust, step, cnfine; kw...
)
