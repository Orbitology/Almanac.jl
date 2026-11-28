# Geometry layer, shared part: the computation method strings (`zzprsmet`), the
# aberration correction without stellar aberration, the lookup of bodies, frames
# and radii with the error kinds of the geometry calls of CSPICE, and the target shape
# (`zzsinutl`: `zzsuelin`, `zzsudski`, `zzraysfx`, `zzraynp`, `zzmaxrad`).
#
# Every function is a pure function of the `KernelSet` snapshot and of its arguments.
# The private names of the geometry layer start with `_geo_`, `_Geo` or `_GEO_`.

# --- method strings ---------------------------------------------------------------------

const _GEO_NEARPOINT = 0x01     # sub-point type "NEAR POINT" (legacy ellipsoid forms)
const _GEO_NADIR = 0x02         # sub-point type "NADIR"
const _GEO_INTERCEPT = 0x03     # sub-point type "INTERCEPT"
const _GEO_TANGENT = 0x01       # limb or terminator point definition "TANGENT"
const _GEO_GUIDED = 0x02        # limb or terminator point definition "GUIDED"
const _GEO_UMBRAL = 0x01        # terminator type "UMBRAL"
const _GEO_PENUMBRAL = 0x02     # terminator type "PENUMBRAL"

"The surface list of a method without the `SURFACES` keyword (never modified)."
const _GEO_NOSURF = Int[]

"""
    _GeoMethod

A parsed method string, the outputs of `zzprsmet`: `dsk` (the shape is `DSK`, else
`ELLIPSOID`), `subtyp` (0 or one of `_GEO_NEARPOINT`, `_GEO_NADIR`, `_GEO_INTERCEPT`),
`pntdef` (0, `_GEO_TANGENT`, `_GEO_GUIDED`), `trmtyp` (0, `_GEO_UMBRAL`,
`_GEO_PENUMBRAL`) and the surface ID list `srflst` (empty: all surfaces).
"""
struct _GeoMethod
    dsk::Bool
    subtyp::UInt8
    pntdef::UInt8
    trmtyp::UInt8
    srflst::Vector{Int}
end

@inline _geo_ell(subtyp::UInt8, pntdef::UInt8 = 0x00, trmtyp::UInt8 = 0x00) =
    _GeoMethod(false, subtyp, pntdef, trmtyp, _GEO_NOSURF)
@inline _geo_dsk(subtyp::UInt8, pntdef::UInt8 = 0x00, trmtyp::UInt8 = 0x00) =
    _GeoMethod(true, subtyp, pntdef, trmtyp, _GEO_NOSURF)

"""
    _geo_method(ks, bodyid, method) -> _GeoMethod

`zzprsmet`: parse the method string `method` for the body `bodyid`. The common
spellings are compared as they are and do not allocate; every other string goes through
`_dsk_parse_method`. Throws `SpiceError` (`:badmethodsyntax`, `:badpriorityspec`,
`:idcodenotfound`, `:toomanysurfaces`, `:syntaxerror`).
"""
@inline function _geo_method(ks::KernelSet, bodyid::Integer, method::AbstractString)
    method == "ELLIPSOID" && return _geo_ell(0x00)
    method == "NEAR POINT/ELLIPSOID" && return _geo_ell(_GEO_NEARPOINT)
    method == "INTERCEPT/ELLIPSOID" && return _geo_ell(_GEO_INTERCEPT)
    method == "ELLIPSOID/NEAR POINT" && return _geo_ell(_GEO_NEARPOINT)
    method == "ELLIPSOID/INTERCEPT" && return _geo_ell(_GEO_INTERCEPT)
    method == "NEAR POINT: ELLIPSOID" && return _geo_ell(_GEO_NEARPOINT)
    method == "INTERCEPT: ELLIPSOID" && return _geo_ell(_GEO_INTERCEPT)
    method == "Near point: ellipsoid" && return _geo_ell(_GEO_NEARPOINT)
    method == "Intercept: ellipsoid" && return _geo_ell(_GEO_INTERCEPT)
    method == "NADIR/ELLIPSOID" && return _geo_ell(_GEO_NADIR)
    method == "TANGENT/ELLIPSOID" && return _geo_ell(0x00, _GEO_TANGENT)
    method == "GUIDED/ELLIPSOID" && return _geo_ell(0x00, _GEO_GUIDED)
    method == "ELLIPSOID/TANGENT" && return _geo_ell(0x00, _GEO_TANGENT)
    method == "ELLIPSOID/GUIDED" && return _geo_ell(0x00, _GEO_GUIDED)
    method == "UMBRAL/TANGENT/ELLIPSOID" && return _geo_ell(0x00, _GEO_TANGENT, _GEO_UMBRAL)
    method == "PENUMBRAL/TANGENT/ELLIPSOID" && return _geo_ell(0x00, _GEO_TANGENT, _GEO_PENUMBRAL)
    method == "UMBRAL/GUIDED/ELLIPSOID" && return _geo_ell(0x00, _GEO_GUIDED, _GEO_UMBRAL)
    method == "PENUMBRAL/GUIDED/ELLIPSOID" && return _geo_ell(0x00, _GEO_GUIDED, _GEO_PENUMBRAL)
    method == "DSK/UNPRIORITIZED" && return _geo_dsk(0x00)
    method == "NADIR/DSK/UNPRIORITIZED" && return _geo_dsk(_GEO_NADIR)
    method == "INTERCEPT/DSK/UNPRIORITIZED" && return _geo_dsk(_GEO_INTERCEPT)
    method == "TANGENT/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_TANGENT)
    method == "GUIDED/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_GUIDED)
    method == "UMBRAL/TANGENT/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_TANGENT, _GEO_UMBRAL)
    method == "PENUMBRAL/TANGENT/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_TANGENT, _GEO_PENUMBRAL)
    method == "UMBRAL/GUIDED/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_GUIDED, _GEO_UMBRAL)
    method == "PENUMBRAL/GUIDED/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_GUIDED, _GEO_PENUMBRAL)
    return _geo_method_slow(ks, bodyid, method)
end

@noinline function _geo_method_slow(ks::KernelSet, bodyid::Integer, method::AbstractString)
    shape, subtyp, _, srflst = _dsk_parse_method(ks, bodyid, method)
    sub = subtyp == "" ? 0x00 : subtyp == "NEAR POINT" ? _GEO_NEARPOINT :
        subtyp == "NADIR" ? _GEO_NADIR : _GEO_INTERCEPT
    pntdef = 0x00
    trmtyp = 0x00
    # The parser has accepted the string, so each of these keywords occurs one time at
    # most, and `zzprsmet` takes a keyword only where it follows a '/'. The legacy
    # ellipsoid forms return before the keywords are read.
    legacy = shape === :ellipsoid && (
        sub == _GEO_NEARPOINT || (
            sub == 0x00 &&
                _eqstr(method, "ELLIPSOID")
        ) || (
            sub == _GEO_INTERCEPT && (
                _eqstr(method, "INTERCEPT:ELLIPSOID") || _eqstr(method, "ELLIPSOID:INTERCEPT") ||
                    _eqstr(method, "INTERCEPT/ELLIPSOID") || _eqstr(method, "ELLIPSOID/INTERCEPT")
            )
        )
    )
    if !legacy
        tok = pushfirst!(_dsk_lex_method(normalize_name(method)), "/")
        for i in 2:length(tok)           # invariant: the keywords of tok[2:i-1] are recorded
            tok[i - 1] == "/" || continue
            t = tok[i]
            t == "TANGENT" && (pntdef = _GEO_TANGENT)
            t == "GUIDED" && (pntdef = _GEO_GUIDED)
            t == "UMBRAL" && (trmtyp = _GEO_UMBRAL)
            t == "PENUMBRAL" && (trmtyp = _GEO_PENUMBRAL)
        end
    end
    return _GeoMethod(shape === :dsk, sub, pntdef, trmtyp, isempty(srflst) ? _GEO_NOSURF : srflst)
end

# --- aberration corrections, bodies, frames, radii --------------------------------------

"The correction `c` without stellar aberration (`loccor` of `zzsfxcor`, the correction of `zzgfocin`)."
@inline _geo_lt_only(c::Aberration) = _aberration(c.flags & ~_AB_S)

"`s` of the geometry calls: -1 for reception, +1 for transmission, 0 without light time."
@inline _geo_sign(c::Aberration) = c.uselt ? (c.transmit ? 1.0 : -1.0) : 0.0

"`frinfo` without the class ID: `(center, class, found)`."
@inline function _geo_frinfo(ks::KernelSet, id::Int32)
    f = _frame_def_or_nothing(ks, id)
    f === nothing && return Int32(0), Int32(0), false
    return f.center, f.class, true
end

"The body-fixed frame of the geometry calls: `_dsk_fixfid` with `:noframe` for an unknown frame."
@inline _geo_fixfid(ks::KernelSet, fixref::NameOrID, trgcde::Int32, target) =
    _dsk_fixfid(ks, fixref, trgcde, target, :noframe, :noframe)

# --- positions --------------------------------------------------------------------------

"`spkezp_`: position of `targ` relative to `obs` in frame `ref` and the light time."
@inline function _geo_spkezp(
        ks::KernelSet, targ::Int32, et, ref::Int32, corr::Aberration,
        obs::Int32
    )
    if corr.geometric || targ == obs
        (r,) = spk_jet(ks, targ, et, ref, obs, Val(0))
        return r, _m_norm(r) / CLIGHT
    end
    return abcorr_position(ks, targ, et, ref, corr, obs)
end

# --- the target shape -------------------------------------------------------------------

"""
    _GeoShape

The target shape of a geometry call, the saved values of `zzsinutl`: `dsk` (DSK data,
else the ellipsoid with the radii `a`, `b`, `c`), the target `trgcde`, its body-fixed
frame `fixfid`, the radius `maxrad` of the outer bounding sphere, the surface list and
the DSK tolerances. For DSK data, `segs` is the list of loaded segments of the target
(`ks.dsk[trgcde]` of the snapshot that built the shape), and `segok` tells whether the
checks of `_dsk_body_segments` pass on it; the ray intercepts of `_geo_rayx` use both,
so the segments are looked up once per shape, not once per ray.
"""
struct _GeoShape
    dsk::Bool
    trgcde::Int32
    fixfid::Int32
    a::Float64
    b::Float64
    c::Float64
    maxrad::Float64
    srflst::Vector{Int}
    tols::DSKTolerances
    segs::Vector{DSKSegment}
    segok::Bool
end

"The segment list of a shape without DSK data (never modified)."
const _GEO_NOSEGS = DSKSegment[]

# A DSK shape built without its segments looks them up at each ray (`segok` false).
_GeoShape(dsk, trgcde, fixfid, a, b, c, maxrad, srflst, tols) =
    _GeoShape(dsk, trgcde, fixfid, a, b, c, maxrad, srflst, tols, _GEO_NOSEGS, !dsk)

"`zzsuelin`: the ellipsoid of the body `trgcde`."
@inline function _geo_ellshape(ks::KernelSet, trgcde::Int32, fixfid::Int32)
    a, b, c = _positive_radii(ks, trgcde)
    return _GeoShape(
        false, trgcde, fixfid, a, b, c, max(max(a, b), c),
        _GEO_NOSURF, _DSK_DEFAULT_TOLS
    )
end

"`zzsudski`: the DSK surface of the body `trgcde` with the surfaces `srflst`."
function _geo_dskshape(ks::KernelSet, tols::DSKTolerances, trgcde::Int32, fixfid::Int32, srflst)
    _, maxrad = _dsk_body_sph(ks, trgcde, srflst)
    # `_dsk_body_sph` has thrown unless a DSK file is loaded and the target has segments.
    segs = ks.dsk[trgcde]
    sh = _GeoShape(
        true, trgcde, fixfid, 0.0, 0.0, 0.0, maxrad, srflst, tols, segs,
        _dsk_segments_ok(ks, segs)
    )
    @assert !isempty(sh.segs) && sh.maxrad >= 0.0
    return sh
end

@inline _geo_shape(ks::KernelSet, tols::DSKTolerances, m::_GeoMethod, trgcde::Int32, fixfid::Int32) =
    m.dsk ? _geo_dskshape(ks, tols, trgcde, fixfid, m.srflst) : _geo_ellshape(ks, trgcde, fixfid)

@noinline _geo_dsk_epoch() = throw(
    SpiceError(
        :notsupported,
        "a DSK shape needs a plain (Float64) epoch: the plate model is not differentiated"
    )
)

@inline _geo_plain(v::Vec3{Float64}) = v
@inline _geo_plain(x::Float64) = x
_geo_plain(::Any) = _geo_dsk_epoch()

"""
    _geo_rayx(ks, sh, vertex, raydir, et) -> (found, spoint)

`zzraysfx`: the first point where the ray from `vertex` along `raydir` (body-fixed
frame of the shape, epoch `et` for the DSK segment selection) meets the surface.
Generic over the number type for an ellipsoid; no allocation. `ks` must be the
snapshot that built `sh`. The result is that of `_dsk_bux` for the same ray.
"""
@inline function _geo_rayx(
        ks::KernelSet, sh::_GeoShape, vertex::SVector{3, T},
        raydir::SVector{3, T}, et
    ) where {T}
    if sh.dsk
        # The errors of `zzdsksba` come with the first ray, as in `_dsk_bux`.
        sh.segok || _dsk_body_segments(ks, sh.trgcde)
        xpt, _, _, found = _dsk_bux_ray(
            ks, sh.tols, sh.segs, sh.maxrad, sh.trgcde, sh.srflst, _geo_plain(et),
            sh.fixfid, _geo_plain(vertex), _geo_plain(raydir)
        )
        return found, SVector{3, T}(xpt)
    end
    return _m_surfpt(vertex, raydir, sh.a, sh.b, sh.c)
end

"`zzraynp`: the point of the ray's line nearest to the ellipsoid (or to the outer sphere of DSK data) and the distance."
@inline function _geo_raynp(sh::_GeoShape, vertex::SVector{3}, raydir::SVector{3})
    sh.dsk && return _m_npedln(sh.maxrad, sh.maxrad, sh.maxrad, vertex, raydir, false)
    return _m_npedln(sh.a, sh.b, sh.c, vertex, raydir, false)
end
