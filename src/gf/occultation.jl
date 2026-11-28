# Geometry finder, occultations and transits: `gfoclt` and `gfocce` of CSPICE N0067,
# with the initializer `zzgfocin` and the state `zzgfocst` of `zzgfocu`.
#
# An occultation search has no scalar quantity: it is the search for the intervals on
# which a binary state is true, so it is one call of the step-and-bracket search
# (`zzgfsolv`) for each interval of the confinement window, as `gfudb` is. The state at
# an epoch is the classification of `occult` for the roles "front" and "back"
# (`_geo_occ_state` of `geometry/occultation.jl`).

"Names of the occultation types, and the code of `_geo_occ_state` that each one accepts (0: any code but 0)."
const _GF_OCCTYPES = ("ANNULAR" => 2, "ANY" => 0, "PARTIAL" => 1, "FULL" => 3)

"""
    _GFOccultation

The checked inputs of an occultation search, the saved values of `zzgfocin`: the two
targets with their shapes and their extents (radii, bounding spheres and frames, of
`_geo_occ_extent`), the observer, the correction without stellar aberration, the code
`want` of the occultation type and the DSK tolerances.
"""
struct _GFOccultation
    ks::KernelSet
    tols::DSKTolerances
    front::_GeoOccBody
    fext::_GeoOccExtent
    back::_GeoOccBody
    bext::_GeoOccExtent
    obs::Int32
    corr::Aberration
    want::Int
end

"""
    _gf_occultation(ks, tols, occtyp, front, fshape, fframe, back, bshape, bframe,
                    abcorr, obsrvr) -> _GFOccultation

The checks of `gfocce` and `zzgfocin`, in their order: the two shapes are not both
`POINT` (`:invalidshapecombo`); the three bodies are known (`:idcodenotfound`) and
distinct (`:bodiesnotdistinct`); each shape is `POINT`, `ELLIPSOID` or a DSK method
string (`:invalidshape`, and the errors of the method string); a DSK shape needs a
point as the other target (`:invalidshapecombo`); `occtyp` is `ANY`, `FULL`, `ANNULAR`
or `PARTIAL` (`:invalidocctype`), and `ANY` if a target is a point
(`:badtypeshapecombo`); `abcorr` is a correction of `spkezr` (`:invalidoption`); an
extended target has radii or DSK data, and a body-fixed frame that is known and
centered on it (`:invalidframe`).
"""
function _gf_occultation(
        ks::KernelSet, tols::DSKTolerances, occtyp::AbstractString,
        front::NameOrID, fshape::AbstractString, fframe::NameOrID, back::NameOrID,
        bshape::AbstractString, bframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID
    )
    (_geo_occ_upper(fshape) == "POINT" && _geo_occ_upper(bshape) == "POINT") &&
        throw(
        SpiceError(
            :invalidshapecombo, "The front and back target shape " *
                "specifications are both POINT; at least one of these targets must be an " *
                "extended object."
        )
    )
    idf = _body_id(ks, front, "front target object")
    idb = _body_id(ks, back, "back target object")
    obs = _body_id(ks, obsrvr, "observer")
    (idf == idb || idf == obs || idb == obs) && throw(
        SpiceError(
            :bodiesnotdistinct,
            "The observer and both targets must be distinct objects, but are not: " *
                "OBSRVR = $obsrvr; FRONT = $front; BACK = $back."
        )
    )
    pf = _geo_occ_parse(ks, idf, fshape, "front target")
    pf.kind == 0x00 && _geo_occ_bad_shape(fshape)
    pb = _geo_occ_parse(ks, idb, bshape, "back target")
    pb.kind == 0x00 && _geo_occ_bad_shape(bshape)
    fb = _GeoOccBody(idf, pf.kind, pf.m)
    bb = _GeoOccBody(idb, pb.kind, pb.m)
    _geo_occ_dsk_combo(fb, bb) && throw(
        SpiceError(
            :invalidshapecombo, "Front target shape from FSHAPE string was <$fshape>; back " *
                "target shape from BSHAPE was <$bshape>. When one shape is DSK, the other must " *
                "be POINT."
        )
    )
    typ = _ljucase(occtyp)
    k = findfirst(p -> p.first == typ, _GF_OCCTYPES)
    k === nothing && throw(
        SpiceError(
            :invalidocctype, "The occultation type $occtyp is " *
                "not recognized.  Supported types are: ANNULAR, ANY, PARTIAL,  FULL."
        )
    )
    ((fb.kind == _GEO_OCC_POINT || bb.kind == _GEO_OCC_POINT) && typ != "ANY") &&
        throw(
        SpiceError(
            :badtypeshapecombo, "Occultation type $occtyp is not allowed " *
                "when either target body is modeled as a point. Set OCCTYP to ANY for use " *
                "with point targets."
        )
    )
    corr = _geo_lt_only(Aberration(abcorr))
    # Radii, frames and bounding spheres, the front target first. CSPICE translates the
    # surface names of the back target with the ID of the front target; the surfaces
    # here are those of the back target.
    fext = _geo_occ_extent(ks, fb, fframe, fb.m.srflst)
    bext = _geo_occ_extent(ks, bb, bframe, bb.m.srflst)
    return _GFOccultation(ks, tols, fb, fext, bb, bext, obs, corr, _GF_OCCTYPES[k].second)
end

"""
The state of `zzgfocst`: true if at `et` the front target occults the back target in
the way that the occultation type asks for.
"""
@inline function _gf_occulted(o::_GFOccultation, et::Float64)
    code = _geo_occ_state(o.ks, o.tols, o.front, o.fext, o.back, o.bext, o.obs, o.corr, et)
    return o.want == 0 ? code != 0 : code == o.want
end

"""
    _gf_oclt(ks, tols, occtyp, front, fshape, fframe, back, bshape, bframe, abcorr,
             obsrvr, cnfine, run) -> Window

`gfocce` after the checks of the step and of the tolerance: the checks of
`_gf_occultation`, then the search of each interval of `cnfine` for the intervals on
which the state `_gf_occulted` is true.
"""
function _gf_oclt(
        ks::KernelSet, tols::DSKTolerances, occtyp::AbstractString,
        front::NameOrID, fshape::AbstractString, fframe::NameOrID, back::NameOrID,
        bshape::AbstractString, bframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
        cnfine::Window, run::_GFRun
    )
    o = _gf_occultation(
        ks, tols, occtyp, front, fshape, fframe, back, bshape, bframe,
        abcorr, obsrvr
    )
    _gf_repi(run, cnfine, "Occultation/transit search ", "done.")
    result = _gf_scan!(Window(), et -> _gf_occulted(o, et), cnfine, run, false)
    _gf_repf(run)
    return result
end

"The CHKFSTR checks of `gfoclt_c` and `gfocce_c`, which check the same strings; `who` names the call."
function _gf_oclt_chkfstr(who, occtyp, front, fshape, back, bshape, abcorr, obsrvr)
    _chkfstr(who, "string `occtyp`", occtyp)       # CHKFSTR of gfoclt_c, gfocce_c
    _chkfstr(who, "string `front`", front)
    _chkfstr(who, "string `fshape`", fshape)
    _chkfstr(who, "string `back`", back)
    _chkfstr(who, "string `bshape`", bshape)
    _chkfstr(who, "string `abcorr`", abcorr)
    _chkfstr(who, "string `obsrvr`", obsrvr)
    return nothing
end

"""
    gfoclt([ctx,] occtyp, front, fshape, fframe, back, bshape, bframe, abcorr, obsrvr,
           step, cnfine; tol=gfstol(ctx), bail=nothing, report=nothing,
           threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which the body `front`
occults the body `back` as seen from `obsrvr`, as
`gfoclt_c`. A transit of `front`
across `back` is the same event.

The result is a new `Window` that is a subset of the confinement window `cnfine`;
`cnfine` is not modified. The test at each epoch is the one of [`occult`](@ref).

# Arguments
- `occtyp`: `"FULL"`, `"ANNULAR"`, `"PARTIAL"` or `"ANY"`, in any case; only `"ANY"`
  with a point target.
- `front`, `back`, `obsrvr`: names or IDs of three different bodies.
- `fshape`, `bshape`: `"ELLIPSOID"`, `"POINT"` or
  `"DSK/UNPRIORITIZED[/SURFACES = <list>]"`. At least one target must be extended, and
  a target with a DSK shape needs a point as the other target. Surface names in
  `bshape` are surfaces of `back`; `gfoclt_c` of CSPICE N0067 looks them up as
  surfaces of `front`.
- `fframe`, `bframe`: body-fixed frames of the extended targets. The frame of a point
  target is not read; a blank string is enough.
- `abcorr`: `"NONE"` or a light time correction. Stellar aberration is not applied, as
  in CSPICE: `"LT+S"` acts as `"LT"`.
- `step`: step of the search in seconds, positive; it must be shorter than the shortest
  occultation and the shortest gap between two of them.

The keywords are those of [`gfdist`](@ref). A DSK shape uses the tolerances that
[`dskstl`](@ref) set on the context. An observer on the surface of a target, such as a
ground station under `front = "EARTH"`, gives `:notdisjoint` or a wrong window; for the
rising and setting of a body, search the elevation of [`azlcpo`](@ref) with
[`gfuds`](@ref).

Throws `SpiceError(:invalidstep)` if `step <= 0` or a step does not advance the time
(CSPICE then loops for ever), `SpiceError(:invalidtolerance)` if
`tol <= 0`, and
- `SpiceError(:invalidshapecombo)` if both shapes are `"POINT"`, or if one shape is DSK
  and the other is not a point;
- `SpiceError(:idcodenotfound)` if a body name is not known, and
  `SpiceError(:bodiesnotdistinct)` if two of the three bodies are the same;
- `SpiceError(:invalidshape)` if a shape is not one of the three;
- `SpiceError(:invalidocctype)` if `occtyp` is not known, and
  `SpiceError(:badtypeshapecombo)` if it is not `"ANY"` with a point target;
- `SpiceError(:invalidoption)` if `abcorr` is not a correction;
- `SpiceError(:invalidframe)` if the frame of an extended target is not known or not
  centered on it;
- `SpiceError(:notdisjoint)` if, at an epoch of the search, the observer is inside a
  target or the targets overlap.

Also throws the kinds of [`srfrec`](@ref) if the radii of an ellipsoid are missing or
not valid, and the kinds of [`occult`](@ref) for a DSK shape. The errors of
[`spkpos`](@ref) and [`sincpt`](@ref) come from the epochs that the search samples.

# Examples
The example needs a leap seconds kernel, an SPK, and a PCK with the radii and the
orientation of the Sun and the Moon, loaded with [`furnsh`](@ref).
```julia
cnfine = Window([str2et("2025 JAN 01"), str2et("2026 JAN 01")])

# Solar eclipses of 2025 as seen from the center of the Earth, with a step of 3 minutes.
w = gfoclt("ANY", "MOON", "ELLIPSOID", "IAU_MOON", "SUN", "ELLIPSOID", "IAU_SUN", "LT",
    "EARTH", 180.0, cnfine)
for i in 1:wncard(w)
    start, stop = wnfetd(w, i)
    println(et2utc(start, "C", 0), " to ", et2utc(stop, "C", 0))
end
```

See also [`gfocce`](@ref), [`occult`](@ref), [`gfdist`](@ref).

# References

- [NAIF documentation of `gfoclt_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gfoclt_c.html)
"""
function gfoclt(
        ks::KernelSet, occtyp::AbstractString, front::NameOrID,
        fshape::AbstractString, fframe::NameOrID, back::NameOrID, bshape::AbstractString,
        bframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window;
        tol::Real = _GF_CNVTOL, bail = nothing, report = nothing, threads::Bool = false
    )
    _gf_oclt_chkfstr("gfoclt", occtyp, front, fshape, back, bshape, abcorr, obsrvr)
    return _gf_oclt(
        ks, _DSK_DEFAULT_TOLS, occtyp, front, fshape, fframe, back, bshape,
        bframe, abcorr, obsrvr, cnfine, _gf_run(step, tol; bail, report, threads)
    )
end
function gfoclt(
        ctx::SpiceContext, occtyp::AbstractString, front::NameOrID, fshape::AbstractString,
        fframe::NameOrID, back::NameOrID, bshape::AbstractString, bframe::NameOrID,
        abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window;
        tol::Real = gfstol(ctx), bail = nothing, report = nothing, threads::Bool = false
    )
    _gf_oclt_chkfstr("gfoclt", occtyp, front, fshape, back, bshape, abcorr, obsrvr)
    return _gf_oclt(
        kernelset(ctx), dsktolerances(ctx), occtyp, front, fshape, fframe, back,
        bshape, bframe, abcorr, obsrvr, cnfine, _gf_run(step, tol; bail, report, threads)
    )
end
gfoclt(
    occtyp::AbstractString, front::NameOrID, fshape::AbstractString, fframe::NameOrID,
    back::NameOrID, bshape::AbstractString, bframe::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, step::Real, cnfine::Window; kw...
) =
    gfoclt(
    default_context(), occtyp, front, fshape, fframe, back, bshape, bframe, abcorr,
    obsrvr, step, cnfine; kw...
)

"""
    gfocce([ctx,] occtyp, front, fshape, fframe, back, bshape, bframe, abcorr, obsrvr,
           tol, udstep, udrefn, cnfine; bail=nothing, report=nothing, threads=false) -> Window

Return the occultation window of [`gfoclt`](@ref), found with a step and a refinement
function that the caller supplies, as
`gfocce_c`.

`tol`, `udstep` and `udrefn` are those of [`gfevnt`](@ref). The other arguments, the
keywords `bail`, `report` and `threads`, the result and the errors are those of
[`gfoclt`](@ref). `gfocce(..., tol, step, gfrefn, cnfine)` gives the window of
`gfoclt(..., step, cnfine; tol)`.

See also [`gfoclt`](@ref), [`gfevnt`](@ref), [`gffove`](@ref).

# References

- [NAIF documentation of `gfocce_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gfocce_c.html)
"""
function gfocce(
        ks::KernelSet, occtyp::AbstractString, front::NameOrID,
        fshape::AbstractString, fframe::NameOrID, back::NameOrID, bshape::AbstractString,
        bframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, tol::Real, udstep, udrefn,
        cnfine::Window; bail = nothing, report = nothing, threads::Bool = false
    )
    _gf_oclt_chkfstr("gfocce", occtyp, front, fshape, back, bshape, abcorr, obsrvr)
    return _gf_oclt(
        ks, _DSK_DEFAULT_TOLS, occtyp, front, fshape, fframe, back, bshape,
        bframe, abcorr, obsrvr, cnfine,
        _gf_run(udstep, tol; refine = udrefn, bail, report, threads)
    )
end
function gfocce(
        ctx::SpiceContext, occtyp::AbstractString, front::NameOrID, fshape::AbstractString,
        fframe::NameOrID, back::NameOrID, bshape::AbstractString, bframe::NameOrID,
        abcorr::AbcorrArg, obsrvr::NameOrID, tol::Real, udstep, udrefn, cnfine::Window;
        bail = nothing, report = nothing, threads::Bool = false
    )
    _gf_oclt_chkfstr("gfocce", occtyp, front, fshape, back, bshape, abcorr, obsrvr)
    return _gf_oclt(
        kernelset(ctx), dsktolerances(ctx), occtyp, front, fshape, fframe, back,
        bshape, bframe, abcorr, obsrvr, cnfine,
        _gf_run(udstep, tol; refine = udrefn, bail, report, threads)
    )
end
gfocce(
    occtyp::AbstractString, front::NameOrID, fshape::AbstractString, fframe::NameOrID,
    back::NameOrID, bshape::AbstractString, bframe::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, tol::Real, udstep, udrefn, cnfine::Window; kw...
) =
    gfocce(
    default_context(), occtyp, front, fshape, fframe, back, bshape, bframe, abcorr,
    obsrvr, tol, udstep, udrefn, cnfine; kw...
)
