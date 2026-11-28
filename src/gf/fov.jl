# Geometry finder, visibility in a field of view: `gftfov`, `gfrfov` and `gffove` of
# CSPICE N0067, with the initializer `zzgffvin` and the state `zzgffvst` of `zzgffvu`.
#
# A visibility search has no scalar quantity: it is the search for the intervals on
# which a binary state is true, as `gfoclt` is. The two parts of `zzgffvu` are in
# `geometry/fov.jl`, shared with `fovtrg` and `fovray`: `_geo_fov_setup` reads the field
# of view and makes every check, once for a search, and `_geo_fov_visible` is the test
# at an epoch.

"""
    _gf_fove(ks, inst, tshape, raydir, target, tframe, abcorr, obsrvr, cnfine, run) -> Window

`gffove` after the checks of the step and of the tolerance: the checks of `zzgffvin`
(the observer, then `tshape`, which is `POINT`, `ELLIPSOID` or `RAY`
(`:invalidshape`), then those of `_geo_fov_setup`, which projects the vertices of a
polygon here), then the search of each interval of `cnfine` for the intervals on which
the state `_geo_fov_visible` is true.
"""
function _gf_fove(
        ks::KernelSet, inst::NameOrID, tshape::AbstractString,
        raydir::Vec3{Float64}, target::NameOrID, tframe::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, cnfine::Window, run::_GFRun
    )
    _body_id(ks, obsrvr, "observer")
    shp = _geo_fov_tshape(tshape)
    v = _geo_fov_setup(ks, inst, shp, raydir, target, tframe, abcorr, obsrvr, true)
    _gf_repi(run, cnfine, "Target visibility search ", "done.")
    result = _gf_scan!(Window(), et -> _geo_fov_visible(v, et), cnfine, run, false)
    _gf_repf(run)
    return result
end

@noinline _gf_ray_in_gftfov() = throw(
    SpiceError(
        :invalidoption, "The target shape RAY " *
            "is not supported by this routine. Use the routine GFRFOV instead."
    )
)

"""
    gftfov([ctx,] inst, target, tshape, tframe, abcorr, obsrvr, step, cnfine;
           tol=gfstol(ctx), bail=nothing, report=nothing, threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which the body `target` is in
the field of view of the instrument `inst`, as seen from `obsrvr`, as
`gftfov_c`.

The result is a new `Window` that is a subset of the confinement window `cnfine`;
`cnfine` is not modified. The test at each epoch is the one of [`fovtrg`](@ref).

# Arguments
- `inst`: name or ID of an instrument whose field of view is in the kernel pool (see
  [`getfov`](@ref)). The instrument frame is evaluated at the epoch under test, as in
  [`fovtrg`](@ref).
- `tshape`: `"ELLIPSOID"` or `"POINT"`, in any case. For an ellipsoid, any part of the
  reference ellipsoid in the field of view counts.
- `tframe`: body-fixed frame of an ellipsoid, centered on `target`. It is not read for a
  point.
- `abcorr`: any aberration correction of [`spkezr`](@ref).
- `step`: step of the search in seconds, positive. The search samples the visibility
  from the start of each interval of `cnfine` at this step, so `step` must be shorter
  than the shortest interval on which the target is visible, and than the shortest gap
  between two such intervals. Each change is then bisected until the bracket is not
  longer than `tol`.

The keywords are those of [`gfdist`](@ref).

Throws `SpiceError(:invalidstep)` if `step <= 0` or a step does not advance the time
(CSPICE then loops for ever), `SpiceError(:invalidtolerance)` if
`tol <= 0`, and
- `SpiceError(:invalidoption)` if `tshape` is `"RAY"` (use [`gfrfov`](@ref)) or `abcorr`
  is not a correction;
- `SpiceError(:invalidshape)` if `tshape` is not known;
- `SpiceError(:idcodenotfound)` if the observer, the target or the instrument is not
  known, and `SpiceError(:bodiesnotdistinct)` if the target is the observer;
- `SpiceError(:invalidframe)` if `tframe` is blank, not known or not centered on the
  target, and `SpiceError(:frameinfonotfound)` if the frame has no data (ellipsoid only);
- `SpiceError(:fovtoowide)` if the angular radius of the field of view about its axis
  is not below 90 degrees minus 1e-6 rad;
- `SpiceError(:unknownframe)` if the frame of the instrument is not known. This check is
  made at the first epoch that is searched.

Also throws the kinds of [`srfrec`](@ref) if the radii of an ellipsoid are missing or
not valid, the errors of [`getfov`](@ref), and the kinds of [`fovtrg`](@ref) for the
shape of the field of view and for an observer inside the target. The errors of
[`spkpos`](@ref) come from the epochs that the search samples.

See also [`gfrfov`](@ref), [`gffove`](@ref), [`fovtrg`](@ref).

# References

- [NAIF documentation of `gftfov_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gftfov_c.html)
"""
function gftfov(
        ks::KernelSet, inst::NameOrID, target::NameOrID, tshape::AbstractString,
        tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window;
        tol::Real = _GF_CNVTOL, bail = nothing, report = nothing, threads::Bool = false
    )
    _chkfstr("gftfov", "string `inst`", inst)       # CHKFSTR of gftfov_c
    _chkfstr("gftfov", "string `target`", target)
    _chkfstr("gftfov", "string `tshape`", tshape)
    _chkfstr("gftfov", "string `abcorr`", abcorr)
    _chkfstr("gftfov", "string `obsrvr`", obsrvr)
    _eqstr(tshape, "RAY") && _gf_ray_in_gftfov()
    run = _gf_run(step, tol; bail, report, threads)
    return _gf_fove(ks, inst, tshape, _GF_ZERO3, target, tframe, abcorr, obsrvr, cnfine, run)
end
gftfov(
    ctx::SpiceContext, inst::NameOrID, target::NameOrID, tshape::AbstractString,
    tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window;
    tol::Real = gfstol(ctx), kw...
) =
    gftfov(
    kernelset(ctx), inst, target, tshape, tframe, abcorr, obsrvr, step, cnfine;
    tol, kw...
)
gftfov(
    inst::NameOrID, target::NameOrID, tshape::AbstractString, tframe::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window; kw...
) =
    gftfov(
    default_context(), inst, target, tshape, tframe, abcorr, obsrvr, step, cnfine;
    kw...
)

"""
    gfrfov([ctx,] inst, raydir, rframe, abcorr, obsrvr, step, cnfine;
           tol=gfstol(ctx), bail=nothing, report=nothing, threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which the ray from `obsrvr`
along `raydir` is in the field of view of the instrument `inst`, as
`gfrfov_c`.

`raydir` is a vector of 3 elements in the frame `rframe`. `abcorr` is `"NONE"`, `"S"`
(the direction is corrected for the stellar aberration of the observer, as for light
that is received) or `"XS"` (as for light that is sent). The test at each epoch is the
one of [`fovray`](@ref). `inst`, `step`, `cnfine`, the keywords and the result are those
of [`gftfov`](@ref).

Throws the errors of [`gftfov`](@ref) for `step`, `tol`, the observer, the instrument
and its field of view, and `SpiceError` with the kind
- `:badarg` if `raydir` does not have 3 elements, `:zerovector` if it is zero;
- `:invalidframe` if `rframe` is blank or not known;
- `:invalidoption` if `abcorr` is not one of the three, such as a correction with light
  time.

See also [`gftfov`](@ref), [`gffove`](@ref), [`fovray`](@ref).

# References

- [NAIF documentation of `gfrfov_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gfrfov_c.html)
"""
function gfrfov(
        ks::KernelSet, inst::NameOrID, raydir::AbstractVector, rframe::NameOrID,
        abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window;
        tol::Real = _GF_CNVTOL, bail = nothing, report = nothing, threads::Bool = false
    )
    _chkfstr("gfrfov", "string `inst`", inst)       # CHKFSTR of gfrfov_c
    _chkfstr("gfrfov", "string `rframe`", rframe)
    _chkfstr("gfrfov", "string `abcorr`", abcorr)
    run = _gf_run(step, tol; bail, report, threads)
    return _gf_fove(ks, inst, "RAY", _dsk_vec3(raydir), 0, rframe, abcorr, obsrvr, cnfine, run)
end
gfrfov(
    ctx::SpiceContext, inst::NameOrID, raydir::AbstractVector, rframe::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window; tol::Real = gfstol(ctx),
    kw...
) =
    gfrfov(kernelset(ctx), inst, raydir, rframe, abcorr, obsrvr, step, cnfine; tol, kw...)
gfrfov(
    inst::NameOrID, raydir::AbstractVector, rframe::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, step::Real, cnfine::Window; kw...
) =
    gfrfov(default_context(), inst, raydir, rframe, abcorr, obsrvr, step, cnfine; kw...)

"""
    gffove([ctx,] inst, tshape, raydir, target, tframe, abcorr, obsrvr, tol, udstep,
           udrefn, cnfine; bail=nothing, report=nothing, threads=false) -> Window

Return the visibility window of [`gftfov`](@ref) or [`gfrfov`](@ref), found with a
step and a refinement function that the caller supplies, as
`gffove_c`.

`tshape` is `"ELLIPSOID"`, `"POINT"` or `"RAY"`. For a ray, `raydir` is the direction,
`tframe` its frame, and `target` is not read. For a body, `raydir` is not read, but it
must have 3 elements. `tol`, `udstep` and `udrefn` are those of [`gfevnt`](@ref).

The other arguments, the keywords `bail`, `report` and `threads`, the result and the
errors are those of [`gftfov`](@ref) and [`gfrfov`](@ref).
`gffove(inst, tshape, raydir, target, tframe, abcorr, obsrvr, tol, step, gfrefn, cnfine)`
gives the window of `gftfov` or `gfrfov` with the step `step` and the tolerance `tol`.

See also [`gftfov`](@ref), [`gfrfov`](@ref), [`gfevnt`](@ref).

# References

- [NAIF documentation of `gffove_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/gffove_c.html)
"""
function gffove(
        ks::KernelSet, inst::NameOrID, tshape::AbstractString,
        raydir::AbstractVector, target::NameOrID, tframe::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, tol::Real, udstep, udrefn, cnfine::Window; bail = nothing,
        report = nothing, threads::Bool = false
    )
    _chkfstr("gffove", "string `inst`", inst)       # CHKFSTR of gffove_c
    _chkfstr("gffove", "string `tshape`", tshape)
    _chkfstr("gffove", "string `abcorr`", abcorr)
    run = _gf_run(udstep, tol; refine = udrefn, bail, report, threads)
    return _gf_fove(
        ks, inst, tshape, _dsk_vec3(raydir), target, tframe, abcorr, obsrvr,
        cnfine, run
    )
end
gffove(
    ctx::SpiceContext, inst::NameOrID, tshape::AbstractString, raydir::AbstractVector,
    target::NameOrID, tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, tol::Real,
    udstep, udrefn, cnfine::Window; kw...
) =
    gffove(
    kernelset(ctx), inst, tshape, raydir, target, tframe, abcorr, obsrvr, tol,
    udstep, udrefn, cnfine; kw...
)
gffove(
    inst::NameOrID, tshape::AbstractString, raydir::AbstractVector, target::NameOrID,
    tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, tol::Real, udstep, udrefn,
    cnfine::Window; kw...
) =
    gffove(
    default_context(), inst, tshape, raydir, target, tframe, abcorr, obsrvr, tol,
    udstep, udrefn, cnfine; kw...
)
