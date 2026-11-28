# Illumination angles at a surface point, the phase angle between bodies and the angular
# separation of two bodies: `ilumin`, `illumg`, `illumf`, `phaseq` and `trgsep` of CSPICE
# N0067 (with `zzgfpaq`, and `zzspin` and `zzsepq` for `trgsep`).
#
# The functions are generic over the epoch number type (a dual number passes through;
# the angles need the arcsine, which the `Taylor` type does not have).

"""
    _geo_illum_states(ks, corr, trgcde, obscde, srccde, fixfid, spoint, et)
        -> (trgepc, srfvec, tipos)

The two states of `illumg`: the surface point as seen from the observer (`spkcpt` with
the locus `TARGET`), which gives `srfvec` and the epoch `trgepc` at the point, and the
illumination source as seen from the surface point at `trgepc` (`spkcpo` with the
locus `OBSERVER`), which gives `tipos`.
"""
@inline function _geo_illum_states(
        ks::KernelSet, corr::Aberration, trgcde::Int32,
        obscde::Int32, srccde::Int32, fixfid::Int32, spoint::Vec3{Float64}, et
    )
    f = _GeoFixed(spoint, zero(Vec3{Float64}), 0.0, trgcde, fixfid)
    opstat, lt = _geo_cv(ks, f, et, fixfid, _GEO_AT_TARGET, corr, obscde)
    trgepc = et + _geo_sign(corr) * lt
    tistat, _ = _geo_cv(ks, srccde, trgepc, fixfid, _GEO_AT_OBSERVER, corr, f)
    return trgepc, SVector(opstat[1], opstat[2], opstat[3]), SVector(tistat[1], tistat[2], tistat[3])
end

"""
    _geo_illum_prep(ks, method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
        -> (m, trgcde, fixfid, sp, e, trgepc, srfvec, tipos)

The checks that `illumg` and `illumf` share, in the order of CSPICE, the checked epoch
`e`, and the states of `_geo_illum_states`.
"""
@inline function _geo_illum_prep(
        ks::KernelSet, method::AbstractString, target::NameOrID, ilusrc::NameOrID,
        et::Real, fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
        spoint::AbstractVector
    )
    corr = Aberration(abcorr)
    trgcde = _body_id(ks, target, "target")
    m = _geo_method(ks, trgcde, method)
    m.subtyp == 0x00 || _geo_sub_in_method(method)
    obscde = _body_id(ks, obsrvr, "observer")
    obscde == trgcde && _geo_same_bodies("illumination angles", obsrvr)
    fixfid = _geo_fixfid(ks, fixref, trgcde, target)
    sp = _dsk_vec3(spoint)
    e = check_epoch(et)
    srccde = _body_id(ks, ilusrc, "illumination source")
    trgepc, srfvec, tipos = _geo_illum_states(ks, corr, trgcde, obscde, srccde, fixfid, sp, e)
    return m, trgcde, fixfid, sp, e, trgepc, srfvec, tipos
end

"The phase, incidence and emission angles for the surface normal `normal`."
@inline _geo_illum_angles(normal, obspos, tipos) =
    (vsep(obspos, tipos), vsep(normal, tipos), vsep(normal, obspos))

@inline function _geo_ellipsoid_normal(ks::KernelSet, trgcde::Int32, sp::Vec3{Float64})
    a, b, c = _positive_radii(ks, trgcde)
    return _m_surfnm(a, b, c, sp)
end

function _geo_illumg(
        ks::KernelSet, tols::DSKTolerances, method::AbstractString,
        target::NameOrID, ilusrc::NameOrID, et::Real, fixref::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, spoint::AbstractVector
    )
    _chkfstr("illumg/ilumin", "string `method`", method)       # CHKFSTR of illumg_c and ilumin_c
    _chkfstr("illumg/ilumin", "string `target`", target)
    _chkfstr("illumg/ilumin", "string `ilusrc`", ilusrc)
    _chkfstr("illumg/ilumin", "string `fixref`", fixref)
    _chkfstr("illumg/ilumin", "string `abcorr`", abcorr)
    _chkfstr("illumg/ilumin", "string `obsrvr`", obsrvr)
    m, trgcde, fixfid, sp, e, trgepc, srfvec, tipos =
        _geo_illum_prep(ks, method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
    normal = if m.dsk
        _dsk_bun(ks, tols, trgcde, m.srflst, _geo_plain(e), fixfid, sp)
    else
        _geo_ellipsoid_normal(ks, trgcde, sp)
    end
    return trgepc, srfvec, _geo_illum_angles(normal, -srfvec, tipos)...
end

function _geo_illumf(
        ks::KernelSet, tols::DSKTolerances, method::AbstractString,
        target::NameOrID, ilusrc::NameOrID, et::Real, fixref::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, spoint::AbstractVector
    )
    _chkfstr("illumf", "string `method`", method)       # CHKFSTR of illumf_c
    _chkfstr("illumf", "string `target`", target)
    _chkfstr("illumf", "string `ilusrc`", ilusrc)
    _chkfstr("illumf", "string `fixref`", fixref)
    _chkfstr("illumf", "string `abcorr`", abcorr)
    _chkfstr("illumf", "string `obsrvr`", obsrvr)
    m, trgcde, fixfid, sp, e, trgepc, srfvec, tipos =
        _geo_illum_prep(ks, method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
    obspos = -srfvec
    if !m.dsk
        phase, incdnc, emissn = _geo_illum_angles(_geo_ellipsoid_normal(ks, trgcde, sp), obspos, tipos)
        return trgepc, srfvec, phase, incdnc, emissn, primal(emissn) <= _HALFPI,
            primal(incdnc) <= _HALFPI
    end
    # As in CSPICE, the DSK branch of `illumf` uses every surface of the target (the
    # surface list of the method string is not used) and the epoch at the point.
    sh = _geo_dskshape(ks, tols, trgcde, fixfid, _GEO_NOSURF)
    tp = _geo_plain(trgepc)
    normal = _m_hat(_dsk_bun(ks, tols, trgcde, _GEO_NOSURF, tp, fixfid, sp))
    phase, incdnc, emissn = _geo_illum_angles(normal, obspos, tipos)
    visibl = emissn <= _HALFPI
    lit = incdnc <= _HALFPI
    if lit || visibl
        # A point is hidden or in shadow if the ray from just above it meets the surface.
        vertex = sp + (sh.maxrad * 1.0e-10) * normal
        visibl && (visibl = !_geo_rayx(ks, sh, vertex, obspos, tp)[1])
        lit && (lit = !_geo_rayx(ks, sh, vertex, tipos, tp)[1])
    end
    return trgepc, srfvec, phase, incdnc, emissn, visibl, lit
end

"""
    illumg([ctx,] method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
        -> (trgepc, srfvec, phase, incdnc, emissn)

Return the illumination angles at the surface point `spoint` of `target`, for the
illumination source `ilusrc` and the observer `obsrvr` at `et` (seconds past J2000
TDB), as
`illumg_c` of CSPICE.

- `trgepc`: the epoch at the point (seconds past J2000 TDB). With a light time
  correction it is `et` minus (reception) or plus (transmission) the light time between
  the observer and the point. Without one it is `et`.
- `srfvec`: the vector from the observer to the point, in km, in the body-fixed frame
  `fixref` at `trgepc`.
- `phase`: the angle at the point between the directions to the source and to the
  observer (rad).
- `incdnc`: the angle between the outward surface normal and the direction to the
  source (rad).
- `emissn`: the angle between the outward surface normal and the direction to the
  observer (rad).

`spoint` is in km in `fixref` (3 elements). The call does not check that it is on the
surface. `method` selects the surface normal: `"ELLIPSOID"` or `"DSK/UNPRIORITIZED"`;
the DSK form can end with `"/SURFACES = <list>"`. `abcorr` is a correction string of
[`spkezr`](@ref) or an [`Aberration`](@ref). With a correction, the source is taken
where it appears from the point at `trgepc`.

With `"ELLIPSOID"`, `et` can be a `ForwardDiff.Dual`, and the call does not allocate
when the bodies and the frame are integer IDs. A DSK method uses the tolerances that
[`dskstl`](@ref) set on the context.

Throws `SpiceError` with the kind
- `:invalidoption` if `abcorr` is not a valid correction;
- `:idcodenotfound` if a body name is not known;
- `:badmethodsyntax` if `method` is malformed, `:badpriorityspec` if a DSK method does
  not have `UNPRIORITIZED`, `:invalidmethod` if it holds a sub-point type;
- `:bodiesnotdistinct` if `target` and `obsrvr` are the same body; `ilusrc` can be
  either of them (`ilusrc == target` puts the source at the target center, so
  `incdnc > π/2` on the ellipsoid);
- `:noframe` if `fixref` is not known, `:invalidframe` if it is not centered on `target`;
- `:badarg` if `spoint` does not have 3 elements;
- the kinds of [`srfrec`](@ref) if the radii of `target` are missing or not valid
  (ellipsoid).

Also throws the errors of [`spkcpt`](@ref) and, with a DSK method, of [`srfnrm`](@ref).

See also [`ilumin`](@ref), [`illumf`](@ref), [`gfilum`](@ref).

# References

- [NAIF documentation of `illumg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/illumg_c.html)
"""
function illumg end

"""
    illumf([ctx,] method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
        -> (trgepc, srfvec, phase, incdnc, emissn, visibl, lit)

Return the outputs of [`illumg`](@ref) and two flags, as
`illumf_c` of CSPICE: `visibl`
tells whether the observer sees the point, and `lit` whether the source lights it.

With `"ELLIPSOID"`, `visibl` is `emissn <= π/2` and `lit` is `incdnc <= π/2`. With a
DSK method, a flag is also false if the surface blocks the ray from the point toward
the observer or toward the source.

With a DSK method, the angles and the flags use all loaded surfaces of the target: the
`SURFACES` list of `method` is ignored. The DSK segments are those that cover
`trgepc`; [`illumg`](@ref) takes those that cover `et`. Both are as in CSPICE.

The arguments, the notes on `et`, allocation and tolerances, and the errors are those of
[`illumg`](@ref), except `:dskdatanotfound` if DSK files are loaded but none has data
for `target`.

See also [`ilumin`](@ref), [`illumg`](@ref).

# References

- [NAIF documentation of `illumf_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/illumf_c.html)
"""
function illumf end

for (f, g) in ((:illumg, :_geo_illumg), (:illumf, :_geo_illumf))
    @eval begin
        $f(
            ks::KernelSet, method::AbstractString, target::NameOrID, ilusrc::NameOrID, et::Real,
            fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, spoint::AbstractVector
        ) =
            $g(ks, _DSK_DEFAULT_TOLS, method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
        $f(
            ctx::SpiceContext, method::AbstractString, target::NameOrID, ilusrc::NameOrID,
            et::Real, fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
            spoint::AbstractVector
        ) = $g(
            kernelset(ctx), dsktolerances(ctx), method, target,
            ilusrc, et, fixref, abcorr, obsrvr, spoint
        )
        $f(
            method::AbstractString, target::NameOrID, ilusrc::NameOrID, et::Real,
            fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, spoint::AbstractVector
        ) =
            $f(default_context(), method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
    end
end

"""
    ilumin([ctx,] method, target, et, fixref, abcorr, obsrvr, spoint)
        -> (trgepc, srfvec, phase, incdnc, emissn)

Return the illumination angles at the surface point `spoint` of `target`, with the Sun
as the illumination source, as
`ilumin_c` of CSPICE.

`phase`, `incdnc` and `emissn` are the phase angle, the solar incidence angle and the
emission angle (rad). `trgepc` is the epoch at the point (seconds past J2000 TDB) and
`srfvec` the vector from the observer to the point (km, body-fixed frame `fixref`).
This is [`illumg`](@ref)
with `ilusrc = 10`; see there for the arguments, the notes and the errors.

# Examples
```julia
# Needs an SPK, a PCK with the radii and the rotation of Mars, and a leap seconds kernel.
et = str2et("2025 MAR 01 12:00:00 UTC")
spoint, = subpnt("NEAR POINT/ELLIPSOID", "MARS", et, "IAU_MARS", "LT+S", "EARTH")
trgepc, srfvec, phase, incdnc, emissn = ilumin("ELLIPSOID", "MARS", et, "IAU_MARS", "LT+S", "EARTH", spoint)
daylight = incdnc <= halfpi()          # the Sun is above the horizon of the point
```

See also [`illumg`](@ref), [`illumf`](@ref), [`phaseq`](@ref), [`subslr`](@ref).

# References

- [NAIF documentation of `ilumin_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ilumin_c.html)
"""
ilumin(
    ks::KernelSet, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, spoint::AbstractVector
) =
    _geo_illumg(ks, _DSK_DEFAULT_TOLS, method, target, 10, et, fixref, abcorr, obsrvr, spoint)
ilumin(
    ctx::SpiceContext, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, spoint::AbstractVector
) =
    _geo_illumg(
    kernelset(ctx), dsktolerances(ctx), method, target, 10, et, fixref, abcorr,
    obsrvr, spoint
)
ilumin(
    method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, spoint::AbstractVector
) =
    ilumin(default_context(), method, target, et, fixref, abcorr, obsrvr, spoint)

"""
    phaseq([ctx,] et, target, illmn, obsrvr, abcorr) -> phase

Return the phase angle at `target` at `et` (seconds past J2000 TDB), in radians in
[0, π], as
`phaseq_c` of CSPICE. It is the angle between the directions from the target
to the source `illmn` and to the observer `obsrvr`.

With a correction, the target is taken as seen by the observer, and the source as seen
by the target at `et - lt`, where `lt` is the light time between the observer and the
target. Bodies are names or integer IDs. `et` can be a `ForwardDiff.Dual`. The call
does not allocate when the bodies are integer IDs.

Throws `SpiceError(:invalidoption)` if `abcorr` is not a valid correction or is a
transmission correction, `SpiceError(:bodiesnotdistinct)` if two of the three bodies
are the same, and `SpiceError(:idcodenotfound)` if a body name is not known. Also
throws the errors of [`spkpos`](@ref).

See also [`ilumin`](@ref), [`gfpa`](@ref).

# References

- [NAIF documentation of `phaseq_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/phaseq_c.html)
"""
function phaseq(
        ks::KernelSet, et::Real, target::NameOrID, illmn::NameOrID, obsrvr::NameOrID,
        abcorr::AbcorrArg
    )
    _chkfstr("phaseq", "string `target`", target)       # CHKFSTR of phaseq_c
    _chkfstr("phaseq", "string `illmn`", illmn)
    _chkfstr("phaseq", "string `obsrvr`", obsrvr)
    _chkfstr("phaseq", "string `abcorr`", abcorr)
    targ = _body_id(ks, target, "target")
    illum = _body_id(ks, illmn, "illumination source")
    obs = _body_id(ks, obsrvr, "observer")
    corr = Aberration(abcorr)
    corr.transmit && throw(
        SpiceError(
            :invalidoption,
            "the aberration correction $(string(corr)) calls for transmission; phaseq " *
                "accepts reception corrections only"
        )
    )
    (targ == obs || targ == illum || obs == illum) && throw(
        SpiceError(
            :bodiesnotdistinct,
            "the target $target, the illumination source $illmn and the observer $obsrvr " *
                "must be three different bodies"
        )
    )
    return _geo_phase(ks, targ, illum, obs, corr, check_epoch(et))
end

"""
    _geo_phase(ks, targ, illum, obs, corr, et)

The phase angle as `zzgfpaq` computes it, for three distinct bodies: pi minus the angle
between the position of `targ` seen from `obs` and the position of `illum` seen from
`targ` at the epoch `et - lt` (at `et` without correction), in J2000.
"""
@inline function _geo_phase(
        ks::KernelSet, targ::Int32, illum::Int32, obs::Int32, corr::Aberration, et
    )
    pv1, lt = _geo_spkezp(ks, targ, et, _J2000_FRAME, corr, obs)
    pv2, _ = _geo_spkezp(ks, illum, corr.geometric ? et : et - lt, _J2000_FRAME, corr, targ)
    return pi - vsep(pv1, pv2)
end
phaseq(
    ctx::SpiceContext, et::Real, target::NameOrID, illmn::NameOrID, obsrvr::NameOrID,
    abcorr::AbcorrArg
) = phaseq(kernelset(ctx), et, target, illmn, obsrvr, abcorr)
phaseq(et::Real, target::NameOrID, illmn::NameOrID, obsrvr::NameOrID, abcorr::AbcorrArg) =
    phaseq(default_context(), et, target, illmn, obsrvr, abcorr)

"Throws `SpiceError(:emptystring)` for an empty string argument, as the C interface of `trgsep_c`."
@inline _geo_nonempty(s::AbstractString, what::String) =
    isempty(s) && throw(SpiceError(:emptystring, "trgsep: the $what string is empty"))
@inline _geo_nonempty(::Any, ::String) = false

"""
    _geo_sepshape(s) -> UInt8

The shape of `zzspin`: 1 for `POINT`, 2 for `SPHERE`, 0 for another string. As
`ljucrs` and `isrchc`, leading blanks are removed, the case is ignored and only the
first six characters count, so `"SPHERES"` is `SPHERE`. Does not allocate.
"""
function _geo_sepshape(s::AbstractString)
    point = true
    sphere = true
    k = 0
    for c in s     # invariant: point and sphere tell whether the k characters match
        (k == 0 && c == ' ') && continue
        k += 1
        k > 6 && break
        u = uppercase(c)
        point &= k <= 5 ? u == "POINT"[k] : u == ' '
        sphere &= u == "SPHERE"[k]
    end
    point && k >= 5 && return 0x01
    sphere && k >= 6 && return 0x02
    return 0x00
end

"Radius of the shape `shape` of `body` for `trgsep`: 0 for a point, the largest radius for a sphere (`zzspin`)."
function _geo_sepradius(ks::KernelSet, body::Int32, shape::AbstractString)
    sh = _geo_sepshape(shape)
    sh == 0x01 && return 0.0
    sh == 0x02 && return maximum(_positive_radii(ks, body))
    throw(
        SpiceError(
            :notrecognized,
            "the body shape $shape is not recognized; the supported shapes are POINT and SPHERE"
        )
    )
end

"Angular radius of a sphere of radius `r` at the distance `range` (`zzsepq`)."
@inline function _geo_angrad(r::Float64, range, obs, body, which)
    primal(range) > r || throw(
        SpiceError(
            :insidebody, "the observer $obs is inside the " *
                "surface of target $which, $body: range $(primal(range)) km, radius $r km"
        )
    )
    return r == 0.0 ? zero(range) : asin(r / range)
end

"""
    trgsep([ctx,] et, targ1, shape1, frame1, targ2, shape2, frame2, obsrvr, abcorr) -> angle

Return the angular separation (rad) of two targets as seen by the observer `obsrvr` at
`et` (seconds past J2000 TDB), as
`trgsep_c` of CSPICE.

Each target is a point (`shape = "POINT"`) or a sphere (`"SPHERE"`) whose radius is the
largest of the three radii of the body in the kernel pool. For spheres the angle is
the separation of the limbs: the angle between the centers minus the two angular radii,
so it is negative when the disks overlap. The frames `frame1` and `frame2` are not used
by these shapes, as in CSPICE. Bodies are names or integer IDs. `abcorr` takes the
values of [`spkezr`](@ref), and each position is corrected as seen by the observer.
`et` can be a `ForwardDiff.Dual`.

Throws `SpiceError`:
- `:emptystring` if a string argument is empty;
- `:idcodenotfound` if a body name is not known;
- `:bodiesnotdistinct` if two of the three bodies are the same;
- `:invalidoption` if `abcorr` is not a valid correction;
- `:notrecognized` if a shape is not `POINT` or `SPHERE`;
- `:kernelvarnotfound`, `:typemismatch`, `:invalidcount` or `:badaxislength` if the
  radii of a sphere are missing or not valid;
- `:insidebody` if the observer is not outside a sphere;
- `:badarg` if `et` is not finite;
- the errors of [`spkpos`](@ref).

# Examples
```julia
# Needs an SPK with the Sun, the Earth and the Moon, and a PCK with their radii.
et = str2et("2024 APR 08 18:18:00 UTC")
sep = trgsep(et, "MOON", "SPHERE", " ", "SUN", "SPHERE", " ", "EARTH", "LT+S")
overlap = sep < 0           # the disks of the Moon and of the Sun overlap
```

See also [`phaseq`](@ref), [`gfsep`](@ref).

# References

- [NAIF documentation of `trgsep_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/trgsep_c.html)
"""
function trgsep(
        ks::KernelSet, et::Real, targ1::NameOrID, shape1::AbstractString,
        frame1::NameOrID, targ2::NameOrID, shape2::AbstractString, frame2::NameOrID,
        obsrvr::NameOrID, abcorr::AbcorrArg
    )
    _geo_nonempty(targ1, "target 1")
    _geo_nonempty(shape1, "shape 1")
    _geo_nonempty(frame1, "frame 1")
    _geo_nonempty(targ2, "target 2")
    _geo_nonempty(shape2, "shape 2")
    _geo_nonempty(frame2, "frame 2")
    _geo_nonempty(obsrvr, "observer")
    _geo_nonempty(abcorr, "aberration correction")
    b1 = _body_id(ks, targ1, "target 1")
    b2 = _body_id(ks, targ2, "target 2")
    obs = _body_id(ks, obsrvr, "observer")
    (obs == b1 || obs == b2 || b1 == b2) && throw(
        SpiceError(
            :bodiesnotdistinct,
            "the targets $targ1 and $targ2 and the observer $obsrvr must be three different bodies"
        )
    )
    corr = Aberration(abcorr)
    r1 = _geo_sepradius(ks, b1, shape1)
    r2 = _geo_sepradius(ks, b2, shape2)
    e = check_epoch(et)
    pv1, _ = _geo_spkezp(ks, b1, e, Int32(1), corr, obs)
    pv2, _ = _geo_spkezp(ks, b2, e, Int32(1), corr, obs)
    ang1 = _geo_angrad(r1, vnorm(pv1), obsrvr, targ1, 1)
    ang2 = _geo_angrad(r2, vnorm(pv2), obsrvr, targ2, 2)
    return vsep(pv1, pv2) - ang1 - ang2
end
trgsep(
    ctx::SpiceContext, et::Real, targ1::NameOrID, shape1::AbstractString, frame1::NameOrID,
    targ2::NameOrID, shape2::AbstractString, frame2::NameOrID, obsrvr::NameOrID,
    abcorr::AbcorrArg
) =
    trgsep(kernelset(ctx), et, targ1, shape1, frame1, targ2, shape2, frame2, obsrvr, abcorr)
trgsep(
    et::Real, targ1::NameOrID, shape1::AbstractString, frame1::NameOrID, targ2::NameOrID,
    shape2::AbstractString, frame2::NameOrID, obsrvr::NameOrID, abcorr::AbcorrArg
) =
    trgsep(default_context(), et, targ1, shape1, frame1, targ2, shape2, frame2, obsrvr, abcorr)
