# Azimuth and elevation: `azlrec`, `recazl`, `dazldr`, `drdazl` and `azlcpo` of CSPICE
# N0067. The coordinate calls are generic over the number type and do not allocate.

"""
    azlrec(range, az, el, azccw, elplsz) -> SVector{3}
    azlrec(range, az, el; azccw=true, elplsz=true) -> SVector{3}

Return the rectangular coordinates of the point at distance `range`, azimuth `az` and
elevation `el` (radians), as
`azlrec_c` of CSPICE.

The azimuth is zero along +X. `azccw` tells whether it grows counterclockwise about +Z,
seen from +Z. `elplsz` tells whether the elevation grows toward +Z. The ranges of the
inputs are not checked. The result has the unit of `range`.

See also [`recazl`](@ref), [`drdazl`](@ref).

# References

- [NAIF documentation of `azlrec_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/azlrec_c.html)
"""
function azlrec(range::Real, az::Real, el::Real, azccw::Bool, elplsz::Bool)
    r, a, e = promote(float(range), float(az), float(el))
    # As in CSPICE, a zero angle keeps its sign.
    return _m_latrec(r, (azccw || primal(a) == 0.0) ? a : -a, (elplsz || primal(e) == 0.0) ? e : -e)
end
azlrec(range::Real, az::Real, el::Real; azccw::Bool = true, elplsz::Bool = true) =
    azlrec(range, az, el, azccw, elplsz)

"""
    recazl(rectan, azccw, elplsz) -> (range, az, el)
    recazl(rectan; azccw=true, elplsz=true) -> (range, az, el)

Return the distance, the azimuth in [0, 2π) and the elevation in [-π/2, π/2] (radians)
of the point `rectan`, as
`recazl_c` of CSPICE.

`rectan` has 3 elements. The flags are those of [`azlrec`](@ref), which is the inverse
conversion.

See also [`dazldr`](@ref), [`azlcpo`](@ref).

# References

- [NAIF documentation of `recazl_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recazl_c.html)
"""
function recazl(rectan::AbstractVector, azccw::Bool, elplsz::Bool)
    # `recazl_` calls the Fortran `recrad_` and `reclat_`, whose longitude is the
    # arc tangent of the coordinates as they are (`reclat_c`, and so `recrad_c`, scale
    # them first).
    range, az, el = _m_reclat(_m_v3(rectan), false)
    primal(az) < 0.0 && (az += _TWOPI)
    if !azccw && primal(az) > 0.0
        az = 2.0 * Float64(pi) - az
        primal(az) < 0.0 && (az = zero(az))
    end
    (elplsz || primal(el) == 0.0) || (el = -el)
    return range, az, el
end
recazl(rectan::AbstractVector; azccw::Bool = true, elplsz::Bool = true) =
    recazl(rectan, azccw, elplsz)

"""
    dazldr(x, y, z, azccw, elplsz) -> SMatrix{3,3}
    dazldr(x, y, z; azccw=true, elplsz=true) -> SMatrix{3,3}

Return the Jacobian of the conversion from rectangular coordinates to (range, azimuth,
elevation) at the point `(x, y, z)`, as
`dazldr_c` of CSPICE.

Row `i` holds the partial derivatives of coordinate `i` of (range, azimuth, elevation)
with respect to `x`, `y` and `z`. The flags are those of [`azlrec`](@ref).

Throws `SpiceError(:pointonzaxis)` if `x` and `y` are both zero.

See also [`drdazl`](@ref), [`recazl`](@ref).

# References

- [NAIF documentation of `dazldr_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/dazldr_c.html)
"""
function dazldr(x::Real, y::Real, z::Real, azccw::Bool, elplsz::Bool)
    (primal(x) == 0.0 && primal(y) == 0.0) && throw(
        SpiceError(
            :pointonzaxis,
            "the Jacobian is not defined for a point on the Z axis"
        )
    )
    j = dlatdr(x, y, z)
    sa = azccw ? 1.0 : -1.0
    se = elplsz ? 1.0 : -1.0
    return SMatrix{3, 3}(
        j[1, 1], sa * j[2, 1], se * j[3, 1], j[1, 2], sa * j[2, 2],
        se * j[3, 2], j[1, 3], sa * j[2, 3], se * j[3, 3]
    )
end
dazldr(x::Real, y::Real, z::Real; azccw::Bool = true, elplsz::Bool = true) =
    dazldr(x, y, z, azccw, elplsz)

"""
    drdazl(range, az, el, azccw, elplsz) -> SMatrix{3,3}
    drdazl(range, az, el; azccw=true, elplsz=true) -> SMatrix{3,3}

Return the Jacobian of the conversion from (range, azimuth, elevation) to rectangular
coordinates, as
`drdazl_c` of CSPICE.

Column `j` holds the partial derivatives of `x`, `y` and `z` with respect to coordinate
`j` of (range, azimuth, elevation). The angles are in radians. The flags are those of
[`azlrec`](@ref).

Throws `SpiceError(:valueoutofrange)` if `range` is negative.

See also [`dazldr`](@ref), [`azlrec`](@ref).

# References

- [NAIF documentation of `drdazl_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/drdazl_c.html)
"""
function drdazl(range::Real, az::Real, el::Real, azccw::Bool, elplsz::Bool)
    primal(range) < 0.0 && throw(
        SpiceError(
            :valueoutofrange,
            "the range must not be negative; it is $(primal(range))"
        )
    )
    sa = azccw ? 1.0 : -1.0
    se = elplsz ? 1.0 : -1.0
    j = drdlat(range, sa * az, se * el)
    return SMatrix{3, 3}(
        j[1, 1], j[2, 1], j[3, 1], sa * j[1, 2], sa * j[2, 2], sa * j[3, 2],
        se * j[1, 3], se * j[2, 3], se * j[3, 3]
    )
end
drdazl(range::Real, az::Real, el::Real; azccw::Bool = true, elplsz::Bool = true) =
    drdazl(range, az, el, azccw, elplsz)

"""
    azlcpo([ctx,] method, target, et, abcorr, azccw, elplsz, obspos, obsctr, obsref)
        -> (azlsta, lt)

Return the range, azimuth and elevation of `target` and their rates, seen from an
observer at a constant position on the body `obsctr`, and the one-way light time `lt`
(s), as
`azlcpo_c` of CSPICE.

`azlsta` is an `SVector{6}`: range (km), azimuth and elevation (rad), then their rates
(km/s, rad/s). `obspos` is the observer position relative to the center of `obsctr`
(km, 3 elements), in the body-fixed frame `obsref`, which must be centered on `obsctr`.
`et` is in seconds past J2000 TDB.

The angles are measured in the topocentric frame at the observer. Its +Z axis is the
outward normal of the reference ellipsoid of `obsctr` at the point nearest to the
observer, and its +X axis points north. For an observer on the polar axis, +X is the +X
axis of `obsref`. `azccw` and `elplsz` are the flags of [`azlrec`](@ref).

`method` must be `"ELLIPSOID"`. `abcorr` is a correction string of [`spkezr`](@ref) or
an [`Aberration`](@ref). `et` can be a `ForwardDiff.Dual`. The call does not allocate
when the bodies and the frame are integer IDs.

Throws `SpiceError` with the kind
- `:invalidmethod` if `method` is not `"ELLIPSOID"`;
- `:invalidoption` if `abcorr` is not a valid correction;
- `:idcodenotfound` if a body name is not known;
- `:unknownframe` if `obsref` is not known, `:invalidframe` if it is not centered on
  `obsctr`;
- `:badarg` if `obspos` does not have 3 elements;
- `:pointonzaxis` if the target is on the local vertical;
- the kinds of [`srfrec`](@ref) if the radii of `obsctr` are missing or not valid.

Also throws the errors of [`spkcpo`](@ref).

# Examples
```julia
# Needs an SPK, a PCK with the radii and the rotation of the Earth, and a leap seconds kernel.
et = str2et("2025 MAR 01 12:00:00 UTC")
station = georec(-116.79 * rpd(), 35.25 * rpd(), 1.0, 6378.1366, 1 / 298.257)   # km
azlsta, lt = azlcpo("ELLIPSOID", "MOON", et, "LT+S", false, true, station, "EARTH", "IAU_EARTH")
az, el = azlsta[2] * dpr(), azlsta[3] * dpr()    # degrees; azimuth clockwise from north
```

See also [`spkcpo`](@ref), [`recazl`](@ref), [`dazldr`](@ref).

# References

- [NAIF documentation of `azlcpo_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/azlcpo_c.html)
"""
function azlcpo(
        ks::KernelSet, method::AbstractString, target::NameOrID, et::Real,
        abcorr::AbcorrArg, azccw::Bool, elplsz::Bool, obspos::AbstractVector, obsctr::NameOrID,
        obsref::NameOrID
    )
    _chkfstr("azlcpo", "string `method`", method)       # CHKFSTR of azlcpo_c
    _chkfstr("azlcpo", "string `target`", target)
    _chkfstr("azlcpo", "string `abcorr`", abcorr)
    _chkfstr("azlcpo", "string `obsctr`", obsctr)
    _chkfstr("azlcpo", "string `obsref`", obsref)
    obscde = _body_id(ks, obsctr, "observer center")
    fxfcde = _frame_code(ks, obsref)
    center, _, found = _geo_frinfo(ks, fxfcde)
    found || _geo_unknown_frame(obsref)
    center == obscde || _dsk_bad_frame(:invalidframe, obsref, obsctr, center)
    (method == "ELLIPSOID" || _eqstr(method, "ELLIPSOID")) || throw(
        SpiceError(
            :invalidmethod, "the method <$method> is not supported; use ELLIPSOID"
        )
    )
    pos = _dsk_vec3(obspos)
    xftopo = if pos[1] == 0.0 && pos[2] == 0.0
        # On the polar axis the local +X axis is the +X axis of the body-fixed frame.
        pos[3] < 0.0 ? _m_rotmat(one(Mat3{Float64}), Float64(pi), 1) : one(Mat3{Float64})
    else
        a, b, c = _positive_radii(ks, obscde)
        obsspt, _ = _m_nearpt(pos, a, b, c)
        twovec(_m_surfnm(a, b, c, obsspt), 3, SVector(0.0, 0.0, 1.0), 1)
    end
    corr = Aberration(abcorr)
    trgcde = _body_id(ks, target, "target")
    f = _GeoFixed(pos, zero(Vec3{Float64}), 0.0, obscde, fxfcde)
    state, lt = _geo_cv(ks, trgcde, check_epoch(et), fxfcde, _GEO_AT_OBSERVER, corr, f)
    r, v = _m_pv(state)
    lr = _m_mxv(xftopo, r)
    lv = _m_mxv(xftopo, v)
    range, az, el = recazl(lr, azccw, elplsz)
    rates = _m_mxv(dazldr(lr[1], lr[2], lr[3], azccw, elplsz), lv)
    return SVector(range, az, el, rates[1], rates[2], rates[3]), lt
end
azlcpo(
    ctx::SpiceContext, method::AbstractString, target::NameOrID, et::Real,
    abcorr::AbcorrArg, azccw::Bool, elplsz::Bool, obspos::AbstractVector, obsctr::NameOrID,
    obsref::NameOrID
) =
    azlcpo(kernelset(ctx), method, target, et, abcorr, azccw, elplsz, obspos, obsctr, obsref)
azlcpo(
    method::AbstractString, target::NameOrID, et::Real, abcorr::AbcorrArg, azccw::Bool,
    elplsz::Bool, obspos::AbstractVector, obsctr::NameOrID, obsref::NameOrID
) =
    azlcpo(default_context(), method, target, et, abcorr, azccw, elplsz, obspos, obsctr, obsref)
