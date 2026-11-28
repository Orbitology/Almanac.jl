# Local solar time (`et2lst`) and planetocentric longitude of the Sun (`lspcn`).

"`HH:MM:SS` of the numbers `0 <= h, m, s <= 99`, then `suffix`, in one buffer."
function _lst_clock(h::Int, m::Int, s::Int, suffix::String)
    buf = Base.StringVector(8 + ncodeunits(suffix))
    for (k, v) in ((1, h), (4, m), (7, s))
        buf[k] = UInt8('0') + v ÷ 10
        buf[k + 1] = UInt8('0') + v % 10
    end
    buf[3] = buf[6] = UInt8(':')
    copyto!(buf, 9, codeunits(suffix), 1, ncodeunits(suffix))
    return String(buf)
end

"""
    et2lst([ctx,] et, body, lon, type) -> (hr, mn, sc, time, ampm)

Return the local solar time at `et` for the longitude `lon` on the body with ID
`body`, as
`et2lst_c`.

- `et`: seconds past J2000 TDB, a plain real number (not a `ForwardDiff.Dual`).
- `lon`: longitude in radians.
- `type`: `"PLANETOCENTRIC"` (`lon` positive east) or `"PLANETOGRAPHIC"` (`lon` with
  the sense of [`recpgr`](@ref) for the body: positive west for a prograde body such as
  Mars), in any letter case.

`hr`, `mn` and `sc` are the integer hour, minute and second on a 24 hour clock. `time`
is the string `HH:MM:SS`, and `ampm` is the string `HH:MM:SS A.M.` or `HH:MM:SS P.M.`
on a 12 hour clock.

The local solar time is the angle between the meridian of `lon` and the meridian of
the Sun, with 12:00:00 at the sub-solar meridian and 24 "hours" in one rotation; an
"hour" is not 3600 s. The position of the Sun has the correction `"LT+S"` and is in
the body-fixed frame of [`cidfrm`](@ref). The sign of the second value of
`BODY<body>_PM` gives the direction of the rotation. For the Sun (`body = 10`) the
result is always 12:00:00. The call needs SPK data for the body and the Sun and PCK
data for the body.

Throws `SpiceError`, with the checks in this order, as in CSPICE after the first:
- `:badarg` if `et` is a `ForwardDiff.Dual`;
- `:emptystring` if `type` is the empty string;
- `:invalidoption`, for `"PLANETOGRAPHIC"`, if `BODY<body>_PGR_POSITIVE_LON` is not
  `EAST` or `WEST`;
- `:missingdata`, for `"PLANETOGRAPHIC"`, if `BODY<body>_PGR_POSITIVE_LON` and
  `BODY<body>_PM` are both absent and the body is not the Earth;
- `:unknownsystem` if `type` is not one of the two;
- `:cantfindframe` if the body has no body-fixed frame;
- the kinds of [`spkezr`](@ref);
- `:cantgetrotationtype` if `BODY<body>_PM` is absent or is not a numeric variable of
  two or more values.

# Examples
```julia
furnsh("naif0012.tls"); furnsh("pck00011.tpc"); furnsh("de440s.bsp")
et = str2et("2025-03-01T12:00:00")
hr, mn, sc, time, ampm = et2lst(et, 399, deg2rad(-118.0), "PLANETOCENTRIC")
```

# References

- [NAIF documentation of `et2lst_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/et2lst_c.html)
"""
function et2lst(ks::KernelSet, et::Real, body::Integer, lon::Real, type::AbstractString)
    _plain_epoch(et, "et2lst")
    _chkfstr("et2lst", "coordinate system", type)          # CHKFSTR of et2lst_c
    system = strip(type)
    if _eqci(system, "PLANETOGRAPHIC", true)
        g = _m_pgrsense(ks, body) * Float64(lon)    # the sense of `pgrrec_`
        y, x = sincos(g)
        mylong = (x == 0.0 && y == 0.0) ? 0.0 : atan(y, x)
    elseif _eqci(system, "PLANETOCENTRIC", true)
        mylong = Float64(lon)
    else
        throw(
            SpiceError(
                :unknownsystem, "the coordinate system '$type' is not a " *
                    "recognized system of longitude; the recognized systems are " *
                    "'PLANETOCENTRIC' and 'PLANETOGRAPHIC'"
            )
        )
    end
    body == 10 && return 12, 0, 0, "12:00:00", "12:00:00 P.M."
    fix = cidfrm(ks, body)
    fix === nothing && throw(
        SpiceError(
            :cantfindframe, "the body-fixed frame " *
                "associated with body $body could not be determined"
        )
    )
    state, _ = spkezr(ks, 10, et, fix[1], Aberration("LT+S"), body)
    slong = (state[1] == 0.0 && state[2] == 0.0) ? 0.0 : atan(state[2], state[1])
    _, angle = _rmaind(mylong - slong, _TWOPI)
    angle > pi && (angle -= _TWOPI)
    pm = get(ks.pool, _body_key(body, "PM"), nothing)
    (pm isa Vector{Float64} && length(pm) >= 2) || throw(
        SpiceError(
            :cantgetrotationtype,
            "the rotation type for the body $body could not be determined because the " *
                "BODY$(body)_PM keyword was not found in the pool or was not of the expected " *
                "type or dimension"
        )
    )
    pm[2] < 0.0 && (angle = -angle)
    secnds = clamp(angle * 86400.0 / _TWOPI, -43200.0, 43200.0)
    hours, rest = _rmaind(secnds, 3600.0)
    mins, secnds = _rmaind(rest, 60.0)
    hr, mn, sc = trunc(Int, hours) + 12, trunc(Int, mins), trunc(Int, secnds)
    hr == 24 && (hr = 0)
    h12 = hr == 0 ? 12 : hr > 12 ? hr - 12 : hr
    return hr, mn, sc, _lst_clock(hr, mn, sc, ""), _lst_clock(h12, mn, sc, hr >= 12 ? " P.M." : " A.M.")
end
et2lst(ctx::SpiceContext, et::Real, body::Integer, lon::Real, type::AbstractString) =
    et2lst(kernelset(ctx), et, body, lon, type)
et2lst(et::Real, body::Integer, lon::Real, type::AbstractString) =
    et2lst(default_context(), et, body, lon, type)

"""
    lspcn([ctx,] body, et, abcorr) -> Float64

Return the planetocentric longitude of the Sun (`Ls`) for `body` at `et`, in radians
in `[0, 2pi)`, as
`lspcn_c`.

`Ls` is the season angle of the body: 0 at its spring equinox of the north, `pi/2` at
the summer solstice of the north. It is the right ascension of the Sun as seen from
the body, in this frame:

- the Z axis is the orbital angular momentum of the body about the Sun;
- the X axis is the cross product of the north pole of the body and the Z axis.

`body` is a name or a NAIF ID, `et` is in seconds past J2000 TDB (a plain real number,
not a `ForwardDiff.Dual`), and `abcorr` is as for [`spkezr`](@ref). The call needs PCK data for the pole and SPK data for the body
and the Sun.

Throws `SpiceError(:badarg)` if `et` carries derivatives, `SpiceError(:dependentvectors)`
if the pole is parallel to the orbital angular momentum, and the errors of
[`spkezr`](@ref) and [`tipbod`](@ref).

# Examples
```julia
furnsh("naif0012.tls"); furnsh("pck00011.tpc"); furnsh("de440s.bsp")
ls = rad2deg(lspcn("EARTH", str2et("2025-03-01"), "LT+S"))    # degrees
```

# References

- [NAIF documentation of `lspcn_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/lspcn_c.html)
"""
function lspcn(ks::KernelSet, body::NameOrID, et::Real, abcorr::AbcorrArg)
    _chkfstr("lspcn", "string `body`", body)       # CHKFSTR of lspcn_c
    _plain_epoch(et, "lspcn")
    id = body_code(ks, body)
    tipm = tipbod(ks, 1, id, et)
    npole = Vec3{Float64}(tipm[3, 1], tipm[3, 2], tipm[3, 3])
    r, v = spk_jet(ks, id, check_epoch(et), _J2000_FRAME, Int32(10), Val(1))
    # The arithmetic of `lspcn_`: `ucrss`, `twovec` (Z along the angular momentum, Y in
    # the plane of the pole), `mxv` and `recrad` (with `reclat_`).
    trans = twovec(_m_ucrss(r, v), 3, npole, 2)
    s, _ = spkezr(ks, 10, et, 1, Aberration(abcorr), id)
    pos = _m_mxv(trans, Vec3{Float64}(s[1], s[2], s[3]))
    _, ra, _ = _m_reclat(pos, false)
    return ra < 0.0 ? ra + _TWOPI : ra
end
lspcn(ctx::SpiceContext, body::NameOrID, et::Real, abcorr::AbcorrArg) =
    lspcn(kernelset(ctx), body, et, abcorr)
lspcn(body::NameOrID, et::Real, abcorr::AbcorrArg) = lspcn(default_context(), body, et, abcorr)
