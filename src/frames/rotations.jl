# Frame orientation: text-PCK and binary-PCK body frames (class 2), fixed-offset
# frames (class 4), built-in inertial frames (class 1), and the chain logic that
# joins them. Everything that depends on the epoch is generic over the epoch type, so
# an epoch of type `Taylor{N,Float64}` gives the time derivatives to order N.

const _SEC_PER_CENTURY = _SPD * 36525.0          # seconds per Julian century

"""
Number of frames of the chain of the first frame that `frame_rotation` keeps, as the
arrays of `zzrefch0_` hold them: further links are folded into the last one, so a chain
can be of any length, but a frame of the second chain meets only those frames.
"""
const MAX_FRAME_CHAIN = 10

"""
Number of links after which a frame chain is taken to be a cycle of frame definitions
and throws `SpiceError(:framechain)`; `zzrefch0_` has no bound and never returns.
"""
const _FRAME_CHAIN_LIMIT = 1000

# --- small rotation helpers -----------------------------------------------------------

"""
    _eul313(w, delta, phi) -> Mat3

`[w]_3 * [delta]_1 * [phi]_3`, the same matrix as `eul2m_(w, delta, phi, 3, 1, 3)` with
the arithmetic of `rotate_` and `rotmat_` (the products with the constant 0 and 1
elements are left out, which does not change a result). Generic over the angle type.
"""
@inline function _eul313(w::T, delta::T, phi::T) where {T <: Real}
    sp, cp = sincos(phi)
    sd, cd = sincos(delta)
    sw, cw = sincos(w)
    # rows 2 and 3 after the rotation about axis 1
    m21 = -(cd * sp)
    m22 = cd * cp
    m31 = sd * sp
    m32 = -(sd * cp)
    return Mat3{T}(
        cw * cp + sw * m21, cw * m21 - sw * cp, m31,
        cw * sp + sw * m22, cw * m22 - sw * sp, m32,
        sw * sd, cw * sd, cd
    )
end

# --- text PCK constants ---------------------------------------------------------------

"""
    _bodbry(body) -> Int

Barycenter code of the planetary system of `body` (planets and satellites `PNN` and
`PXNNN` give `P`), or `body` itself for any other body, as `zzbodbry_` of CSPICE.
"""
_bodbry(body::Integer) =
    100 <= body <= 999 ? Int(body) ÷ 100 : 10000 <= body <= 99999 ? Int(body) ÷ 10000 : Int(body)

_pool_numbers_or_empty(pool, name) = something(_pool_numbers(pool, name), Float64[])

"First value of the first of the numeric pool variables `names` that exists, else `default`."
function _pool_first(pool, names, default::Float64)
    for name in names
        v = _pool_numbers(pool, name)
        v !== nothing && !isempty(v) && return v[1]
    end
    return default
end

"""
    build_bodyrot(pool) -> Dict{Int32,BodyRot}

Orientation constants of every body that has numeric `BODY<id>_PM`, `_POLE_RA` and
`_POLE_DEC` in `pool`, read as `tisbod_` reads them: the nutation-precession angles,
the reference frame and the reference epoch come from the barycenter of the body
system (`_bodbry`), the amplitudes from the body itself. `nut_prec_angles` is the
stored vector; its layout is `phase_degree + 1` coefficients per angle, with
`phase_degree` from `BODY<sys>_MAX_PHASE_DEGREE` (default 1).
"""
function build_bodyrot(pool::AbstractDict)
    out = Dict{Int32, BodyRot}()
    for (name, value) in pool
        value isa Vector{Float64} || continue
        # The test of the ends spares the regular expression for most names.
        (startswith(name, "BODY") && endswith(name, "_PM")) || continue
        m = match(r"^BODY(-?\d+)_PM$", name)
        m === nothing && continue
        id = tryparse(Int32, m.captures[1])
        id === nothing && continue
        # The names are joined from strings: `*` of strings is faster than interpolation.
        bid = "BODY" * string(id)
        ra = _pool_numbers(pool, bid * "_POLE_RA")
        dec = _pool_numbers(pool, bid * "_POLE_DEC")
        (ra === nothing || dec === nothing) && continue
        sys = _bodbry(id)
        bs = "BODY" * string(sys)
        epoch1, epoch2 = bs * "_CONSTANTS_JED_EPOCH", bs * "_CONSTS_JED_EPOCH"
        ref1, ref2 = bs * "_CONSTANTS_REF_FRAME", bs * "_CONSTS_REF_FRAME"
        epoch = _pool_first(pool, (epoch1, epoch2), _JD2000)
        ref = _pool_first(pool, (ref1, ref2), 1.0)
        both(a, b) = _pool_numbers(pool, a) !== nothing && _pool_numbers(pool, b) !== nothing
        # `tisbod_` reads the frame with `gipool_` after the epoch check and before the
        # check of the two frame forms.
        refid = _gipool_value(ref)
        specerr = both(epoch1, epoch2) ? :competingepochspec :
            refid === nothing ? :intoutofrange :
            both(ref1, ref2) ? :competingframespec : :none
        out[id] = BodyRot(
            ra, dec, value,
            _pool_numbers_or_empty(pool, bid * "_NUT_PREC_RA"),
            _pool_numbers_or_empty(pool, bid * "_NUT_PREC_DEC"),
            _pool_numbers_or_empty(pool, bid * "_NUT_PREC_PM"),
            _pool_numbers_or_empty(pool, bs * "_NUT_PREC_ANGLES"),
            something(refid, _J2000_FRAME), epoch, _phase_degree(pool, sys), specerr
        )
    end
    return out
end

"""
    _phase_degree(pool, sys) -> Int32

Value of `BODY<sys>_MAX_PHASE_DEGREE` rounded as CSPICE rounds it, or 1 if the variable
is absent. The value is not checked here; `_phase_ncoef` checks it at evaluation.
"""
function _phase_degree(pool::AbstractDict, sys::Integer)
    v = _pool_numbers(pool, "BODY" * string(sys) * "_MAX_PHASE_DEGREE")
    (v === nothing || isempty(v)) && return Int32(1)
    return Int32(_nint(v[1]))                  # `i_dnnt` of `bodvcd_`'s value
end

"""
    _bodyrot_check(br, body)

The checks of `tisbod_` before an evaluation: both forms of the epoch or the frame
variable (`br.specerr`), then the degree of the phase angles whenever the system has
phase angles, whether or not the body has terms that use them.
"""
@inline function _bodyrot_check(br::BodyRot, body::Integer)
    sys = _bodbry(body)
    br.specerr === :competingepochspec && throw(
        SpiceError(
            :competingepochspec, "both BODY$(sys)_CONSTANTS_JED_EPOCH and " *
                "BODY$(sys)_CONSTS_JED_EPOCH are in the kernel pool; at most one may be"
        )
    )
    br.specerr === :intoutofrange && throw(
        SpiceError(
            :intoutofrange, "the value of BODY$(sys)_CONSTANTS_REF_FRAME or " *
                "BODY$(sys)_CONSTS_REF_FRAME is outside the range of integers"
        )
    )
    br.specerr === :competingframespec && throw(
        SpiceError(
            :competingframespec, "both BODY$(sys)_CONSTANTS_REF_FRAME and " *
                "BODY$(sys)_CONSTS_REF_FRAME are in the kernel pool; at most one may be"
        )
    )
    isempty(br.nut_prec_angles) || _phase_ncoef(br, body)
    return nothing
end

"""
    _phase_ncoef(br, body) -> Int

Number of polynomial coefficients per nutation-precession angle: `br.phase_degree + 1`.
Throws `SpiceError(:degreeoutofrange)` unless the degree is in 1:3.
"""
@inline function _phase_ncoef(br::BodyRot, body::Integer)
    deg = Int(br.phase_degree)
    1 <= deg <= 3 || throw(
        SpiceError(
            :degreeoutofrange, "maximum phase angle degree " *
                "for body $(_bodbry(body)) must be in the range 1:3 but was $deg"
        )
    )
    return deg + 1
end

"""
    _phase_ncoef(br, body, nterm) -> Int

`_phase_ncoef(br, body)` after the check of `tisbod_` that the phase angle list holds
the `nterm` angles that the longest amplitude list needs: throws
`SpiceError(:insufficientangles)` otherwise.
"""
@inline function _phase_ncoef(br::BodyRot, body::Integer, nterm::Int)
    nco = _phase_ncoef(br, body)
    nphase = length(br.nut_prec_angles) ÷ nco
    nterm <= nphase || throw(
        SpiceError(
            :insufficientangles,
            "insufficient number of nutation/precession angles for body $body: " *
                "number of angles is $nphase, number required is $nterm"
        )
    )
    return nco
end

"Value at `x` of the polynomial with coefficients `c` (lowest order first) by Horner's rule."
@inline function _horner(c::Vector{Float64}, x::T) where {T <: Real}
    n = length(c)
    n == 0 && return zero(T)
    acc = zero(T) + @inbounds c[n]
    for k in (n - 1):-1:1
        acc = @inbounds(c[k]) + x * acc
    end
    return acc
end

"Rotation from the inertial frame `req` to the inertial frame `ref`, as `irfrot_`. Requires `req != ref`."
@inline _irf_rot(req::Integer, ref::Integer) =
    req == _J2000_FRAME ? inertial_rotation(ref) :
    ref == _J2000_FRAME ? transpose(inertial_rotation(req)) :
    _m_mxmt(inertial_rotation(ref), inertial_rotation(req))

"""
    _text_pck_rotation(br, body, et::T, req = J2000) -> Mat3{T}

Matrix from the inertial frame `req` to the body-fixed frame of `body` from the
text-PCK constants `br`, as `tisbod_` of CSPICE computes it (same order of operations):
RA, DEC and W polynomials, trigonometric terms, `eul2m(W, pi/2 - DEC, pi/2 + RA, 3, 1,
3)`, then, if the reference frame of the constants is not `req`, the rotation from `req`
to that frame (`irfrot_`). Throws
`SpiceError(:insufficientangles)` if an amplitude vector is longer than the angle list.
"""
function _text_pck_rotation(
        br::BodyRot, body::Integer, et::T, req::Integer = _J2000_FRAME
    ) where {T <: Real}
    _bodyrot_check(br, body)
    epoch = et - _SPD * (br.jed_epoch - _JD2000)
    td = epoch / _SPD
    tc = epoch / _SEC_PER_CENTURY
    ra = _horner(br.pole_ra, tc)
    dec = _horner(br.pole_dec, tc)
    w = _horner(br.pm, td)

    ac, dc, wc, th = br.nut_prec_ra, br.nut_prec_dec, br.nut_prec_pm, br.nut_prec_angles
    na, nd, nw = length(ac), length(dc), length(wc)
    nterm = max(na, nd, nw)
    if nterm > 0
        nco = _phase_ncoef(br, body, nterm)
        sa = zero(T)
        sd = zero(T)
        sw = zero(T)
        @inbounds for i in 1:nterm    # invariant: sa, sd, sw hold the sums of terms 1:i-1
            k = nco * (i - 1)
            if nco == 2
                theta = (th[k + 1] + tc * th[k + 2]) * _RPD
            else
                theta = zero(T) + th[k + 1]
                p = tc
                for j in 2:nco
                    theta += p * th[k + j]
                    p = tc * p
                end
                theta *= _RPD
            end
            s, c = sincos(theta)
            i <= na && (sa += ac[i] * s)
            i <= nd && (sd += dc[i] * c)
            i <= nw && (sw += wc[i] * s)
        end
        ra += sa
        dec += sd
        w += sw
    end

    ra *= _RPD
    dec *= _RPD
    w = f2c_dmod(w * _RPD, _TWOPI)
    tipm = _eul313(w, _HALFPI - dec, ra + _HALFPI)
    br.ref == req && return tipm
    return _m_mxm(tipm, _irf_rot(req, br.ref))
end

"An empty segment list: the binary PCK segments of a body that has none. Never mutated."
const _NO_SEGMENTS = SegmentList(Segment[])

"""
    _body_state(ks, classid, et::Float64, req = J2000) -> Mat3{Taylor{1,Float64}}

Return the matrix of `_body_rotation` with its time derivative as `tisbod_` gives
them: from a binary PCK segment if one covers `et` (`_binary_pck_state`), else from the
text PCK constants (`_text_pck_state`). Throws as `_body_rotation`.
"""
_body_state(ks::KernelSet, classid::Int32, et::Float64, req::Integer = _J2000_FRAME) =
    _body_state(ks, get(ks.pck, classid, _NO_SEGMENTS), classid, et, req)

# `segs` are the binary PCK segments of `classid` in `ks`, `_NO_SEGMENTS` for none.
function _body_state(
        ks::KernelSet, segs::SegmentList{Segment}, classid::Int32, et::Float64,
        req::Integer = _J2000_FRAME
    )
    i = find_segment(segs, et)
    if i > 0
        seg = @inbounds segs[i]
        return _binary_pck_state(ks.daf[seg.file], seg, et, req)
    end
    br = get(ks.bodyrot, classid, nothing)
    br === nothing && return _body_rotation(ks, segs, classid, taylor_var(et, Val(1)), req)::Mat3{Taylor{1, Float64}}
    return _text_pck_state(br, classid, et, req)
end

"`x^n` for `n >= 0` by the repeated squaring of `pow_di` of the f2c library."
@inline function _powdi(x::Float64, n::Int)
    p = 1.0
    while true     # invariant: p * x^n is the result
        isodd(n) && (p *= x)
        n >>= 1
        n == 0 && break
        x *= x
    end
    return p
end

"""
    _text_pck_state(br, body, et::Float64, req = J2000) -> Mat3{Taylor{1,Float64}}

Return the matrix of `_text_pck_rotation` with its time derivative, with the
arithmetic of `tisbod_` of CSPICE: the rates of RA, DEC and W from the derivatives of
their polynomials and trigonometric terms, written out as `tisbod_` writes them, and the
derivative of the matrix from these angles and rates by `eul2xf_`. The derivative of
the `Taylor` evaluation of `_text_pck_rotation` is the same function, rounded
otherwise. One exception: the rate of a cubic nutation-precession angle is the true
one, where `tisbod_` of N0067 divides its cubic term by the seconds per century once
too often. Throws as `_text_pck_rotation`.
"""
function _text_pck_state(br::BodyRot, body::Integer, et::Float64, req::Integer = _J2000_FRAME)
    _bodyrot_check(br, body)
    d = _SPD
    t = _SEC_PER_CENTURY
    epoch = et - d * (br.jed_epoch - _JD2000)
    td = epoch / d
    tc = epoch / t
    coef(c, i) = i <= length(c) ? @inbounds(c[i]) : 0.0
    rc, dc0, wc0 = br.pole_ra, br.pole_dec, br.pm
    # `tisbod_` holds three coefficients of each polynomial, the missing ones zero
    ra = coef(rc, 1) + tc * (coef(rc, 2) + tc * coef(rc, 3))
    dec = coef(dc0, 1) + tc * (coef(dc0, 2) + tc * coef(dc0, 3))
    w = coef(wc0, 1) + td * (coef(wc0, 2) + td * coef(wc0, 3))
    dra = (coef(rc, 2) + 2.0 * tc * coef(rc, 3)) / t
    ddec = (coef(dc0, 2) + 2.0 * tc * coef(dc0, 3)) / t
    dw = (coef(wc0, 2) + 2.0 * td * coef(wc0, 3)) / d
    ac, dc, wc, th = br.nut_prec_ra, br.nut_prec_dec, br.nut_prec_pm, br.nut_prec_angles
    na, nd, nw = length(ac), length(dc), length(wc)
    nterm = max(na, nd, nw)
    if nterm > 0
        nco = _phase_ncoef(br, body, nterm)
        sa = sd = sw = 0.0
        dsa = dsd = dsw = 0.0
        @inbounds for i in 1:nterm    # invariant: the sums hold the terms 1:i-1
            k = nco * (i - 1)
            if nco == 2
                theta = (th[k + 1] + tc * th[k + 2]) * _RPD
                dtheta = (th[k + 2] / t) * _RPD
            else
                theta = 0.0
                for j in 1:nco
                    theta += _powdi(tc, j - 1) * th[k + j]
                end
                theta *= _RPD
                dtheta = th[k + 2] / t
                for j in 3:nco
                    # `tisbod_` of N0067 divides by `t^(l-1)`, which is wrong for a cubic
                    # term (l = 3); the rate of the angle is divided by `t` once.
                    l = j - 1
                    dtheta += (l * _powdi(tc, l - 1) * th[k + j]) / t
                end
                dtheta *= _RPD
            end
            s, c = sincos(theta)
            ds = c * dtheta
            dcs = -s * dtheta
            i <= na && (sa += ac[i] * s; dsa += ac[i] * ds)
            i <= nd && (sd += dc[i] * c; dsd += dc[i] * dcs)
            i <= nw && (sw += wc[i] * s; dsw += wc[i] * ds)
        end
        ra += sa
        dec += sd
        w += sw
        dra += dsa
        ddec += dsd
        dw += dsw
    end
    ra *= _RPD
    dec *= _RPD
    w = f2c_dmod(w * _RPD, _TWOPI)
    dra *= _RPD
    ddec *= _RPD
    dw *= _RPD
    x = eul2xf(SVector(w, _HALFPI - dec, ra + _HALFPI, dw, -ddec, dra), 3, 1, 3)
    r = _eul313(w, _HALFPI - dec, ra + _HALFPI)
    dr = _m_blocks(x)[2]
    if br.ref != req
        c = _irf_rot(req, br.ref)
        r, dr = _m_mxm(r, c), _m_mxm(dr, c)
    end
    return map((a, b) -> Taylor((a, b)), r, dr)
end

"""
    _binary_pck_rotation(data, seg, et::T, req = J2000) -> Mat3{T}

Matrix from the inertial frame `req` to the body frame of the binary PCK segment `seg`
at `et`, as `pcke02_`, `pcke03_`, `pckmat_` and `tisbod_` give it: 3-1-3 Euler angles of
the segment (W reduced by whole turns for types 2 and 20), then, if the reference frame
of the segment is not `req`, the rotation from `req` to that frame. Requires `seg.t0 <= primal(et) <= seg.t1`.
"""
@inline function _binary_pck_rotation(
        data::Vector{Float64}, seg::Segment, et::T, req::Integer = _J2000_FRAME
    ) where {T <: Real}
    ang = pck_eval(data, seg, et, Val(0))[1]
    # `pcke02_` and `pcke20_` reduce W with `d_mod`, and `pcke03_` does not; the same
    # arithmetic keeps the value of CSPICE.
    w = seg.dtype == 3 ? ang[3] : f2c_dmod(ang[3], _TWOPI)
    tipm = _eul313(w, ang[2], ang[1])
    seg.frame == req && return tipm
    return _m_mxm(tipm, _irf_rot(req, seg.frame))
end

"""
    _binary_pck_state(data, seg, et::Float64, req = J2000) -> Mat3{Taylor{1,Float64}}

Return the matrix of `_binary_pck_rotation` with its time derivative, with the
arithmetic of `pckmat_` and `tisbod_`: the angles and their rates of the segment
(`chbint_` for type 2, `spke20_` for type 20; W reduced by whole turns as `pcke02_`
and `pcke20_` do) go to `eul2xf_`, and both blocks are multiplied on the right by the
rotation from `req` to the frame of the segment if the two differ. The derivative of the `Taylor`
evaluation of `_binary_pck_rotation` is the same function, rounded otherwise. Type 3
takes the derivative from the angular velocity series of the segment (`_pck03_state`),
which is not the derivative of the `Taylor` evaluation.
"""
function _binary_pck_state(data::Vector{Float64}, seg::Segment, et::Float64, req::Integer = _J2000_FRAME)
    if seg.dtype == 3
        r, dr = _pck03_state(data, seg, et)
    else
        ang, rate = pck_eval(data, seg, et, Val(1))
        x = eul2xf(
            SVector(f2c_dmod(ang[3], _TWOPI), ang[2], ang[1], rate[3], rate[2], rate[1]), 3, 1, 3
        )
        r, dr = _m_blocks(x)
    end
    if seg.frame != req
        c = _irf_rot(req, seg.frame)
        r, dr = _m_mxm(r, c), _m_mxm(dr, c)
    end
    return map((a, b) -> Taylor((a, b)), r, dr)
end

"""
    _missing_pole(ks, body)

With `BODY<body>_PM` in the pool and a pole angle absent, throw
`SpiceError(:kernelvarnotfound)` for the first absent one, as `bodvcd_` in `tisbod_`
does; return otherwise.
"""
@noinline function _missing_pole(ks::KernelSet, body::Int32)
    haskey(ks.pool, "BODY$(body)_PM") || return nothing
    for item in ("POLE_RA", "POLE_DEC")
        name = "BODY$(body)_$item"
        _pool_numbers(ks.pool, name) === nothing && throw(
            SpiceError(
                :kernelvarnotfound,
                "the variable $name could not be found in the kernel pool"
            )
        )
    end
    return nothing
end

"""
    _body_rotation(ks, classid, et::T, req = J2000) -> Mat3{T}

Matrix from the inertial frame `req` to the class 2 frame with class ID `classid`: from a binary PCK
segment if one covers `primal(et)`, otherwise from the text-PCK constants, as
`tisbod_` does. Throws `SpiceError(:framedatanotfound)` if neither exists, and
`SpiceError(:kernelvarnotfound)` if `BODY<id>_PM` exists without a pole angle.
"""
_body_rotation(ks::KernelSet, classid::Int32, et::Real, req::Integer = _J2000_FRAME) =
    _body_rotation(ks, get(ks.pck, classid, _NO_SEGMENTS), classid, et, req)

# `segs` are the binary PCK segments of `classid` in `ks`, `_NO_SEGMENTS` for none.
function _body_rotation(
        ks::KernelSet, segs::SegmentList{Segment}, classid::Int32, et::T,
        req::Integer = _J2000_FRAME
    ) where {T <: Real}
    i = find_segment(segs, primal(et))
    if i > 0
        seg = @inbounds segs[i]
        return _binary_pck_rotation(ks.daf[seg.file], seg, et, req)::Mat3{T}
    end
    br = get(ks.bodyrot, classid, nothing)
    br === nothing && _missing_pole(ks, classid)
    br === nothing && throw(
        SpiceError(
            :framedatanotfound,
            "PCK data required to compute the orientation of the body-fixed frame with " *
                "class ID $classid at ET $(primal(et)) were not found (no binary PCK segment " *
                "covers the epoch and no text PCK constants are loaded)"
        )
    )
    return _text_pck_rotation(br, classid, et, req)::Mat3{T}
end

# --- fixed-offset (TK) frames ---------------------------------------------------------

"""
    _q2m(q) -> Mat3{Float64}

Rotation matrix of the SPICE-style quaternion `q` (scalar first), as `q2m_` of CSPICE
after the normalization that `tkfram_` applies (`vhatg_`). Requires `q` of length 4.
"""
function _q2m(q::AbstractVector{Float64})
    # `vnormg_`: the norm scaled by the largest component
    qmax = max(abs(q[1]), abs(q[2]), abs(q[3]), abs(q[4]))
    n = qmax == 0.0 ? 0.0 :
        qmax * sqrt((q[1] / qmax)^2 + (q[2] / qmax)^2 + (q[3] / qmax)^2 + (q[4] / qmax)^2)
    q0, q1, q2, q3 = n > 0 ? (q[1] / n, q[2] / n, q[3] / n, q[4] / n) : (q[1], q[2], q[3], q[4])
    return _m_q2m(q0, q1, q2, q3, false)
end

"""
    _tk_key(pool, idkey, namekey, suffix) -> String

Name of the TK frame variable `suffix` that `tkfram_` reads: the form keyed by frame
name if it is in `pool`, the form keyed by ID otherwise (present or not). Throws
`SpiceError(:competingframespec)` if both forms are present.
"""
function _tk_key(
        pool::AbstractDict, idkey::AbstractString, namekey::AbstractString,
        suffix::AbstractString
    )
    a = idkey * suffix
    b = namekey * suffix
    namekey == idkey && return a
    haskey(pool, a) && haskey(pool, b) && throw(
        SpiceError(
            :competingframespec,
            "frame definition keywords '$a' and '$b' are both present in the pool"
        )
    )
    return haskey(pool, b) ? b : a
end

"""
    _tk_relative_id(pool, frames, name) -> Union{Int32,Nothing}

ID of the frame called `name`: a frame of `frames` (defined in the pool), a
`FRAME_<name>` pool variable, a built-in frame of `frames.jl`, or `nothing`.
"""
function _tk_relative_id(pool::AbstractDict, frames::AbstractDict, name::AbstractString)
    key = normalize_name(name)
    for fd in values(frames)
        normalize_name(fd.name) == key && return fd.id
    end
    # `namfrm_`: `gipool_(FRAME_<name>, 1, 8, ...)`, then the test of one value
    k = "FRAME_" * key
    v = _gipool_first(get(pool, k, nothing), k; room = 8, single = true)
    v === nothing || return v
    return get(BUILTIN_FRAMENAME2ID, key, nothing)
end

"""
    _tk_matrix(pool, idkey, namekey, spec) -> Mat3{Float64}

Constant matrix of one TK frame (from the frame to its relative frame) for the
specification `spec` (`MATRIX`, `ANGLES` or `QUATERNION`), as `tkfram_` of CSPICE builds
it. Throws the `SpiceError` of `tkfram_` for a missing or malformed variable.
"""
function _tk_matrix(
        pool::AbstractDict, idkey::AbstractString, namekey::AbstractString,
        spec::AbstractString
    )
    kind = uppercase(rstrip(spec, ' '))       # ucase_, then s_cmp: trailing blanks only
    if kind == "MATRIX"
        key = _tk_key(pool, idkey, namekey, "_MATRIX")
        _badkpv(pool, "TKFRAM", key, "=", 9, 1, 'N')
        m = Mat3{Float64}(pool[key])            # stored column by column
        r = _m_sharpr(m)
        c2 = dot(r[:, 2], m[:, 2]) < 0 ? -r[:, 2] : r[:, 2]
        c3 = dot(r[:, 3], m[:, 3]) < 0 ? -r[:, 3] : r[:, 3]
        return hcat(r[:, 1], c2, c3)
    elseif kind == "ANGLES"
        kang = _tk_key(pool, idkey, namekey, "_ANGLES")
        kaxes = _tk_key(pool, idkey, namekey, "_AXES")
        units = get(pool, _tk_key(pool, idkey, namekey, "_UNITS"), nothing)
        _badkpv(pool, "TKFRAM", kang, "=", 3, 1, 'N')
        _badkpv(pool, "TKFRAM", kaxes, "=", 3, 1, 'N')
        ang, axes = pool[kang], pool[kaxes]
        # The axes as `gipool_` reads them, before the angles are converted.
        ax = map((1, 2, 3)) do i
            a = _gipool_value(axes[i])
            a === nothing ? _gipool_range_error(kaxes, i, axes[i]) : Int(a)
        end
        # Each angle converted by `convrt_` (to degrees, then over the degrees per
        # radian), from RADIANS without a unit, as `tkfram_` does.
        u = units isa Vector{String} && !isempty(units) ? units[1] : "RADIANS"
        rad = map(x -> _convrt(x, u, "RADIANS"), (ang[1], ang[2], ang[3]))
        return _m_eul2m(Float64(rad[1]), Float64(rad[2]), Float64(rad[3]), ax[1], ax[2], ax[3])
    elseif kind == "QUATERNION"
        key = _tk_key(pool, idkey, namekey, "_Q")
        _badkpv(pool, "TKFRAM", key, "=", 4, 1, 'N')
        return _q2m(pool[key])
    end
    throw(
        SpiceError(
            :unknownframespec, "the frame specification \"$(idkey)_SPEC = '$spec'\" " *
                "is not MATRIX, ANGLES or QUATERNION"
        )
    )
end

"""
    _tk_parse(pool, frames, frcode) -> TKFrame

Relative frame and constant matrix of the TK frame with the class ID `frcode`, from the
`TKFRAME_*` variables of `pool`, as `tkfram_` of CSPICE. The variables keyed by name use
the name of the frame with the ID `frcode`, as in CSPICE. Throws the `SpiceError` of
`tkfram_`, in its order: `:zeroframeid`, `:incompleteframe` (no frame has the ID
`frcode`), `:competingframespec`, the kinds of `badkpv` for `_SPEC` and `_RELATIVE`,
`:badframespec` (unknown relative frame), `:badframespec2` (relative to itself), and
the kinds of `_tk_matrix`.
"""
function _tk_parse(pool::AbstractDict, frames::AbstractDict, frcode::Integer)
    frcode == 0 && throw(
        SpiceError(
            :zeroframeid,
            "frame identification codes are required to be non-zero"
        )
    )
    name = frame_name(frames, frcode)
    name === nothing && throw(
        SpiceError(
            :incompleteframe,
            "the Text Kernel (TK) frame with ID code $frcode does not have a recognized name"
        )
    )
    idkey = "TKFRAME_$(frcode)"
    namekey = "TKFRAME_" * name
    kspec = _tk_key(pool, idkey, namekey, "_SPEC")
    krel = _tk_key(pool, idkey, namekey, "_RELATIVE")
    _badkpv(pool, "TKFRAM", kspec, "=", 1, 1, 'C')
    _badkpv(pool, "TKFRAM", krel, "=", 1, 1, 'C')
    rel = pool[krel][1]
    relid = _tk_relative_id(pool, frames, rel)
    (relid === nothing || relid == 0) && throw(
        SpiceError(
            :badframespec, "the frame to " *
                "which frame $frcode is relatively defined is not recognized; the kernel pool " *
                "specification of the relative frame is '$rel'"
        )
    )
    relid == frcode && throw(
        SpiceError(
            :badframespec2, "bad fixed offset frame " *
                "specification: the frame '$rel' (frame ID $frcode) is defined relative to itself"
        )
    )
    return TKFrame(relid, _tk_matrix(pool, idkey, namekey, pool[kspec][1]))
end

"""
    build_tk(pool, frames) -> Dict{Int32,TKFrame}

Fixed-offset data of every class 4 frame of `frames` and of the built-in table, keyed
by frame ID, from `_tk_parse` of its class ID. A frame whose definition is missing or
malformed is left out, so that loading a kernel never fails; a query of the frame
(`tkfram`, `frame_rotation`) then parses it again and throws the error of CSPICE.
"""
function build_tk(pool::AbstractDict, frames::AbstractDict)
    out = Dict{Int32, TKFrame}()
    builtin = (fd for fd in values(BUILTIN_FRAMES) if !haskey(frames, fd.id))
    for fd in Iterators.flatten((values(frames), builtin))
        fd.class == 4 || continue
        try
            tk = _tk_parse(pool, frames, fd.classid)
            tk.relative == fd.id && continue      # relative to itself by frame ID
            out[fd.id] = tk
        catch err
            err isa SpiceError || rethrow()
        end
    end
    return out
end

"""
    _tk_missing(ks, fd)

Throw the error of the class 4 frame `fd` that `build_tk` left out: the error of
`_tk_parse`, or `SpiceError(:framedatanotfound)` if the frame is relative to itself.
"""
@noinline function _tk_missing(ks::KernelSet, fd::FrameDef)
    _tk_parse(ks.pool, ks.frames, fd.classid)
    throw(
        SpiceError(
            :framedatanotfound,
            "no valid TKFRAME_* definition is loaded for the fixed-offset frame $(fd.id)"
        )
    )
end

# --- frame chains -----------------------------------------------------------------------

"""
    _frame_step(ks, id, et::T) -> (parent::Int32, rot::Mat3{T})

One link of the frame chain: the frame that `id` is defined relative to and the matrix
from frame `id` to that parent at `et`, for a frame of class 1, 2, 4 or 7; `_chain_step`
links the frames of classes 3, 5 and 6. Requires `id != 1`. A class 4 frame without
data throws the error of `_tk_missing`.
"""
function _frame_step(ks::KernelSet, id::Int32, et::T) where {T <: Real}
    fd = frame_def(ks, id)
    class = fd.class
    if class == 1
        return _J2000_FRAME, Mat3{T}(transpose(inertial_rotation(id)))
    elseif class == 2
        return _J2000_FRAME, Mat3{T}(transpose(_body_rotation(ks, fd.classid, et)))
    elseif class == 4
        tk = get(ks.tk, id, nothing)
        tk === nothing && _tk_missing(ks, fd)
        return tk.relative, Mat3{T}(tk.rot)
    elseif class == FRAME_PROVIDER_CLASS
        p = frame_provider(ks, fd)
        p === nothing || return p.relative, provider_rotation(p, et)
    end
    throw(
        SpiceError(
            :notsupported,
            "frame $id has class $class; only classes 1 to 6 are supported"
        )
    )
end

"""
    _frame_parent(ks, id) -> Int32

The frame that `id != 1` is defined relative to, without any evaluation, for a frame
of class 1, 4 or 7 (`_ck_in_chain` stops at the other classes): the relative frame of a
class 4 frame, J2000 for class 1. Throws as `_frame_step` does.
"""
function _frame_parent(ks::KernelSet, id::Int32)
    fd = frame_def(ks, id)
    class = fd.class
    class == 1 && return _J2000_FRAME
    if class == 4
        tk = get(ks.tk, id, nothing)
        tk === nothing && _tk_missing(ks, fd)
        return tk.relative
    end
    if class == FRAME_PROVIDER_CLASS
        p = frame_provider(ks, fd)
        p === nothing || return p.relative
    end
    throw(
        SpiceError(
            :notsupported,
            "frame $id has class $class; only classes 1 to 6 are supported"
        )
    )
end

@noinline _chain_too_long(id) = throw(
    SpiceError(
        :framechain,
        "the chain of frame $id does not reach J2000 in $_FRAME_CHAIN_LIMIT links; the frame definitions may form a cycle"
    )
)

"""
    frame_rotation(ks, from, to, et::T) -> Mat3{T}

Return the matrix that rotates vectors from frame `from` to frame `to` at `et` (ET
seconds), generic over the epoch type. The frame chains are followed as `zzrefch0_` of
CSPICE follows them; `from == to` gives the identity. A CK frame is evaluated as
`ckfrot_` does; with a `Taylor` epoch the result holds the derivatives of the
interpolated rotation.
"""
function frame_rotation(ks::KernelSet, from::Int32, to::Int32, et::T) where {T <: Real}
    from == to && return one(Mat3{T})
    _both_inertial(from, to) && return _inertial_pair(T, from, to)
    return _chain_rotation(_CKExact(), ks, from, to, et)::Mat3{T}
end

frame_rotation(ks::KernelSet, from::Integer, to::Integer, et::Real) =
    frame_rotation(ks, Int32(from), Int32(to), et)

"True if `from` and `to` are built-in inertial frames (IDs 1 to 21), which no kernel redefines."
@inline _both_inertial(from::Int32, to::Int32) =
    1 <= from <= N_INERTIAL_FRAMES && 1 <= to <= N_INERTIAL_FRAMES

"""
    _inertial_pair(T, from, to) -> Mat3{T}

Return `frame_rotation` from the built-in inertial frame `from` to the built-in
inertial frame `to != from`, without the lookups of the chains: both chains are one link
to J2000, so the operations are those of the general path (`_frame_step` gives the
transposed `inertial_rotation`, and the two links meet at J2000).
"""
@inline function _inertial_pair(::Type{T}, from::Int32, to::Int32) where {T}
    from == _J2000_FRAME && return Mat3{T}(inertial_rotation(to))
    up = Mat3{T}(transpose(inertial_rotation(from)))
    to == _J2000_FRAME && return up
    return Mat3{T}(_m_mtxm(Mat3{T}(transpose(inertial_rotation(to))), up))
end

"""
    rotation_jet(ks, from, to, et::Real, Val(N)) -> NTuple{N+1,Mat3}

The matrix from frame `from` to frame `to` at `et` and its ordinary time derivatives
of orders 1 to N (element `k+1` is the derivative of order `k`, per second^k).
"""
@inline rotation_jet(ks::KernelSet, from::Integer, to::Integer, et::Real, ::Val{0}) =
    (frame_rotation(ks, Int32(from), Int32(to), et),)

function rotation_jet(ks::KernelSet, from::Integer, to::Integer, et::Real, v::Val{N}) where {N}
    N isa Integer && N >= 0 || throw(SpiceError(:badarg, "the derivative order must be an integer >= 0"))
    # Two built-in inertial frames: a constant matrix, whose derivatives are zero.
    et isa Float64 && from != to && _both_inertial(Int32(from), Int32(to)) &&
        return (_inertial_pair(Float64, Int32(from), Int32(to)), ntuple(_ -> zero(Mat3{Float64}), Val(N))...)
    r = frame_rotation(ks, Int32(from), Int32(to), taylor_var(et, v))
    d = map(derivatives, r)
    return ntuple(k -> map(x -> x[k], d), Val(N + 1))
end
