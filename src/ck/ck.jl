# Orientation kernels (CK): evaluation of segments of data types 1 to 6, the segment
# search, and the CK-based reference frames (frame class 3).
#
# A CK segment gives the C-matrix of an instrument or structure: the matrix that rotates
# vectors from the reference frame of the segment to the instrument frame, as a function
# of the encoded spacecraft clock (ticks). The functions here follow the readers
# `ckrNN` and the evaluators `ckeNN` of CSPICE, but read the memory-mapped data array
# directly and keep no state, so every result is a pure function of the `KernelSet`.
#
# Everything that depends on the clock value is generic over its number type: the
# record, the interval and every other discrete choice is made with `primal`, and the
# interpolation is carried out in the given type, so a `Taylor` clock value gives the
# derivatives of the interpolated rotation.

# --- searches in a sorted part of the data array ----------------------------------------

# `_lstled` and `_lstltd` are in `spk_more.jl`.

"""
    _lstcld(data, base, n, x) -> Int

Index in `1:n` of the element of `data[base+1:base+n]` (ascending) closest to `x`, as
`lstcld_`: of two equally close the later one, and of a run of equal elements the last
one. Requires `n >= 1`.
"""
@inline function _lstcld(data::Vector{Float64}, base::Int, n::Int, x::Float64)
    @inbounds begin
        i = if x <= data[base + 1]
            1
        elseif data[base + n] <= x
            n
        else
            j = _lstltd(data, base, n, x)      # data[j] < x <= data[j+1], 1 <= j < n
            x - data[base + j] < data[base + j + 1] - x ? j : j + 1
        end
        while i < n && data[base + i] == data[base + i + 1]   # invariant: data[i] is closest
            i += 1
        end
        return i
    end
end

# --- small rotation helpers, as the CSPICE routines of the same names --------------------

"Matrix of the quaternion stored at `data[a:a+3]`, as `q2m_` (`_m_q2m` with the scaling of a non-`Float64` length)."
@inline _ck_q2m(data::Vector{Float64}, a::Int) =
    @inbounds _m_q2m(data[a], data[a + 1], data[a + 2], data[a + 3], true)

@inline _ck_none(::Type{T}) where {T} = (false, zero(Mat3{T}), zero(Vec3{T}), 0.0)

"Angular velocity stored at `data[a:a+2]` as `Vec3{T}`, or zero if `has` is false."
@inline _ck_av(::Type{T}, data::Vector{Float64}, a::Int, has::Bool) where {T} =
    has ? @inbounds(Vec3{T}(data[a], data[a + 1], data[a + 2])) : zero(Vec3{T})

# --- type 1: discrete pointing ----------------------------------------------------------

"""
    _ck_eval01(data, seg, sclk::T, tol, needav) -> (found, cmat::Mat3{T}, av::Vec3{T}, clkout)

Type 1 segment, as `ckr01_` and `cke01_`: the pointing instance closest to `sclk` among
the hundred instances of its directory group, if it is within `tol` ticks. The result is
constant in `sclk`.
"""
function _ck_eval01(
        data::Vector{Float64}, seg::CKSegment, sclk::T, tol::Float64,
        needav::Bool
    ) where {T <: Real}
    p = primal(sclk)
    beg = seg.baddr
    psiz = seg.avflag ? 7 : 4
    nrec = _cint(data[seg.eaddr])          # `(integer)` of `ckr01_`
    (nrec >= 1 && beg + (psiz + 1) * nrec - 1 <= seg.eaddr) ||
        _throw_spice(:baddaf, "CK type 1 segment with a wrong record count")
    ndir = (nrec - 1) ÷ 100
    tbase = beg + nrec * psiz - 1               # times are data[tbase+1:tbase+nrec]
    group = ndir == 0 ? 1 : _lstled(data, tbase + nrec, ndir, p) + 1
    skip = (group - 1) * 100
    n = min(100, nrec - skip)
    i = skip + _lstcld(data, tbase + skip, n, p)
    @inbounds t = data[tbase + i]
    abs(p - t) > tol && return _ck_none(T)
    a = beg + psiz * (i - 1)
    return true, Mat3{T}(_ck_q2m(data, a)), _ck_av(T, data, a + 4, seg.avflag), t
end

# --- type 2: intervals of constant angular velocity -------------------------------------

"""
    _ck_eval02(data, seg, sclk::T, tol, needav) -> (found, cmat::Mat3{T}, av::Vec3{T}, clkout)

Type 2 segment, as `ckr02_` and `cke02_`: the stored C-matrix of an interval, rotated
about the constant angular velocity by the angle `|av| * rate * (clkout - start)`.
Inside an interval `clkout` is `sclk` and the rotation carries its number type; within
`tol` of an interval the nearer end of the nearer interval is used.
"""
function _ck_eval02(
        data::Vector{Float64}, seg::CKSegment, sclk::T, tol::Float64,
        needav::Bool
    ) where {T <: Real}
    p = primal(sclk)
    beg = seg.baddr
    nrec = round(Int, ((seg.eaddr - beg + 1) * 100.0 + 1.0) / 1001.0, RoundNearestTiesAway)
    nrec >= 1 || _throw_spice(:baddaf, "CK type 2 segment with a wrong record count")
    ndir = (nrec - 1) ÷ 100
    group = ndir == 0 ? 1 : _lstled(data, beg + nrec * 10 - 1, ndir, p) + 1
    skip = (group - 1) * 100
    sb = beg + 8 * nrec - 1 + skip              # interval starts are data[sb+1:sb+n]
    n = min(100, nrec - skip)
    i = _lstled(data, sb, n, p)
    start = 0.0
    index = 0
    clkout = sclk
    @inbounds if i == 0
        p + tol >= data[sb + 1] || return _ck_none(T)
        start = data[sb + 1]
        clkout = convert(T, start)
        index = 1
    else
        stopi = data[beg + nrec * 9 + skip + i - 1]
        if p <= stopi
            start = data[sb + i]
            index = i
        elseif i == n
            p - tol <= stopi || return _ck_none(T)
            start = data[sb + i]
            clkout = convert(T, stopi)
            index = i
        else
            diff1 = p - stopi
            diff2 = data[sb + i + 1] - p
            min(diff1, diff2) <= tol || return _ck_none(T)
            if diff2 <= diff1
                start = data[sb + i + 1]
                clkout = convert(T, start)
                index = i + 1
            else
                start = data[sb + i]
                clkout = convert(T, stopi)
                index = i
            end
        end
    end
    a = beg + 8 * (skip + index - 1)            # quaternion, angular velocity, rate
    @inbounds av = Vec3{Float64}(data[a + 4], data[a + 5], data[a + 6])
    @inbounds rate = data[a + 7]
    angle = ((clkout - start) * rate) * _m_norm(av)
    # `axisar_` gives the identity for a zero axis; a series epoch keeps the zero rates
    # of `one`.
    rot = _m_iszero(av) ? one(Mat3{T}) : Mat3{T}(_m_axisar(av, angle))
    cmat = Mat3{T}(_m_mxmt(_ck_q2m(data, a), rot))
    return true, cmat, Vec3{T}(av), primal(clkout)
end

# --- type 3: linear interpolation between pointing instances ----------------------------

"Pointing instance `i` of the type 3 segment as a constant result (`t1 == t2` in `cke03_`)."
@inline function _ck_single03(
        ::Type{T}, data::Vector{Float64}, seg::CKSegment, psiz::Int,
        tbase::Int, i::Int
    ) where {T}
    a = seg.baddr + psiz * (i - 1)
    return true, Mat3{T}(_ck_q2m(data, a)), _ck_av(T, data, a + 4, seg.avflag),
        @inbounds(data[tbase + i])
end

"""
    _ck_eval03(data, seg, sclk::T, tol, needav) -> (found, cmat::Mat3{T}, av::Vec3{T}, clkout)

Type 3 segment, as `ckr03_` and `cke03_`. Between two pointing instances of one
interpolation interval the C-matrix is `C1 * R(axis, frac * angle)ᵀ`, where `R(axis,
angle) = C2ᵀ C1` and `frac` is the fraction of the clock interval; the angular velocity
is the linear combination of the two stored vectors. `frac` carries the number type of
`sclk`. Outside the interpolation intervals the nearer instance is used if it is within
`tol` ticks, as a constant.
"""
function _ck_eval03(
        data::Vector{Float64}, seg::CKSegment, sclk::T, tol::Float64,
        needav::Bool
    ) where {T <: Real}
    p = primal(sclk)
    beg, en = seg.baddr, seg.eaddr
    psiz = seg.avflag ? 7 : 4
    numint = round(Int, data[en - 1])
    numrec = round(Int, data[en])
    (numrec >= 1 && numint >= 1 && beg + (psiz + 1) * numrec + numint + 1 <= en) ||
        _throw_spice(:baddaf, "CK type 3 segment with wrong record or interval counts")
    nrdir = (numrec - 1) ÷ 100
    tbase = beg + psiz * numrec - 1             # times are data[tbase+1:tbase+numrec]
    group = nrdir == 0 ? 1 : _lstltd(data, tbase + numrec, nrdir, p) + 1
    skip = (group - 1) * 100
    n = min(100, numrec - skip)
    i = _lstled(data, tbase + skip, n, p)
    l = 0
    @inbounds if i == 0
        if group == 1
            data[tbase + 1] - p <= tol || return _ck_none(T)
            return _ck_single03(T, data, seg, psiz, tbase, 1)
        end
        l = skip
    elseif i == n
        p - data[tbase + skip + n] <= tol || return _ck_none(T)
        return _ck_single03(T, data, seg, psiz, tbase, skip + n)
    else
        l = skip + i
    end
    r = l + 1
    @inbounds lsclk, rsclk = data[tbase + l], data[tbase + r]
    # The first interval start after `p`; the two instances are joined only if the
    # right one is before it.
    sb = tbase + numrec + nrdir                 # starts are data[sb+1:sb+numint]
    j = _lstled(data, sb, numint, p)
    nstart = j < numint ? @inbounds(data[sb + j + 1]) : floatmax(Float64)
    if !(rsclk < nstart)
        ldiff = p - lsclk
        rdiff = rsclk - p
        (ldiff <= tol || rdiff <= tol) || return _ck_none(T)
        return _ck_single03(T, data, seg, psiz, tbase, ldiff < rdiff ? l : r)
    end
    a1 = beg + psiz * (l - 1)
    a2 = a1 + psiz
    frac = (sclk - lsclk) / (rsclk - lsclk)
    cmat1 = _ck_q2m(data, a1)
    cmat2 = _ck_q2m(data, a2)
    axis, angle = _m_raxisa_q(_m_m2q_unchecked(_m_mtxm(cmat2, cmat1)))
    delta = Mat3{T}(_m_axisar(axis, angle * frac))
    cmat = Mat3{T}(_m_mxmt(cmat1, delta))
    av = zero(Vec3{T})
    if seg.avflag
        @inbounds av1 = Vec3{Float64}(data[a1 + 4], data[a1 + 5], data[a1 + 6])
        @inbounds av2 = Vec3{Float64}(data[a2 + 4], data[a2 + 5], data[a2 + 6])
        av = Vec3{T}((1.0 - frac) * av1 + frac * av2)
    end
    return true, cmat, av, p
end

# --- type 4: Chebyshev polynomials over intervals of variable length ---------------------

"Value at `x` of the Chebyshev expansion with the `ncf >= 1` coefficients `data[cb+1:cb+ncf]` on the interval `mid ± rad`, as `chbval_`."
@inline function _ck_chbval(
        data::Vector{Float64}, cb::Int, ncf::Int, mid::Float64,
        rad::Float64, x::T
    ) where {T <: Real}
    s = (x - mid) / rad
    s2 = 2.0 * s
    b1 = b2 = (zero(T),)
    for j in ncf:-1:2       # invariant: b1, b2 are the Clenshaw terms of the coefficients j+1:ncf
        b1, b2 = _clenshaw_step(@inbounds(data[cb + j]), s2, b1, b2), b1
    end
    return _clenshaw_end(@inbounds(data[cb + 1]), s, b1, b2)[1]
end

"""
    _ck_meta04(data, seg) -> (refbas, nref, pdrbas, pdrtyp, pktbas, npkt, pktsz, pktoff)

Meta data of the generic segment `seg` that a type 4 CK segment needs, as `sgmeta_`
(base addresses are absolute, the item after a base is the first item). Throws
`SpiceError` unless the reference index type is 3 (explicit values, last value `<=` the
request), which is the type that the CK type 4 writer uses.
"""
@inline function _ck_meta04(data::Vector{Float64}, seg::CKSegment)
    en = seg.eaddr
    metasz = round(Int, data[en])
    metasz >= 15 || _throw_spice(:baddaf, "generic segment with fewer than 15 meta data items")
    b, hasoff = _sg_meta_base(metasz, en)       # item k is data[b+k]
    b >= seg.baddr || _throw_spice(:baddaf, "generic segment shorter than its meta data")
    item(k) = round(Int, @inbounds data[b + k])
    begm1 = seg.baddr - 1
    item(5) == 3 || _throw_spice(
        :notsupported,
        "CK type 4 segment with a reference directory type other than 3"
    )
    return item(6) + begm1, item(7), item(8) + begm1, item(10), item(11) + begm1,
        item(12), item(15), (hasoff ? item(16) : 0)
end

"Address of the first value (the interval midpoint) of the fixed-size packet `k` of a generic segment, as `sgfpkt_`."
@inline _ck04_fixed_packet(pktbas::Int, pktsz::Int, pktoff::Int, k::Int) =
    pktbas + (k - 1) * (pktsz + pktoff) + 1 + pktoff

"""
    _ck04_packet(data, meta, k) -> Int

Address of the first value of packet `k` of a type 4 segment with the meta data `meta`
of `_ck_meta04`, as `sgfpkt_`.
"""
@inline function _ck04_packet(data::Vector{Float64}, meta::NTuple{8, Int}, k::Int)
    _, _, pdrbas, pdrtyp, pktbas, _, pktsz, pktoff = meta
    return pdrtyp == 0 ? _ck04_fixed_packet(pktbas, pktsz, pktoff, k) :
        pktbas + trunc(Int, data[pdrbas + k])
end

"""
    _ck04_counts(packed, toint) -> NTuple{7,Int}

The seven coefficient counts of a type 4 packet, packed in one double base 128, as
`zzck4d2i_`; `toint` turns the quotient into an integer.
"""
@inline function _ck04_counts(packed::Float64, toint)
    x = 128.0^6
    ncf = ntuple(_ -> 0, Val(7))
    for k in 7:-1:1         # invariant: counts k+1:7 are decoded, `packed` holds the rest
        v = toint(packed / x)
        ncf = Base.setindex(ncf, v, k)
        packed -= v * x
        x /= 128.0
    end
    return ncf
end

"""
    _ck_eval04(data, seg, sclk::T, tol, needav) -> (found, cmat::Mat3{T}, av::Vec3{T}, clkout)

Type 4 segment, as `ckr04_` and `cke04_`: Chebyshev polynomials for the four quaternion
components (and the three angular velocity components) over the interval that holds
`sclk`; the quaternion is normalized. In a gap the nearer end of the nearer interval is
used if it is within `tol` ticks.
"""
function _ck_eval04(
        data::Vector{Float64}, seg::CKSegment, sclk::T, tol::Float64,
        needav::Bool
    ) where {T <: Real}
    p = primal(sclk)
    meta = _ck_meta04(data, seg)
    refbas, nref, npkt = meta[1], meta[2], meta[6]
    (npkt >= 1 && nref >= 1) || _throw_spice(:baddaf, "CK type 4 segment without records")
    packet(k) = _ck04_packet(data, meta, k)
    idx = _lstled(data, refbas, nref, p)
    clkout = sclk
    a = 0
    @inbounds if idx == 0
        a = packet(1)
        left = data[a] - data[a + 1]
        p >= left - tol || return _ck_none(T)
        clkout = convert(T, left)
    elseif idx >= npkt
        a = packet(npkt)
        right = data[a] + data[a + 1]
        p <= right + tol || return _ck_none(T)
        p >= right && (clkout = convert(T, right))
    else
        a = packet(idx)
        rbnd1 = data[a] + data[a + 1]
        if p > rbnd1
            a2 = packet(idx + 1)
            lbnd2 = data[a2] - data[a2 + 1]
            if p - rbnd1 <= lbnd2 - p
                p <= rbnd1 + tol || return _ck_none(T)
                clkout = convert(T, rbnd1)
            else
                p >= lbnd2 - tol || return _ck_none(T)
                a = a2
                clkout = convert(T, lbnd2)
            end
        end
    end
    @inbounds mid, rad, packed = data[a], data[a + 1], data[a + 2]
    ncf = _ck04_counts(packed, x -> trunc(Int, x))
    cb = a + 2                                  # coefficient j of a set is data[cb+j]
    q0 = _ck_chbval(data, cb, ncf[1], mid, rad, clkout); cb += ncf[1]
    q1 = _ck_chbval(data, cb, ncf[2], mid, rad, clkout); cb += ncf[2]
    q2 = _ck_chbval(data, cb, ncf[3], mid, rad, clkout); cb += ncf[3]
    q3 = _ck_chbval(data, cb, ncf[4], mid, rad, clkout); cb += ncf[4]
    mag = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
    cmat = primal(mag) > 0.0 ? _m_q2m(q0 / mag, q1 / mag, q2 / mag, q3 / mag, true) :
        one(Mat3{T})
    av = zero(Vec3{T})
    if seg.avflag && needav
        a1 = _ck_chbval(data, cb, ncf[5], mid, rad, clkout); cb += ncf[5]
        a2 = _ck_chbval(data, cb, ncf[6], mid, rad, clkout); cb += ncf[6]
        a3 = _ck_chbval(data, cb, ncf[7], mid, rad, clkout)
        av = Vec3{T}(a1, a2, a3)
    end
    return true, cmat, av, primal(clkout)
end

# --- types 5 and 6: Hermite and Lagrange interpolation of quaternions --------------------

"Packet sizes of the subtypes 0 to 3 of CK types 5 and 6."
const _CK5_PACKSZ = (8, 4, 14, 7)
"Largest window sizes (packets) of the subtypes 0 to 3."
const _CK5_MAXWND = (12, 24, 12, 24)
"Largest number of values that one interpolation uses (24 packets, or 12 with derivatives)."
const _CK5_MAXN = 24

"No check of the abscissae in the CK kernels: `cke05_` does not check them either."
@inline _ck_noequal(_, _) = nothing

"""
    _ck_lagrange!(w, dw, data, xb, yb, stride, n, sg, x::T) -> (p, dp)

Value and derivative at `x` of the polynomial through the `n` points with abscissae
`data[xb+j]` and ordinates `sg[j] * data[yb+(j-1)*stride]`, `j = 1:n`, by Neville's
scheme as `lgrind_` (`_neville_d!`). `w` and `dw` are work vectors of length `>= n`.
Requires `1 <= n` and distinct abscissae.
"""
@inline function _ck_lagrange!(
        w::MVector{_CK5_MAXN, T}, dw::MVector{_CK5_MAXN, T},
        data::Vector{Float64}, xb::Int, yb::Int, stride::Int, n::Int,
        sg::MVector{_CK5_MAXN, Float64}, x::T
    ) where {T <: Real}
    @inbounds for j in 1:n
        w[j] = convert(T, sg[j] * data[yb + (j - 1) * stride])
    end
    return _neville_d!(w, dw, n, i -> @inbounds(data[xb + i]), x, _ck_noequal)
end

"""
    _ck_hermite!(w, dw, data, xb, yb, doff, stride, n, rate, x::T) -> (f, df)

Value and derivative at `x` of the Hermite polynomial that has the values
`data[yb+(j-1)*stride]` and the derivatives `rate * data[yb+doff+(j-1)*stride]` at the
abscissae `data[xb+j]`, `j = 1:n`, as `hrmint_` (`_hermite_table!`). `w` and `dw` are
work vectors of length `>= 2n`. Requires `1 <= n <= 12` and distinct abscissae.
"""
@inline function _ck_hermite!(
        w::MVector{_CK5_MAXN, T}, dw::MVector{_CK5_MAXN, T},
        data::Vector{Float64}, xb::Int, yb::Int, doff::Int, stride::Int, n::Int,
        rate::Float64, x::T
    ) where {T <: Real}
    @inbounds for j in 1:n
        w[2j - 1] = convert(T, data[yb + (j - 1) * stride])
        w[2j] = convert(T, data[yb + doff + (j - 1) * stride] * rate)
    end
    return _hermite_table!(w, dw, n, i -> @inbounds(data[xb + i]), x, _ck_noequal)
end

@noinline _ck_badsign() = throw(
    SpiceError(
        :badquatsign,
        "quaternion sign error: a quaternion of a Hermite CK record is farther from its " *
            "predecessor than its negative"
    )
)
@noinline _ck_zeroquat() = throw(
    SpiceError(
        :dividebyzero,
        "the interpolated quaternion of a CK segment has magnitude zero"
    )
)

"Norm of `(s1, s2, s3, s4)` as `vnormg_`: the components are divided by the largest magnitude before they are squared."
@inline function _ck_vnormg(s1::T, s2::T, s3::T, s4::T) where {T <: Real}
    smax = max(abs(primal(s1)), abs(primal(s2)), abs(primal(s3)), abs(primal(s4)))
    smax == 0.0 && return zero(T)
    a1, a2, a3, a4 = s1 / smax, s2 / smax, s3 / smax, s4 / smax
    return smax * sqrt(a1 * a1 + a2 * a2 + a3 * a3 + a4 * a4)
end

"""
    _ck_qav(q, mags, ds, rate) -> Vec3

Angular velocity (radians per second) from the unit quaternion `q = s / mags`, where
`s` is the interpolated quaternion and `ds` its derivative with respect to ticks, as
`cke05_` computes it with `qdq2av_`; `rate` is the nominal seconds per tick.
"""
@inline function _ck_qav(q::NTuple{4, T}, mags::T, ds::NTuple{4, T}, rate::Float64) where {T <: Real}
    r = 1.0 / mags
    qds = (q[1] * ds[1] + q[2] * ds[2] + q[3] * ds[3] + q[4] * ds[4]) / mags
    dq = ntuple(i -> r * ds[i] - qds * q[i], Val(4))
    # `qdq2av_` normalizes q again (`vhatg_`); qtemp = conj(q) * dq, and the angular
    # velocity is -2 times its vector part.
    m = _ck_vnormg(q[1], q[2], q[3], q[4])
    s = q[1] / m
    v = Vec3{T}(-(q[2] / m), -(q[3] / m), -(q[4] / m))
    dv = Vec3{T}(dq[2], dq[3], dq[4])
    vec = s * dv + dq[1] * v + cross(v, dv)
    return (-2.0 / rate) * vec
end

"""
    _ck_eval05(data, pb, xb, n, subtyp, rate, t::T, needav) -> (cmat::Mat3{T}, av::Vec3{T})

Evaluate the window of `n` packets of subtype `subtyp` (0 Hermite of quaternion and its
derivative, 1 Lagrange of quaternion, 2 Hermite with angular velocity, 3 Lagrange with
angular velocity), as `cke05_`. Element `k` of packet `j` is
`data[pb+(j-1)*packsz+k]` and its epoch is `data[xb+j]`. The polynomials are evaluated
in the number type of `t` and the quaternion is normalized. `av` is zero unless
`needav`. Requires `1 <= n <= 24` (12 for the Hermite subtypes).
"""
@inline function _ck_eval05(
        data::Vector{Float64}, pb::Int, xb::Int, n::Int, subtyp::Int,
        rate::Float64, t::T, needav::Bool
    ) where {T <: Real}
    packsz = _CK5_PACKSZ[subtyp + 1]
    lagrange = subtyp == 1 || subtyp == 3
    w = MVector{_CK5_MAXN, T}(undef)
    dw = MVector{_CK5_MAXN, T}(undef)
    sg = MVector{_CK5_MAXN, Float64}(undef)
    # Signs that make each quaternion the closer of ±q to its predecessor. The
    # Lagrange subtypes apply them; the Hermite subtypes require all to be +1.
    @inbounds sg[1] = 1.0
    @inbounds for j in 2:n      # invariant: sg[1:j-1] are final
        a = pb + (j - 2) * packsz
        b = a + packsz
        dneg = 0.0
        dpos = 0.0
        for k in 1:4
            prev = sg[j - 1] * data[a + k]
            dneg += (prev + data[b + k])^2
            dpos += (prev - data[b + k])^2
        end
        flip = dneg < dpos
        flip && !lagrange && _ck_badsign()
        sg[j] = flip ? -1.0 : 1.0
    end
    if lagrange
        s1, d1 = _ck_lagrange!(w, dw, data, xb, pb + 1, packsz, n, sg, t)
        s2, d2 = _ck_lagrange!(w, dw, data, xb, pb + 2, packsz, n, sg, t)
        s3, d3 = _ck_lagrange!(w, dw, data, xb, pb + 3, packsz, n, sg, t)
        s4, d4 = _ck_lagrange!(w, dw, data, xb, pb + 4, packsz, n, sg, t)
    else
        s1, d1 = _ck_hermite!(w, dw, data, xb, pb + 1, 4, packsz, n, rate, t)
        s2, d2 = _ck_hermite!(w, dw, data, xb, pb + 2, 4, packsz, n, rate, t)
        s3, d3 = _ck_hermite!(w, dw, data, xb, pb + 3, 4, packsz, n, rate, t)
        s4, d4 = _ck_hermite!(w, dw, data, xb, pb + 4, 4, packsz, n, rate, t)
    end
    mags = _ck_vnormg(s1, s2, s3, s4)
    if subtyp == 3
        # `vhatg_`: a division by the norm, and the zero quaternion for a norm of zero,
        # which `q2m_` turns into the identity matrix
        z = zero(s1)
        q = primal(mags) == 0.0 ? (z, z, z, z) : (s1 / mags, s2 / mags, s3 / mags, s4 / mags)
    else
        # `vsclg_` by the reciprocal of the norm, which must not be zero
        primal(mags) == 0.0 && _ck_zeroquat()
        r = 1.0 / mags
        q = (r * s1, r * s2, r * s3, r * s4)
    end
    cmat = _m_q2m(q[1], q[2], q[3], q[4], true)
    av = zero(Vec3{T})
    if needav
        if subtyp == 0 || subtyp == 1
            av = _ck_qav(q, mags, (d1, d2, d3, d4), rate)
        elseif subtyp == 3
            @inbounds for k in 1:n
                sg[k] = 1.0
            end
            a1, _ = _ck_lagrange!(w, dw, data, xb, pb + 5, packsz, n, sg, t)
            a2, _ = _ck_lagrange!(w, dw, data, xb, pb + 6, packsz, n, sg, t)
            a3, _ = _ck_lagrange!(w, dw, data, xb, pb + 7, packsz, n, sg, t)
            av = Vec3{T}(a1, a2, a3)
        else
            a1, _ = _ck_hermite!(w, dw, data, xb, pb + 9, 3, packsz, n, rate, t)
            a2, _ = _ck_hermite!(w, dw, data, xb, pb + 10, 3, packsz, n, rate, t)
            a3, _ = _ck_hermite!(w, dw, data, xb, pb + 11, 3, packsz, n, rate, t)
            av = Vec3{T}(a1, a2, a3)
        end
    end
    return cmat, av
end

"""
    _ck_check05(subtyp, wndsiz, kind=:notsupported) -> Nothing

Throw `SpiceError` unless `subtyp` is in 0:3 and `wndsiz` is a valid window size for it.
A bad subtype is `kind`: `:notsupported` in `ckr05_`, `:invalidsubtype` in `ckr06_`.
"""
@inline function _ck_check05(subtyp::Int, wndsiz::Int, kind::Symbol = :notsupported)
    0 <= subtyp <= 3 || _throw_spice(kind, "unexpected subtype of a CK type 5 or 6 segment")
    (2 <= wndsiz <= _CK5_MAXWND[subtyp + 1] && iseven(wndsiz)) || _throw_spice(
        :invalidvalue,
        "window size of a CK type 5 or 6 segment is odd or outside the allowed range"
    )
    return nothing
end

"""
    _ck_eval05seg(data, seg, sclk::T, tol, needav) -> (found, cmat::Mat3{T}, av::Vec3{T}, clkout)

Type 5 segment, as `ckr05_` and `cke05_`: the window of packets around `sclk`, cut at
the interpolation interval that holds it, interpolated by `_ck_eval05`. A request
outside the packets or in a gap between two interpolation intervals gives the nearer
packet if it is within `tol` ticks, as a constant.
"""
function _ck_eval05seg(
        data::Vector{Float64}, seg::CKSegment, sclk::T, tol::Float64,
        needav::Bool
    ) where {T <: Real}
    p = primal(sclk)
    (p + tol < seg.t0 || p - tol > seg.t1) && return _ck_none(T)
    t = clamp(p, seg.t0, seg.t1)
    beg, en = seg.baddr, seg.eaddr
    @inbounds rate = data[en - 4]
    @inbounds subtyp, wndsiz = round(Int, data[en - 3]), round(Int, data[en - 2])
    @inbounds nints, n = round(Int, data[en - 1]), round(Int, data[en])
    _ck_check05(subtyp, wndsiz)
    packsz = _CK5_PACKSZ[subtyp + 1]
    eb = beg + n * packsz - 1                   # epochs are data[eb+1:eb+n]
    nidir = (nints - 1) ÷ 100
    sbase = en - 5 - nidir - nints              # interval starts are data[sbase+1:sbase+nints]
    (n >= 1 && nints >= 1 && sbase >= eb + n) ||
        _throw_spice(:baddaf, "CK type 5 segment with wrong packet or interval counts")
    @inbounds begin
        if t < data[eb + 1]
            p + tol < data[eb + 1] && return _ck_none(T)
            t = data[eb + 1]
        elseif t > data[eb + n]
            p - tol > data[eb + n] && return _ck_none(T)
            t = data[eb + n]
        end
        i = _lstltd(data, eb, n, t)
        low = i == 0 ? 1 : i
        high = i == 0 ? min(2, n) : low + 1
        lepoch, hepoch = data[eb + low], data[eb + high]
        i = _lstled(data, sbase, nints, t)
        i == 0 && return _ck_none(T)
        start = data[sbase + i]
        nstart = i < nints ? data[sbase + i + 1] : floatmax(Float64)
        nnstrt = i + 1 < nints ? data[sbase + i + 2] : floatmax(Float64)
        if hepoch == nstart
            # `low` ends an interpolation interval and `high` starts the next one.
            if abs(t - lepoch) <= abs(hepoch - t)
                abs(t - lepoch) > tol && return _ck_none(T)
                t = lepoch
                high = low
            else
                abs(hepoch - t) > tol && return _ck_none(T)
                t = hepoch
                low = high
                start = nstart
                nstart = nnstrt
            end
        end
        lsize = min(wndsiz ÷ 2, low)
        rsize = min(wndsiz ÷ 2, n - high + 1)
        first = low - lsize + 1
        i = _lstltd(data, eb + first - 1, lsize + rsize, start)
        j = _lstltd(data, eb + first - 1, lsize + rsize, nstart)
        i == j && return _ck_none(T)
        first += i
    end
    tt = t == p ? sclk : convert(T, t)
    cmat, av = _ck_eval05(
        data, beg - 1 + (first - 1) * packsz, eb + first - 1, j - i,
        subtyp, rate, tt, needav && seg.avflag
    )
    return true, cmat, av, t
end

"""
    _ck_eval06seg(data, seg, sclk::T, tol, needav) -> (found, cmat::Mat3{T}, av::Vec3{T}, clkout)

Type 6 segment, as `ckr06_` and `cke06_`: the mini-segment whose interval holds `sclk`
(the earlier one at a common bound, or the later one if the segment says so), and in it
the window of packets around `sclk`, interpolated by `_ck_eval05`. After the last
packet of a mini-segment the last packet, or the start of the next mini-segment if it
is nearer, is used if it is within `tol` ticks.

Deviation from CSPICE: when `ckr06_` moves to the next mini-segment it rounds the clock
rate of that mini-segment to a whole number; the stored rate is used here.
"""
function _ck_eval06seg(
        data::Vector{Float64}, seg::CKSegment, sclk::T, tol::Float64,
        needav::Bool
    ) where {T <: Real}
    p = primal(sclk)
    (p + tol < seg.t0 || p - tol > seg.t1) && return _ck_none(T)
    t = clamp(p, seg.t0, seg.t1)
    beg, en = seg.baddr, seg.eaddr
    @inbounds ivlsel = round(Int, data[en - 1]) == 1
    @inbounds n = round(Int, data[en])
    ndir = n ÷ 100
    ptrbas = en - 2 - (n + 1)                   # pointers are data[ptrbas+1:ptrbas+n+1]
    ivbas = ptrbas - ndir - (n + 1)             # bounds are data[ivbas+1:ivbas+n+1]
    (n >= 1 && ivbas >= beg) ||
        _throw_spice(:baddaf, "CK type 6 segment with a wrong mini-segment count")
    miniix = ivlsel ? min(_lstled(data, ivbas, n + 1, t) + 1, n + 1) - 1 :
        max(_lstltd(data, ivbas, n + 1, t), 1)
    1 <= miniix <= n || return _ck_none(T)
    while true
        # At most two passes: the second one is for the next mini-segment.
        @inbounds minib = round(Int, data[ptrbas + miniix]) + beg - 1
        @inbounds minie = round(Int, data[ptrbas + miniix + 1]) + beg - 2
        (beg <= minib && minib + 4 <= minie <= ivbas) ||
            _throw_spice(:baddaf, "CK type 6 segment with a wrong mini-segment pointer")
        @inbounds rate = data[minie - 3]
        @inbounds subtyp, wndsiz = round(Int, data[minie - 2]), round(Int, data[minie - 1])
        @inbounds npkt = round(Int, data[minie])
        npkdir = (npkt - 1) ÷ 100
        @inbounds lstepc = data[minie - 4 - npkdir]
        if t > lstepc
            tol == 0.0 && return _ck_none(T)
            @inbounds b1 = data[ivbas + miniix + 1]
            (t - lstepc > tol && b1 - t > tol) && return _ck_none(T)
            if miniix == n
                t - lstepc <= tol || return _ck_none(T)
                t = lstepc
            elseif t - lstepc <= (b1 - lstepc) / 2
                t = lstepc
            else
                miniix += 1
                t = b1
                continue
            end
        end
        _ck_check05(subtyp, wndsiz, :invalidsubtype)
        packsz = _CK5_PACKSZ[subtyp + 1]
        eb = minib - 1 + npkt * packsz          # epochs are data[eb+1:eb+npkt]
        (npkt >= 1 && eb + npkt <= minie) ||
            _throw_spice(:baddaf, "CK type 6 mini-segment with a wrong packet count")
        i = _lstltd(data, eb, npkt, t)
        low = i == 0 ? 1 : i
        high = low + 1
        lsize = min(wndsiz ÷ 2, low)
        rsize = min(wndsiz ÷ 2, npkt - high + 1)
        first = low - lsize + 1
        tt = t == p ? sclk : convert(T, t)
        cmat, av = _ck_eval05(
            data, minib - 1 + (first - 1) * packsz, eb + first - 1,
            lsize + rsize, subtyp, rate, tt, needav && seg.avflag
        )
        return true, cmat, av, t
    end
    return
end

# --- one segment, and the segment search -------------------------------------------------

"""
    _ck_eval(data, seg, sclk::T, tol, needav) -> (found, cmat::Mat3{T}, av::Vec3{T}, clkout::Float64)

Pointing of the CK segment `seg` (data array `data`) for the encoded SCLK `sclk` with
the tolerance `tol >= 0` ticks, as `ckpfs_`: `cmat` rotates vectors from the frame
`seg.frame` to the instrument frame, `av` is the angular velocity of the instrument
frame in `seg.frame` (radians per second, zero if the segment has none or `needav` is
false and the type computes it separately), and `clkout` is the clock time of the
result. `found` is false if the segment has no pointing within `tol` of `sclk`.
`cmat` and `av` are functions of `sclk` in its number type where the data type
interpolates, and constants where it takes a stored instance. Throws `SpiceError(:ckunknowndatatype)` for a data type outside 1:6.
"""
@inline function _ck_eval(
        data::Vector{Float64}, seg::CKSegment, sclk::T, tol::Float64,
        needav::Bool
    ) where {T <: Real}
    dtype = seg.dtype
    dtype == 3 && return _ck_eval03(data, seg, sclk, tol, needav)
    dtype == 5 && return _ck_eval05seg(data, seg, sclk, tol, needav)
    dtype == 6 && return _ck_eval06seg(data, seg, sclk, tol, needav)
    dtype == 2 && return _ck_eval02(data, seg, sclk, tol, needav)
    dtype == 1 && return _ck_eval01(data, seg, sclk, tol, needav)
    dtype == 4 && return _ck_eval04(data, seg, sclk, tol, needav)
    return _throw_spice(:ckunknowndatatype, "CK data type $dtype is not supported")
end

"""
    _ck_search(ks, inst, sclk::T, tol, needav) -> (found, cmat, av, clkout, ref::Int32)

Pointing of the instrument `inst` for the encoded SCLK `sclk`, from the first segment in
priority order that gives one, as the search of `ckbss_`/`cksns_` with `ckpfs_`: a
segment is tried if its coverage meets `[sclk - tol, sclk + tol]` and, if `needav`, it
has angular velocity. `ref` is the reference frame of that segment (0 if not found).
Requires `tol >= 0`.
"""
function _ck_search(
        ks::KernelSet, inst::Int32, sclk::T, tol::Float64,
        needav::Bool
    ) where {T <: Real}
    segs = get(ks.ck, inst, nothing)
    if segs !== nothing
        p = primal(sclk)
        alpha = p - tol
        omega = p + tol
        # With `tol == 0` a segment is tried only if it covers `p`, so the search starts at
        # the first one that does (0 for none: no segment is tried).
        start = tol == 0.0 ? find_segment(segs, p) : 1
        for k in (start == 0 ? (1:0) : (start:length(segs)))
            # invariant: no segment before `k` gives a pointing
            seg = @inbounds segs[k]
            (seg.t0 <= seg.t1 && omega >= seg.t0 && alpha <= seg.t1) || continue
            (!needav || seg.avflag) || continue
            found, cmat, av, clkout = _ck_eval(ks.daf[seg.file], seg, sclk, tol, needav)
            found && return true, cmat, av, clkout, seg.frame
        end
    end
    return false, zero(Mat3{T}), zero(Vec3{T}), 0.0, Int32(0)
end

# --- CK meta data -------------------------------------------------------------------------

@inline _ckmeta_default(ckid::Integer) = ckid <= -1000 ? Int32(ckid ÷ 1000) : Int32(0)

"""
    build_ckmeta(pool) -> Dict{Int32,NTuple{2,Int32}}

The pair (clock ID, SPK ID) for each CK ID that has a numeric variable `CK_<id>_SCLK`
or `CK_<id>_SPK` in the pool (its first value, rounded by `gipool_`); the one of the two
that is absent holds the default of `ckmeta_`.
"""
function build_ckmeta(pool::AbstractDict)
    out = Dict{Int32, NTuple{2, Int32}}()
    for (name, value) in pool
        m = match(r"^CK_(-?\d+)_(SCLK|SPK)$", name)
        m === nothing && continue
        id = tryparse(Int32, m.captures[1])
        id === nothing && continue
        # `gipool_(name, 1, 1, ...)`: the first value, rounded
        v = _gipool_first(value, name; strict = false)
        v === nothing && continue
        old = get(out, id, (_ckmeta_default(id), _ckmeta_default(id)))
        out[id] = m.captures[2] == "SCLK" ? (v, old[2]) : (old[1], v)
    end
    return out
end

"Clock ID of the CK ID `ckid` (`ckmeta_` with `SCLK`)."
@inline function _ck_sclkid(ks::KernelSet, ckid::Int32)
    m = get(ks.ckmeta, ckid, nothing)
    return m === nothing ? _ckmeta_default(ckid) : m[1]
end

# --- frame class 3 ------------------------------------------------------------------------

"Mode of a CK frame evaluation: the rotation in the number type of the epoch."
struct _CKExact end
"Mode of a CK frame evaluation: rotation and CSPICE angular velocity at a `Float64` epoch."
struct _CKState end

@noinline _ck_noconnect(fd::FrameDef, et, why::String) = throw(
    SpiceError(
        :noframeconnect,
        "the CK frame $(fd.id) (CK ID $(fd.classid)) has no orientation at ET $(primal(et)): $why"
    )
)

"Clock data of the CK frame `fd`; throws `SpiceError(:noframeconnect)` if no SCLK kernel is loaded for it."
@inline function _ck_frame_clock(ks::KernelSet, fd::FrameDef, et)
    clock = get(ks.sclk, _ck_sclkid(ks, fd.classid), nothing)
    clock === nothing && _ck_noconnect(
        fd, et,
        "no type 1 SCLK data are loaded for the clock $(_ck_sclkid(ks, fd.classid))"
    )
    return clock
end

"""
    _ck_frame_step(mode, ks, fd, et) -> (ref::Int32, rot)

One link of a frame chain for the class 3 frame `fd`: the reference frame `ref` of the
CK segment that gives the orientation at `et`, and the matrix from frame `fd` to `ref`.

With `_CKExact()` this is `ckfrot_`: tolerance zero, no angular velocity needed, and
`rot::Mat3{T}` is the interpolated rotation as a function of `et::T`, through the
piecewise linear clock conversion `sce2c`.

With `_CKState()` this is `ckfxfm_`: tolerance zero, only segments with angular
velocity, `et::Float64`, and `rot::Mat3{Taylor{1,Float64}}` holds the rotation `R` and
the derivative `[av]ₓ R` that CSPICE builds from the angular velocity `av` of the
segment (`rav2xf_`, `invstm_`).

Throws `SpiceError(:noframeconnect)` if no segment gives a pointing.
"""
@inline function _ck_frame_step(mode, ks::KernelSet, fd::FrameDef, et::Real)
    found, ref, rot = _ck_frame_step_found(mode, ks, fd, et)
    found || _ck_noconnect(
        fd, et, mode isa _CKState ? "no CK segment with angular velocity covers this epoch" :
            "no CK segment covers this epoch"
    )
    return ref, rot
end

"`_ck_frame_step` with `found` first: false (and the identity) where it throws."
@inline function _ck_frame_step_found(::_CKExact, ks::KernelSet, fd::FrameDef, et::T) where {T <: Real}
    ticks = _sce2c(ks, _ck_frame_clock(ks, fd, et), et)
    found, cmat, _, _, ref = _ck_search(ks, fd.classid, ticks, 0.0, false)
    found || return false, fd.id, one(Mat3{T})
    return true, ref, Mat3{T}(transpose(cmat))
end

@inline function _ck_frame_step_found(::_CKState, ks::KernelSet, fd::FrameDef, et::Float64)
    ticks = _sce2c(ks, _ck_frame_clock(ks, fd, et), et)
    found, cmat, av, _, ref = _ck_search(ks, fd.classid, ticks, 0.0, true)
    found || return false, fd.id, one(Mat3{Taylor{1, Float64}})
    r = Mat3{Float64}(transpose(cmat))
    dr = _m_mxm(_ck_skew(av), r)
    return true, ref, map((a, b) -> Taylor((a, b)), r, dr)
end

"Cross-product matrix `[v]ₓ`, so that `[v]ₓ w = v × w`."
@inline _ck_skew(v::Vec3{T}) where {T} =
    Mat3{T}(zero(T), v[3], -v[2], -v[3], zero(T), v[1], v[2], -v[1], zero(T))

"""
    _sx_mxm(a, b) -> Mat3{Taylor{1,Float64}}

The product `a * b` of two rotations with their derivatives, with the
arithmetic of the 6 by 6 product of their state transformations in `zzmsxf_` and
`mxmg_`: the value is the three products of `a.v` and `b.v` summed from left to right;
the derivative is `a.d * b.v` and then `a.v * b.d`, the six products summed from left to
right. The product of `Taylor` numbers sums the terms in another order, which changes
the last bit.
"""
@inline function _sx_mxm(a::Mat3{Taylor{1, Float64}}, b::Mat3{Taylor{1, Float64}})
    return Mat3{Taylor{1, Float64}}(
        ntuple(Val(9)) do k
            i, j = _m_rowcol(k)
            v = a[i, 1].v * b[1, j].v + a[i, 2].v * b[2, j].v + a[i, 3].v * b[3, j].v
            d = a[i, 1].d[1] * b[1, j].v + a[i, 2].d[1] * b[2, j].v + a[i, 3].d[1] * b[3, j].v +
                a[i, 1].v * b[1, j].d[1] + a[i, 2].v * b[2, j].d[1] + a[i, 3].v * b[3, j].d[1]
            Taylor((v, d))
        end
    )
end

"""
    _sx_mxv(m, s) -> Vec3{Taylor{1,Float64}}

The product of the rotation `m` with its derivative and the vector `s` with
its derivative, with the arithmetic of `mxvg_` applied to the 6 by 6 state
transformation and the state: the value is the three products of `m.v` and `s.v`, the
derivative `m.d * s.v` and then `m.v * s.d`, each sum from left to right.
"""
@inline function _sx_mxv(m::Mat3{Taylor{1, Float64}}, s::Vec3{Taylor{1, Float64}})
    return Vec3{Taylor{1, Float64}}(
        ntuple(Val(3)) do i
            v = m[i, 1].v * s[1].v + m[i, 2].v * s[2].v + m[i, 3].v * s[3].v
            d = m[i, 1].d[1] * s[1].v + m[i, 2].d[1] * s[2].v + m[i, 3].d[1] * s[3].v +
                m[i, 1].v * s[1].d[1] + m[i, 2].v * s[2].d[1] + m[i, 3].v * s[3].d[1]
            Taylor((v, d))
        end
    )
end

"`_sx_mxm(transpose(a), b)`: the product with the inverse of `a`, as `invstm_` and `zzmsxf_`."
@inline _sx_mtxm(a::Mat3{Taylor{1, Float64}}, b::Mat3{Taylor{1, Float64}}) =
    _sx_mxm(Mat3{Taylor{1, Float64}}(transpose(a)), b)

"Matrix product of two links of a chain: `_sx_mxm` in the state mode, `_m_mxm` in the exact mode."
@inline _ck_mxm(::_CKState, a, b) = _sx_mxm(a, b)
@inline _ck_mxm(::_CKExact, a, b) = _m_mxm(a, b)

@inline _ck_epoch(::_CKExact, et) = et
@inline _ck_epoch(::_CKState, et::Float64) = taylor_var(et, Val(1))
@inline _ck_numtype(::_CKExact, ::Type{T}) where {T} = T
@inline _ck_numtype(::_CKState, ::Type{Float64}) = Taylor{1, Float64}

"""
    _ck_pck_step(mode, ks, fd, et) -> (J2000, Mat3)

One link of a frame chain for the PCK frame `fd`: the matrix to J2000, in the state mode
with the derivative of `tisbod_` (`_body_state`), in the exact mode as `_frame_step`.
"""
@inline _ck_pck_step(::_CKExact, ks::KernelSet, fd::FrameDef, et) = _frame_step(ks, fd.id, et)
@inline _ck_pck_step(::_CKState, ks::KernelSet, fd::FrameDef, et::Float64) =
    (_J2000_FRAME, Mat3{Taylor{1, Float64}}(transpose(_body_state(ks, fd.classid, et))))

"""
    _chain_step(mode, ks, id, et) -> (parent, Mat3)

One link of a frame chain as `zzrotgt0_` (exact mode) and `zzfrmgt0_` (state mode) give
it: a CK frame is linked to the reference frame of its segment for `et`, a switch frame
as `_sw_frame_step` links it, every other frame as `_frame_step` links it.
"""
@inline function _chain_step(mode, ks::KernelSet, id::Int32, et::T) where {T <: Real}
    fd = frame_def(ks, id)
    if fd.class == 3
        return _ck_frame_step(mode, ks, fd, et)
    elseif fd.class == 5 || fd.class == 6     # dynamic and switch frames, `dynframes.jl`
        return _dyn_chain_step(mode, ks, fd, et)
    elseif fd.class == 2
        return _ck_pck_step(mode, ks, fd, et)
    end
    return _frame_step(ks, id, _ck_epoch(mode, et))
end

"""
    _chain_step_found(mode, ks, id, et) -> (found, parent, Mat3)

`_chain_step`, but a CK frame without a segment for `et` gives `found` false instead of
the error, as `zzrotgt0_` and `zzfrmgt0_` report it (`id` and the identity).
"""
@inline function _chain_step_found(mode, ks::KernelSet, id::Int32, et::T) where {T <: Real}
    fd = frame_def(ks, id)
    fd.class == 3 && return _ck_frame_step_found(mode, ks, fd, et)
    parent, rot = _chain_step(mode, ks, id, et)
    return true, parent, rot
end

"`transpose(a) * b` in the mode: `_sx_mtxm` in the state mode, `_m_mtxm` in the exact mode."
@inline _ck_mtxm(::_CKState, a, b) = _sx_mtxm(a, b)
@inline _ck_mtxm(::_CKExact, a, b) = _m_mtxm(a, b)

"Throw `SpiceError(:unknownframe)` unless `from`, then `to`, is a known frame, as `zzrefch0_` checks them before the first link."
@inline function _check_frame_pair(ks::KernelSet, from::Int32, to::Int32)
    frame_def(ks, from)
    frame_def(ks, to)
    return nothing
end

"`(k, rots[k])` for the first `k <= len` with `nodes[k] == id`, or `(0, rots[1])`."
@inline function _chain_find(nodes::NTuple{M, Int32}, rots, len::Int, id::Int32) where {M}
    for k in 1:len
        @inbounds nodes[k] == id && return k, rots[k]
    end
    return 0, rots[1]
end

"""
    _chain_rotation(mode, ks, from, to, et) -> Mat3

Matrix from frame `from` to frame `to` at `et` as `zzrefch0_` (exact mode) and
`zzfrmch0_` (state mode) compute it: the chain of `from` is evaluated link by link, so
that a CK or switch frame is followed through the frame its segment or base gives for
`et`, until it reaches `to` or J2000; the chain of `to` is then evaluated until it meets
a frame of that chain, and the two parts are joined there. The first frame met can be
other than J2000 only when it is the frame of a CK segment or of a switch base (the
static parents of those frames are J2000). As in `zzrefch0_`, the chain of `from` keeps
`MAX_FRAME_CHAIN` frames and folds the links after them into the last one, and the chain
of `to` meets only the frames kept; each chain throws `SpiceError(:framechain)` after
`_FRAME_CHAIN_LIMIT` links, where CSPICE loops forever.
"""
function _chain_rotation(mode, ks::KernelSet, from::Int32, to::Int32, et::T) where {T <: Real}
    S = _ck_numtype(mode, T)
    from == to && return one(Mat3{S})
    _check_frame_pair(ks, from, to)
    nodes = ntuple(_ -> Int32(0), Val(MAX_FRAME_CHAIN))
    rots = ntuple(_ -> one(Mat3{S}), Val(MAX_FRAME_CHAIN))
    nodes = Base.setindex(nodes, from, 1)
    len = 1
    links = 0
    cur = from
    rot = one(Mat3{S})
    tail = rot
    found = true
    # invariant: nodes[1:len] is the chain of `from` to cur, with the frames after
    # nodes[MAX_FRAME_CHAIN - 1] but the last folded away; rots[k] is the matrix from
    # `from` to nodes[k], rot = rots[len], and once len == MAX_FRAME_CHAIN, tail is the
    # matrix from nodes[len - 1] to cur, which `zzrefch0_` folds into its last link and
    # multiplies last. The chain ends at a CK frame without a segment for `et` (found
    # false), as in `zzrefch0_`.
    while cur != _J2000_FRAME && cur != to
        links += 1
        links > _FRAME_CHAIN_LIMIT && _chain_too_long(from)
        found, cur, step = _chain_step_found(mode, ks, cur, et)
        found || break
        if len < MAX_FRAME_CHAIN
            rot = len == 1 ? step : _ck_mxm(mode, step, rot)
            len += 1
            tail = step
        else
            tail = _ck_mxm(mode, step, tail)
            rot = _ck_mxm(mode, tail, rots[len - 1])
        end
        nodes = Base.setindex(nodes, cur, len)
        rots = Base.setindex(rots, rot, len)
    end
    cur == to && return rot::Mat3{S}
    down = one(Mat3{S})
    cur = to
    n = 0
    k, up = _chain_find(nodes, rots, len, cur)
    # invariant: down is the matrix from `to` to cur, and no frame of the chain of `to`
    # before cur is on the chain of `from`
    while k == 0
        # no common frame: the error of the link where the chain of `from` ended
        cur == _J2000_FRAME && _chain_step(mode, ks, nodes[len], et)
        n += 1
        n > _FRAME_CHAIN_LIMIT && _chain_too_long(to)
        cur, step = _chain_step(mode, ks, cur, et)
        down = n == 1 ? step : _ck_mxm(mode, step, down)
        k, up = _chain_find(nodes, rots, len, cur)
    end
    k == 1 && return Mat3{S}(transpose(down))
    return Mat3{S}(_ck_mtxm(mode, down, up))
end

"""
    _ck_in_chain(ks, id) -> Bool

True if the chain of frame `id` towards J2000 holds a frame of class 3, 5 or 6 (the
classes whose CSPICE state transformation is not the time derivative of the rotation)
or of class 2 (whose derivative `tisbod_` computes with other arithmetic), before any
other link needs an evaluation. Throws as `_frame_parent` does.
"""
function _ck_in_chain(ks::KernelSet, id::Int32)
    cur = id
    n = 0
    while cur != _J2000_FRAME       # invariant: no frame from id to before cur is class 2, 3, 5 or 6
        n += 1
        n > _FRAME_CHAIN_LIMIT && _chain_too_long(id)
        class = frame_def(ks, cur).class
        (class == 2 || class == 3 || class == 5 || class == 6) && return true
        cur = _frame_parent(ks, cur)
    end
    return false
end

"""
    _ck_state_rotation(ks, from, to, et::Float64) -> Mat3{Taylor{1,Float64}}

Matrix from frame `from` to frame `to` at `et` and its derivative as `frmchg_` of
CSPICE gives them (the two blocks of `sxform_c`): each CK link contributes the angular
velocity of its segment, and needs a segment that has one. The two chains are followed
to their first common frame, as in `frame_rotation`.
"""
function _ck_state_rotation(ks::KernelSet, from::Int32, to::Int32, et::Float64)
    return _chain_rotation(_CKState(), ks, from, to, et)
end
