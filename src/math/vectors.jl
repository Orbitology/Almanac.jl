# Vector and matrix calls of CSPICE: three-dimensional forms, general-dimension (`g`)
# forms, derivatives of vector functions of a state, and frame construction.
#
# Each public function has the CSPICE name, returns its results instead of filling
# output arguments, and follows the arithmetic of the CSPICE source (the `_c` file where
# CSPICE has native C code, the f2c file otherwise), so that a `Float64` input gives the
# CSPICE result. Fixed-size results are `StaticArrays` values. The functions are generic
# over the number type; every discrete decision is made on `primal`.
#
# Every private helper of the files in `src/math/` has the prefix `_m_`.

# --- three-dimensional vectors -----------------------------------------------------

"""
    vadd(v1, v2) -> SVector{3}

Return the sum `v1 + v2` of two vectors of 3 elements, as
`vadd_c` of CSPICE.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

# Examples
```jldoctest
julia> vadd([1, 2, 3], [4, 5, 6]) == [5, 7, 9]
true
```

See also [`vsub`](@ref), [`vlcom`](@ref), [`vaddg`](@ref).

# References

- [NAIF documentation of `vadd_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vadd_c.html)
"""
vadd(v1::AbstractVector, v2::AbstractVector) = _m_v3(v1) + _m_v3(v2)

"""
    vsub(v1, v2) -> SVector{3}

Return the difference `v1 - v2` of two vectors of 3 elements, as
`vsub_c` of CSPICE.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vadd`](@ref), [`vsubg`](@ref).

# References

- [NAIF documentation of `vsub_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vsub_c.html)
"""
vsub(v1::AbstractVector, v2::AbstractVector) = _m_v3(v1) - _m_v3(v2)

"""
    vscl(s, v1) -> SVector{3}

Return the product `s * v1` of a scalar and a vector of 3 elements, as
`vscl_c` of
CSPICE.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

See also [`vlcom`](@ref), [`vsclg`](@ref).

# References

- [NAIF documentation of `vscl_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vscl_c.html)
"""
vscl(s::Real, v1::AbstractVector) = _m_f(s) * _m_v3(v1)

"""
    vminus(v1) -> SVector{3}

Return `-v1` for a vector of 3 elements, as
`vminus_c` of CSPICE.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

See also [`vminug`](@ref).

# References

- [NAIF documentation of `vminus_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vminus_c.html)
"""
vminus(v1::AbstractVector) = -_m_v3(v1)

"""
    vdot(v1, v2)

Return the dot product of two vectors of 3 elements, as
`vdot_c` of CSPICE.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vcrss`](@ref), [`dvdot`](@ref), [`vdotg`](@ref).

# References

- [NAIF documentation of `vdot_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vdot_c.html)
"""
vdot(v1::AbstractVector, v2::AbstractVector) = _m_dot(_m_v3(v1), _m_v3(v2))

"""
    vcrss(v1, v2) -> SVector{3}

Return the cross product `v1 × v2` of two vectors of 3 elements, as
`vcrss_c` of CSPICE.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

# Examples
```jldoctest
julia> vcrss([1, 0, 0], [0, 1, 0]) == [0, 0, 1]
true
```

See also [`ucrss`](@ref), [`dvcrss`](@ref), [`vdot`](@ref).

# References

- [NAIF documentation of `vcrss_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vcrss_c.html)
"""
vcrss(v1::AbstractVector, v2::AbstractVector) = _m_cross(_m_v3(v1), _m_v3(v2))

"""
    vnorm(v1)

Return the Euclidean norm of a vector of 3 elements, as
`vnorm_c` of CSPICE.

The components are divided by the one largest in magnitude before they are squared, so
the squares do not overflow or underflow. The result is zero only for the zero vector.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

# Examples
```jldoctest
julia> vnorm([3, 4, 0])
5.0

julia> isfinite(vnorm([3e200, 4e200, 0]))
true
```

See also [`vhat`](@ref), [`unorm`](@ref), [`vdist`](@ref), [`vnormg`](@ref).

# References

- [NAIF documentation of `vnorm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vnorm_c.html)
"""
vnorm(v1::AbstractVector) = _m_norm(_m_v3(v1))

"""
    vhat(v1) -> SVector{3}

Return the unit vector along `v1`, as
`vhat_c` of CSPICE. The zero vector gives the
zero vector.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

See also [`unorm`](@ref), [`dvhat`](@ref), [`vhatg`](@ref).

# References

- [NAIF documentation of `vhat_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vhat_c.html)
"""
vhat(v1::AbstractVector) = _m_hat(_m_v3(v1))

"""
    unorm(v1) -> (vout::SVector{3}, vmag)

Return the unit vector along `v1` and the norm of `v1`, as
`unorm_c` of CSPICE. The
zero vector gives the zero vector and the norm 0.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

# Examples
```jldoctest
julia> unorm([3, 4, 0])
([0.6, 0.8, 0.0], 5.0)
```

See also [`vhat`](@ref), [`vnorm`](@ref), [`unormg`](@ref).

# References

- [NAIF documentation of `unorm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/unorm_c.html)
"""
unorm(v1::AbstractVector) = _m_unorm(_m_v3(v1))

"""
    vdist(v1, v2)

Return the distance `vnorm(v1 - v2)` between two points, as
`vdist_c` of CSPICE.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vnorm`](@ref), [`vrel`](@ref), [`vdistg`](@ref).

# References

- [NAIF documentation of `vdist_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vdist_c.html)
"""
vdist(v1::AbstractVector, v2::AbstractVector) = _m_norm(_m_v3(v1) - _m_v3(v2))

"""
    vsep(v1, v2)

Return the angle between two vectors of 3 elements, in radians in [0, pi], as
`vsep_c`
of CSPICE.

The method keeps the precision near 0 and pi. The result is 0 if a vector is zero.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

# Examples
```jldoctest
julia> vsep([1, 0, 0], [0, 3, 0])
1.5707963267948966
```

See also [`dvsep`](@ref), [`vsepg`](@ref).

# References

- [NAIF documentation of `vsep_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vsep_c.html)
"""
vsep(v1::AbstractVector, v2::AbstractVector) = _m_vsep(_m_v3(v1), _m_v3(v2))

"""
    vrel(v1, v2)

Return the relative difference `|v1 - v2| / max(|v1|, |v2|)` of two vectors of 3
elements, as
`vrel_c` of CSPICE. The result is 0 if the vectors are equal, also if both
are zero.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vdist`](@ref), [`vrelg`](@ref).

# References

- [NAIF documentation of `vrel_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vrel_c.html)
"""
function vrel(v1::AbstractVector, v2::AbstractVector)
    a = _m_v3(v1)
    b = _m_v3(v2)
    nunorm = _m_norm(a - b)
    primal(nunorm) == 0.0 && return nunorm
    n1 = _m_norm(a)
    n2 = _m_norm(b)
    return nunorm / (primal(n1) >= primal(n2) ? n1 : n2)
end

"""
    vproj(a, b) -> SVector{3}

Return the projection of the vector `a` on the vector `b`, as
`vproj_c` of CSPICE. The
result is the zero vector if `a` or `b` is zero.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

# Examples
```jldoctest
julia> vproj([1, 1, 0], [2, 0, 0]) == [1, 0, 0]
true

julia> vperp([1, 1, 0], [2, 0, 0]) == [0, 1, 0]
true
```

See also [`vperp`](@ref), [`vprjp`](@ref), [`vprojg`](@ref).

# References

- [NAIF documentation of `vproj_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vproj_c.html)
"""
vproj(a::AbstractVector, b::AbstractVector) = _m_proj(_m_v3(a), _m_v3(b), true)

"""
    vperp(a, b) -> SVector{3}

Return the component of the vector `a` that is perpendicular to the vector `b`, as
`vperp_c` of CSPICE.

!!! warning "Difference from CSPICE"
    For a zero `b` this call returns `a`, as the header of `vperp_c` says; the code of
    `vperp_c` in CSPICE version N0067 returns the zero vector.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vproj`](@ref), which has an example.

# References

- [NAIF documentation of `vperp_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vperp_c.html)
"""
vperp(a::AbstractVector, b::AbstractVector) = _m_perp(_m_v3(a), _m_v3(b), true)

"""
    vlcom(a, v1, b, v2) -> SVector{3}

Return the linear combination `a*v1 + b*v2` of two vectors of 3 elements, as
`vlcom_c`
of CSPICE.

Throws `SpiceError(:badarg)` if a vector does not have 3 elements.

See also [`vlcom3`](@ref), [`vlcomg`](@ref).

# References

- [NAIF documentation of `vlcom_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vlcom_c.html)
"""
vlcom(a::Real, v1::AbstractVector, b::Real, v2::AbstractVector) =
    _m_f(a) * _m_v3(v1) + _m_f(b) * _m_v3(v2)

"""
    vlcom3(a, v1, b, v2, c, v3) -> SVector{3}

Return the linear combination `a*v1 + b*v2 + c*v3` of three vectors of 3 elements, as
`vlcom3_c` of CSPICE.

Throws `SpiceError(:badarg)` if a vector does not have 3 elements.

See also [`vlcom`](@ref).

# References

- [NAIF documentation of `vlcom3_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vlcom3_c.html)
"""
vlcom3(a::Real, v1::AbstractVector, b::Real, v2::AbstractVector, c::Real, v3::AbstractVector) =
    _m_f(a) * _m_v3(v1) + _m_f(b) * _m_v3(v2) + _m_f(c) * _m_v3(v3)

"""
    vrotv(v, axis, theta) -> SVector{3}

Return the vector `v` rotated about the vector `axis` by `theta` radians, as
`vrotv_c`
of CSPICE.

The rotation follows the right-hand rule. `axis` need not have unit length; a zero
`axis` gives `v`.

Throws `SpiceError(:badarg)` if a vector does not have 3 elements.

# Examples
```jldoctest
julia> vrotv([1, 0, 0], [0, 0, 2], pi / 2) ≈ [0, 1, 0]
true
```

See also [`axisar`](@ref), [`rotvec`](@ref).

# References

- [NAIF documentation of `vrotv_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vrotv_c.html)
"""
vrotv(v::AbstractVector, axis::AbstractVector, theta::Real) =
    _m_vrotv(_m_v3(v), _m_v3(axis), _m_f(theta))

"""
    vzero(v) -> Bool

Return `true` if all 3 components of `v` are zero, as
`vzero_c` of CSPICE.

Throws `SpiceError(:badarg)` if `v` does not have 3 elements.

See also [`vzerog`](@ref).

# References

- [NAIF documentation of `vzero_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vzero_c.html)
"""
vzero(v::AbstractVector) = _m_iszero(_m_v3(v))

"""
    vequ(vin) -> SVector{3}

Return the vector `vin` of 3 elements as an `SVector{3}`, as
`vequ_c` of CSPICE.

Throws `SpiceError(:badarg)` if `vin` does not have 3 elements.

See also [`vpack`](@ref), [`vequg`](@ref).

# References

- [NAIF documentation of `vequ_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vequ_c.html)
"""
vequ(vin::AbstractVector) = _m_v3(vin)

"""
    vpack(x, y, z) -> SVector{3}

Return the vector with the components `x`, `y`, `z`, as
`vpack_c` of CSPICE.

See also [`vupack`](@ref).

# References

- [NAIF documentation of `vpack_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vpack_c.html)
"""
vpack(x::Real, y::Real, z::Real) = SVector(_m_f(x), _m_f(y), _m_f(z))

"""
    vupack(v) -> (x, y, z)

Return the three components of the vector `v` as a tuple, as
`vupack_c` of CSPICE.

Throws `SpiceError(:badarg)` if `v` does not have 3 elements.

See also [`vpack`](@ref).

# References

- [NAIF documentation of `vupack_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vupack_c.html)
"""
vupack(v::AbstractVector) = Tuple(_m_v3(v))

"""
    ucrss(v1, v2) -> SVector{3}

Return the unit vector along the cross product `v1 × v2`, as
`ucrss_c` of CSPICE.

Each vector is scaled by its component largest in magnitude before the product, so the
product does not overflow. The result is the zero vector if the vectors are parallel or
one of them is zero.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vcrss`](@ref), [`ducrss`](@ref).

# References

- [NAIF documentation of `ucrss_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ucrss_c.html)
"""
ucrss(v1::AbstractVector, v2::AbstractVector) = _m_ucrss(_m_v3(v1), _m_v3(v2))

"""
    vtmv(v1, matrix, v2)

Return the bilinear form `v1' * matrix * v2` of two vectors of 3 elements and a 3 by 3
matrix, as
`vtmv_c` of CSPICE.

Throws `SpiceError(:badarg)` if a vector does not have 3 elements or the matrix is not
3 by 3.

See also [`vtmvg`](@ref).

# References

- [NAIF documentation of `vtmv_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vtmv_c.html)
"""
function vtmv(v1::AbstractVector, matrix::AbstractMatrix, v2::AbstractVector)
    a = _m_v3(v1)
    m = _m_m3(matrix)
    b = _m_v3(v2)
    val = zero(a[1] * m[1, 1] * b[1])
    for k in 1:3, l in 1:3        # invariant: val is the sum of the terms before (k, l)
        val += a[k] * m[k, l] * b[l]
    end
    return val
end

# --- three-dimensional matrices ------------------------------------------------------

"""
    mxv(m, vin) -> SVector{3}

Return the product `m * vin` of a 3 by 3 matrix and a vector of 3 elements, as
`mxv_c`
of CSPICE.

Throws `SpiceError(:badarg)` if `m` is not 3 by 3 or `vin` does not have 3 elements.

# Examples
```jldoctest
julia> m = [0 1 0; -1 0 0; 0 0 1];

julia> mxv(m, [1, 2, 3]) == [2, -1, 3]
true

julia> mtxv(m, [1, 2, 3]) == [-2, 1, 3]
true
```

See also [`mtxv`](@ref), [`mxm`](@ref), [`mxvg`](@ref).

# References

- [NAIF documentation of `mxv_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mxv_c.html)
"""
mxv(m::AbstractMatrix, vin::AbstractVector) = _m_mxv(_m_m3(m), _m_v3(vin))

"""
    mtxv(m, vin) -> SVector{3}

Return the product of the transpose of the 3 by 3 matrix `m` and the vector `vin`, as
`mtxv_c` of CSPICE. For a rotation matrix this is the inverse rotation.

Throws `SpiceError(:badarg)` if `m` is not 3 by 3 or `vin` does not have 3 elements.

See also [`mxv`](@ref), which has an example, and [`mtxvg`](@ref).

# References

- [NAIF documentation of `mtxv_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mtxv_c.html)
"""
mtxv(m::AbstractMatrix, vin::AbstractVector) = _m_mtxv(_m_m3(m), _m_v3(vin))

"""
    mxm(m1, m2) -> SMatrix{3,3}

Return the product `m1 * m2` of two 3 by 3 matrices, as
`mxm_c` of CSPICE.

Throws `SpiceError(:badarg)` if a matrix is not 3 by 3.

See also [`mtxm`](@ref), [`mxmt`](@ref), [`mxv`](@ref), [`mxmg`](@ref).

# References

- [NAIF documentation of `mxm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mxm_c.html)
"""
mxm(m1::AbstractMatrix, m2::AbstractMatrix) = _m_mxm(_m_m3(m1), _m_m3(m2))

"""
    mtxm(m1, m2) -> SMatrix{3,3}

Return the product of the transpose of `m1` and `m2`, for 3 by 3 matrices, as
`mtxm_c`
of CSPICE.

Throws `SpiceError(:badarg)` if a matrix is not 3 by 3.

See also [`mxm`](@ref), [`mxmt`](@ref), [`mtxmg`](@ref).

# References

- [NAIF documentation of `mtxm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mtxm_c.html)
"""
mtxm(m1::AbstractMatrix, m2::AbstractMatrix) = _m_mtxm(_m_m3(m1), _m_m3(m2))

"""
    mxmt(m1, m2) -> SMatrix{3,3}

Return the product of `m1` and the transpose of `m2`, for 3 by 3 matrices, as
`mxmt_c`
of CSPICE.

Throws `SpiceError(:badarg)` if a matrix is not 3 by 3.

See also [`mxm`](@ref), [`mtxm`](@ref), [`mxmtg`](@ref).

# References

- [NAIF documentation of `mxmt_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mxmt_c.html)
"""
mxmt(m1::AbstractMatrix, m2::AbstractMatrix) = _m_mxmt(_m_m3(m1), _m_m3(m2))

"""
    xpose(m1) -> SMatrix{3,3}

Return the transpose of a 3 by 3 matrix, as
`xpose_c` of CSPICE.

Throws `SpiceError(:badarg)` if `m1` is not 3 by 3.

See also [`xpose6`](@ref), [`xposeg`](@ref).

# References

- [NAIF documentation of `xpose_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/xpose_c.html)
"""
xpose(m1::AbstractMatrix) = _m_xpose(_m_m3(m1))

"""
    xpose6(m1) -> SMatrix{6,6}

Return the transpose of a 6 by 6 matrix, as
`xpose6_c` of CSPICE.

Throws `SpiceError(:badarg)` if `m1` is not 6 by 6.

See also [`xpose`](@ref), [`xposeg`](@ref).

# References

- [NAIF documentation of `xpose6_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/xpose6_c.html)
"""
function xpose6(m1::AbstractMatrix)
    m = _m_m6(m1)
    return SMatrix{6, 6}(ntuple(k -> m[(k - 1) ÷ 6 + 1, (k - 1) % 6 + 1], Val(36)))
end

"""
    invstm(mat) -> SMatrix{6,6}

Return the inverse of a state transformation matrix, as
`invstm_c` of CSPICE.

A state transformation `[R 0; dR R]`, with a rotation `R` and its time derivative `dR`,
as [`sxform`](@ref) returns it, has the inverse
`[R' 0; dR' R']`. The call transposes each of the four 3 by 3 blocks of `mat` in its
place. It does not check that `mat` is a state transformation.

Throws `SpiceError(:badarg)` if `mat` is not 6 by 6.

# Examples
```jldoctest
julia> x = eul2xf([0.1, 0.2, 0.3, 1e-3, 2e-3, 3e-3], 3, 1, 3);

julia> invstm(x) * x ≈ [i == j for i in 1:6, j in 1:6]
true
```

See also [`xpose6`](@ref), [`invert`](@ref).

# References

- [NAIF documentation of `invstm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/invstm_c.html)
"""
function invstm(mat::AbstractMatrix)
    m = _m_m6(mat)
    # element k in column order is (i, j); its source is the transposed place in its block
    return SMatrix{6, 6}(
        ntuple(Val(36)) do k
            i, j = (k - 1) % 6, (k - 1) ÷ 6
            bi, bj = 3 * (i ÷ 3), 3 * (j ÷ 3)
            m[bi + j - bj + 1, bj + i - bi + 1]
        end
    )
end

"""
    det(m1)

Return the determinant of a 3 by 3 matrix, as
`det_c` (expansion along the first row).
Not exported, because `LinearAlgebra` exports the same name: call `Almanac.det`.

# References

- [NAIF documentation of `det_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/det_c.html)
"""
det(m1::AbstractMatrix) = _m_det(_m_m3(m1))

"""
    trace(matrix)

Return the trace of a 3 by 3 matrix, as
`trace_c` of CSPICE.

Throws `SpiceError(:badarg)` if `matrix` is not 3 by 3.

# References

- [NAIF documentation of `trace_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/trace_c.html)
"""
function trace(matrix::AbstractMatrix)
    m = _m_m3(matrix)
    return m[1, 1] + m[2, 2] + m[3, 3]
end

"""
    invert(m) -> SMatrix{3,3}

Return the inverse of a 3 by 3 matrix, as
`invert_c` of CSPICE.

A matrix whose determinant is below `1e-16` in absolute value counts as singular: the
result is then the zero matrix, and no error is thrown.

Throws `SpiceError(:badarg)` if `m` is not 3 by 3.

# Examples
```jldoctest
julia> invert([2 0 0; 0 4 0; 0 0 1]) == [0.5 0 0; 0 0.25 0; 0 0 1]
true

julia> iszero(invert([1 2 3; 2 4 6; 0 0 1]))
true
```

See also [`invort`](@ref), [`xpose`](@ref).

# References

- [NAIF documentation of `invert_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/invert_c.html)
"""
function invert(m1::AbstractMatrix)
    m = _m_m3(m1)
    mdet = _m_det(m)
    abs(primal(mdet)) < 1.0e-16 && return zero(m) * zero(mdet)
    invdet = 1.0 / mdet
    return invdet * SMatrix{3, 3}(
        (m[2, 2] * m[3, 3] - m[3, 2] * m[2, 3]),
        -(m[2, 1] * m[3, 3] - m[3, 1] * m[2, 3]),
        (m[2, 1] * m[3, 2] - m[3, 1] * m[2, 2]),
        -(m[1, 2] * m[3, 3] - m[3, 2] * m[1, 3]),
        (m[1, 1] * m[3, 3] - m[3, 1] * m[1, 3]),
        -(m[1, 1] * m[3, 2] - m[3, 1] * m[1, 2]),
        (m[1, 2] * m[2, 3] - m[2, 2] * m[1, 3]),
        -(m[1, 1] * m[2, 3] - m[2, 1] * m[1, 3]),
        (m[1, 1] * m[2, 2] - m[2, 1] * m[1, 2])
    )
end

"""
    invort(m) -> SMatrix{3,3}

Return the inverse of a 3 by 3 matrix whose columns are orthogonal and have any
non-zero length, as
`invort_c` of CSPICE.

The result is the transpose with each row divided by the squared length of the
matching column. The columns are not checked to be orthogonal.

Throws `SpiceError(:zerolengthcolumn)` if a column is zero. Throws
`SpiceError(:columntoosmall)` if the length of a column has no finite reciprocal.
Throws `SpiceError(:badarg)` if `m` is not 3 by 3.

See also [`invert`](@ref).

# References

- [NAIF documentation of `invort_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/invort_c.html)
"""
invort(m1::AbstractMatrix) = _m_invort(_m_m3(m1))

"""
    isrot(m, ntol, dtol) -> Bool

Return `true` if the 3 by 3 matrix `m` is a rotation matrix within two tolerances, as
`isrot_c` of CSPICE.

Each column must have a norm in `1 ± ntol`, and the matrix of the unit columns must
have a determinant in `1 ± dtol`.

Throws `SpiceError(:valueoutofrange)` if `ntol` or `dtol` is negative, and
`SpiceError(:badarg)` if `m` is not 3 by 3.

See also [`sharpr`](@ref).

# References

- [NAIF documentation of `isrot_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/isrot_c.html)
"""
function isrot(m::AbstractMatrix, ntol::Real, dtol::Real)
    ntol >= 0 || throw(SpiceError(:valueoutofrange, "ntol should be non-negative; it is $ntol"))
    dtol >= 0 || throw(SpiceError(:valueoutofrange, "dtol should be non-negative; it is $dtol"))
    return _m_isrot(_m_m3(m), Float64(ntol), Float64(dtol))
end

"""
    ident() -> SMatrix{3,3,Float64}

Return the 3 by 3 identity matrix, as
`ident_c` of CSPICE.

# References

- [NAIF documentation of `ident_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ident_c.html)
"""
ident() = SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

"""
    rotate(angle, iaxis) -> SMatrix{3,3}

Return the matrix that rotates a coordinate system by `angle` radians about the axis
`iaxis` (1, 2, 3 for x, y, z), as
`rotate_c` of CSPICE.

The product of the result and a vector gives the components of that vector in the
rotated system. The vector itself is not rotated: [`axisar`](@ref) and [`vrotv`](@ref)
do that, and their matrix is the transpose. As in CSPICE, any integer `iaxis` is taken
modulo 3, so 4 is x, 0 is z and -1 is y.

# Examples
```jldoctest
julia> Matrix(round.(rotate(pi / 2, 3); digits=12))
3×3 Matrix{Float64}:
  0.0  1.0  0.0
 -1.0  0.0  0.0
  0.0  0.0  1.0
```

See also [`rotmat`](@ref), [`rotvec`](@ref), [`eul2m`](@ref), [`axisar`](@ref).

# References

- [NAIF documentation of `rotate_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/rotate_c.html)
"""
rotate(angle::Real, iaxis::Integer) = _m_rotate(_m_f(angle), iaxis)

"""
    rotmat(m1, angle, iaxis) -> SMatrix{3,3}

Return `rotate(angle, iaxis) * m1` for a 3 by 3 matrix `m1`, as
`rotmat_c` of CSPICE.
This is the transformation `m1` followed by a rotation of the coordinate system by
`angle` radians about the axis `iaxis`.

Throws `SpiceError(:badarg)` if `m1` is not 3 by 3.

See also [`rotate`](@ref), [`rotvec`](@ref).

# References

- [NAIF documentation of `rotmat_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/rotmat_c.html)
"""
rotmat(m1::AbstractMatrix, angle::Real, iaxis::Integer) = _m_rotmat(_m_m3(m1), _m_f(angle), iaxis)

"""
    rotvec(v1, angle, iaxis) -> SVector{3}

Return `rotate(angle, iaxis) * v1`, as
`rotvec_c` of CSPICE: the components of the
vector `v1` in a coordinate system rotated by `angle` radians about the axis `iaxis`.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

# Examples
```jldoctest
julia> rotvec([1, 0, 0], pi / 2, 3) ≈ [0, -1, 0]
true
```

See also [`rotate`](@ref), [`rotmat`](@ref), [`vrotv`](@ref).

# References

- [NAIF documentation of `rotvec_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/rotvec_c.html)
"""
function rotvec(v1::AbstractVector, angle::Real, iaxis::Integer)
    v = _m_v3(v1)
    s, c = sincos(_m_f(angle))
    i1, i2, i3 = _m_axes(iaxis)
    a = v[i1] - zero(c)
    b = c * v[i2] + s * v[i3]
    d = -s * v[i2] + c * v[i3]
    return i1 == 1 ? SVector(a, b, d) : i1 == 2 ? SVector(d, a, b) : SVector(b, d, a)
end

@noinline function _m_check_twovec(indexa::Integer, indexp::Integer)
    (max(indexp, indexa) > 3 || min(indexp, indexa) < 1) && throw(
        SpiceError(
            :badindex,
            "the definition indexes must be in 1:3; they are $indexa and $indexp"
        )
    )
    indexa == indexp && throw(
        SpiceError(
            :undefinedframe,
            "the two definition indexes are both $indexa; they must be different"
        )
    )
    return nothing
end

"""
    twovec(axdef, indexa, plndef, indexp) -> SMatrix{3,3}

Return the rotation matrix from the input frame to the right-handed frame that two
vectors define, as
`twovec_c` of CSPICE.

Axis `indexa` of the new frame is along `axdef`. Axis `indexp` is in the plane of
`axdef` and `plndef`, on the side of `plndef`. The rows of the result are the axes of
the new frame, so the product with a vector gives its components in the new frame.

Throws:
- `SpiceError(:badindex)` if an index is not 1, 2 or 3.
- `SpiceError(:undefinedframe)` if the indexes are equal.
- `SpiceError(:dependentvectors)` if the vectors are parallel or one is zero.
- `SpiceError(:badarg)` if a vector does not have 3 elements.

# Examples
The frame with z along the input x and with x toward the input y:
```jldoctest
julia> Matrix(twovec([2, 0, 0], 3, [1, 1, 0], 1))
3×3 Matrix{Float64}:
 0.0  1.0  0.0
 0.0  0.0  1.0
 1.0  0.0  0.0
```

See also [`twovxf`](@ref).

# References

- [NAIF documentation of `twovec_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/twovec_c.html)
"""
function twovec(axdef::AbstractVector, indexa::Integer, plndef::AbstractVector, indexp::Integer)
    _m_check_twovec(indexa, indexp)
    a = _m_v3(axdef)
    p = _m_v3(plndef)
    i1 = Int(indexa)
    i2 = i1 % 3 + 1
    r1 = _m_hat(a) - zero(p)
    if indexp == i2
        r3 = _m_ucrss(a, p)
        r2 = _m_ucrss(r3, a)
    else
        r2 = _m_ucrss(p, a)
        r3 = _m_ucrss(a, r2)
    end
    _m_iszero(r2) && throw(
        SpiceError(
            :dependentvectors,
            "the vectors axdef and plndef are linearly dependent"
        )
    )
    return i1 == 1 ? _m_rows(r1, r2, r3) : i1 == 2 ? _m_rows(r3, r1, r2) : _m_rows(r2, r3, r1)
end

"""
    twovxf(axdef, indexa, plndef, indexp) -> SMatrix{6,6}

Return the 6 by 6 state transformation to the frame that two states define, as
`twovxf_c` of CSPICE.

The position parts of the states `axdef` and `plndef` define the frame as in
[`twovec`](@ref); the velocity parts give its time derivative. The upper left and lower
right 3 by 3 blocks of the result are the `twovec` matrix. The lower left block is its
time derivative. The upper right block is zero.

Throws `SpiceError(:badindex)`, `SpiceError(:undefinedframe)` and
`SpiceError(:dependentvectors)` in the cases of [`twovec`](@ref), and
`SpiceError(:badarg)` if a state does not have 6 elements.

# References

- [NAIF documentation of `twovxf_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/twovxf_c.html)
"""
function twovxf(axdef::AbstractVector, indexa::Integer, plndef::AbstractVector, indexp::Integer)
    _m_check_twovec(indexa, indexp)
    a = _m_v6(axdef)
    p = _m_v6(plndef)
    i1 = Int(indexa)
    i2 = i1 % 3 + 1
    s1 = _m_dvhat(a, false) - zero(p)
    if indexp == i2
        s3 = _m_ducrss(a, p)
        s2 = _m_ducrss(s3, a)
    else
        s2 = _m_ducrss(p, a)
        s3 = _m_ducrss(a, s2)
    end
    _m_iszero(_m_pv(s2)[1]) && throw(
        SpiceError(
            :dependentvectors,
            "the direction vectors of the states axdef and plndef are linearly dependent"
        )
    )
    x, y, z = i1 == 1 ? (s1, s2, s3) : i1 == 2 ? (s3, s1, s2) : (s2, s3, s1)
    r = _m_rows(_m_pv(x)[1], _m_pv(y)[1], _m_pv(z)[1])
    d = _m_rows(_m_pv(x)[2], _m_pv(y)[2], _m_pv(z)[2])
    return _m_xform(r, d)
end

@inline function _m_xform(r::SMatrix{3, 3}, d::SMatrix{3, 3})
    z = zero(r[1, 1] + d[1, 1])
    return SMatrix{6, 6}(
        ntuple(Val(36)) do k
            i = (k - 1) % 6 + 1
            j = (k - 1) ÷ 6 + 1
            j <= 3 ? (i <= 3 ? r[i, j] + z : d[i - 3, j] + z) : (i <= 3 ? z : r[i - 3, j - 3] + z)
        end
    )
end

"""
    frame(x) -> (x̂, y, z)

Return a right-handed orthonormal frame whose first axis is along `x`, as
`frame_c` of
CSPICE.

`x̂` is the unit vector along `x`, and `y` and `z` complete the frame; the three are
`SVector{3}` values. `y[i] == 0` for the first index `i` of the smallest
`abs(x̂[i])`. The zero vector gives the axes of the input frame. CSPICE overwrites its
argument `x` with `x̂`; this call does not change `x`.

Throws `SpiceError(:badarg)` if `x` does not have 3 elements.

# Examples
```jldoctest
julia> frame([0, 0, 2]) == ([0, 0, 1], [0, -1, 0], [1, 0, 0])
true
```

See also [`twovec`](@ref).

# References

- [NAIF documentation of `frame_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/frame_c.html)
"""
frame(x::AbstractVector) = _m_frame(_m_v3(x))

"Right-handed orthonormal frame `(x̂, y, z)` of `frame_c`, for a static vector."
@inline function _m_frame(x0::SVector{3})
    x = _m_hat(x0)
    a = x[1] * x[1]
    b = x[2] * x[2]
    c = x[3] * x[3]
    z0 = zero(a)
    o = one(a)
    primal(a + b + c) == 0.0 &&
        return SVector(o, z0, z0), SVector(z0, o, z0), SVector(z0, z0, o)
    # The smallest component (s0) is zero in y; the other two (s1, s2) are swapped, with
    # the sign of the first changed: y[s1] = -x[s2]/f, y[s2] = x[s1]/f, z[s0] = f,
    # z[s1] = -x[s0] y[s2], z[s2] = x[s0] y[s1]. The cyclic order of (s0, s1, s2) makes
    # the frame right-handed.
    a0, b0, c0 = primal(a), primal(b), primal(c)
    if a0 <= b0 && a0 <= c0                     # (s0, s1, s2) = (1, 2, 3)
        f = sqrt(b + c)
        y2 = -x[3] / f
        y3 = x[2] / f
        return x, SVector(z0, y2, y3), SVector(f, -x[1] * y3, x[1] * y2)
    elseif b0 <= a0 && b0 <= c0                 # (2, 3, 1)
        f = sqrt(a + c)
        y3 = -x[1] / f
        y1 = x[3] / f
        return x, SVector(y1, z0, y3), SVector(x[2] * y3, f, -x[2] * y1)
    end
    f = sqrt(a + b)                             # (3, 1, 2)
    y1 = -x[2] / f
    y2 = x[1] / f
    return x, SVector(y1, y2, z0), SVector(-x[3] * y2, x[3] * y1, f)
end

# --- derivatives of vector functions of a state -------------------------------------

"""
    dvhat(s1) -> SVector{6}

Return the unit vector along the position of the state `s1`, followed by the time
derivative of that unit vector, as
`dvhat_c` of CSPICE.

A state has 6 elements: a position and a velocity. A zero position gives a zero unit
vector and the velocity unchanged.

Throws `SpiceError(:badarg)` if `s1` does not have 6 elements.

# Examples
```jldoctest
julia> dvhat([2, 0, 0, 1, 2, 3]) == [1, 0, 0, 0, 1, 1.5]
true
```

See also [`vhat`](@ref), [`dvnorm`](@ref), [`ducrss`](@ref).

# References

- [NAIF documentation of `dvhat_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/dvhat_c.html)
"""
dvhat(s1::AbstractVector) = _m_dvhat(_m_v6(s1))

"""
    dvdot(s1, s2)

Return the time derivative of the dot product of the positions of two states, as
`dvdot_c` of CSPICE.

Throws `SpiceError(:badarg)` if a state does not have 6 elements.

See also [`vdot`](@ref), [`dvcrss`](@ref).

# References

- [NAIF documentation of `dvdot_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/dvdot_c.html)
"""
function dvdot(s1::AbstractVector, s2::AbstractVector)
    a = _m_v6(s1)
    b = _m_v6(s2)
    return a[1] * b[4] + a[2] * b[5] + a[3] * b[6] + a[4] * b[1] + a[5] * b[2] + a[6] * b[3]
end

"""
    dvcrss(s1, s2) -> SVector{6}

Return the cross product of the positions of two states, followed by its time
derivative, as
`dvcrss_c` of CSPICE.

Throws `SpiceError(:badarg)` if a state does not have 6 elements.

See also [`vcrss`](@ref), [`ducrss`](@ref).

# References

- [NAIF documentation of `dvcrss_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/dvcrss_c.html)
"""
dvcrss(s1::AbstractVector, s2::AbstractVector) = _m_dvcrss(_m_v6(s1), _m_v6(s2))

"""
    ducrss(s1, s2) -> SVector{6}

Return the unit vector along the cross product of the positions of two states, followed
by its time derivative, as
`ducrss_c` of CSPICE.

Each state is scaled by its largest position component first. Parallel positions give
a zero unit vector. If a scaled velocity or a product of the scaled states is not
finite, the derivative part is NaN; CSPICE returns numbers without meaning there.

Throws `SpiceError(:badarg)` if a state does not have 6 elements.

See also [`ucrss`](@ref), [`dvcrss`](@ref), [`dvhat`](@ref).

# References

- [NAIF documentation of `ducrss_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ducrss_c.html)
"""
ducrss(s1::AbstractVector, s2::AbstractVector) = _m_ducrss(_m_v6(s1), _m_v6(s2))

"""
    dvnorm(state)

Return the time derivative of the norm of the position of `state`, as
`dvnorm_c` of
CSPICE. A zero position gives 0.

Throws `SpiceError(:badarg)` if `state` does not have 6 elements.

See also [`vnorm`](@ref), [`dvhat`](@ref).

# References

- [NAIF documentation of `dvnorm_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/dvnorm_c.html)
"""
function dvnorm(state::AbstractVector)
    p, v = _m_pv(_m_v6(state))
    u, len = _m_unorm(p)
    primal(len) == 0.0 && return zero(_m_dot(v, u))
    return _m_dot(v, u)
end

"""
    dvsep(s1, s2)

Return the time derivative of the angle between the positions of two states, in rad/s
for velocities per second, as
`dvsep_c` of CSPICE. Parallel positions, or a zero position,
give 0.

Throws `SpiceError(:numericoverflow)` if the size of the result is above 1e308, and
`SpiceError(:badarg)` if a state does not have 6 elements.

See also [`vsep`](@ref).

# References

- [NAIF documentation of `dvsep_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/dvsep_c.html)
"""
function dvsep(s1::AbstractVector, s2::AbstractVector)
    u1, d1 = _m_pv(_m_dvhat(_m_v6(s1), false))
    u2, d2 = _m_pv(_m_dvhat(_m_v6(s2), false))
    pcross = _m_cross(u1, u2)
    numr = _m_dot(u1, d2) + _m_dot(d1, u2)
    _m_iszero(pcross) && return zero(numr)
    return _m_div(-numr, _m_norm(pcross))
end

# --- general-dimension forms ---------------------------------------------------------
#
# A "g" call gives a static result when every array argument is static, and a `Vector`
# or a `Matrix` that counts from 1 otherwise. The two forms do the same arithmetic.

"Elementwise `f` over `v`: a static vector if `v` and every array of `others` are static, a `Vector` otherwise."
@inline _m_gmap(f::F, v::StaticVector, others::StaticArray...) where {F} = map(f, v)
_m_gmap(f::F, v::AbstractVector, others::AbstractArray...) where {F} =
    [f(v[begin + i - 1]) for i in 1:length(v)]

"Elementwise `f` over `a` and `b` of equal length, with the result container of `_m_gmap`."
@inline _m_gmap2(f::F, a::StaticVector, b::StaticVector) where {F} = map(f, a, b)
_m_gmap2(f::F, a::AbstractVector, b::AbstractVector) where {F} =
    [f(a[begin + i - 1], b[begin + i - 1]) for i in 1:length(a)]

@inline function _m_samelen(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) ||
        throw(SpiceError(:badarg, "the vectors have $(length(a)) and $(length(b)) elements"))
    return nothing
end

"""
    vaddg(v1, v2)

Return the sum `v1 + v2` of two vectors of equal length, as
`vaddg_c` of CSPICE.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vadd`](@ref), [`vsubg`](@ref), [`vlcomg`](@ref).

# References

- [NAIF documentation of `vaddg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vaddg_c.html)
"""
vaddg(v1::AbstractVector, v2::AbstractVector) =
    (_m_samelen(v1, v2); _m_gmap2((x, y) -> _m_f(x) + _m_f(y), v1, v2))

"""
    vsubg(v1, v2)

Return the difference `v1 - v2` of two vectors of equal length, as
`vsubg_c` of
CSPICE.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vsub`](@ref), [`vaddg`](@ref).

# References

- [NAIF documentation of `vsubg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vsubg_c.html)
"""
vsubg(v1::AbstractVector, v2::AbstractVector) =
    (_m_samelen(v1, v2); _m_gmap2((x, y) -> _m_f(x) - _m_f(y), v1, v2))

"""
    vsclg(s, v1)

Return the product `s * v1` of a scalar and a vector of any length, as
`vsclg_c` of
CSPICE.

See also [`vscl`](@ref).

# References

- [NAIF documentation of `vsclg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vsclg_c.html)
"""
function vsclg(s::Real, v1::AbstractVector)
    sf = _m_f(s)
    return _m_gmap(x -> sf * _m_f(x), v1)
end

"""
    vminug(vin)

Return `-vin` for a vector of any length, as
`vminug_c` of CSPICE.

See also [`vminus`](@ref).

# References

- [NAIF documentation of `vminug_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vminug_c.html)
"""
vminug(vin::AbstractVector) = _m_gmap(x -> -_m_f(x), vin)

"""
    vequg(vin)

Return a copy of the vector `vin` of any length, as
`vequg_c` of CSPICE.

See also [`vequ`](@ref).

# References

- [NAIF documentation of `vequg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vequg_c.html)
"""
vequg(vin::AbstractVector) = _m_gmap(_m_f, vin)

"""
    vdotg(v1, v2)

Return the dot product of two vectors of equal length, as
`vdotg_c` of CSPICE.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vdot`](@ref).

# References

- [NAIF documentation of `vdotg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vdotg_c.html)
"""
function vdotg(v1::AbstractVector, v2::AbstractVector)
    _m_samelen(v1, v2)
    dot = zero(_m_f(zero(eltype(v1))) * _m_f(zero(eltype(v2))))
    @inbounds for i in eachindex(v1, v2)       # invariant: dot is the sum of terms before i
        dot += _m_f(v1[i]) * _m_f(v2[i])
    end
    return dot
end

@inline function _m_maxabsg(v::AbstractVector)
    scale = 0.0
    @inbounds for x in v                       # invariant: scale bounds the elements seen
        scale = max(scale, abs(primal(x)))
    end
    return scale
end

"""
    vnormg(v1)

Return the Euclidean norm of a vector of any length, as
`vnormg_c` of CSPICE. The
elements are divided by the one largest in magnitude before they are squared, as in
[`vnorm`](@ref).

See also [`unormg`](@ref), [`vdistg`](@ref).

# References

- [NAIF documentation of `vnormg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vnormg_c.html)
"""
function vnormg(v1::AbstractVector)
    scale = _m_maxabsg(v1)
    norm = zero(_m_f(zero(eltype(v1)))) / 1.0
    scale == 0.0 && return norm
    @inbounds for x in v1                      # invariant: norm is the sum of squares seen
        t = _m_f(x) / scale
        norm += t * t
    end
    return scale * sqrt(norm)
end

"""
    vhatg(v1)

Return the unit vector along a vector of any length, as
`vhatg_c` of CSPICE. The zero
vector gives the zero vector.

See also [`vhat`](@ref), [`unormg`](@ref).

# References

- [NAIF documentation of `vhatg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vhatg_c.html)
"""
vhatg(v1::AbstractVector) = unormg(v1)[1]

"""
    unormg(v1) -> (vout, vmag)

Return the unit vector along a vector of any length and the norm of that vector, as
`unormg_c` of CSPICE. The zero vector gives the zero vector and the norm 0.

See also [`unorm`](@ref), [`vhatg`](@ref).

# References

- [NAIF documentation of `unormg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/unormg_c.html)
"""
function unormg(v1::AbstractVector)
    vmag = vnormg(v1)
    if primal(vmag) > 0.0
        return _m_gmap(x -> _m_f(x) / vmag, v1), vmag
    end
    z = zero(vmag)
    return _m_gmap(x -> _m_f(x) * z, v1), vmag
end

"""
    vdistg(v1, v2)

Return the distance between two points of equal dimension, as
`vdistg_c` of CSPICE.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vdist`](@ref).

# References

- [NAIF documentation of `vdistg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vdistg_c.html)
"""
function vdistg(v1::AbstractVector, v2::AbstractVector)
    _m_samelen(v1, v2)
    return _m_vdistg(v1, v2)
end

"""
    vsepg(v1, v2)

Return the angle between two vectors of equal length, in radians in [0, pi], with the
method of [`vsep`](@ref). The result is 0 if a vector is zero.

`vsepg_c` of CSPICE N0067 can lose up to half of the digits near 0 and pi, so the
results of the two differ there.

Throws `SpiceError(:badarg)` if the lengths differ.

# References

- [NAIF documentation of `vsepg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vsepg_c.html)
"""
function vsepg(v1::AbstractVector, v2::AbstractVector)
    _m_samelen(v1, v2)
    mag1 = vnormg(v1)
    mag2 = vnormg(v2)
    (primal(mag1) == 0.0 || primal(mag2) == 0.0) && return zero(asin(zero(mag1 * mag2)))
    d = primal(vdotg(v1, v2))
    r1 = 1.0 / mag1
    r2 = 1.0 / mag2
    d == 0.0 && return Float64(pi) / 2 - zero(asin(zero(mag1 * mag2)))
    sgn = d > 0.0 ? -1.0 : 1.0
    acc = zero(mag1 * mag2)
    @inbounds for i in eachindex(v1, v2)       # invariant: acc is the sum of squares seen
        t = _m_f(v1[i]) * r1 + sgn * (_m_f(v2[i]) * r2)
        acc += t * t
    end
    half = 2.0 * asin(0.5 * sqrt(acc))
    return d > 0.0 ? half : Float64(pi) - half
end

"""
    vrelg(v1, v2)

Return the relative difference `|v1 - v2| / max(|v1|, |v2|)` of two vectors of equal
length, as
`vrelg_c` of CSPICE. The result is 0 for equal vectors.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vrel`](@ref).

# References

- [NAIF documentation of `vrelg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vrelg_c.html)
"""
function vrelg(v1::AbstractVector, v2::AbstractVector)
    nunorm = vdistg(v1, v2)
    primal(nunorm) == 0.0 && return nunorm
    n1 = vnormg(v1)
    n2 = vnormg(v2)
    return nunorm / (primal(n1) >= primal(n2) ? n1 : n2)
end

"""
    vprojg(a, b)

Return the projection of the vector `a` on the vector `b` of equal length, as
`vprojg_c` of CSPICE. A zero `b` gives the zero vector.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vproj`](@ref).

# References

- [NAIF documentation of `vprojg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vprojg_c.html)
"""
function vprojg(a::AbstractVector, b::AbstractVector)
    adotb = vdotg(a, b)
    bdotb = vdotg(b, b)
    scale = primal(bdotb) == 0.0 ? zero(adotb / one(bdotb)) : adotb / bdotb
    return _m_gmap(x -> scale * _m_f(x), b, a)
end

"""
    vlcomg(a, v1, b, v2)

Return the linear combination `a*v1 + b*v2` of two vectors of equal length, as
`vlcomg_c` of CSPICE.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vlcom`](@ref).

# References

- [NAIF documentation of `vlcomg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vlcomg_c.html)
"""
function vlcomg(a::Real, v1::AbstractVector, b::Real, v2::AbstractVector)
    _m_samelen(v1, v2)
    af, bf = _m_f(a), _m_f(b)
    return _m_gmap2((x, y) -> af * _m_f(x) + bf * _m_f(y), v1, v2)
end

"""
    vzerog(v) -> Bool

Return `true` if every element of the vector `v` of any length is zero, as
`vzerog_c`
of CSPICE. A vector without elements is not the zero vector, and gives `false`.

See also [`vzero`](@ref).

# References

- [NAIF documentation of `vzerog_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vzerog_c.html)
"""
vzerog(v::AbstractVector) = !isempty(v) && all(x -> primal(x) == 0.0, v)

"""
    vtmvg(v1, matrix, v2)

Return the bilinear form `v1' * matrix * v2` for any dimensions that fit, as
`vtmvg_c`
of CSPICE.

Throws `SpiceError(:badarg)` if `size(matrix) != (length(v1), length(v2))`.

See also [`vtmv`](@ref).

# References

- [NAIF documentation of `vtmvg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/vtmvg_c.html)
"""
function vtmvg(v1::AbstractVector, matrix::AbstractMatrix, v2::AbstractVector)
    size(matrix) == (length(v1), length(v2)) ||
        throw(SpiceError(:badarg, "the matrix size does not fit the vectors"))
    val = zero(_m_f(zero(eltype(v1))) * _m_f(zero(eltype(matrix))) * _m_f(zero(eltype(v2))))
    @inbounds for k in 1:length(v1), l in 1:length(v2)   # invariant: sum of terms before (k, l)
        val += _m_f(v1[begin + k - 1]) * _m_f(matrix[begin + k - 1, begin + l - 1]) * _m_f(v2[begin + l - 1])
    end
    return val
end

@inline _m_eltype(a, b) = typeof(_m_f(zero(eltype(a))) * _m_f(zero(eltype(b))))

"""
    mxvg(m1, v2)

Return the product `m1 * v2` of a matrix and a vector of any sizes that fit, as
`mxvg_c` of CSPICE.

Throws `SpiceError(:badarg)` if `size(m1, 2) != length(v2)`.

# Examples
```jldoctest
julia> mxvg([1 2 3; 4 5 6], [1, 1, 1])
2-element Vector{Float64}:
  6.0
 15.0
```

See also [`mxv`](@ref), [`mtxvg`](@ref), [`mxmg`](@ref).

# References

- [NAIF documentation of `mxvg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mxvg_c.html)
"""
function mxvg(m1::AbstractMatrix, v2::AbstractVector)
    size(m1, 2) == length(v2) || throw(SpiceError(:badarg, "the matrix and the vector do not fit"))
    out = zeros(_m_eltype(m1, v2), size(m1, 1))
    @inbounds for i in 1:size(m1, 1), k in 1:size(m1, 2)   # invariant: out[i] sums terms before k
        out[i] += _m_f(m1[begin + i - 1, begin + k - 1]) * _m_f(v2[begin + k - 1])
    end
    return out
end
"Sum of `g(1)`, ..., `g(M)` added to zero in that order: the static forms give the bits of the array forms."
@inline function _m_ordsum(g::G, ::Val{M}) where {G, M}
    t = ntuple(g, Val(M))
    return foldl(+, t; init = zero(t[1]))
end
@inline _m_ordsum(g::G, ::Val{0}) where {G} = 0.0

mxvg(m1::StaticMatrix{N, M}, v2::StaticVector{M}) where {N, M} =
    SVector{N}(ntuple(i -> _m_ordsum(k -> _m_f(m1[i, k]) * _m_f(v2[k]), Val(M)), Val(N)))

"""
    mtxvg(m1, v2)

Return the product of the transpose of `m1` and the vector `v2`, for any sizes that
fit, as
`mtxvg_c` of CSPICE.

Throws `SpiceError(:badarg)` if `size(m1, 1) != length(v2)`.

See also [`mtxv`](@ref), [`mxvg`](@ref).

# References

- [NAIF documentation of `mtxvg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mtxvg_c.html)
"""
function mtxvg(m1::AbstractMatrix, v2::AbstractVector)
    size(m1, 1) == length(v2) || throw(SpiceError(:badarg, "the matrix and the vector do not fit"))
    out = zeros(_m_eltype(m1, v2), size(m1, 2))
    @inbounds for i in 1:size(m1, 2), k in 1:size(m1, 1)   # invariant: out[i] sums terms before k
        out[i] += _m_f(m1[begin + k - 1, begin + i - 1]) * _m_f(v2[begin + k - 1])
    end
    return out
end
mtxvg(m1::StaticMatrix{N, M}, v2::StaticVector{N}) where {N, M} =
    SVector{M}(ntuple(i -> _m_ordsum(k -> _m_f(m1[k, i]) * _m_f(v2[k]), Val(N)), Val(M)))

"""
    mxmg(m1, m2)

Return the product `m1 * m2` of two matrices of any sizes that fit, as
`mxmg_c` of
CSPICE.

Throws `SpiceError(:badarg)` if `size(m1, 2) != size(m2, 1)`.

See also [`mxm`](@ref), [`mtxmg`](@ref), [`mxmtg`](@ref).

# References

- [NAIF documentation of `mxmg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mxmg_c.html)
"""
function mxmg(m1::AbstractMatrix, m2::AbstractMatrix)
    size(m1, 2) == size(m2, 1) || throw(SpiceError(:badarg, "the matrix sizes do not fit"))
    out = zeros(_m_eltype(m1, m2), size(m1, 1), size(m2, 2))
    @inbounds for i in 1:size(m1, 1), j in 1:size(m2, 2), k in 1:size(m1, 2)
        # invariant: out[i, j] is the sum of the terms before k
        out[i, j] += _m_f(m1[begin + i - 1, begin + k - 1]) * _m_f(m2[begin + k - 1, begin + j - 1])
    end
    return out
end
mxmg(m1::StaticMatrix{N, M}, m2::StaticMatrix{M, P}) where {N, M, P} =
    SMatrix{N, P}(ntuple(q -> _m_ordsum(k -> _m_f(m1[(q - 1) % N + 1, k]) * _m_f(m2[k, (q - 1) ÷ N + 1]), Val(M)), Val(N * P)))

"""
    mtxmg(m1, m2)

Return the product of the transpose of `m1` and `m2`, for any sizes that fit, as
`mtxmg_c` of CSPICE.

Throws `SpiceError(:badarg)` if `size(m1, 1) != size(m2, 1)`.

See also [`mtxm`](@ref), [`mxmg`](@ref).

# References

- [NAIF documentation of `mtxmg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mtxmg_c.html)
"""
function mtxmg(m1::AbstractMatrix, m2::AbstractMatrix)
    size(m1, 1) == size(m2, 1) || throw(SpiceError(:badarg, "the matrix sizes do not fit"))
    out = zeros(_m_eltype(m1, m2), size(m1, 2), size(m2, 2))
    @inbounds for i in 1:size(m1, 2), j in 1:size(m2, 2), k in 1:size(m1, 1)
        # invariant: out[i, j] is the sum of the terms before k
        out[i, j] += _m_f(m1[begin + k - 1, begin + i - 1]) * _m_f(m2[begin + k - 1, begin + j - 1])
    end
    return out
end
mtxmg(m1::StaticMatrix{N, M}, m2::StaticMatrix{N, P}) where {N, M, P} =
    SMatrix{M, P}(ntuple(q -> _m_ordsum(k -> _m_f(m1[k, (q - 1) % M + 1]) * _m_f(m2[k, (q - 1) ÷ M + 1]), Val(N)), Val(M * P)))

"""
    mxmtg(m1, m2)

Return the product of `m1` and the transpose of `m2`, for any sizes that fit, as
`mxmtg_c` of CSPICE.

Throws `SpiceError(:badarg)` if `size(m1, 2) != size(m2, 2)`.

See also [`mxmt`](@ref), [`mxmg`](@ref).

# References

- [NAIF documentation of `mxmtg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/mxmtg_c.html)
"""
function mxmtg(m1::AbstractMatrix, m2::AbstractMatrix)
    size(m1, 2) == size(m2, 2) || throw(SpiceError(:badarg, "the matrix sizes do not fit"))
    out = zeros(_m_eltype(m1, m2), size(m1, 1), size(m2, 1))
    @inbounds for i in 1:size(m1, 1), j in 1:size(m2, 1), k in 1:size(m1, 2)
        # invariant: out[i, j] is the sum of the terms before k
        out[i, j] += _m_f(m1[begin + i - 1, begin + k - 1]) * _m_f(m2[begin + j - 1, begin + k - 1])
    end
    return out
end
mxmtg(m1::StaticMatrix{N, M}, m2::StaticMatrix{P, M}) where {N, M, P} =
    SMatrix{N, P}(ntuple(q -> _m_ordsum(k -> _m_f(m1[(q - 1) % N + 1, k]) * _m_f(m2[(q - 1) ÷ N + 1, k]), Val(M)), Val(N * P)))

"""
    xposeg(matrix)

Return the transpose of a matrix of any size, as
`xposeg_c` of CSPICE.

See also [`xpose`](@ref), [`xpose6`](@ref).

# References

- [NAIF documentation of `xposeg_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/xposeg_c.html)
"""
xposeg(matrix::AbstractMatrix) =
    [_m_f(matrix[begin + j - 1, begin + i - 1]) for i in 1:size(matrix, 2), j in 1:size(matrix, 1)]
xposeg(matrix::StaticMatrix) = SMatrix(transpose(_m_f.(matrix)))
