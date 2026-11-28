# Validation

Almanac.jl follows the arithmetic of CSPICE N0067 operation by operation, in the same
order of summation, and its tests call each routine and its CSPICE counterpart, through
SPICE.jl, on the same input. The inputs include NAIF kernels and kernels written by the
CSPICE writers, so that every data type is exercised. With `de440s.bsp` and all nine
aberration corrections, states and light times agree bit for bit in inertial frames and
to 4e-15 relative in body-fixed frames, and a kernel written by Almanac.jl holds the same
data, word for word, as the one that CSPICE writes from the same input. NAIF's test suite
of N0067 (tspice_c, the SPICELIB tests translated to C) is also replayed against
Almanac.jl with the correctly rounded functions of CORE-MATH on both sides, and the
replay fails on any difference from SPICELIB in the last bit that is not in its list of
known differences and that CSPICE does not share.

## Known limits of accuracy

The elementary functions (`sin`, `cos`, `exp`, ...) are those of Julia, whereas CSPICE
calls the C library of the platform; for about 3 % of arguments, `sin` and `cos` differ
in the last bit, and so can a result. The difference grows where a computation carries
it over a long span, as in a type 5 SPK segment propagated over six years (3e-8 km at
2e8 km from the Sun), and where the geometry amplifies it, as at a limb point or a
grazing intercept on a DSK shape: up to 2.5e-3 km in NAIF's test suite replayed with the
functions of Julia, at a limb point that a change of the epoch by one unit in the last
place moves by 3e-3 km in CSPICE itself. With the functions of CORE-MATH on both sides,
these results are SPICELIB's to the last bit.

Three limits come from the size of a number or from the method, and hold in CSPICE as
well:

- The rotation angle of `IAU_EARTH` reaches about 50,000 rad this century, and its
  rounding alone moves the Sun by some meters in that frame.
- An aberration-corrected position has millimeter steps at the boundaries between SPK
  [records](@ref "Smoothness at record boundaries"), which target and observer cross at
  different epochs.
- `spkezr` and `spkpos` differ by up to 2e-11 relative with `CN` or `CN+S`, because they
  stop iterating by different rules.

Where CSPICE approximates a velocity or a rate, the derivative of Almanac.jl differs from
it; the tutorial section
[Differences from `spkezr` and `sxform`](@ref "Differences from spkezr and sxform") gives
each case.

## Other checks

- **Errors:** an input that CSPICE refuses gives the kind of the CSPICE error, and a check
  that CSPICE does not make throws an error kind of Almanac.jl, such as `:badarg`.
- **Derivatives:** each order of a jet is compared with automatic differentiation of the
  order below, and every AD backend with ForwardDiff.
- **Threads:** queries from many threads, while kernels are loaded and unloaded, give the
  answers of one thread.
- **Allocation:** AllocCheck.jl proves that the vector, matrix, coordinate and table calls
  do not allocate on `Float64` static-array inputs, apart from their error paths, and the
  tests measure the allocations of [`spkezr`](@ref), [`spkpos`](@ref), [`spkder`](@ref),
  [`pxform`](@ref), [`sxform`](@ref) and [`pxder`](@ref) call by call.
- **Package hygiene:** Aqua.jl and ExplicitImports.jl.
