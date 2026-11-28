# The elementary functions of CORE-MATH in place of Julia's, for the replay with
# TSPICE_MATH=core (the default of replay.jl). The recording and the oracle take sin,
# cos, exp, ... from CORE-MATH, whose results are correctly rounded; Julia's own differ
# from them in the last bit for some arguments, and so do the outputs of Almanac that
# depend on them. Each method below is the Float64 method of a function of Base that
# calls the function of CORE-MATH in TSPICE_CRMATH_LIB (OUTDIR/oracle/lib/libcrmath.so,
# which cspice/build_cspice.sh builds), and throws a DomainError where Julia's method
# throws one (CORE-MATH returns NaN there), so that only the last bit of a result
# changes. replay.jl includes this file before `using Almanac`, and Almanac is compiled
# with these methods.
#
# The load checks that the library exports each function, that each method is the one of
# this file, and that `sin` gives the value of CORE-MATH at an argument where Julia's own
# `sin` does not; replay.jl checks the same of a function of Almanac after loading it.

using Libdl: Libdl

"The library of the functions of CORE-MATH (TSPICE_CRMATH_LIB)."
const CRMATH_LIB = get(ENV, "TSPICE_CRMATH_LIB", "")

let h = isfile(CRMATH_LIB) ? Libdl.dlopen(CRMATH_LIB; throw_error = false) : nothing
    h === nothing &&
        error("TSPICE_MATH=core, but TSPICE_CRMATH_LIB=\"$CRMATH_LIB\" is not a library: give replay.jl the environment of math_env in oracle.jl, as run.jl does, or set TSPICE_MATH=julia")
    names = [:cr_sin, :cr_cos, :cr_tan, :cr_asin, :cr_acos, :cr_atan, :cr_atan2, :cr_exp, :cr_log, :cr_pow, :cr_sinh, :cr_cosh, :cr_tanh, :cr_sincos]
    absent = filter(n -> Libdl.dlsym(h, n; throw_error = false) === nothing, names)
    isempty(absent) || error("$CRMATH_LIB does not export $(join(absent, ", "))")
end

# Julia's sin(2.5) is 0.5984721441039564, 0.503 units in the last place below the exact
# value; CORE-MATH gives the nearest double, 0.5984721441039565 (Julia 1.10 to 1.13)
const CRMATH_X = 2.5
const CRMATH_JULIA_SIN = sin(CRMATH_X)

@noinline crmath_domain(f, x) =
    throw(DomainError(x, "$f is not defined at this argument (the DomainError of Julia's method)"))

@eval begin
    Base.sin(x::Float64) =
        isinf(x) ? crmath_domain(:sin, x) : ccall((:cr_sin, $CRMATH_LIB), Float64, (Float64,), x)
    Base.cos(x::Float64) =
        isinf(x) ? crmath_domain(:cos, x) : ccall((:cr_cos, $CRMATH_LIB), Float64, (Float64,), x)
    Base.tan(x::Float64) =
        isinf(x) ? crmath_domain(:tan, x) : ccall((:cr_tan, $CRMATH_LIB), Float64, (Float64,), x)
    Base.asin(x::Float64) =
        abs(x) > 1 ? crmath_domain(:asin, x) : ccall((:cr_asin, $CRMATH_LIB), Float64, (Float64,), x)
    Base.acos(x::Float64) =
        abs(x) > 1 ? crmath_domain(:acos, x) : ccall((:cr_acos, $CRMATH_LIB), Float64, (Float64,), x)
    Base.log(x::Float64) =
        x < 0 ? crmath_domain(:log, x) : ccall((:cr_log, $CRMATH_LIB), Float64, (Float64,), x)
    Base.atan(x::Float64) = ccall((:cr_atan, $CRMATH_LIB), Float64, (Float64,), x)
    Base.exp(x::Float64) = ccall((:cr_exp, $CRMATH_LIB), Float64, (Float64,), x)
    Base.sinh(x::Float64) = ccall((:cr_sinh, $CRMATH_LIB), Float64, (Float64,), x)
    Base.cosh(x::Float64) = ccall((:cr_cosh, $CRMATH_LIB), Float64, (Float64,), x)
    Base.tanh(x::Float64) = ccall((:cr_tanh, $CRMATH_LIB), Float64, (Float64,), x)
    Base.atan(y::Float64, x::Float64) =
        ccall((:cr_atan2, $CRMATH_LIB), Float64, (Float64, Float64), y, x)
    function Base.:^(x::Float64, y::Float64)
        x < 0 && isfinite(y) && !isinteger(y) && crmath_domain(:^, x)
        return ccall((:cr_pow, $CRMATH_LIB), Float64, (Float64, Float64), x, y)
    end
    function Base.sincos(x::Float64)
        isinf(x) && crmath_domain(:sincos, x)
        s, c = Ref(0.0), Ref(0.0)
        ccall((:cr_sincos, $CRMATH_LIB), Cvoid, (Float64, Ref{Float64}, Ref{Float64}), x, s, c)
        return (s[], c[])
    end
end

let here = @__FILE__
    overridden = [
        [(f, (Float64,)) for f in (sin, cos, tan, asin, acos, atan, exp, log, sinh, cosh, tanh, sincos)];
        [(atan, (Float64, Float64)), (^, (Float64, Float64))]
    ]
    for (f, types) in overridden
        m = which(f, types)
        String(m.file) == here || error("$f$types is the method $m, not the one of $here")
    end
    cr = ccall((:cr_sin, CRMATH_LIB), Float64, (Float64,), CRMATH_X)
    CRMATH_JULIA_SIN != cr ||
        error("Julia's own sin($CRMATH_X) is the value of CORE-MATH: the check needs another argument, at which they differ")
    sin(CRMATH_X) == cr || error("sin($CRMATH_X) gives $(sin(CRMATH_X)), not the value of CORE-MATH $cr")
end
