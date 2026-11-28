# Group 9: aberration corrections against CSPICE.
#
# Error measure: the Euclidean norm of the difference from the CSPICE vector divided by
# the norm of the CSPICE vector (position and velocity apart), and the relative
# difference of the light time. The table that is printed has the largest value of each
# over all pairs and epochs, per frame and correction, with `NONE` as the baseline.
#
# Why the agreement is near 1e-13 and not 1e-16. A light-time corrected state is the
# difference of two barycentric vectors (target at `et ± lt`, observer at `et`), as in
# `spkltc`. One unit in the last place of a barycentric coordinate (1.5e8 km for the
# Earth and the Moon) is 3e-8 km, which is 8e-14 of the Earth-Moon distance. The two
# libraries sum the Chebyshev series in a different order, so they differ by such units.
#
# Derivatives. `spkder` with a correction returns exact derivatives of the position
# function that CSPICE computes, so it is checked against the central difference of the
# order below (tolerances as in group 3) and against nested `ForwardDiff` of `spkpos`.
# The difference between order 1 and the CSPICE velocity is measured and printed.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Almanac keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALMANAC_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff

const AB_CORRECTIONS = ("LT", "LT+S", "CN", "CN+S", "XLT", "XLT+S", "XCN", "XCN+S")

"Thresholds on the largest relative difference of position, velocity and light time."
# Velocity with `+S`: `spkacs` takes the observer acceleration from the difference of two
# velocities 2 s apart. One unit in the last place of a velocity of 30 km/s moves that
# acceleration by 3e-15 km/s^2, and the derivative of the stellar aberration multiplies
# it by `|r| / c` (15000 s for Neptune): 5e-11 km/s, or 2e-12 of the velocity.
const AB_TOL_INERTIAL = (pos = 2.0e-12, vel = 5.0e-11, lt = 2.0e-12)
# Body-fixed frames: the rounding of the prime-meridian angle far from J2000 (see
# `TOL_BODYFIXED` of group 1) is on top.
const AB_TOL_BODYFIXED = (pos = 1.0e-9, vel = 1.0e-9, lt = 2.0e-12)

const AB_FD_INERTIAL = (h = 2000.0, rtol = 2.0e-4)
const AB_FD_FAST = (h = 20.0, rtol = 2.0e-5)
"Tolerance of `spkder` against nested ForwardDiff of `spkpos`, orders 0 to 3."
const AB_AD_RTOL = (1.0e-13, 1.0e-11, 1.0e-10, 1.0e-9)

"Relative difference `|a - b| / |b|`; for `b == 0` the absolute difference."
function ab_reldiff(a, b)
    nb = norm(b)
    return nb == 0 ? norm(a - b) : norm(a - b) / nb
end

"""
    ab_compare(ctx, pairs, epochs, ref, abcorr) -> NamedTuple

Largest relative differences from CSPICE of `spkezr` (position, velocity, light time)
and of `spkpos` (position, light time) over `pairs` and `epochs`.
Precondition: CSPICE has the same kernels loaded as `ctx`.
"""
function ab_compare(ctx, pairs, epochs, ref::AbstractString, abcorr::AbstractString)
    @assert !isempty(pairs) && !isempty(epochs)
    pos = vel = lt = ppos = plt = 0.0
    for (targ, obs) in pairs, et in epochs
        cs, cl = SPICE.spkezr(string(targ), et, ref, abcorr, string(obs))
        s, l = spkezr(ctx, targ, et, ref, abcorr, obs)
        pos = max(pos, ab_reldiff(s[SOneTo(3)], cs[1:3]))
        vel = max(vel, ab_reldiff(s[SVector(4, 5, 6)], cs[4:6]))
        lt = max(lt, abs(l - cl) / cl)
        cp, cpl = SPICE.spkpos(string(targ), et, ref, abcorr, string(obs))
        p, pl = spkpos(ctx, targ, et, ref, abcorr, obs)
        ppos = max(ppos, ab_reldiff(p, cp))
        plt = max(plt, abs(pl - cpl) / cpl)
    end
    return (; pos, vel, lt, ppos, plt)
end

"Nested forward-mode derivatives of `f` at `t`, orders 0 to 3."
function ab_ad_jet(f, t::Float64)
    d1(s) = ForwardDiff.derivative(f, s)
    d2(s) = ForwardDiff.derivative(d1, s)
    d3(s) = ForwardDiff.derivative(d2, s)
    return (f(t), d1(t), d2(t), d3(t))
end

ab_alloc_ezr(ctx, et, ref, ab) = @allocated spkezr(ctx, 301, et, ref, ab, 10)
ab_alloc_pos(ctx, et, ref, ab) = @allocated spkpos(ctx, 301, et, ref, ab, 10)
ab_alloc_der(ctx, et, ref, ab, ::Val{N}) where {N} =
    @allocated spkder(ctx, 301, et, ref, ab, 10, Val(N))
ab_alloc_parse(s) = @allocated Aberration(s)
ab_alloc_ssb(ctx, et) = @allocated spkssb(ctx, 399, et, 1)
ab_alloc_ltc(ctx, et, ab, st) = @allocated spkltc(ctx, 301, et, 1, ab, st)
ab_alloc_aps(ctx, et, ab, st, acc) = @allocated spkaps(ctx, 301, et, 1, ab, st, acc)
ab_alloc_acs(ctx, et, ab) = @allocated spkacs(ctx, 301, et, 1, ab, 399)
ab_alloc_apo(ctx, et, ab, st) = @allocated spkapo(ctx, 301, et, 1, st, ab)
ab_alloc_ltime(ctx, et) = @allocated ltime(ctx, et, 399, "->", 5)
ab_alloc_stelab(p, v) = @allocated stelab(p, v)

"Allocation of the second call of `f(args...)` (the first call compiles)."
function ab_warm(f, args...)
    f(args...)
    return f(args...)
end

group("9. aberration corrections", Kernels.FULL_SET_FILES...) do
    rng = MersenneTwister(SEED + 9)
    set = Kernels.full_set()
    oracle_load(set)
    ctx = neo_context(set)
    # The light time moves the target epoch by up to 5 hours (Neptune), and the stellar
    # aberration needs the observer at ± 1 s: stay one day inside the coverage.
    span = (DE440S_COVERAGE[1] + 86400.0, DE440S_COVERAGE[2] - 86400.0)

    @apart @testset "correction strings" begin
        for (s, uselt, conv, stel, xmit) in (
                ("NONE", false, false, false, false),
                ("LT", true, false, false, false), ("LT+S", true, false, true, false),
                ("CN", true, true, false, false), ("CN+S", true, true, true, false),
                ("XLT", true, false, false, true), ("XLT+S", true, false, true, true),
                ("XCN", true, true, false, true), ("XCN+S", true, true, true, true),
            )
            c = Aberration(s)
            @test (c.uselt, c.converged, c.stellar, c.transmit) == (uselt, conv, stel, xmit)
            @test c.geometric == (s == "NONE")
            @test string(c) == s
            @test isbits(c)
            # Case and blanks are ignored, as `ljucrs` of CSPICE does.
            spaced = " " * join(collect(lowercase(s)), " ") * "  "
            @test Aberration(lowercase(s)) == c
            @test Aberration(spaced) == c
            @test Aberration(SubString(s * "!", 1, length(s))) == c
            @test Aberration(c) === c
            @test Aberration(; uselt, converged = conv, stellar = stel, transmit = xmit) == c
            @test ab_warm(ab_alloc_parse, s) == 0
            @test ab_warm(ab_alloc_parse, spaced) == 0
        end
        @test (@inferred Aberration("cn + s")) == Aberration("CN+S")
        # `propertynames` lists the documented properties only
        @test propertynames(Aberration("LT")) == (:uselt, :converged, :stellar, :transmit, :geometric)
        # Values that CSPICE recognizes and rejects, and values that it does not know.
        for bad in (
                "RL", "RL+S", "XRL", "XRL+S", "S", "XS", "", "   ", "LT+", "LTS", "L",
                "NONE+S", "XNONE", "LT+S+S", "CN+SS", "λT", "LT\t",
            )
            @test throws_kind(() -> Aberration(bad), :invalidoption)
            @test_throws Exception SPICE.spkezr("301", 0.0, "J2000", bad, "399")
        end
        # The calls parse the string with `Aberration`: one bad value through each of them.
        @test throws_kind(() -> spkezr(ctx, 301, 0.0, "J2000", "RL+S", 399), :invalidoption)
        @test throws_kind(() -> spkpos(ctx, 301, 0.0, "J2000", "RL+S", 399), :spkinvalidoption)
        @test throws_kind(() -> spkder(ctx, 301, 0.0, "J2000", "RL+S", 399, Val(1)), :invalidoption)
        # Regression: the string forms checked the correction first. CSPICE checks that no
        # string is empty, then the body names, then, in `spkezr`, the correction;
        # `spkpos` checks the frame and the SPK data of the observer before the correction,
        # and the correction before the data of the target.
        # The error of `spkezr_c` or `spkpos_c`; the exception of SPICE.jl holds the long
        # message only.
        function ckind(call, targ, ref, ab, obs)
            out, lt = zeros(6), Ref(0.0)
            if call === :spkezr
                ccall(
                    (:spkezr_c, SPICE.libcspice), Cvoid,
                    (Cstring, Cdouble, Cstring, Cstring, Cstring, Ptr{Cdouble}, Ref{Cdouble}),
                    targ, 0.0, ref, ab, obs, out, lt
                )
            else
                ccall(
                    (:spkpos_c, SPICE.libcspice), Cvoid,
                    (Cstring, Cdouble, Cstring, Cstring, Cstring, Ptr{Cdouble}, Ref{Cdouble}),
                    targ, 0.0, ref, ab, obs, out, lt
                )
            end
            return cspice_kind()
        end
        for (call, targ, ref, ab, obs, kind) in (
                (:spkezr, "", "J2000", "ABC", "EARTH", :emptystring),
                (:spkezr, "NOBODY", "J2000", "ABC", "EARTH", :idcodenotfound),
                (:spkpos, "MOON", "J2000", "", "EARTH", :emptystring),
                (:spkpos, "MOON", "NO_SUCH_FRAME", "ABC", "EARTH", :unknownframe),
                (:spkpos, "MOON", "J2000", "ABC", "-999", :spkinsuffdata),
                (:spkpos, "-999", "J2000", "ABC", "EARTH", :spkinvalidoption),
            )
            f = getfield(Almanac, call)
            @test throws_kind(() -> f(ctx, targ, 0.0, ref, ab, obs), kind)
            @test ckind(call, targ, ref, ab, obs) == kind
        end
        # CSPICE keeps the first five non-blank characters only, so it reads "XLT+SX"
        # as "XLT+S". Almanac does not truncate: the value is an error.
        @test SPICE.spkezr("301", 0.0, "J2000", "XLT+SX", "399")[2] > 0
        @test throws_kind(() -> Aberration("XLT+SX"), :invalidoption)
    end

    pairs = [
        (301, 399), (399, 301), (10, 399), (399, 10), (5, 301), (4, 0), (3, 199),
        (0, 301), (301, 10), (299, 5), (8, 399), (10, 0),
    ]
    frames = Tuple{String, Any}[
        ("J2000", AB_TOL_INERTIAL), ("ECLIPJ2000", AB_TOL_INERTIAL),
        ("IAU_EARTH", AB_TOL_BODYFIXED), ("IAU_MOON", AB_TOL_BODYFIXED),
        ("MOON_PA_DE440", AB_TOL_BODYFIXED), ("MOON_ME", AB_TOL_BODYFIXED),
    ]

    @apart @testset "states and light time against CSPICE" begin
        epochs = random_epochs(rng, 40, span...)
        measured(
            @sprintf(
                "  %-14s %-6s %-9s %-9s %-9s %-9s %-9s", "frame", "abcorr",
                "ezr pos", "ezr vel", "ezr lt", "pos pos", "pos lt"
            )
        )
        for (ref, tol) in frames, ab in ("NONE", AB_CORRECTIONS...)
            d = ab_compare(ctx, pairs, epochs, ref, ab)
            measured(
                @sprintf(
                    "  %-14s %-6s %.2e  %.2e  %.2e  %.2e  %.2e", ref, ab, d.pos,
                    d.vel, d.lt, d.ppos, d.plt
                )
            )
            @test d.pos <= tol.pos
            @test d.vel <= tol.vel
            @test d.lt <= tol.lt
            @test d.ppos <= tol.pos
            @test d.plt <= tol.lt
        end
    end

    # Regression: the light time is `vnorm_` of the position divided by `clight_`
    # (`spkgeo_`, `spkltc_`, `spkapo_`), the norm scaled by the largest component; the
    # plain norm gave another last bit in about one call in three. The states of
    # de440s.bsp in these inertial frames need no sine or cosine, so the light times are
    # equal to the last bit. `zzstelab_` (f2c routine, no sine or cosine) takes its
    # norms, projections and unit vectors scaled in the same way: equal to the last bit
    # on random input. Own random stream.
    @apart @testset "light time and zzstelab_ to the last bit" begin
        r2 = MersenneTwister(SEED + 901)
        fails = Failures()
        for et in random_epochs(r2, 10, span...), ref in ("J2000", "ECLIPJ2000"),
                ab in ("NONE", "LT", "CN+S", "XLT+S"), (targ, obs) in ((301, 399), (10, 301), (5, 399))
            l1 = spkezr(ctx, targ, et, ref, ab, obs)[2]
            l2 = SPICE.spkezr(string(targ), et, ref, ab, string(obs))[2]
            check!(fails, l1 == l2, "spkezr light time", et, ref, ab, targ, obs)
            l1 = spkpos(ctx, targ, et, ref, ab, obs)[2]
            l2 = SPICE.spkpos(string(targ), et, ref, ab, string(obs))[2]
            check!(fails, l1 == l2, "spkpos light time", et, ref, ab, targ, obs)
        end
        for k in 1:500
            r = SVector{3}(randn(r2, 3) .* 10.0^rand(r2, 3:9))
            v = SVector{3}(randn(r2, 3) .* 10.0^rand(r2, -3:2))
            # an observer speed of 1e-4 km/s takes the numerical derivative of `zzstelab_`
            vo = SVector{3}(randn(r2, 3) .* (iseven(k) ? 30.0 : 1.0e-4))
            acc = SVector{3}(randn(r2, 3) .* 1.0e-5)
            xmit = isodd(k ÷ 2)
            sc, dsc = Almanac._zzstelab(xmit, acc, vo, r, v)
            a = zeros(3); b = zeros(3)
            ccall(
                (:zzstelab_, SPICE.libcspice), Cint,
                (Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
                xmit, collect(acc), collect(vo), vcat(r, v), a, b
            )
            check!(fails, sc == a && dsc == b, "zzstelab_", r, v, vo, acc, xmit)
        end
        @test isempty(fails)
    end

    @apart @testset "forms of the call" begin
        # The forms do not depend on the epoch: one epoch in the default run.
        fails = Failures()
        for et in random_epochs(rng, nsample(1, 5), span...), ab in AB_CORRECTIONS, ref in ("J2000", "IAU_EARTH")
            byid = spkezr(ctx, 301, et, ref, ab, 10)
            check!(fails, byid[1] isa SVector{6, Float64} && byid[2] isa Float64, "spkezr types", et, ab, ref)
            check!(fails, spkezr(ctx, "MOON", et, ref, ab, "SUN") == byid, "spkezr by name", et, ab, ref)
            check!(fails, spkezr(ctx, 301, et, ref, Aberration(ab), 10) == byid, "spkezr Aberration", et, ab, ref)
            check!(fails, spkezr(ctx, 301, et, ref, " " * lowercase(ab), 10) == byid, "spkezr lower case", et, ab, ref)
            check!(fails, spkezr(Almanac.kernelset(ctx), 301, et, ref, ab, 10) == byid, "spkezr KernelSet", et, ab, ref)
            pos, lt = spkpos(ctx, 301, et, ref, ab, 10)
            check!(fails, pos isa SVector{3, Float64}, "spkpos type", et, ab, ref)
            check!(fails, spkpos(ctx, "MOON", et, ref, Aberration(ab), "SUN") == (pos, lt), "spkpos by name", et, ab, ref)
            # spkpos follows `spkapo` and spkezr follows `spkltc`: rounding apart.
            check!(fails, isapprox(pos, byid[1][SOneTo(3)]; rtol = 1.0e-12), "spkpos ≈ spkezr", et, ab, ref)
            check!(fails, isapprox(lt, byid[2]; rtol = 1.0e-12), "spkpos lt ≈ spkezr lt", et, ab, ref)
            # The light time is the norm of the position without stellar aberration,
            # and stellar aberration keeps the norm.
            check!(fails, isapprox(lt, norm(pos) / 299792.458; rtol = 1.0e-13), "lt = |r| / c", et, ab, ref)
        end
        # A body relative to itself: zero state and zero light time, as CSPICE.
        for ab in AB_CORRECTIONS, ref in ("J2000", "IAU_MOON")
            state, lt = spkezr(ctx, 399, 1.0e8, ref, ab, 399)
            check!(fails, all(iszero, state) && lt == 0, "spkezr self", ab, ref)
            check!(fails, all(iszero, SPICE.spkezr("399", 1.0e8, ref, ab, "399")[1]), "CSPICE spkezr self", ab, ref)
            check!(fails, all(iszero, spkpos(ctx, 399, 1.0e8, ref, ab, 399)[1]), "spkpos self", ab, ref)
            check!(fails, all(all(iszero, d) for d in spkder(ctx, 399, 1.0e8, ref, ab, 399, Val(2))), "spkder self", ab, ref)
        end
        @test isempty(fails)
        # Reception and transmission differ, and the correction has the expected size:
        # the Sun seen from the Earth moves by about 20 arc seconds of aberration.
        geo = spkpos(ctx, 10, 0.0, "J2000", "NONE", 399)[1]
        app = spkpos(ctx, 10, 0.0, "J2000", "LT+S", 399)[1]
        xmt = spkpos(ctx, 10, 0.0, "J2000", "XLT+S", 399)[1]
        @test 15 < acos(clamp(dot(geo, app) / (norm(geo) * norm(app)), -1, 1)) * 206264.8 < 25
        @test norm(app - xmt) > 1.0e3
    end

    @apart @testset "default context" begin
        kclear()
        foreach(furnsh, set)
        for ab in ("LT+S", "XCN"), ref in ("J2000", "IAU_MOON")
            @test spkezr("MOON", 1.0e7, ref, ab, "EARTH") == spkezr(ctx, "MOON", 1.0e7, ref, ab, "EARTH")
            @test spkpos(301, 1.0e7, ref, Aberration(ab), 399) == spkpos(ctx, 301, 1.0e7, ref, ab, 399)
            @test spkder(301, 1.0e7, ref, ab, 399, Val(2)) == spkder(ctx, 301, 1.0e7, ref, ab, 399, Val(2))
        end
        @test spkssb(399, 1.0e7, "J2000") == spkssb(ctx, 399, 1.0e7, "J2000")
        @test spkacs(301, 1.0e7, "J2000", "CN+S", 399) == spkacs(ctx, 301, 1.0e7, "J2000", "CN+S", 399)
        @test ltime(1.0e7, 399, "->", 5) == ltime(ctx, 1.0e7, 399, "->", 5)
        kclear()
    end

    @apart @testset "low-level calls against CSPICE" begin
        # Each epoch runs every correction through every low-level call for two frames and
        # four pairs of bodies (about 1100 checks, which `fails` collects): three epochs in
        # the default run.
        epochs = random_epochs(rng, nsample(3, 15), span...)
        worst = Dict{String, Float64}()
        fails = Failures()
        function note(name, x, tol, inputs...)
            worst[name] = max(get(worst, name, 0.0), x)
            return within!(fails, x, tol, name, inputs...)
        end
        tol = 1.0e-11
        for et in epochs, ref in ("J2000", "ECLIPJ2000"), (targ, obs) in ((301, 399), (5, 301), (399, 10), (10, 0))
            cssb = SPICE.spkssb(obs, et, ref)
            ssb = spkssb(ctx, obs, et, ref)
            check!(fails, ssb isa SVector{6, Float64}, "spkssb type", et, ref, obs)
            note("spkssb", ab_reldiff(ssb, cssb), 1.0e-13, et, ref, obs)
            check!(fails, ssb == spkgeo(ctx, obs, et, ref, 0)[1], "spkssb == spkgeo", et, ref, obs)
            acc = SVector(1.0e-6, -2.0e-6, 5.0e-7)
            for ab in ("NONE", AB_CORRECTIONS...)
                cs, cl, cd = SPICE.spkltc(targ, et, ref, ab, cssb)
                s, l, d = spkltc(ctx, targ, et, ref, ab, cssb)
                note("spkltc state", ab_reldiff(s[SOneTo(3)], cs[1:3]) + ab_reldiff(s[SVector(4, 5, 6)], cs[4:6]), tol, et, ref, targ, obs, ab)
                note("spkltc lt", abs(l - cl) / cl, tol, et, ref, targ, obs, ab)
                note("spkltc dlt", abs(d - cd) / max(abs(cd), 1.0e-9), 1.0e-8, et, ref, targ, obs, ab)
                check!(fails, spkltc(ctx, targ, et, ref, Aberration(ab), ssb) == spkltc(ctx, targ, et, ref, ab, ssb), "spkltc Aberration form", et, ref, targ, obs, ab)

                cs, cl, cd = SPICE.spkaps(targ, et, ref, ab, cssb, vcat(acc, zeros(3)))   # SPICE.jl wants 6
                s, l, d = spkaps(ctx, targ, et, ref, ab, cssb, acc)
                note("spkaps state", ab_reldiff(s[SOneTo(3)], cs[1:3]) + ab_reldiff(s[SVector(4, 5, 6)], cs[4:6]), tol, et, ref, targ, obs, ab)
                note("spkaps lt", abs(l - cl) / cl, tol, et, ref, targ, obs, ab)
                note("spkaps dlt", abs(d - cd) / max(abs(cd), 1.0e-9), 1.0e-8, et, ref, targ, obs, ab)

                cs, cl, cd = SPICE.spkacs(targ, et, ref, ab, obs)
                s, l, d = spkacs(ctx, targ, et, ref, ab, obs)
                note("spkacs state", ab_reldiff(s[SOneTo(3)], cs[1:3]) + ab_reldiff(s[SVector(4, 5, 6)], cs[4:6]), tol, et, ref, targ, obs, ab)
                note("spkacs lt", abs(l - cl) / cl, tol, et, ref, targ, obs, ab)
                note("spkacs dlt", abs(d - cd) / max(abs(cd), 1.0e-9), 1.0e-8, et, ref, targ, obs, ab)
                # spkacs in an inertial frame is the state of spkezr. (With NONE, spkacs
                # takes the difference of two barycentric states and spkezr does not.)
                ab == "NONE" || check!(fails, (s, l) == spkezr(ctx, targ, et, ref, ab, obs), "spkacs == spkezr", et, ref, targ, obs, ab)

                cp, cpl = SPICE.spkapo(targ, et, ref, cssb, ab)
                p, pl = spkapo(ctx, targ, et, ref, cssb, ab)
                note("spkapo position", ab_reldiff(p, cp), tol, et, ref, targ, obs, ab)
                note("spkapo lt", abs(pl - cpl) / cpl, tol, et, ref, targ, obs, ab)
                if ab != "NONE"
                    q, ql = spkpos(ctx, targ, et, ref, ab, obs)
                    check!(fails, isapprox(p, q; rtol = 1.0e-12) && isapprox(pl, ql; rtol = 1.0e-12), "spkapo ≈ spkpos", et, ref, targ, obs, ab)
                end
            end
            for dir in ("->", "<-")
                cet, cel = SPICE.ltime(et, obs, dir, targ)
                tt, el = ltime(ctx, et, obs, dir, targ)
                note("ltime elapsed", abs(el - cel) / cel, tol, et, ref, targ, obs, dir)
                check!(fails, abs(tt - cet) <= 4 * eps(abs(et)) + tol * cel, "ltime target epoch", et, ref, targ, obs, dir)
            end
            # Stellar aberration of a position, with the true observer velocity and with
            # a velocity near the limit.
            p0 = spkpos(ctx, targ, et, ref, "LT", obs)[1]
            for v in (ssb[SVector(4, 5, 6)], SVector(1.0e5, -2.0e5, 3.0e4), SVector(0.0, 0.0, 0.0))
                capp = SPICE.stelab(collect(p0), collect(v))
                app = stelab(p0, v)
                check!(fails, app isa SVector{3, Float64}, "stelab type", et, ref, targ, obs, v)
                note("stelab", ab_reldiff(app, capp), 1.0e-14, et, ref, targ, obs, v)
                check!(fails, isapprox(norm(app), norm(p0); rtol = 1.0e-14), "stelab keeps the length", et, ref, targ, obs, v)
                note("stlabx", ab_reldiff(stlabx(p0, v), SPICE.stelab(collect(p0), collect(-v))), 1.0e-14, et, ref, targ, obs, v)
            end
        end
        @test isempty(fails)
        for name in sort!(collect(keys(worst)))
            measured(@sprintf("  %-58s max=%.3e", name * ", relative", worst[name]))
        end
        # Errors of the low-level calls, with the short messages of CSPICE N0067 (checked
        # with SPICE.jl): a frame that is not inertial is BADFRAME, and a bad correction
        # is reported before it.
        st = spkssb(ctx, 399, 0.0, "J2000")
        @test throws_kind(() -> spkltc(ctx, 301, 0.0, "IAU_EARTH", "LT", st), :badframe)
        @test throws_kind(() -> spkaps(ctx, 301, 0.0, "IAU_EARTH", "LT", st, zeros(3)), :badframe)
        @test throws_kind(() -> spkacs(ctx, 301, 0.0, "IAU_EARTH", "LT", 399), :badframe)
        @test throws_kind(() -> spkapo(ctx, 301, 0.0, "IAU_EARTH", st, "LT"), :badframe)
        @test throws_kind(() -> spkltc(ctx, 301, 0.0, "IAU_EARTH", "XYZ", st), :invalidoption)
        @test throws_kind(() -> spkaps(ctx, 301, 0.0, "IAU_EARTH", "S", st, zeros(3)), :invalidoption)
        @test throws_kind(() -> spkapo(ctx, 301, 0.0, "IAU_EARTH", st, "XYZ"), :spkinvalidoption)
        @test throws_kind(() -> spkltc(ctx, 301, 0.0, "J2000", "LT", zeros(5)), :badarg)
        @test throws_kind(() -> spkaps(ctx, 301, 0.0, "J2000", "LT", st, zeros(2)), :badarg)
        @test throws_kind(() -> spkltc(ctx, 301, 0.0, "J2000", "S", st), :invalidoption)
        @test throws_kind(() -> ltime(ctx, 0.0, 399, "<>", 301), :baddirection)
        @test throws_kind(() -> stelab(SVector(1.0, 0.0, 0.0), SVector(3.0e5, 0.0, 0.0)), :valueoutofrange)
        @test throws_kind(() -> stlabx(SVector(1.0, 0.0, 0.0), SVector(3.0e5, 0.0, 0.0)), :valueoutofrange)
        @test throws_kind(() -> stelab([1.0, 0.0], [0.0, 0.0, 0.0]), :badarg)
        # `ltime` reads the first two characters of the direction (`ltime.c`), and an
        # empty direction is EMPTYSTRING (`ltime_c`)
        @test ltime(ctx, 0.0, 399, "-> ", 301) == ltime(ctx, 0.0, 399, "->", 301)
        @test SPICE.ltime(0.0, 399, "-> ", 301)[2] == SPICE.ltime(0.0, 399, "->", 301)[2]
        @test ltime(ctx, 0.0, 399, "<-x", 301) == ltime(ctx, 0.0, 399, "<-", 301)
        @test throws_kind(() -> ltime(ctx, 0.0, 399, "-", 301), :baddirection)
        @test throws_kind(() -> ltime(ctx, 0.0, 399, "", 301), :emptystring)
        ccall(
            (:ltime_c, SPICE.libcspice), Cvoid,
            (Cdouble, Cint, Cstring, Cint, Ref{Cdouble}, Ref{Cdouble}),
            0.0, 399, "", 301, Ref(0.0), Ref(0.0)
        )
        @test cspice_kind() === :emptystring
        # `zzstelab` signals DIVIDEBYZERO when the part of the observer velocity normal
        # to the line of sight is not below the speed of light
        fast = [st[1], st[2], st[3], 0.0, 0.0, 2 * clight()]
        @test throws_kind(() -> spkaps(ctx, 301, 0.0, "J2000", "LT+S", fast, zeros(3)), :dividebyzero)
        ccall(
            (:spkaps_c, SPICE.libcspice), Cvoid,
            (
                Cint, Cdouble, Cstring, Cstring, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                Ref{Cdouble}, Ref{Cdouble},
            ),
            301, 0.0, "J2000", "LT+S", fast, zeros(3), zeros(6), Ref(0.0), Ref(0.0)
        )
        @test cspice_kind() === :dividebyzero
        # `spkltc` signals DIVIDEBYZERO for a range rate near the speed of light: a
        # target in memory that comes towards the barycenter at (1 - 1e-11) c, which
        # CSPICE cannot load (the kind is the one of `spkltc.c`)
        fctx = SpiceContext()
        vfast = (1 - 1.0e-11) * clight()
        load_function!(
            fctx, "fast", t -> [1.0e8 - vfast * t, 0.0, 0.0, -vfast, 0.0, 0.0], -1000, 0,
            "J2000", -1000.0, 1000.0
        )
        @test throws_kind(() -> spkltc(fctx, -1000, 0.0, "J2000", "LT", zeros(6)), :dividebyzero)
        # extreme magnitudes: the unit vector is scaled as in `vhat` (audit of 2026-11-16)
        for mag in (1.0e-300, 1.0e-200, 1.0e200, 1.0e300), d in ((1.0, 0.0, 0.0), (0.6, 0.8, 0.0))
            p = mag .* collect(d)
            @test stelab(p, [0.0, 30.0, 0.0]) == SPICE.stelab(p, [0.0, 30.0, 0.0])
        end
        # An epoch whose light-time corrected epoch leaves the coverage.
        @test throws_kind(() -> spkezr(ctx, 8, DE440S_COVERAGE[1] - 86000.0, "J2000", "LT", 399), :spkinsuffdata)
    end

    cases = [
        (301, 399, "J2000", AB_FD_INERTIAL), (399, 10, "J2000", AB_FD_INERTIAL),
        (10, 301, "ECLIPJ2000", AB_FD_INERTIAL), (5, 399, "J2000", AB_FD_INERTIAL),
        (301, 399, "IAU_EARTH", AB_FD_FAST), (10, 399, "IAU_EARTH", AB_FD_FAST),
        (399, 301, "MOON_PA_DE440", AB_FD_INERTIAL), (10, 301, "IAU_MOON", AB_FD_INERTIAL),
        (4, 399, "IAU_MOON", AB_FD_INERTIAL),
    ]
    depochs = record_safe_epochs(rng, 6)

    @apart @testset "spkder: order 0 and the CSPICE velocity" begin
        measured("  spkder order 1 against the velocity of CSPICE spkezr, largest relative difference")
        fails = Failures()
        for ab in AB_CORRECTIONS
            w0 = 0.0
            w1i = 0.0
            w1b = 0.0
            for (targ, obs, ref, _) in cases, et in depochs[1:nsample(2, 6)]
                jet = spkder(ctx, targ, et, ref, ab, obs, Val(1))
                check!(fails, jet isa NTuple{2, SVector{3, Float64}}, "spkder type", ab, targ, obs, ref, et)
                state, _ = spkezr(ctx, targ, et, ref, ab, obs)
                cstate, _ = SPICE.spkezr(string(targ), et, ref, ab, string(obs))
                w0 = max(w0, ab_reldiff(jet[1], state[SOneTo(3)]))
                d = ab_reldiff(jet[2], cstate[4:6])
                ref in ("J2000", "ECLIPJ2000") ? (w1i = max(w1i, d)) : (w1b = max(w1b, d))
                # A lower order is a prefix of a higher order.
                full = spkder(ctx, targ, et, ref, ab, obs, Val(3))
                check!(fails, all(isapprox(jet[k], full[k]; rtol = 1.0e-12) for k in 1:2), "order 1 prefix of 3", ab, targ, obs, ref, et)
                check!(fails, isapprox(spkder(ctx, targ, et, ref, ab, obs, Val(0))[1], full[1]; rtol = 1.0e-13), "order 0 prefix of 3", ab, targ, obs, ref, et)
            end
            measured(
                @sprintf(
                    "    %-6s order 0 vs spkezr %.2e; order 1 vs CSPICE velocity: inertial %.2e, body-fixed %.2e",
                    ab, w0, w1i, w1b
                )
            )
            @test w0 <= 1.0e-12
            # `CN` forms: rounding. `LT` forms: CSPICE takes the light-time rate of the
            # converged equation for a position that is one iteration from the geometric
            # light time, which is a relative difference of order (v/c) * (range rate/c).
            @test w1i <= (occursin("CN", ab) ? 5.0e-12 : 1.0e-6)
            @test w1b <= (occursin("CN", ab) ? 5.0e-11 : 1.0e-5)
        end
        @test isempty(fails)
        # NONE through the corrected form is the geometric form.
        @test spkder(ctx, 301, depochs[1], "IAU_EARTH", "NONE", 399, Val(3)) ==
            spkder(ctx, 301, depochs[1], "IAU_EARTH", 399, Val(3))
        @test throws_kind(() -> spkder(ctx, 301, 0.0, "J2000", "LT", 399, Val(-1)), :badarg)
    end

    @apart @testset "spkder: finite differences" begin
        for (targ, obs, ref, fd) in cases, ab in AB_CORRECTIONS
            worst = zeros(3)
            wpos = 0.0
            for et in depochs
                jet(t) = spkder(ctx, targ, t, ref, ab, obs, Val(3))
                j0 = jet(et); jp = jet(et + fd.h); jm = jet(et - fd.h)
                worst .= max.(worst, [ab_reldiff((jp[k] - jm[k]) / (2 * fd.h), j0[k + 1]) for k in 1:3])
                # Order 1 against the central difference of the position of spkpos.
                pp = spkpos(ctx, targ, et + fd.h, ref, ab, obs)[1]
                pm = spkpos(ctx, targ, et - fd.h, ref, ab, obs)[1]
                wpos = max(wpos, ab_reldiff((pp - pm) / (2 * fd.h), j0[2]))
            end
            ab in ("LT", "XCN+S") && measured(
                @sprintf(
                    "  FD  %3d rel %3d %-13s %-6s orders 1..3: %s; order 1 vs spkpos %.2e (tol %.0e)",
                    targ, obs, ref, ab, join((@sprintf("%.2e", w) for w in worst), " "), wpos, fd.rtol
                )
            )
            @test all(<=(fd.rtol), worst)
            @test wpos <= fd.rtol
        end
    end

    @apart @testset "spkder: nested ForwardDiff of spkpos (dual-number epoch)" begin
        probe, plt = spkpos(ctx, 301, ForwardDiff.Dual(depochs[1], 1.0), "IAU_EARTH", "CN+S", 399)
        @test eltype(probe) <: ForwardDiff.Dual && plt isa ForwardDiff.Dual
        @test ForwardDiff.value.(probe) == spkpos(ctx, 301, depochs[1], "IAU_EARTH", "CN+S", 399)[1]
        for ab in AB_CORRECTIONS
            worst = zeros(4)
            for (targ, obs, ref, _) in cases, et in depochs[1:3]
                pos(t) = spkpos(ctx, targ, t, ref, ab, obs)[1]
                ad = ab_ad_jet(pos, et)
                jet = spkder(ctx, targ, et, ref, ab, obs, Val(3))
                # Order 0 is the position of spkezr, not the one of spkpos. For the
                # converged forms the two differ in CSPICE as well: spkltc (spkezr)
                # stops when the light time changes by less than 1e-17 |epoch|, which
                # is 4e-9 s at 4.5e8 s and can be after one iteration, and spkapo
                # (spkpos) always makes three. The difference is then of the order of
                # the square of speed over light speed; a C reproducer measured up to
                # 2.0e-11 relative between `spkezr_c` and `spkpos_c`.
                ezr = spkezr(ctx, targ, et, ref, ab, obs)[1][SOneTo(3)]
                worst[1] = max(worst[1], ab_reldiff(jet[1], ezr))
                worst[2:4] .= max.(worst[2:4], [ab_reldiff(jet[k], ad[k]) for k in 2:4])
            end
            measured(
                @sprintf(
                    "  AD  %-6s orders 0..3: %s", ab,
                    join((@sprintf("%.2e", w) for w in worst), " ")
                )
            )
            @test all(worst[k] <= AB_AD_RTOL[k] for k in 1:4)
        end
        # A dual-number epoch through spkezr and through spkder itself.
        et = depochs[2]
        for ab in ("LT+S", "XCN"), ref in ("J2000", "IAU_MOON")
            jet = spkder(ctx, 10, et, ref, ab, 301, Val(3))
            dstate = spkezr(ctx, 10, ForwardDiff.Dual(et, 1.0), ref, ab, 301)[1]
            @test isapprox(ForwardDiff.value.(dstate), spkezr(ctx, 10, et, ref, ab, 301)[1]; rtol = 1.0e-14)
            @test isapprox(ForwardDiff.partials.(dstate, 1)[SOneTo(3)], jet[2]; rtol = 1.0e-6)
            d3 = ForwardDiff.derivative(t -> spkder(ctx, 10, t, ref, ab, 301, Val(2))[3], et)
            @test isapprox(d3, jet[4]; rtol = 1.0e-10)
        end
    end

    @apart @testset "zero allocation and inferred types" begin
        et = 1.2345e8
        refs = (
            1, 17, Int(namfrm(ctx, "IAU_EARTH")), Int(namfrm(ctx, "MOON_PA_DE440")),
            Int(namfrm(ctx, "MOON_ME")),
        )
        fails = Failures()
        for ref in refs, ab in ("NONE", AB_CORRECTIONS...)
            check!(fails, ab_warm(ab_alloc_ezr, ctx, et, ref, ab) == 0, "spkezr", ref, ab)
            check!(fails, ab_warm(ab_alloc_pos, ctx, et, ref, ab) == 0, "spkpos", ref, ab)
            check!(fails, ab_warm(ab_alloc_ezr, ctx, et, ref, Aberration(ab)) == 0, "spkezr Aberration", ref, ab)
            for N in 0:3
                check!(fails, ab_warm(ab_alloc_der, ctx, et, ref, ab, Val(N)) == 0, "spkder", ref, ab, N)
            end
        end
        @test isempty(fails)
        st = spkssb(ctx, 399, et, 1)
        acc = SVector(1.0e-6, 2.0e-6, 3.0e-6)
        @test ab_warm(ab_alloc_ssb, ctx, et) == 0
        @test ab_warm(ab_alloc_ltime, ctx, et) == 0
        @test ab_warm(ab_alloc_stelab, st[SOneTo(3)], st[SVector(4, 5, 6)]) == 0
        for ab in ("LT", "CN+S", "XLT+S")
            @test ab_warm(ab_alloc_ltc, ctx, et, ab, st) == 0
            @test ab_warm(ab_alloc_aps, ctx, et, ab, st, acc) == 0
            @test ab_warm(ab_alloc_acs, ctx, et, ab) == 0
            @test ab_warm(ab_alloc_apo, ctx, et, ab, st) == 0
        end
        ks = Almanac.kernelset(ctx)
        for ref in (1, refs[3]), ab in ("LT", "XCN+S")
            @test (@inferred spkezr(ks, 301, et, ref, ab, 399)) isa Tuple{SVector{6, Float64}, Float64}
            @test (@inferred spkpos(ks, 301, et, ref, ab, 399)) isa Tuple{SVector{3, Float64}, Float64}
            @test (@inferred spkder(ks, 301, et, ref, ab, 399, Val(3))) isa NTuple{4, SVector{3, Float64}}
            @test (@inferred spkezr(ctx, 301, et, ref, Aberration(ab), 399)) isa Tuple{SVector{6, Float64}, Float64}
        end
        @test (@inferred spkssb(ks, 399, et, 1)) isa SVector{6, Float64}
        @test (@inferred spkltc(ks, 301, et, 1, "CN", st)) isa Tuple{SVector{6, Float64}, Float64, Float64}
        @test (@inferred spkaps(ks, 301, et, 1, "CN+S", st, acc)) isa Tuple{SVector{6, Float64}, Float64, Float64}
        @test (@inferred spkacs(ks, 301, et, 1, "CN+S", 399)) isa Tuple{SVector{6, Float64}, Float64, Float64}
        @test (@inferred spkapo(ks, 301, et, 1, st, "CN+S")) isa Tuple{SVector{3, Float64}, Float64}
        @test (@inferred ltime(ks, et, 399, "<-", 301)) isa Tuple{Float64, Float64}
        @test (@inferred stelab(st[SOneTo(3)], st[SVector(4, 5, 6)])) isa SVector{3, Float64}
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
