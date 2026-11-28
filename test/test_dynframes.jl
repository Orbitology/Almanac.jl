# Group 17: dynamic frames (class 5) and switch frames (class 6) against CSPICE.
#
# The frame kernels are written by `test/dynfixtures.jl`: every frame family, every kind
# of defining vector, frozen frames, frames with an inertial rotation state, dynamic
# frames that are defined relative to dynamic frames, and switch frames.
#
# `pxform` is compared element by element. The derivative block of `sxform` and the
# velocity of `spkezr` are compared relative to their size. CSPICE takes the derivative
# of a velocity vector by a central difference over about one second, so its derivative
# block has a rounding noise of about 1e-16 divided by the rotation rate of the frame;
# Almanac makes the same difference, and the two results agree to that noise only. The
# thresholds below are therefore per frame, each about ten times the measured value.
#
# Sizes. Each order of `pxder` above zero is a new number type, and the first evaluation
# of a dynamic frame with a new number type compiles the whole evaluation again (about
# 4 s for one order; about 9 s more for `DT_TWOVEC_IN_DYN`, whose vectors are given in
# dynamic frames, so the number types nest one level deeper). The default run therefore
# uses one order above zero, `DF_ORDER` = 2, and does not call `pxder` for
# `DT_TWOVEC_IN_DYN`; `pxform`, `sxform` and `spkezr` are compared for every frame in
# both runs. The exhaustive run (`ALMANAC_TEST_FULL=1`) uses the orders 1, 2 and 3 for
# every frame.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Almanac keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALMANAC_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

isdefined(@__MODULE__, :CKFix) || include("ckfixtures.jl")
isdefined(@__MODULE__, :DynFix) || include("dynfixtures.jl")

"""
    df_oracle_kind(call, from, to, et) -> Symbol

Short message of the error that CSPICE signals for `pxform_c` (`call = :pxform_c`) or
`sxform_c` from frame `from` to frame `to` at `et`, as a lower-case symbol without
`SPICE(...)`; `:none` without an error. The C functions are called directly, because
SPICE.jl keeps the long message only.
"""
function df_oracle_kind(call::Symbol, from::AbstractString, to::AbstractString, et::Float64)
    out = zeros(36)
    if call === :pxform_c
        ccall(
            (:pxform_c, SPICE.libcspice), Cvoid, (Cstring, Cstring, Cdouble, Ptr{Float64}),
            from, to, et, out
        )
    else
        ccall(
            (:sxform_c, SPICE.libcspice), Cvoid, (Cstring, Cstring, Cdouble, Ptr{Float64}),
            from, to, et, out
        )
    end
    ccall((:failed_c, SPICE.libcspice), Cint, ()) == 0 && return :none
    buf = zeros(UInt8, 64)
    ccall((:getmsg_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Ptr{UInt8}), "SHORT", 64, buf)
    ccall((:reset_c, SPICE.libcspice), Cvoid, ())
    short = unsafe_string(pointer(buf))
    m = match(r"SPICE\(([A-Z0-9_]+)\)", short)
    return m === nothing ? :unknown : Symbol(lowercase(m.captures[1]))
end

"Kind of the `SpiceError` that Almanac throws for `f()`; `:none` without an error."
function df_neo_kind(f)
    try
        f()
    catch e
        e isa Almanac.SpiceError || rethrow()
        return e.kind
    end
    return :none
end

"""
    df_rel(a, b, floor=0.0)

Largest element of `|a - b|`, divided by the largest element of `|b|` or by `floor` if
that is larger. If both arrays are below 1e-17 the result is zero: the derivative block
of a frame that does not rotate is the rounding noise of a cancellation in both engines.
"""
function df_rel(a, b, floor = 0.0)
    s = maximum(abs, b)
    (s <= 1.0e-17 && maximum(abs, a) <= 1.0e-17) && return 0.0
    return maximum(abs.(a .- b)) / max(s, floor, floatmin())
end

df_alloc_px(ks, from, to, et) = @allocated pxform(ks, from, to, et)
df_alloc_sx(ks, from, to, et) = @allocated sxform(ks, from, to, et)
df_alloc_der(ks, from, to, et, ::Val{N}) where {N} = @allocated pxder(ks, from, to, et, Val(N))
df_alloc_ezr(ks, et, ref, ab) = @allocated spkezr(ks, 301, et, ref, ab, 399)
df_alloc_pos(ks, et, ref, ab) = @allocated spkpos(ks, 301, et, ref, ab, 399)

"Order of `pxder` in the checks of the second derivative and of the allocation."
const DF_ORDER = 2
"Order of `pxder` where the first derivative is enough: 1 in the exhaustive run."
const DF_LOW = FULL ? 1 : DF_ORDER
"Order of `pxder` in the thread check: 3 in the exhaustive run."
const DF_THREAD = FULL ? 3 : DF_ORDER

"True if `pxder` is called for the frame `name` in this run (see \"Sizes\" above)."
df_has_der(name::AbstractString) = FULL || name != "DT_TWOVEC_IN_DYN"

"Five-point central difference of the matrix function `f` at `t` with the step `h`."
df_fd(f, t, h) = (8 .* (f(t + h) .- f(t - h)) .- (f(t + 2h) .- f(t - 2h))) ./ (12h)

"""
Thresholds per frame, `(pxform, sxform derivative block, spkezr position, spkezr
velocity)`: the default, and the frames whose derivative comes from a numerical
difference in CSPICE.
"""
const DF_DEFAULT_TOL = (5.0e-14, 2.0e-12, 5.0e-14, 2.0e-12)
const DF_TOL = Dict(
    "DT_GSE" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-10),
    "DT_GSE_LTS" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-10),
    "DT_GSE_CN" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-10),
    "DT_NADIR" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-10),
    "DT_NADIR_LT" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-10),
    "DT_VEL_PCK" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-7),
    "DT_CONST_ABCORR" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-10),
    "DT_EULER_ON_GSE" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-10),
    "DT_SWITCH_TIME" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-10),
    "DT_SWITCH_NEST" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-10),
    "DT_PRODUCT" => (5.0e-14, 5.0e-9, 5.0e-14, 5.0e-10),
    # A velocity that is given in a dynamic frame: CSPICE differentiates a quantity
    # that already holds a numerical derivative. Measured (2026-10-12): the residual of
    # a cubic fit of the `pxform` elements over 20 s is 1.1e-11 in both engines (7e-16
    # for `DT_GSE`), and the two engines differ by up to 8e-12, which is that noise.
    "DT_TWOVEC_IN_DYN" => (5.0e-10, 5.0e-4, 5.0e-10, 5.0e-4),
)

"""
    df_thread_check(ctx, rng)

Testset: many tasks get the same answers as one task, for `pxform`, `sxform`, `spkezr`
and `pxder` of order `DF_THREAD` on random frames and epochs. The exhaustive run calls
it after the other checks, as before. The default run calls it first, so that the order
`DF_THREAD` was not used before and the tasks race to make its function pointer.
"""
function df_thread_check(ctx, rng)
    if Threads.nthreads() == 1
        skip_notice("17. threads", "one thread only")
        return nothing
    end
    @testset "threads: the same answers from many tasks" begin
        names = first.(DynFix.FRAMES)
        jobs = [(rand(rng, names), (2 * rand(rng) - 1) * 2.9e8) for _ in 1:nsample(200, 400)]
        serial = [
            (
                pxform(ctx, "J2000", n, et), sxform(ctx, "J2000", n, et),
                spkezr(ctx, 301, et, n, "LT+S", 399)[1],
            ) for (n, et) in jobs
        ]
        results = Vector{Any}(undef, length(jobs))
        der = Vector{Any}(undef, length(jobs))
        Threads.@threads for i in eachindex(jobs)
            n, et = jobs[i]
            der[i] = df_has_der(n) ? pxder(ctx, "J2000", n, et, Val(DF_THREAD)) : nothing
            results[i] = (
                pxform(ctx, "J2000", n, et), sxform(ctx, "J2000", n, et),
                spkezr(ctx, 301, et, n, "LT+S", 399)[1],
            )
        end
        @test all(i -> results[i] == serial[i], eachindex(jobs))
        @test all(
            i -> !df_has_der(jobs[i][1]) ||
                der[i] == pxder(ctx, "J2000", jobs[i][1], jobs[i][2], Val(DF_THREAD)), eachindex(jobs)
        )
    end
    return nothing
end

group(
    "17. dynamic frames and switch frames against CSPICE", Kernels.LSK, Kernels.PCK11,
    Kernels.DE440S
) do
    rng = MersenneTwister(SEED + 17)
    set = [Kernels.LSK, Kernels.PCK11, Kernels.DE440S, DynFix.frame_kernel()]
    oracle_load(set)
    ctx = neo_context(set)
    ks = Almanac.kernelset(ctx)
    span = 5.5e8
    nepoch = nsample(30, 150)

    @apart @testset "frame definitions are read" begin
        fails = Failures()
        for (name, family, center) in DynFix.FRAMES
            id = namfrm(ctx, name)
            check!(fails, id == SPICE.namfrm(name) != 0, "namfrm", name)
            check!(fails, frinfo(ctx, id) == SPICE.frinfo(id), "frinfo", name)
            check!(fails, frinfo(ctx, id)[1] == center, "center", name)
            class = frinfo(ctx, id)[2]
            check!(fails, class in (5, 6), "class", name, class)
            check!(fails, haskey(class == 5 ? ks.dyn : ks.switch, Int32(id)), "in the snapshot", name)
        end
        @test isempty(fails)
        @test length(ks.dyn) + length(ks.switch) == length(DynFix.FRAMES) + length(DynFix.BASE_FRAMES)
        @test all(d -> d.errkind === :none, values(ks.dyn))
        @test all(d -> d.errkind === :none, values(ks.switch))
        # A snapshot without a dynamic frame has no inner snapshot.
        @test Almanac.kernelset(neo_context([Kernels.LSK, Kernels.PCK11])).inner === nothing
        @test ks.inner !== nothing && ks.inner.ks.inner !== nothing
        @test ks.inner.ks.inner.ks.inner === nothing
    end

    FULL || df_thread_check(ctx, rng)

    rows = String[]
    @apart @testset "pxform, sxform and spkezr for each frame" begin
        fails = Failures()
        for (name, family, center) in DynFix.FRAMES
            tol = get(DF_TOL, name, DF_DEFAULT_TOL)
            dp = 0.0; ds = 0.0; dr = 0.0; dv = 0.0; dra = 0.0; dva = 0.0
            for et in random_epochs(rng, nepoch, -span, span)
                a = pxform(ctx, "J2000", name, et)
                b = SPICE.pxform("J2000", name, et)
                dp = max(dp, maximum(abs.(a .- b)))
                a6 = sxform(ctx, "J2000", name, et)
                b6 = SPICE.sxform("J2000", name, et)
                dp = max(dp, maximum(abs.(a6[1:3, 1:3] .- b6[1:3, 1:3])))
                ds = max(ds, df_rel(a6[4:6, 1:3], b6[4:6, 1:3]))
                for (targ, obs) in (("MOON", "EARTH"), ("SUN", "MOON"))
                    s, lt = spkezr(ctx, targ, et, name, "NONE", obs)
                    t, lu = SPICE.spkezr(targ, et, name, "NONE", obs)
                    dr = max(dr, norm(s[1:3] .- t[1:3]) / norm(t[1:3]))
                    dv = max(dv, norm(s[4:6] .- t[4:6]) / norm(t[4:6]))
                    s, lt = spkezr(ctx, targ, et, name, "LT+S", obs)
                    t, lu = SPICE.spkezr(targ, et, name, "LT+S", obs)
                    dra = max(dra, norm(s[1:3] .- t[1:3]) / norm(t[1:3]))
                    dva = max(dva, norm(s[4:6] .- t[4:6]) / norm(t[4:6]))
                    check!(fails, abs(lt - lu) <= 1.0e-12 * lu, "LT+S light time", name, targ, obs, et)
                end
            end
            push!(
                rows, @sprintf(
                    "  %-17s px=%.1e sx=%.1e ezr NONE r=%.1e v=%.1e  LT+S r=%.1e v=%.1e  %s",
                    name, dp, ds, dr, dv, dra, dva, family
                )
            )
            ok = dp <= tol[1] && ds <= tol[2] && max(dr, dra) <= tol[3] && max(dv, dva) <= tol[4]
            ok || @info "above the threshold" name tol dp ds dr dv dra dva
            @test dp <= tol[1]
            @test ds <= tol[2]
            @test max(dr, dra) <= tol[3]
            @test max(dv, dva) <= tol[4]
        end
        @test isempty(fails)
        measured(
            "  largest difference from CSPICE, $nepoch epochs per frame (px: element of the " *
                "matrix; sx: derivative block, relative; ezr: relative)"
        )
        foreach(measured, rows)
    end

    # Regression: `zzdynfrm_` builds a two-vector frame with `vhat_`, `dvhat_` and
    # `ducrss_`, `frmchg_` joins the links with `zzmsxf_`, and `spkgeo_` applies the
    # state transformation with `mxvg_`. With the same arithmetic the results are those
    # of CSPICE to the last bit for these frames, whose vectors come from the Chebyshev
    # segments of de440s.bsp and need no sine or cosine (which differ in the last bit
    # between Julia and the C library). The old arithmetic differed by about 1e-16
    # relative in the matrix and by up to 5e-9 in the derivative block. Own random
    # stream, so that the cases after this one keep their data.
    @apart @testset "two-vector frames and spkezr to the last bit" begin
        r2 = MersenneTwister(SEED + 1701)
        fails = Failures()
        for name in ("DT_GSE", "DT_GSE_LTS", "DT_GSE_INERTIAL"), et in random_epochs(r2, 10, -span, span)
            check!(fails, pxform(ctx, name, "J2000", et) == SPICE.pxform(name, "J2000", et), "pxform", name, et)
            check!(fails, sxform(ctx, "J2000", name, et) == SPICE.sxform("J2000", name, et), "sxform", name, et)
            for (targ, obs) in (("MOON", "EARTH"), ("SUN", "MOON"))
                s = spkezr(ctx, targ, et, name, "NONE", obs)[1]
                check!(fails, s == SPICE.spkezr(targ, et, name, "NONE", obs)[1], "spkezr", name, targ, et)
            end
        end
        @test isempty(fails)
    end

    # An of-date frame and a frame with the inertial rotation state whose base is a PCK
    # frame: `zzdynfrm_` builds both relative to J2000 and goes to the base and back
    # with the state transformation of the base. In the derivative block the rotation
    # rate of the Earth (7.3e-5 rad/s) then cancels, in both engines, down to the rate
    # of the precession (8e-12 rad/s) or to zero, so the difference is measured relative
    # to that rate. Measured on 600 epochs for each frame and start frame: matrix
    # elements within 1.4e-15 of CSPICE, derivative block within 5.6e-16 of the rate
    # (4.1e-20 rad/s), from the last bit of a sine or cosine. Own random stream.
    @apart @testset "of-date and inertial-state frames on a PCK base frame" begin
        r4 = MersenneTwister(SEED + 1703)
        rate = 7.292115e-5
        dp = 0.0; ds = 0.0
        for (name, family, center) in DynFix.BASE_FRAMES, from in ("J2000", "IAU_MOON")
            for et in random_epochs(r4, nepoch, -span, span)
                dp = max(dp, maximum(abs.(pxform(ctx, from, name, et) .- SPICE.pxform(from, name, et))))
                a6 = sxform(ctx, from, name, et)
                b6 = SPICE.sxform(from, name, et)
                dp = max(dp, maximum(abs.(a6[1:3, 1:3] .- b6[1:3, 1:3])))
                ds = max(ds, df_rel(a6[4:6, 1:3], b6[4:6, 1:3], rate))
            end
        end
        measured(@sprintf("  frames on a PCK base: px=%.1e sx=%.1e of the rotation rate of the base", dp, ds))
        @test dp <= DF_DEFAULT_TOL[1]
        @test ds <= 5.0e-15
    end

    # Regression: the state transformations of the of-date frames come, as in
    # `zzdynfrm_`, from the angles and their analytic rates of `zzeprc76_`, `zzenut80_`
    # and `zzmobliq_` through `eul2xf_`; they were the derivative of a `Taylor`
    # evaluation, which differed in the last bits of every derivative block. The f2c
    # routines are called directly. `zzmobliq_` has no sine or cosine and is compared
    # to the last bit. The others use `sin` and `cos`, which differ in the last bit
    # between Julia and the C library on some arguments; with the C library swapped in
    # all of them agree to the last bit. Measured with Julia's functions (500 epochs):
    # derivative block equal to the last bit for precession 500 (old arithmetic 1),
    # for nutation 386 (old 0). Own random stream.
    @apart @testset "of-date rates as zzeprc76_, zzenut80_, zzmobliq_" begin
        r3 = MersenneTwister(SEED + 1702)
        n = 200
        nprec = 0
        nnut = 0
        fails = Failures()
        for et in random_epochs(r3, n, -span, span)
            mob = Ref(0.0); dmob = Ref(0.0)
            ccall((:zzmobliq_, SPICE.libcspice), Cint, (Ref{Float64}, Ref{Float64}, Ref{Float64}), et, mob, dmob)
            check!(fails, Almanac._dyn_mobliq_state(et) == (mob[], dmob[]), "zzmobliq_", et)
            for (sym, f) in ((:zzeprc76_, Almanac._dyn_prec_state), (:zzenut80_, Almanac._dyn_enut80_state))
                x = zeros(36)
                if sym === :zzeprc76_
                    ccall((:zzeprc76_, SPICE.libcspice), Cint, (Ref{Float64}, Ptr{Float64}), et, x)
                else
                    ccall((:zzenut80_, SPICE.libcspice), Cint, (Ref{Float64}, Ptr{Float64}), et, x)
                end
                x = reshape(x, 6, 6)       # column major, as Fortran stores it
                m = f(et)
                d = [m[i, j].d[1] for i in 1:3, j in 1:3]
                exact = d == x[4:6, 1:3]
                sym === :zzeprc76_ ? (nprec += exact) : (nnut += exact)
                rel = maximum(abs.(d .- x[4:6, 1:3])) / maximum(abs, x[4:6, 1:3])
                check!(fails, rel <= 1.0e-14, "rate block", sym, et, rel)
            end
        end
        @test isempty(fails)
        @test nprec >= 0.9 * n
        @test nnut >= 0.5 * n
    end

    @apart @testset "other corrections, spkpos, and frame to frame" begin
        d = Float64[]; dv = Float64[]; dq = Float64[]; dx = Float64[]; dxl = Float64[]
        names = first.(DynFix.FRAMES)
        for ab in ("LT", "CN", "CN+S", "XLT", "XLT+S", "XCN", "XCN+S"), name in names
            name == "DT_TWOVEC_IN_DYN" && continue
            for et in random_epochs(rng, 6, -span, span)
                s, lt = spkezr(ctx, 301, et, name, ab, 10)
                t, lu = SPICE.spkezr("301", et, name, ab, "10")
                push!(d, norm(s[1:3] .- t[1:3]) / norm(t[1:3]))
                push!(dv, norm(s[4:6] .- t[4:6]) / norm(t[4:6]))
                p, lp = spkpos(ctx, 301, et, name, ab, 10)
                q, lq = SPICE.spkpos("301", et, name, ab, "10")
                push!(dq, norm(p .- q) / norm(q))
            end
        end
        for _ in 1:400
            a, b = rand(rng, names), rand(rng, names)
            et = (2 * rand(rng) - 1) * 2.9e8
            # `DT_TWOVEC_IN_DYN` has a rounding noise of about 1e-11 in both engines
            # (see `DF_TOL`), and about 8 percent of the pairs have it on one side, so
            # these pairs are a sample of their own with the limit of that frame.
            noisy = (a == "DT_TWOVEC_IN_DYN") != (b == "DT_TWOVEC_IN_DYN")
            push!(noisy ? dxl : dx, maximum(abs.(pxform(ctx, a, b, et) .- SPICE.pxform(a, b, et))))
        end
        s1 = report_quantiles("spkezr, 7 other corrections: position, relative", d)
        s2 = report_quantiles("spkezr, 7 other corrections: velocity, relative", dv)
        s3 = report_quantiles("spkpos, 7 other corrections: position, relative", dq)
        s4 = report_quantiles("pxform between two dynamic or switch frames", dx)
        s5 = report_quantiles("pxform between DT_TWOVEC_IN_DYN and another frame", dxl)
        @test s1.max <= 5.0e-13
        @test s2.max <= 5.0e-7
        @test s3.max <= 5.0e-13
        # The limit of a pair is the sum of the `pxform` limits of its two frames:
        # 5e-14 for each frame without the noise (measured 1.6e-14 on 18000 pairs).
        @test s4.max <= 2 * DF_DEFAULT_TOL[1]
        @test s5.max <= DF_TOL["DT_TWOVEC_IN_DYN"][1] + DF_DEFAULT_TOL[1]
        @test length(dx) > 300 && length(dxl) > 5
    end

    @apart @testset "pxder against finite differences and sxform" begin
        lines = String[]
        fails = Failures()
        for (name, family, center) in DynFix.FRAMES
            id = namfrm(ctx, name)
            def = get(ks.dyn, Int32(id), nothing)
            switch = def === nothing
            frozen = !switch && def.frozen
            inertial = !switch && def.inertial
            loose = name == "DT_TWOVEC_IN_DYN"
            df_has_der(name) || continue
            e1 = 0.0; e2 = 0.0; es = 0.0
            for et in record_safe_epochs(rng, 12; span = 2.8e8)
                # Not near an epoch at which a switch frame changes its base.
                switch && any(e -> abs(et - e) < 200, DynFix.SWITCH_EDGES) && continue
                R, dR, ddR = pxder(ctx, "J2000", name, et, Val(DF_ORDER))
                # Order 0 is the rotation of `pxform`, computed with series arithmetic.
                within!(fails, maximum(abs.(R .- pxform(ctx, "J2000", name, et))), loose ? 1.0e-9 : 1.0e-13, "order 0 vs pxform", name, et)
                if frozen
                    check!(fails, iszero(dR) && iszero(ddR), "frozen: zero derivatives", name, et)
                    check!(fails, iszero(sxform(ctx, "J2000", name, et)[4:6, 1:3]), "frozen: zero sxform rate", name, et)
                    continue
                end
                h = 16.0
                fd1 = df_fd(t -> pxform(ctx, "J2000", name, t), et, h)
                fd2 = df_fd(t -> pxder(ctx, "J2000", name, t, Val(DF_LOW))[2], et, h)
                # The floors are the rounding noise of the differences: a frame that
                # turns slower than this is compared in absolute terms.
                e1 = max(e1, df_rel(dR, fd1, 1.0e-10))
                e2 = max(e2, df_rel(ddR, fd2, 1.0e-15))
                sx = sxform(ctx, "J2000", name, et)[4:6, 1:3]
                if inertial
                    check!(fails, iszero(sx), "inertial: zero sxform rate", name, et)
                    check!(fails, !iszero(dR), "inertial: nonzero derivative", name, et)
                else
                    # The floor is 1e-12 rad/s. The derivative of a true-of-date frame
                    # passes near zero where the rate of the nutation cancels the rate
                    # of the precession (measured: 1e-13 rad/s against 1e-11 typical),
                    # and the difference between `sxform` and `pxder` does not: it is the
                    # change of the nutation coefficients that the rates of `zzwahr_`
                    # omit (measured on 4000 epochs: at most 2.6e-17 rad/s for
                    # `DT_ETE_OF_DATE` and `DT_SWITCH_ALL`). A quotient by the derivative
                    # has no bound there; with the floor it is at most 2.6e-5.
                    es = max(es, df_rel(sx, dR, 1.0e-12))
                end
            end
            frozen && continue
            push!(
                lines, @sprintf(
                    "  %-17s order 1 vs FD=%.1e  order 2 vs FD=%.1e  sxform vs order 1=%.1e",
                    name, e1, e2, es
                )
            )
            @test e1 <= (loose ? 1.0e-3 : 2.0e-6)
            @test e2 <= (loose ? 1.0e-1 : 2.0e-5)
            # The derivative block of CSPICE is not the derivative of its rotation for
            # every definition (see the docstring of `sxform`); the difference is small.
            @test es <= (loose ? 1.0e-3 : 5.0e-5)
        end
        @test isempty(fails)
        measured("  pxder of the rotation from J2000 to each frame (relative to the largest element)")
        foreach(measured, lines)
        # Orders 0 to 3 are consistent with each other (orders 0 and 2 in the default run).
        for name in ("DT_GSE", "DT_ETE_OF_DATE", "DT_EULER", "DT_GSM", "DT_SWITCH_ALL", "DT_NADIR_LT")
            et = record_safe_epochs(rng, 1; span = 2.8e8)[1]
            j0 = pxder(ctx, "J2000", name, et, Val(0))
            if FULL
                j3 = pxder(ctx, "J2000", name, et, Val(3))
                j1 = pxder(ctx, "J2000", name, et, Val(1))
                @test maximum(abs.(j0[1] .- j3[1])) <= 1.0e-13 && maximum(abs.(j1[1] .- j3[1])) <= 1.0e-13
                @test df_rel(j1[2], j3[2]) <= 1.0e-12
            else
                j2 = pxder(ctx, "J2000", name, et, Val(DF_ORDER))
                @test maximum(abs.(j0[1] .- j2[1])) <= 1.0e-13
            end
        end
    end

    @apart @testset "frozen frames and inertial rotation state" begin
        fails = Failures()
        for (frozen, live, t0) in (
                ("DT_GSE_FROZEN", "DT_GSE", "2015-MAR-03 01:02:03 TDB"),
                ("DT_ETE_FROZEN", "DT_ETE_OF_DATE", "1995-JUN-15 00:00:00 TDB"),
            )
            et0 = str2et(ctx, t0)
            fixed = pxform(ctx, "J2000", live, et0)
            for et in random_epochs(rng, 20, -span, span)
                check!(fails, pxform(ctx, "J2000", frozen, et) == fixed, "frozen", frozen, et)
            end
        end
        for (inert, live) in (("DT_ETE_INERTIAL", "DT_ETE_OF_DATE"), ("DT_GSE_INERTIAL", "DT_GSE"))
            for et in random_epochs(rng, 20, -span, span)
                check!(fails, pxform(ctx, "J2000", inert, et) == pxform(ctx, "J2000", live, et), "inertial: pxform", inert, et)
                check!(fails, iszero(sxform(ctx, "J2000", inert, et)[4:6, 1:3]), "inertial: zero rate", inert, et)
                check!(fails, iszero(SPICE.sxform("J2000", inert, et)[4:6, 1:3]), "inertial: zero CSPICE rate", inert, et)
                check!(fails, pxder(ctx, "J2000", inert, et, Val(DF_LOW))[2] == pxder(ctx, "J2000", live, et, Val(DF_LOW))[2], "inertial: pxder", inert, et)
            end
        end
        @test isempty(fails)
    end

    @apart @testset "switch frames at the ends of the intervals" begin
        fails = Failures()
        worst = 0.0
        n = 0
        for edge in DynFix.SWITCH_EDGES, delta in (-1.0, -1.0e-3, 0.0, 1.0e-3, 1.0),
                name in ("DT_SWITCH_TIME", "DT_SWITCH_NEST")
            et = edge + delta
            want = df_oracle_kind(:pxform_c, "J2000", name, et)
            got = df_neo_kind(() -> pxform(ctx, "J2000", name, et))
            check!(fails, got == want, "error kind", name, et, want, got)
            got === :none || continue
            worst = max(worst, maximum(abs.(pxform(ctx, "J2000", name, et) .- SPICE.pxform("J2000", name, et))))
            n += 1
        end
        @test isempty(fails)
        measured(@sprintf("  %-58s n=%-7d max=%.3e", "switch frames at interval ends: pxform element", n, worst))
        @test worst <= 5.0e-14
        @test n >= 40
        # The base that is listed last has priority where the intervals overlap.
        @test pxform(ctx, "J2000", "DT_SWITCH_TIME", 2.0e8) == pxform(ctx, "J2000", "DT_EULER", 2.0e8)
        @test pxform(ctx, "J2000", "DT_SWITCH_TIME", 0.0) == pxform(ctx, "J2000", "ECLIPJ2000", 0.0)
        @test pxform(ctx, "J2000", "DT_SWITCH_TIME", -1.5e8) == pxform(ctx, "J2000", "DT_GSE", -1.5e8)
        @test pxform(ctx, "J2000", "DT_SWITCH_TIME", -5.0e8) == pxform(ctx, "J2000", "IAU_EARTH", -5.0e8)
        @test pxform(ctx, "J2000", "DT_SWITCH_ALL", 1.0e7) == pxform(ctx, "J2000", "DT_ETE_OF_DATE", 1.0e7)
        @test pxform(ctx, "J2000", "DT_SWITCH_NEST", 4.0e8) == pxform(ctx, "J2000", "IAU_MOON", 4.0e8)
        @test throws_kind(() -> pxform(ctx, "J2000", "DT_SWITCH_TIME", 7.0e8), :noframeconnect)
        @test throws_kind(() -> sxform(ctx, "J2000", "DT_SWITCH_TIME", -7.0e8), :noframeconnect)
        @test throws_kind(() -> spkezr(ctx, 301, 7.0e8, "DT_SWITCH_TIME", "NONE", 399), :noframeconnect)
    end

    @apart @testset "no allocation and inferred types" begin
        lts = Aberration("LT+S"); none = Aberration("NONE"); cn = Aberration("XCN+S")
        et = 1.2345e8
        nalloc = 0
        for (name, family, center) in DynFix.FRAMES
            id = namfrm(ctx, name)
            der = df_has_der(name)
            pxform(ks, 1, id, et); sxform(ks, 1, id, et)
            der && (pxder(ks, 1, id, et, Val(DF_ORDER)); pxder(ks, id, 17, et, Val(DF_LOW)))
            for ab in (lts, none, cn)
                spkezr(ks, 301, et, id, ab, 399); spkpos(ks, 301, et, id, ab, 399)
            end
            a = (
                df_alloc_px(ks, 1, id, et), df_alloc_px(ks, id, 10013, et), df_alloc_sx(ks, 1, id, et),
                der ? df_alloc_der(ks, 1, id, et, Val(DF_ORDER)) : 0,
                der ? df_alloc_der(ks, id, 17, et, Val(DF_LOW)) : 0,
                df_alloc_ezr(ks, et, id, lts), df_alloc_ezr(ks, et, id, none), df_alloc_ezr(ks, et, id, cn),
                df_alloc_pos(ks, et, id, lts), df_alloc_pos(ks, et, id, none),
            )
            all(==(0), a) || @info "allocation" name a
            @test all(==(0), a)
            nalloc += length(a)
            @test (@inferred pxform(ks, 1, id, et)) isa SMatrix{3, 3, Float64, 9}
            @test (@inferred sxform(ks, 1, id, et)) isa SMatrix{6, 6, Float64, 36}
            der && @test (@inferred pxder(ks, 1, id, et, Val(DF_ORDER))) isa NTuple{3, SMatrix{3, 3, Float64, 9}}
            @test (@inferred spkezr(ks, 301, et, id, lts, 399)) isa Tuple{SVector{6, Float64}, Float64}
            @test (@inferred spkezr(ks, 301, et, id, none, 399)) isa Tuple{SVector{6, Float64}, Float64}
            @test (@inferred spkpos(ks, 301, et, id, lts, 399)) isa Tuple{SVector{3, Float64}, Float64}
        end
        println("  ", nalloc, " calls with zero allocation")
        # The rotation is generic over the number type of the epoch.
        tv = Almanac.taylor_var(et, Val(2))
        id = namfrm(ctx, "DT_GSE_LTS")
        m = @inferred Almanac.frame_rotation(ks, Int32(1), Int32(id), tv)
        @test m isa SMatrix{3, 3, Almanac.Taylor{2, Float64}, 9}
        @test maximum(abs.(map(x -> x.v, m) .- pxform(ks, 1, id, et))) <= 1.0e-13
    end

    FULL && df_thread_check(ctx, rng)

    @apart @testset "definitions that are not valid, and the recursion limit" begin
        set2 = [set..., DynFix.error_kernel()]
        oracle_load(set2)
        ctx2 = neo_context(set2)     # the load itself does not throw
        et = 1.0e8
        lines = String[]
        for (name, kind) in DynFix.ERROR_FRAMES
            @test namfrm(ctx2, name) == SPICE.namfrm(name) != 0
            want = df_oracle_kind(:pxform_c, "J2000", name, et)
            got = df_neo_kind(() -> pxform(ctx2, "J2000", name, et))
            push!(lines, @sprintf("  %-17s Almanac :%-20s CSPICE :%s", name, got, want))
            got == want || @info "error kind" name got want
            @test got == want
            @test got == kind
            skind = get(DynFix.STATE_KINDS, name, kind)
            @test df_neo_kind(() -> sxform(ctx2, "J2000", name, et)) ==
                df_oracle_kind(:sxform_c, "J2000", name, et) == skind
            @test df_neo_kind(() -> spkezr(ctx2, 301, et, name, "LT+S", 399)) == skind
            if kind === :none
                @test maximum(abs.(pxform(ctx2, "J2000", name, et) .- SPICE.pxform("J2000", name, et))) <= 5.0e-10
            end
        end
        foreach(println, lines)
        # The valid frames of the same context are not affected.
        @test maximum(abs.(pxform(ctx2, "J2000", "DT_GSE", et) .- SPICE.pxform("J2000", "DT_GSE", et))) <= 5.0e-14
        # A class 5 frame without any definition, and a class above 6.
        extra = joinpath(mktempdir(), "empty.tf")
        write(
            extra, "KPL/FK\n\\begindata\n" * DynFix.head("DT_EMPTY", 1850131, 5) *
                DynFix.head("DT_CLASS7", 1850132, 7) * "\\begintext\n"
        )
        SPICE.furnsh(extra); furnsh(ctx2, extra)
        @test df_neo_kind(() -> pxform(ctx2, "J2000", "DT_EMPTY", et)) !== :none
        @test df_oracle_kind(:pxform_c, "J2000", "DT_EMPTY", et) !== :none
        k7 = df_neo_kind(() -> pxform(ctx2, "J2000", "DT_CLASS7", et))
        println(
            "  class 5 frame without a definition: Almanac :",
            df_neo_kind(() -> pxform(ctx2, "J2000", "DT_EMPTY", et)), "; frame of class 7: Almanac :", k7
        )
        @test k7 !== :none
        @test df_oracle_kind(:pxform_c, "J2000", "DT_CLASS7", et) !== :none
    end

    @apart @testset "angle units, a switch frame with frame IDs and time strings, type and size of an item" begin
        # Forms of a definition that the kernel of `DynFix` does not use. Each angle unit
        # has a factor of its own in the code (the kernel has RADIANS and DEGREES only);
        # the bases of a switch frame can be frame IDs and its bounds can be time strings.
        units = (
            ("ARCMINUTES", 600.0), ("ARCSECONDS", 36000.0), ("HOURANGLE", 1.5),
            ("MINUTEANGLE", 90.0), ("SECONDANGLE", 5400.0),
        )      # 10 degrees and 22.5 degrees
        io = IOBuffer()
        print(io, "KPL/FK\n\\begindata\n")
        for (k, (unit, a)) in enumerate(units)
            id = 1850140 + k
            print(
                io, DynFix.head("DT_UNIT_$k", id, 5), """
                       FRAME_$(id)_RELATIVE = 'J2000'
                       FRAME_$(id)_DEF_STYLE = 'PARAMETERIZED'
                       FRAME_$(id)_FAMILY = 'EULER'
                       FRAME_$(id)_EPOCH = @2010-JAN-01/12:00:00
                       FRAME_$(id)_AXES = ( 3, 1, 3 )
                       FRAME_$(id)_UNITS = '$unit'
                       FRAME_$(id)_ANGLE_1_COEFFS = ( $a, $(a * 1.0e-7) )
                       FRAME_$(id)_ANGLE_2_COEFFS = ( $(2a) )
                       FRAME_$(id)_ANGLE_3_COEFFS = ( $(-3a), $(a * 2.0e-8) )
                """
            )
        end
        # IAU_EARTH (10013) from 2000 to 2006, ECLIPJ2000 (17) from 2005 to 2010
        print(
            io, DynFix.head("DT_SWITCH_IDS", 1850150, 6), """
                   FRAME_1850150_ALIGNED_WITH = ( 10013, 17 )
                   FRAME_1850150_START = ( '2000 JAN 01 00:00:00 TDB', '2005 JAN 01 00:00:00 TDB' )
                   FRAME_1850150_STOP = ( '2006 JAN 01 00:00:00 TDB', '2010 JAN 01 00:00:00 TDB' )
            """
        )
        # an item of the wrong type, and an item with more values than it can have
        for (name, id, family) in (("DT_BADTYPE", 1850151, "5"), ("DT_BADSIZE", 1850152, "( 'EULER', 'EULER' )"))
            print(
                io, DynFix.head(name, id, 5), """
                       FRAME_$(id)_RELATIVE = 'J2000'
                       FRAME_$(id)_DEF_STYLE = 'PARAMETERIZED'
                       FRAME_$(id)_FAMILY = $family
                """
            )
        end
        print(io, "\\begintext\n")
        set4 = [set..., joinpath(mktempdir(), "forms.tf")]
        write(set4[end], String(take!(io)))
        oracle_load(set4)
        c = neo_context(set4)
        for k in eachindex(units), et in (-2.0e8, 3.0e8)
            a6 = sxform(c, "J2000", "DT_UNIT_$k", et)
            b6 = SPICE.sxform("J2000", "DT_UNIT_$k", et)
            @test maximum(abs.(a6[1:3, 1:3] .- b6[1:3, 1:3])) <= DF_DEFAULT_TOL[1]
            @test df_rel(a6[4:6, 1:3], b6[4:6, 1:3]) <= DF_DEFAULT_TOL[2]
        end
        # the units are not equal by accident: 10 degrees and 22.5 degrees
        @test pxform(c, "J2000", "DT_UNIT_1", 0.0) != pxform(c, "J2000", "DT_UNIT_3", 0.0)
        y(s) = str2et(c, s)
        for (t, base) in (
                ("2002 JUN 01 TDB", "IAU_EARTH"), ("2005 JUN 01 TDB", "ECLIPJ2000"),
                ("2008 JUN 01 TDB", "ECLIPJ2000"),
            )
            @test pxform(c, "J2000", "DT_SWITCH_IDS", y(t)) == pxform(c, "J2000", base, y(t))
            @test maximum(
                abs.(
                    pxform(c, "J2000", "DT_SWITCH_IDS", y(t)) .-
                        SPICE.pxform("J2000", "DT_SWITCH_IDS", y(t))
                )
            ) <= DF_DEFAULT_TOL[1]
        end
        @test df_neo_kind(() -> pxform(c, "J2000", "DT_SWITCH_IDS", y("2012 JUN 01 TDB"))) ===
            df_oracle_kind(:pxform_c, "J2000", "DT_SWITCH_IDS", y("2012 JUN 01 TDB")) === :noframeconnect
        for (name, kind) in (("DT_BADTYPE", :badvariabletype), ("DT_BADSIZE", :badvariablesize))
            @test df_neo_kind(() -> pxform(c, "J2000", name, 0.0)) ===
                df_oracle_kind(:pxform_c, "J2000", name, 0.0) === kind
        end
    end

    @apart @testset "definitions follow furnsh, unload and pool edits" begin
        c = neo_context([Kernels.LSK, Kernels.PCK11, Kernels.DE440S])
        et = 2.0e8
        @test namfrm(c, "DT_GSE") == 0
        furnsh(c, DynFix.frame_kernel())
        a = pxform(c, "J2000", "DT_EULER", et)
        @test a == pxform(ctx, "J2000", "DT_EULER", et)
        # An edit of a definition gives a new snapshot with the new definition.
        pdpool(c, "FRAME_1850013_ANGLE_2_COEFFS", [0.2])
        SPICE.kclear(); foreach(SPICE.furnsh, set); SPICE.pdpool("FRAME_1850013_ANGLE_2_COEFFS", [0.2])
        b = pxform(c, "J2000", "DT_EULER", et)
        @test b != a
        @test maximum(abs.(b .- SPICE.pxform("J2000", "DT_EULER", et))) <= 5.0e-14
        unload(c, DynFix.frame_kernel())
        @test namfrm(c, "DT_EULER") == 0
        @test Almanac.kernelset(c).inner === nothing
        @test throws_kind(() -> pxform(c, "J2000", "DT_EULER", et), :unknownframe) ||
            throws_kind(() -> pxform(c, "J2000", "DT_EULER", et), :noframe)
    end

    @apart @testset "switch frames over CK frames" begin
        set3 = [
            Kernels.LSK, Kernels.PCK11, Kernels.DE440S, CKFix.all_kernels()...,
            DynFix.switch_ck_kernel(),
        ]
        oracle_load(set3)
        c = neo_context(set3)
        tk = CKFix.ticks_tdb()
        t0 = sct2e(c, CKFix.SC_TDB, tk[1]) - 300.0
        t1 = sct2e(c, CKFix.SC_TDB, tk[end]) + 300.0
        # CSPICE keeps the segment of the last CK search of an instrument, and its answer
        # can then depend on the earlier calls (see `test_ck.jl`). A search far outside
        # the coverage clears that state before each query of the oracle.
        reset() = foreach(
            i -> (@assert SPICE.ckgp(i, 1.0e11, 0.0, "J2000") === nothing),
            (-77001, -77002, -77003, -77013)
        )
        fails = Failures()
        for name in ("DT_SWITCH_CK", "DT_SWITCH_CK_NOAV")
            dp = Float64[]; ds = Float64[]
            bases = name == "DT_SWITCH_CK" ? ("IAU_EARTH", "TESTSC_T3", "TESTSC_T2") :
                ("ECLIPJ2000", "TESTSC_T1", "TESTSC_T3_NOAV")
            used = zeros(Int, 3)
            nerr = 0
            for et in random_epochs(rng, 500, t0, t1)
                reset()
                want = df_oracle_kind(:pxform_c, "J2000", name, et)
                got = df_neo_kind(() -> pxform(c, "J2000", name, et))
                check!(fails, got == want, "pxform error kind", name, et, want, got)
                got === :none || (nerr += 1; continue)
                a = pxform(c, "J2000", name, et)
                reset()
                push!(dp, maximum(abs.(a .- SPICE.pxform("J2000", name, et))))
                for k in 3:-1:1
                    df_neo_kind(() -> pxform(c, "J2000", bases[k], et)) === :none || continue
                    a == pxform(c, "J2000", bases[k], et) && (used[k] += 1; break)
                end
                reset()
                wants = df_oracle_kind(:sxform_c, "J2000", name, et)
                gots = df_neo_kind(() -> sxform(c, "J2000", name, et))
                check!(fails, gots == wants, "sxform error kind", name, et, wants, gots)
                gots === :none || continue
                reset()
                push!(ds, maximum(abs.(sxform(c, "J2000", name, et) .- SPICE.sxform("J2000", name, et))))
            end
            println("  ", name, ": bases used (lowest priority first) ", used, ", epochs with an error: ", nerr)
            s1 = report_quantiles("$name: pxform element", dp)
            s2 = report_quantiles("$name: sxform element", ds)
            @test s1.max <= 1.0e-12
            @test s2.max <= 1.0e-12
            @test count(>(0), used) >= 2
        end
        @test isempty(fails)
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
