# Group 18: kernel-pool edits, the read calls that go with them, and the watchers,
# against CSPICE.
#
# The same sequence of operations runs on CSPICE and on a Almanac context; after each
# step the whole pool is compared (names, types, values). This fixes the rule for
# edits and the load list: an edit stays through `furnsh` and through the `unload` of a
# binary kernel, and is lost with the `unload` of a text kernel and with `kclear`.
#
# Conventions that differ from SPICE.jl: `gnpool` returns the names in ascending order
# (CSPICE: order of its hash table) and an empty vector if nothing matches; `stpool`
# returns `nothing` for an absent string.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Almanac keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALMANAC_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

"`badkpv_c` of CSPICE: its result, or the short error message as a lower-case symbol."
function pw_oracle_badkpv(caller, name, comp, size, divby, type)
    r = ccall(
        (:badkpv_c, SPICE.libcspice), Cint, (Cstring, Cstring, Cstring, Cint, Cint, UInt8),
        caller, name, comp, size, divby, UInt8(type)
    )
    ccall((:failed_c, SPICE.libcspice), Cint, ()) == 0 && return r != 0
    buf = zeros(UInt8, 64)
    ccall((:getmsg_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Ptr{UInt8}), "SHORT", 64, buf)
    ccall((:reset_c, SPICE.libcspice), Cvoid, ())
    return Symbol(lowercase(replace(unsafe_string(pointer(buf)), "SPICE(" => "", ")" => "")))
end

"Pool variable names of CSPICE that match `template`, sorted."
function pw_oracle_names(template::AbstractString)
    out = String[]
    start = 1
    while true
        part = SPICE.gnpool(template, start, 500)
        part === nothing && break
        append!(out, part)
        length(part) < 500 && break
        start += 500
    end
    return sort!(unique!(out))
end

"""
    pw_pool_differences(ctx) -> Vector{String}

Names of the pool variables that differ between CSPICE and the context `ctx`: present
on one side only, or with another type, count or value. Empty if the pools are equal.
"""
function pw_pool_differences(ctx)
    want = pw_oracle_names("*")
    got = gnpool(ctx, "*")
    bad = symdiff(want, got)
    for name in intersect(want, got)
        n, type = SPICE.dtpool(name)
        if dtpool(ctx, name) != (n, type)
            push!(bad, name)
        elseif type === :N
            gdpool(ctx, name) == SPICE.gdpool(name; room = n) || push!(bad, name)
            expool(ctx, name) || push!(bad, name)
        else
            gcpool(ctx, name) == SPICE.gcpool(name; room = n, lenout = 128) || push!(bad, name)
            expool(ctx, name) && push!(bad, name)
        end
    end
    return bad
end

"Write a text kernel with the assignment lines `lines` to a new temporary file."
function pw_kernel(name::AbstractString, lines)
    path = joinpath(mktempdir(), name)
    write(path, "KPL/PCK\n\n\\begindata\n" * join(lines, "\n") * "\n\\begintext\n")
    return path
end

"Run the same operation on CSPICE (`fo`) and on Almanac (`fn`)."
pw_both(fo, fn) = (fo(); fn(); nothing)

group(
    "18. kernel-pool edits and watchers against CSPICE", Kernels.LSK, Kernels.PCK11,
    Kernels.MOON_TF, Kernels.DE440S
) do
    rng = MersenneTwister(SEED + 18)
    ka = pw_kernel("a.tpc", ["PW_A = ( 1.0, 2.0, 3.0 )", "PW_SHARED = 10", "PW_TEXT = 'from a'"])
    kb = pw_kernel(
        "b.tpc", [
            "PW_B = ( 4.0, 5.0 )", "PW_SHARED = 20", "PW_SHARED += 21",
            "PW_WATCHED = 7",
        ]
    )
    kc = pw_kernel("c.tpc", ["PW_C = 'c'", "PW_EDIT_D = -1.0", "PW_DATE = @2000-JAN-01/12:00:00"])

    @apart @testset "edits, furnsh, unload and kclear" begin
        set = [Kernels.LSK, Kernels.PCK11, ka, Kernels.DE440S]
        oracle_load(set)
        ctx = neo_context(set)
        steps = 0
        check(what) = (
            bad = pw_pool_differences(ctx); steps += 1;
            isempty(bad) || @info "pool differs after: $what" first(bad, 8);
            @test isempty(bad); @test ktotal(ctx) == SPICE.ktotal("ALL")
        )
        check("load")
        pw_both(() -> SPICE.pdpool("PW_EDIT_D", [1.5, -2.5, 1.0e300]), () -> pdpool(ctx, "PW_EDIT_D", [1.5, -2.5, 1.0e300]))
        pw_both(() -> SPICE.pipool("PW_EDIT_I", [1, -2, 3]), () -> pipool(ctx, "PW_EDIT_I", [1, -2, 3]))
        pw_both(() -> SPICE.pcpool("PW_EDIT_C", ["one", "two words", "x"]), () -> pcpool(ctx, "PW_EDIT_C", ["one", "two words", "x"]))
        check("pdpool, pipool, pcpool of new variables")
        @test gdpool(ctx, "PW_EDIT_D") == [1.5, -2.5, 1.0e300]
        @test gipool(ctx, "PW_EDIT_I") == [1, -2, 3]
        @test gcpool(ctx, "PW_EDIT_C") == ["one", "two words", "x"]

        pw_both(() -> SPICE.pdpool("BODY399_RADII", [7000.0, 7000.0, 6900.0]), () -> pdpool(ctx, "BODY399_RADII", [7000.0, 7000.0, 6900.0]))
        pw_both(() -> SPICE.pcpool("PW_A", ["now text"]), () -> pcpool(ctx, "PW_A", ["now text"]))
        pw_both(() -> SPICE.pdpool("PW_TEXT", [3.0]), () -> pdpool(ctx, "PW_TEXT", [3.0]))
        check("replace variables of kernels, also with another type")
        @test bodvrd(ctx, "EARTH", "RADII") == SPICE.bodvrd("EARTH", "RADII") == [7000.0, 7000.0, 6900.0]

        pw_both(() -> SPICE.dvpool("PW_EDIT_I"), () -> dvpool(ctx, "PW_EDIT_I"))
        pw_both(() -> SPICE.dvpool("BODY499_PM"), () -> dvpool(ctx, "BODY499_PM"))
        pw_both(() -> SPICE.dvpool("PW_NOT_THERE"), () -> dvpool(ctx, "PW_NOT_THERE"))
        check("dvpool")
        @test !expool(ctx, "PW_EDIT_I") && !SPICE.expool("PW_EDIT_I")

        lines = [
            "PW_L1 = 1", "PW_L2 = ( 1.0, 2.0D0, 3E1 )", "PW_L3 = 'text'", "PW_L2 += 4",
            "PW_L4 = ( 'a', 'b' )", "PW_L4 += 'c'", "PW_L5 = @2001-FEB-03", "PW_EDIT_D = 9.0",
        ]
        pw_both(() -> SPICE.lmpool(lines), () -> lmpool(ctx, lines))
        check("lmpool")
        pw_both(() -> SPICE.ldpool(kc), () -> ldpool(ctx, kc))
        check("ldpool")
        @test gdpool(ctx, "PW_EDIT_D") == [-1.0]

        # An edit stays through a later furnsh; the new kernel assigns on top of it.
        pw_both(() -> SPICE.pdpool("PW_SHARED", [99.0]), () -> pdpool(ctx, "PW_SHARED", [99.0]))
        pw_both(() -> SPICE.pdpool("PW_KEEP", [5.0]), () -> pdpool(ctx, "PW_KEEP", [5.0]))
        pw_both(() -> SPICE.furnsh(kb), () -> furnsh(ctx, kb))
        check("furnsh after edits")
        @test gdpool(ctx, "PW_KEEP") == [5.0]
        @test gdpool(ctx, "PW_SHARED") == [20.0, 21.0]
        pw_both(() -> SPICE.furnsh(Kernels.MOON_TF), () -> furnsh(ctx, Kernels.MOON_TF))
        check("furnsh of a frame kernel after edits")

        # The unload of a binary kernel leaves the pool and the edits as they are.
        pw_both(() -> SPICE.unload(Kernels.DE440S), () -> unload(ctx, Kernels.DE440S))
        check("unload of a binary kernel")
        @test gdpool(ctx, "PW_KEEP") == [5.0]
        @test gdpool(ctx, "BODY399_RADII") == [7000.0, 7000.0, 6900.0]

        # The unload of a text kernel builds the pool again: every edit is lost.
        pw_both(() -> SPICE.unload(kb), () -> unload(ctx, kb))
        check("unload of a text kernel")
        @test !expool(ctx, "PW_KEEP")
        @test !expool(ctx, "PW_L1")
        @test gdpool(ctx, "PW_SHARED") == [10.0]
        @test gdpool(ctx, "BODY399_RADII") == [6378.1366, 6378.1366, 6356.7519]

        # A kernel that is not loaded: nothing changes, also not the edits.
        pw_both(() -> SPICE.pdpool("PW_KEEP", [6.0]), () -> pdpool(ctx, "PW_KEEP", [6.0]))
        pw_both(() -> SPICE.unload(kb), () -> unload(ctx, kb))
        check("unload of a kernel that is not loaded")
        @test gdpool(ctx, "PW_KEEP") == [6.0]

        # clpool empties the pool and keeps the load list.
        pw_both(() -> SPICE.clpool(), () -> clpool(ctx))
        check("clpool")
        @test isempty(gnpool(ctx, "*"))
        @test ktotal(ctx) == 4
        pw_both(() -> SPICE.pdpool("PW_AFTER_CLEAR", [1.0]), () -> pdpool(ctx, "PW_AFTER_CLEAR", [1.0]))
        check("pdpool after clpool")
        pw_both(() -> SPICE.unload(Kernels.MOON_TF), () -> unload(ctx, Kernels.MOON_TF))
        check("unload of a text kernel after clpool")
        @test gdpool(ctx, "PW_A") == [1.0, 2.0, 3.0]
        @test !expool(ctx, "PW_AFTER_CLEAR")

        pw_both(() -> SPICE.pdpool("PW_KEEP", [7.0]), () -> pdpool(ctx, "PW_KEEP", [7.0]))
        pw_both(() -> SPICE.kclear(), () -> kclear(ctx))
        check("kclear")
        @test isempty(gnpool(ctx, "*"))
        @test ktotal(ctx) == 0
        pw_both(() -> SPICE.pdpool("PW_ALONE", [1.0, 2.0]), () -> pdpool(ctx, "PW_ALONE", [1.0, 2.0]))
        check("pdpool into an empty context")
        println("  pool equal to CSPICE after each of ", steps, " steps")
    end

    @apart @testset "random sequence of operations" begin
        SPICE.kclear()
        ctx = SpiceContext()
        files = [
            pw_kernel(
                "r$i.tpc", [
                    "PW_R$(i) = $(i)", "PW_V$(i % 3) = ( $(i).5, 2.0 )",
                    "PW_S$(i % 2) = 'file $i'", "PW_V$(i % 3) += $(i)",
                ]
            ) for i in 1:6
        ]
        loaded = String[]
        vars = ["PW_V0", "PW_V1", "PW_V2", "PW_S0", "PW_S1", "PW_X1", "PW_X2", "PW_X3"]
        nbad = 0
        counts = Dict{Symbol, Int}()
        # 150 operations reach each of the nine kinds with the seed of the suite (the
        # rarest, `clpool` and `kclear`, have a probability of 0.03 each: 4 and 2 times).
        nops = nsample(150, 400)
        for step in 1:nops
            r = rand(rng)
            op = r < 0.2 ? :pdpool : r < 0.35 ? :pcpool : r < 0.45 ? :pipool : r < 0.6 ? :dvpool :
                r < 0.7 ? :lmpool : r < 0.82 ? :furnsh : r < 0.94 ? :unload : r < 0.97 ? :clpool : :kclear
            counts[op] = get(counts, op, 0) + 1
            v = rand(rng, vars)
            if op === :pdpool
                x = randn(rng, rand(rng, 1:5))
                pw_both(() -> SPICE.pdpool(v, x), () -> pdpool(ctx, v, x))
            elseif op === :pcpool
                x = [randstring(rng, rand(rng, 1:20)) for _ in 1:rand(rng, 1:4)]
                pw_both(() -> SPICE.pcpool(v, x), () -> pcpool(ctx, v, x))
            elseif op === :pipool
                x = rand(rng, -1000:1000, rand(rng, 1:5))
                pw_both(() -> SPICE.pipool(v, x), () -> pipool(ctx, v, x))
            elseif op === :dvpool
                pw_both(() -> SPICE.dvpool(v), () -> dvpool(ctx, v))
            elseif op === :lmpool
                x = ["$v = ( $(rand(rng, 1:9)), $(rand(rng, 1:9)) )", "$(rand(rng, vars)) += $(rand(rng, 1:9))"]
                # `+=` on a character variable with a number is an error in both engines.
                ok = try
                    SPICE.lmpool(x); true
                catch e
                    e isa SPICE.SpiceError || rethrow(); false
                end
                okn = try
                    lmpool(ctx, x); true
                catch e
                    e isa Almanac.SpiceError || rethrow(); false
                end
                ok == okn || (nbad += 1)
                (ok && okn) || (SPICE.kclear(); kclear(ctx); empty!(loaded))
            elseif op === :furnsh
                f = rand(rng, files)
                pw_both(() -> SPICE.furnsh(f), () -> furnsh(ctx, f))
                push!(loaded, f)
            elseif op === :unload
                f = rand(rng, files)
                pw_both(() -> SPICE.unload(f), () -> unload(ctx, f))
            elseif op === :clpool
                pw_both(() -> SPICE.clpool(), () -> clpool(ctx))
            else
                pw_both(() -> SPICE.kclear(), () -> kclear(ctx))
                empty!(loaded)
            end
            bad = pw_pool_differences(ctx)
            if !isempty(bad) || ktotal(ctx) != SPICE.ktotal("ALL")
                nbad += 1
                nbad <= 3 && @info "pool differs" step op bad ktotal(ctx) SPICE.ktotal("ALL")
            end
        end
        println(
            "  ", nops, " random operations ", sort!(collect(counts); by = first),
            ", steps with a difference: ", nbad
        )
        @test nbad == 0
    end

    @apart @testset "expool, dtpool, gnpool, szpool, stpool, badkpv" begin
        set = [Kernels.LSK, Kernels.PCK11, Kernels.MOON_TF, ka]
        oracle_load(set)
        ctx = neo_context(set)
        ks = Almanac.kernelset(ctx)
        long = [
            "PW_LONG = ( 'first part //', 'second part //', 'end.', 'alone', 'x //', 'y' )",
            "PW_MARK = ( 'a -', 'b -', 'c', 'd' )", "PW_NUM = ( 1, 2 )",
        ]
        pw_both(() -> SPICE.lmpool(long), () -> lmpool(ctx, long))
        for name in (
                "BODY399_RADII", "FRAME_31008_NAME", "PW_A", "PW_TEXT", "NOT_THERE",
                "body399_radii", "DELTET/DELTA_AT", "PW_LONG", "", "BODY399_RADII ",
            )
            isempty(name) && continue
            @test expool(ctx, name) == SPICE.expool(name)
            @test dtpool(ctx, name) == SPICE.dtpool(name)
        end
        @test dtpool(ctx, "NOT_THERE") == (0, :X)
        @test (@inferred dtpool(ks, "BODY399_RADII")) == (3, :N)
        @test (@inferred expool(ks, "BODY399_RADII")) === true

        for t in (
                "*", "BODY399_*", "BODY%%%_RADII", "*_PM", "FRAME_*_NAME", "body399_r*",
                "PW_*", "BODY3*9_*A*", "%%%%", "*RADII", "BODY399_RADII", "NOTHING*", "*%",
                "TKFRAME_31___*", "DELTET/*",
            )
            want = pw_oracle_names(t)
            got = gnpool(ctx, t)
            got == want || @info "gnpool" t length(got) length(want)
            @test got == want
            @test issorted(got)
        end
        all = gnpool(ctx, "BODY*")
        @test length(all) > 100
        @test gnpool(ctx, "BODY*", 1, 10) == all[1:10]
        @test gnpool(ctx, "BODY*", 11, 5) == all[11:15]
        @test gnpool(ctx, "BODY*", length(all), 5) == all[end:end]
        @test isempty(gnpool(ctx, "BODY*", length(all) + 1, 5))
        @test gnpool(ctx, "BODY*", -3, 2) == all[1:2]
        @test isempty(gnpool(ctx, "NOTHING*")) && SPICE.gnpool("NOTHING*", 1, 10) === nothing
        # The pages of CSPICE have the same sizes (the order inside is its hash order).
        @test length(SPICE.gnpool("BODY*", 11, 5)) == 5
        @test throws_kind(() -> gnpool(ctx, "*", 1, 0), :badarraysize)

        for name in (
                "MAXVAR", "MAXVAL", "MAXLIN", "MAXCHR", "MXNOTE", "MAXLEN", "MAXAGT",
                "maxvar", " MaxLen ", "NOTHING", "MAXVARS",
            )
            @test szpool(name) == SPICE.szpool(name)
        end

        fails = Failures()
        for (item, contin) in (
                    ("PW_LONG", "//"), ("PW_MARK", "-"), ("PW_LONG", "-"),
                    ("PW_MARK", "//"), ("FRAME_31008_NAME", "//"),
                ), nth in 1:7
            got = stpool(ctx, item, nth, contin)
            want = SPICE.stpool(item, nth, contin)
            if got === nothing
                check!(fails, isempty(want), "stpool nothing", item, nth, contin)
            else
                check!(fails, got == want, "stpool", item, nth, contin)
            end
        end
        @test isempty(fails)
        @test stpool(ctx, "PW_LONG", 1, "//") == "first part second part end."
        @test stpool(ctx, "PW_LONG", 2, "//") == "alone"
        @test stpool(ctx, "PW_LONG", 3, "//") == "x y"
        @test stpool(ctx, "PW_LONG", 4, "//") === nothing
        @test stpool(ctx, "PW_NUM", 1, "//") === nothing
        @test stpool(ctx, "NOT_THERE", 1, "//") === nothing
        @test stpool(ctx, "PW_LONG", 0, "//") === nothing

        # badkpv: each comparison, the divisor (negative, zero), the type, and each error
        # kind, against CSPICE. PW_LONG holds 6 strings, BODY399_RADII 3 numbers.
        for args in (
                ("ME", "PW_LONG", "=", 6, 3, 'C'), ("ME", "PW_LONG", "<", 7, 2, 'c'),
                ("ME", "PW_LONG", ">", 5, -3, 'C'), ("ME", "PW_LONG", "<=", 6, 6, 'C'),
                ("ME", "PW_LONG", "=>", 6, 1, 'C'), ("ME", "PW_LONG", "= ", 6, 1, 'C'),
                ("ME", "PW_LONG", "<", 6, 1, 'C'), ("ME", "PW_LONG", "=", 6, 4, 'C'),
                ("ME", "PW_LONG", "=", 6, 0, 'C'), ("ME", "PW_LONG", ">=", 6, 1, 'C'),
                ("ME", "PW_LONG", " =", 6, 1, 'C'), ("ME", "PW_LONG", "=", 6, 1, 'N'),
                ("ME", "BODY399_RADII", "=", 3, 3, 'C'), ("ME", "BODY399_RADII", "=", 3, 3, 'n'),
                ("ME", "BODY399_RADII", "=", 3, 3, 'X'), ("ME", "NOT_THERE", "=", 3, 3, 'N'),
                ("", "BODY399_RADII", "=", 3, 3, 'N'), ("ME", "BODY399_RADII", "", 3, 3, 'N'),
            )
            got = try
                badkpv(ctx, args...)
            catch e
                e isa SpiceError ? e.kind : rethrow()
            end
            want = pw_oracle_badkpv(args...)
            got === want || @info "badkpv" args got want
            @test got === want
        end
    end

    @apart @testset "tables that are derived from the pool follow an edit" begin
        set = [Kernels.LSK, Kernels.PCK11, Kernels.DE440S]
        oracle_load(set)
        ctx = neo_context(set)
        et = 1.234e8
        before = pxform(ctx, "J2000", "IAU_EARTH", et)
        pm = [200.0, 361.0, 0.0]
        pw_both(() -> SPICE.pdpool("BODY399_PM", pm), () -> pdpool(ctx, "BODY399_PM", pm))
        after = pxform(ctx, "J2000", "IAU_EARTH", et)
        @test after != before
        @test maximum(abs.(after .- SPICE.pxform("J2000", "IAU_EARTH", et))) <= 1.0e-13

        frame = [
            "FRAME_PW_TK = 1777001", "FRAME_1777001_NAME = 'PW_TK'", "FRAME_1777001_CLASS = 4",
            "FRAME_1777001_CLASS_ID = 1777001", "FRAME_1777001_CENTER = 399",
            "TKFRAME_1777001_SPEC = 'ANGLES'", "TKFRAME_1777001_RELATIVE = 'IAU_EARTH'",
            "TKFRAME_1777001_ANGLES = ( 10.0, 20.0, 30.0 )", "TKFRAME_1777001_AXES = ( 3, 2, 1 )",
            "TKFRAME_1777001_UNITS = 'DEGREES'",
        ]
        @test namfrm(ctx, "PW_TK") == 0
        pw_both(() -> SPICE.lmpool(frame), () -> lmpool(ctx, frame))
        @test namfrm(ctx, "PW_TK") == SPICE.namfrm("PW_TK") == 1777001
        @test maximum(abs.(pxform(ctx, "PW_TK", "J2000", et) .- SPICE.pxform("PW_TK", "J2000", et))) <= 1.0e-13
        a = spkezr(ctx, "MOON", et, "PW_TK", "LT+S", "EARTH")[1]
        b = SPICE.spkezr("MOON", et, "PW_TK", "LT+S", "EARTH")[1]
        @test norm(a[1:3] .- b[1:3]) <= 1.0e-12 * norm(b[1:3])

        pw_both(() -> SPICE.pcpool("NAIF_BODY_NAME", ["PW_PROBE", "PW_OTHER"]), () -> pcpool(ctx, "NAIF_BODY_NAME", ["PW_PROBE", "PW_OTHER"]))
        pw_both(() -> SPICE.pipool("NAIF_BODY_CODE", [-9901, -9902]), () -> pipool(ctx, "NAIF_BODY_CODE", [-9901, -9902]))
        @test bodn2c(ctx, "pw_probe") == SPICE.bodn2c("pw_probe") == -9901
        @test bodc2n(ctx, -9902) == SPICE.bodc2n(-9902) == "PW_OTHER"
        pw_both(() -> SPICE.dvpool("NAIF_BODY_NAME"), () -> dvpool(ctx, "NAIF_BODY_NAME"))
        pw_both(() -> SPICE.dvpool("NAIF_BODY_CODE"), () -> dvpool(ctx, "NAIF_BODY_CODE"))
        @test bodn2c(ctx, "PW_PROBE") === nothing && SPICE.bodn2c("PW_PROBE") === nothing

        # Another context does not see the edits.
        other = neo_context(set)
        @test pxform(other, "J2000", "IAU_EARTH", et) == before
        @test namfrm(other, "PW_TK") == 0
        @test_throws MethodError pdpool(Almanac.kernelset(ctx), "X", [1.0])
    end

    @apart @testset "boddef" begin
        set = [Kernels.LSK, Kernels.PCK11]
        oracle_load(set)
        ctx = neo_context(set)
        # CSPICE keeps these pairs for the whole process, also through `kclear`, so the
        # codes are ones that no other group uses.
        pairs = (
            ("PW_BODDEF_ONE", 7654321), ("Pw Boddef  Two", -7654322), ("PW_BODDEF_ONE", 7654323),
            ("PW_ALIAS", 7654324), ("PW_ALIAS_AGAIN", 7654324), ("PW_SECOND_NAME", 7654321),
        )
        for (name, code) in pairs
            pw_both(() -> SPICE.boddef(name, code), () -> boddef(ctx, name, code))
        end
        names = (
            "PW_BODDEF_ONE", "pw_boddef_one", "PW BODDEF TWO", "pw boddef two", "PW_ALIAS",
            "PW_ALIAS_AGAIN", "EARTH", "PW_SECOND_NAME", "PW_NOTHING", "MOON",
        )
        codes = (7654321, -7654322, 7654323, 7654324, 399, 301, 5)
        function compare()
            fails = Failures()
            for n in names
                check!(fails, bodn2c(ctx, n) == SPICE.bodn2c(n), "bodn2c", n)
            end
            for c in codes
                check!(fails, bodc2n(ctx, c) == SPICE.bodc2n(c), "bodc2n", c)
            end
            @test isempty(fails)
            return nothing
        end
        compare()
        @test bodn2c(ctx, "PW_BODDEF_ONE") == 7654323
        # The pairs stay through furnsh, unload and kclear, as in CSPICE.
        pw_both(() -> SPICE.furnsh(ka), () -> furnsh(ctx, ka))
        compare()
        pw_both(() -> SPICE.unload(ka), () -> unload(ctx, ka))
        compare()
        pw_both(() -> SPICE.kclear(), () -> kclear(ctx))
        compare()
        @test bodn2c(ctx, "PW_BODDEF_ONE") == 7654323
        # Regression: a blank name was `:badarg`, an empty one was accepted.
        for (name, kind) in (("   ", :blanknameassigned), ("", :emptystring))
            ccall((:boddef_c, SPICE.libcspice), Cvoid, (Cstring, Cint), name, 7654326)
            @test cspice_kind() === kind
            @test throws_kind(() -> boddef(ctx, name, 7654326), kind)
        end
        # Regression: a name over 36 characters counts by its first 36, as in `zzbodtrn`;
        # the key is made from the compressed full name and then cut, the name kept is
        # the first 36 characters as given.
        long = "PW  BODDEF  LONG  NAME  OF  MORE  THAN  36  CHARACTERS"
        # `boddef_c` directly: `SPICE.boddef` refuses a name of more than 35 characters.
        pw_both(
            () -> ccall((:boddef_c, SPICE.libcspice), Cvoid, (Cstring, Cint), long, 7654325),
            () -> boddef(ctx, long, 7654325)
        )
        @test bodc2n(ctx, 7654325) == c_bodc2n(7654325)
        for n in (long, "PW BODDEF LONG NAME OF MORE THAN 36 XYZ", "pw boddef long name of more than 36", first(long, 36))
            @test bodn2c(ctx, n) == SPICE.bodn2c(n)
        end
        @test bodn2c(ctx, long) == 7654325
        @test bodn2c(SpiceContext(), "PW_BODDEF_ONE") === nothing
        # CSPICE N0067 (`zzboddef`) signals TOOMANYPAIRS past 161 pairs; Almanac has no
        # limit. CSPICE is not called: it would keep the pairs for the whole process.
        many = SpiceContext()
        for i in 1:200
            boddef(many, "PW_MANY_$i", 7_600_000 + i)
        end
        @test all(i -> bodn2c(many, "PW_MANY_$i") == 7_600_000 + i, 1:200)
        @test all(i -> bodc2n(many, 7_600_000 + i) == "PW_MANY_$i", 1:200)
    end

    @apart @testset "watchers" begin
        set = [Kernels.LSK, ka]
        oracle_load(set)
        ctx = neo_context(set)
        nsteps = 0
        same(agent) = (
            nsteps += 1; a = cvpool(ctx, agent); b = SPICE.cvpool(agent);
            a == b || @info "cvpool differs" agent nsteps a b; @test a == b; a
        )
        pw_both(() -> SPICE.swpool("PW_AGENT_1", ["PW_WATCHED", "PW_SHARED"]), () -> swpool(ctx, "PW_AGENT_1", ["PW_WATCHED", "PW_SHARED"]))
        pw_both(() -> SPICE.swpool("PW_AGENT_2", ["PW_OTHER"]), () -> swpool(ctx, "PW_AGENT_2", ["PW_OTHER"]))
        @test same("PW_AGENT_1") === true        # true once after swpool
        @test same("PW_AGENT_1") === false
        @test same("PW_AGENT_2") === true
        @test same("PW_AGENT_3") === false       # an agent that never set a watch

        pw_both(() -> SPICE.pdpool("PW_WATCHED", [1.0]), () -> pdpool(ctx, "PW_WATCHED", [1.0]))
        @test same("PW_AGENT_1") === true
        @test same("PW_AGENT_1") === false
        @test same("PW_AGENT_2") === false
        pw_both(() -> SPICE.pdpool("PW_UNWATCHED", [1.0]), () -> pdpool(ctx, "PW_UNWATCHED", [1.0]))
        same("PW_AGENT_1"); same("PW_AGENT_2")
        pw_both(() -> SPICE.pcpool("PW_OTHER", ["x"]), () -> pcpool(ctx, "PW_OTHER", ["x"]))
        same("PW_AGENT_1")
        @test same("PW_AGENT_2") === true
        pw_both(() -> SPICE.pipool("PW_SHARED", [3]), () -> pipool(ctx, "PW_SHARED", [3]))
        @test same("PW_AGENT_1") === true
        pw_both(() -> SPICE.dvpool("PW_WATCHED"), () -> dvpool(ctx, "PW_WATCHED"))
        @test same("PW_AGENT_1") === true
        pw_both(() -> SPICE.dvpool("PW_WATCHED"), () -> dvpool(ctx, "PW_WATCHED"))    # absent now
        same("PW_AGENT_1")
        pw_both(() -> SPICE.lmpool(["PW_OTHER = 'y'", "PW_Q = 1"]), () -> lmpool(ctx, ["PW_OTHER = 'y'", "PW_Q = 1"]))
        same("PW_AGENT_1")
        @test same("PW_AGENT_2") === true
        pw_both(() -> SPICE.furnsh(kb), () -> furnsh(ctx, kb))      # assigns PW_WATCHED and PW_SHARED
        @test same("PW_AGENT_1") === true
        same("PW_AGENT_2")
        pw_both(() -> SPICE.furnsh(kc), () -> furnsh(ctx, kc))      # assigns none of the watched names
        same("PW_AGENT_1"); same("PW_AGENT_2")
        pw_both(() -> SPICE.ldpool(kb), () -> ldpool(ctx, kb))
        @test same("PW_AGENT_1") === true
        same("PW_AGENT_2")
        pw_both(() -> SPICE.unload(kc), () -> unload(ctx, kc))      # a text kernel: the pool is built again
        same("PW_AGENT_1"); same("PW_AGENT_2")
        pw_both(() -> SPICE.clpool(), () -> clpool(ctx))
        same("PW_AGENT_1"); same("PW_AGENT_2")
        pw_both(() -> SPICE.swpool("PW_AGENT_1", ["PW_MORE"]), () -> swpool(ctx, "PW_AGENT_1", ["PW_MORE"]))
        @test same("PW_AGENT_1") === true
        pw_both(() -> SPICE.pdpool("PW_MORE", [1.0]), () -> pdpool(ctx, "PW_MORE", [1.0]))
        @test same("PW_AGENT_1") === true
        pw_both(() -> SPICE.pdpool("PW_WATCHED", [2.0]), () -> pdpool(ctx, "PW_WATCHED", [2.0]))
        same("PW_AGENT_1")
        pw_both(() -> SPICE.kclear(), () -> kclear(ctx))
        same("PW_AGENT_1"); same("PW_AGENT_2")
        pw_both(() -> SPICE.pdpool("PW_WATCHED", [3.0]), () -> pdpool(ctx, "PW_WATCHED", [3.0]))
        same("PW_AGENT_1"); same("PW_AGENT_2")

        # dwpool: an agent with an update that cvpool has not taken is refused; a deleted
        # agent gets no update; an agent that is not known is no error.
        c_dwpool(agent) = (
            ccall((:dwpool_, SPICE.libcspice), Cint, (Cstring, Cint), agent, sizeof(agent));
            cspice_kind()
        )
        pw_both(() -> SPICE.pdpool("PW_SHARED", [4.0]), () -> pdpool(ctx, "PW_SHARED", [4.0]))
        @test throws_kind(() -> dwpool(ctx, "PW_AGENT_1"), c_dwpool("PW_AGENT_1"))
        @test c_dwpool("PW_AGENT_1") === :updatepending
        same("PW_AGENT_1"); same("PW_AGENT_2")
        @test (dwpool(ctx, "PW_AGENT_1"), c_dwpool("PW_AGENT_1")) === (nothing, :none)
        pw_both(() -> SPICE.pdpool("PW_SHARED", [5.0]), () -> pdpool(ctx, "PW_SHARED", [5.0]))
        @test same("PW_AGENT_1") === false
        @test (dwpool(ctx, "PW_AGENT_9"), c_dwpool("PW_AGENT_9")) === (nothing, :none)
        pw_both(() -> SPICE.swpool("PW_AGENT_1", ["PW_SHARED"]), () -> swpool(ctx, "PW_AGENT_1", ["PW_SHARED"]))
        @test same("PW_AGENT_1") === true
        println("  cvpool equal to CSPICE at ", nsteps, " steps")

        # The watches of two contexts are separate.
        c1 = SpiceContext(); c2 = SpiceContext()
        swpool(c1, "A", ["V"]); swpool(c2, "A", ["V"])
        @test cvpool(c1, "A") && cvpool(c2, "A")
        pdpool(c1, "V", [1.0])
        @test cvpool(c1, "A") && !cvpool(c2, "A")

        # The version counter moves with every change of the pool and with nothing else.
        v0 = poolversion(c1)
        pdpool(c1, "W", [1.0]); v1 = poolversion(c1)
        dvpool(c1, "NOT_THERE"); v2 = poolversion(c1)
        dvpool(c1, "W"); v3 = poolversion(c1)
        furnsh(c1, ka); v4 = poolversion(c1)
        furnsh(c1, Kernels.DE440S); v5 = poolversion(c1)
        unload(c1, Kernels.DE440S); v6 = poolversion(c1)
        kclear(c1); v7 = poolversion(c1)
        @test v0 < v1 == v2 < v3 < v4 == v5 == v6 < v7
        @test (@inferred poolversion(c1)) isa Int
        @test poolversion(c2) < v7
    end

    @apart @testset "errors" begin
        ctx = SpiceContext()
        @test throws_kind(() -> pdpool(ctx, "   ", [1.0]), :badvarname)
        # Regression: an empty name is SPICE(EMPTYSTRING) of the C wrapper, and an empty
        # list of values returns before the name check of pdpool_ and pcpool_.
        ccall((:pdpool_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Ptr{Cdouble}), "", 1, [1.0])
        @test throws_kind(() -> pdpool(ctx, "", [1.0]), cspice_kind())
        @test throws_kind(() -> pdpool(ctx, "", Float64[]), :emptystring)
        @test throws_kind(() -> pcpool(ctx, "", String[]), :emptystring)
        ccall((:pdpool_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Ptr{Cdouble}), "X"^33, 0, [1.0])
        @test cspice_kind() === :none
        @test pdpool(ctx, "X"^33, Float64[]) === nothing
        @test pipool(ctx, "X"^33, Int[]) === nothing
        @test pcpool(ctx, "X"^33, String[]) === nothing
        @test throws_kind(() -> pcpool(ctx, "X"^33, ["a"]), :badvarname)
        @test throws_kind(() -> pipool(ctx, "X"^33, [1]), :badvarname)
        pdpool(ctx, "X"^32, [1.0])
        @test gdpool(ctx, "X"^32) == [1.0]
        pdpool(ctx, "PW_EMPTY", Float64[])
        @test !expool(ctx, "PW_EMPTY")
        pcpool(ctx, "PW_CUT", ["a"^100])
        @test gcpool(ctx, "PW_CUT") == ["a"^80]
        # `lmpool` of a string over 80 characters: CSPICE N0067 (`zzrvbf`) writes past the
        # slot of the value, so it is not called. Almanac keeps the first 80 characters,
        # and the next value and the next variable are intact.
        lmpool(ctx, ["PW_LONGSTR = ( '" * "a"^100 * "', 'b' )", "PW_NEXT = 'c'"])
        @test gcpool(ctx, "PW_LONGSTR") == ["a"^80, "b"]
        @test gcpool(ctx, "PW_NEXT") == ["c"]
        # Regression: `lmpool` took `\begintext` and `\begindata` as markers and replaced
        # tabs by blanks. `zzrvbf` has no markers, and a tab is kept inside a string.
        function lmpool_c(lines)
            len = maximum(ncodeunits, lines) + 1
            buf = zeros(UInt8, len, length(lines))
            for (j, l) in enumerate(lines)
                buf[1:ncodeunits(l), j] .= codeunits(l)
            end
            ccall((:lmpool_c, SPICE.libcspice), Cvoid, (Ptr{UInt8}, Cint, Cint), buf, len, length(lines))
            return cspice_kind()
        end
        marked = ["\\begintext", "PW_MARKED = 1"]
        kind = lmpool_c(marked)
        SPICE.dvpool("PW_MARKED")
        @test kind !== :none
        @test throws_kind(() -> lmpool(ctx, marked), kind)
        tabbed = ["PW_TAB = ( 'a\tb'\t'c' )"]
        pw_both(() -> SPICE.lmpool(tabbed), () -> lmpool(ctx, tabbed))
        @test gcpool(ctx, "PW_TAB") == SPICE.gcpool("PW_TAB") == ["a\tb", "c"]
        pw_both(() -> SPICE.dvpool("PW_TAB"), () -> dvpool(ctx, "PW_TAB"))
        @test throws_kind(() -> ldpool(ctx, "/no/such/file.tpc"), :fileopenfailed)
        @test throws_kind(() -> gnpool(ctx, "*", 1, -1), :badarraysize)
    end

    # A malformed text kernel: the kind is the short message of CSPICE (tspice F_POOL,
    # F_BODCOD), for ldpool and for furnsh of a context with and without direct edits.
    @apart @testset "malformed text kernels give the error kinds of CSPICE" begin
        cases = [
            ["PW_BAD"], ["PW_BAD 1 2"], ["PW_BAD = ( ABC )"], ["PW_BAD = ''"],
            ["PW_BAD = ( )"], ["PW_BAD += ( )"], ["PW_BAD = (", " )"], ["PW_BAD = 1, 'a'"],
            ["PW_BAD = 'a', 1"], ["PW_BAD_" * "X"^30 * " = 1"], ["PW_BAD = @"],
            ["PW_BAD = @JUNK"], ["PW_BAD = @2000 JAN 1 PST"],
            ["NAIF_BODY_NAME = 'PW_SPUD'", "NAIF_BODY_CODE = ( 1009, 1010 )"],
            ["NAIF_BODY_NAME = 'PW_SPUD'"], ["NAIF_BODY_NAME = ' '", "NAIF_BODY_CODE = 1009"],
            # Regression: a name that ends in a nonprinting character is NONPRINTINGCHAR.
            ["PW_BAD\x7f = 1"],
            # Regression: a `+=` of a string to a numeric variable of an earlier line is
            # TYPEMISMATCH on its own line, before the error of a later line.
            ["PW_BAD = 1", "PW_BAD += 'x'", "PW_BAD = ( ABC )"],
        ]
        for (i, lines) in enumerate(cases)
            k = pw_kernel("bad$i.tpc", lines)
            ccall((:furnsh_c, SPICE.libcspice), Cvoid, (Cstring,), k)
            kc = cspice_kind()
            ccall((:unload_c, SPICE.libcspice), Cvoid, (Cstring,), k)
            for v in ("PW_BAD", "NAIF_BODY_NAME", "NAIF_BODY_CODE")
                ccall((:dvpool_c, SPICE.libcspice), Cvoid, (Cstring,), v)
            end
            cspice_kind()
            @test kc !== :none
            @test throws_kind(() -> furnsh(SpiceContext(), k), kc)
            edited = SpiceContext()
            pdpool(edited, "PW_EDIT", [1.0])
            @test throws_kind(() -> furnsh(edited, k), kc)
            any(l -> startswith(l, "NAIF"), lines) ||
                @test throws_kind(() -> ldpool(SpiceContext(), k), kc)
        end
    end

    # Regression: NAIF_BODY_NAME and NAIF_BODY_CODE of the wrong types were MISSINGKPV or
    # BADDIMENSIONS. `zzbodker` reads them with `gcpool` and `gipool`, which find neither,
    # and then defines no names and signals nothing.
    @apart @testset "body name variables of the wrong type" begin
        for (i, lines) in enumerate(
                (["NAIF_BODY_NAME = 5", "NAIF_BODY_CODE = 'PW_X'"], ["NAIF_BODY_NAME = ( 5, 6 )"])
            )
            k = pw_kernel("types$i.tpc", lines)
            ctx = SpiceContext()
            try
                pw_both(() -> SPICE.furnsh(k), () -> furnsh(ctx, k))
                @test bodn2c(ctx, "EARTH") == SPICE.bodn2c("EARTH") == 399
                @test bodc2n(ctx, 301) == SPICE.bodc2n(301) == "MOON"
            finally
                SPICE.unload(k)
                cspice_kind()
            end
        end
    end

    if Threads.nthreads() == 1
        skip_notice("18. threads", "one thread only")
    else
        @testset "threads: edits, queries and watchers on one context" begin
            set = [Kernels.LSK, Kernels.PCK11, Kernels.DE440S]
            ctx = neo_context(set)
            ets = random_epochs(rng, 50, -5.0e8, 5.0e8)
            ref = [spkezr(ctx, 301, et, "J2000", "LT+S", 399)[1] for et in ets]
            radii = gdpool(ctx, "BODY399_RADII")
            nwriters, nrounds = 4, 150
            stop = Threads.Atomic{Bool}(false)
            torn = Threads.Atomic{Int}(0)
            wrong = Threads.Atomic{Int}(0)
            reads = Threads.Atomic{Int}(0)
            for w in 1:nwriters
                swpool(ctx, "PW_T$w", ["PW_THREAD_$w"])
                cvpool(ctx, "PW_T$w")
            end
            v0 = poolversion(ctx)
            readers = map(1:4) do r
                Threads.@spawn begin
                    last = 0
                    while !stop[]
                        # A writer always stores equal values: a mixture is a torn read.
                        for w in 1:nwriters
                            ks = Almanac.kernelset(ctx)
                            expool(ks, "PW_THREAD_$w") || continue
                            x = gdpool(ks, "PW_THREAD_$w")
                            all(==(x[1]), x) || Threads.atomic_add!(torn, 1)
                        end
                        gdpool(ctx, "BODY399_RADII") == radii || Threads.atomic_add!(wrong, 1)
                        i = rand(1:length(ets))
                        spkezr(ctx, 301, ets[i], "J2000", "LT+S", 399)[1] == ref[i] ||
                            Threads.atomic_add!(wrong, 1)
                        v = poolversion(ctx)
                        v >= last || Threads.atomic_add!(wrong, 1)
                        last = v
                        Threads.atomic_add!(reads, 1)
                        # The root task stays on thread 1: a reader that never yields there
                        # keeps it from setting `stop` (seen on Julia 1.10).
                        yield()
                    end
                end
            end
            seen = zeros(Int, nwriters)
            writers = map(1:nwriters) do w
                Threads.@spawn begin
                    for k in 1:nrounds
                        pdpool(ctx, "PW_THREAD_$w", fill(Float64(k), 40))
                        cvpool(ctx, "PW_T$w") && (seen[w] += 1)
                        k % 50 == 0 && dvpool(ctx, "PW_SCRATCH_$w")
                        pcpool(ctx, "PW_SCRATCH_$w", ["round $k"])
                    end
                end
            end
            foreach(fetch, writers)
            stop[] = true
            foreach(fetch, readers)
            println(
                "  ", nwriters * nrounds * 2, " concurrent edits, ", reads[], " reader rounds, torn reads: ",
                torn[], ", wrong answers: ", wrong[]
            )
            @test torn[] == 0
            @test wrong[] == 0
            @test reads[] > 0
            for w in 1:nwriters
                @test gdpool(ctx, "PW_THREAD_$w") == fill(Float64(nrounds), 40)
                @test gcpool(ctx, "PW_SCRATCH_$w") == ["round $nrounds"]
                # Only the owner of an agent clears its update, so each edit was seen.
                @test seen[w] == nrounds
                @test !cvpool(ctx, "PW_T$w")
            end
            @test poolversion(ctx) - v0 >= nwriters * nrounds * 2
            @test gdpool(ctx, "BODY399_RADII") == radii
            @test ktotal(ctx) == 3
        end
    end
    @apart @testset "static vectors as values" begin
        ctx = SpiceContext()
        pdpool(ctx, "SV_D", SVector(1.0, 2.0))
        pipool(ctx, "SV_I", SVector(3, 4))
        pcpool(ctx, "SV_C", SVector("a", "b"))
        @test gdpool(ctx, "SV_D") == [1.0, 2.0]
        @test gipool(ctx, "SV_I") == [3, 4]
        @test gcpool(ctx, "SV_C") == ["a", "b"]
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
