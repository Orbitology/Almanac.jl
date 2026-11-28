# Self-test of the gate of run.jl: it must reject what it is meant to reject.
#
#     julia --project=test test/tspice/selftest.jl [DIR]
#
# DIR is a recording of run.jl (default ~/.cache/almanac-tspice) with the families
# F_LTIME, F_POOL, F_DWPOOL and F_TANGPT3, the oracle and the functions of CORE-MATH
# (DIR/oracle, built by cspice/build_cspice.sh). Each case of the gate runs run.jl on
# F_LTIME and F_POOL in a scratch copy of the recording and checks its exit status: 0 for
# the recording as it is, 1 for a mismatch that the baseline lacks and for a fault
# injected into Almanac's results (TSPICE_FAULT: an output off by one unit in the last
# place, which the comparison at tolerance zero finds, or an error of another kind), 2
# for a recording without the oracle, a recording without the functions of CORE-MATH, a
# recording of another source than tspice_c, a family without a log, a family that
# failed to record, and a family whose log cannot be replayed. The baseline files are
# copied, changed and restored; the scratch directories are removed at the end.
#
# One case checks that the oracle keeps no state from one family to the next: replay.jl,
# in one process, gives F_TANGPT3 the same mismatches after F_DWPOOL, which fills the
# watchers of CSPICE's kernel pool, as alone. It calls replay.jl, which replays its
# families in the order of their names, because run.jl gives a process its largest log
# first. It gives replay.jl the oracle and the functions of CORE-MATH as run.jl does
# (oracle.jl).

const HERE = @__DIR__
include(joinpath(HERE, "oracle.jl"))
const FAMS = ["F_LTIME", "F_POOL"]
const ORACLE_FAMS = ["F_DWPOOL", "F_TANGPT3"]

run_status(dir, fams) = run(ignorestatus(`$(Base.julia_cmd()) --startup-file=no --project=$(dirname(HERE)) $(joinpath(HERE, "run.jl")) --dir $dir --jobs 1 $fams`)).exitcode

"""
A scratch recording with the logs of `fams` and the oracle linked from `src`, its
manifest and status.
"""
function scratch(src, fams = FAMS)
    dir = mktempdir()
    symlink(joinpath(src, "oracle"), joinpath(dir, "oracle"))
    logs = mkpath(joinpath(dir, "logs"))
    for f in fams, ext in (".log.gz", ".log.files")
        p = joinpath(src, "logs", f * ext)
        ispath(p) && symlink(p, joinpath(logs, f * ext))
    end
    cp(joinpath(src, "logs", "manifest.tsv"), joinpath(logs, "manifest.tsv"))
    write(joinpath(logs, "status.tsv"), join(["$f\tpassed\n" for f in fams]))
    return dir
end

"""
    replay_ids(dir, fams, fam)

The mismatches of the family `fam`, sorted lines of ids.tsv, when replay.jl replays the
families `fams` of the recording `dir` in one process with the oracle; `nothing` if a
family does not replay to its end.
"""
function replay_ids(dir, fams, fam)
    @assert fam in fams
    out = mktempdir()
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(HERE)) $(joinpath(HERE, "replay.jl")) $(joinpath(dir, "logs")) $fams`
    env = merge(
        ENV, oracle_env(dir, dirname(HERE)), math_env(dir, "core"),
        Dict("TSPICE_OUT" => out, "TSPICE_ORACLE" => "1")
    )
    run(pipeline(ignorestatus(setenv(cmd, env)); stdout = joinpath(out, "out.txt"), stderr = joinpath(out, "err.txt")))
    done = all(l -> endswith(l, "\tdone"), readlines(joinpath(out, "families.tsv")))
    ids = sort(filter(l -> startswith(l, fam * "\t"), readlines(joinpath(out, "ids.tsv"))))
    rm(out; recursive = true)
    return done ? ids : nothing
end

function with_baseline_edit(f, edit)
    path = joinpath(HERE, "baseline_ids.tsv.gz")
    saved = read(path)
    try
        lines = readlines(`gzip -dc $path`)
        open(pipeline(`gzip -9n`; stdout = path), "w") do io
            foreach(l -> println(io, l), edit(lines))
        end
        return f()
    finally
        write(path, saved)
    end
end

function main(args)
    src = isempty(args) ? joinpath(homedir(), ".cache", "almanac-tspice") : abspath(args[1])
    (isfile(oracle_lib(src)) && isfile(crmath_lib(src))) ||
        (println("no oracle or no functions of CORE-MATH in $src: build them with record.sh or cspice/build_cspice.sh"); return 2)
    results = Pair{String, Bool}[]
    check(name, got, want) = (
        push!(results, name => got == want);
        println(got == want ? "ok    " : "FAIL  ", name, ": got ", got, ", expected ", want)
    )

    dir = scratch(src)
    check("the recording as it is", run_status(dir, FAMS), 0)

    # a mismatch of F_POOL removed from the baseline is a new mismatch (a `cspice` call is
    # not a mismatch)
    got = with_baseline_edit(
        ls -> (
            i = findfirst(l -> startswith(l, "F_POOL\t") && split(l, '\t')[3] != "cspice", ls);
            i === nothing ? ls : deleteat!(ls, i)
        )
    ) do
        run_status(dir, FAMS)
    end
    check("a mismatch that the baseline lacks", got, 1)

    # faults of Almanac that the comparison must find: an output off by one unit in the
    # last place, which CSPICE does not share, and an error of another kind
    for (fault, what) in ("value" => "an output off by one ulp", "kind" => "an error of another kind")
        got = withenv(() -> run_status(dir, FAMS), "TSPICE_FAULT" => fault)
        check(what, got, 1)
    end

    rm(joinpath(dir, "oracle"))
    check("a recording without the oracle", run_status(dir, FAMS), 2)
    # the oracle without the functions of CORE-MATH
    lib = mkpath(joinpath(dir, "oracle", "lib"))
    for f in ("libcspice.so", "libcsupport.so")
        symlink(joinpath(src, "oracle", "lib", f), joinpath(lib, f))
    end
    check("a recording without the functions of CORE-MATH", run_status(dir, FAMS), 2)
    rm(joinpath(dir, "oracle"); recursive = true)
    symlink(joinpath(src, "oracle"), joinpath(dir, "oracle"))

    # the baseline is made from the recording of tspice_c: a recording of another
    # source (the older PEPPER recording) is refused
    mf = joinpath(dir, "logs", "manifest.tsv")
    manifest = read(mf, String)
    write(mf, replace(manifest, r"^source\t.*$"m => "source\tpepper"))
    check("a recording of another source", run_status(dir, FAMS), 2)
    write(mf, manifest)

    rm(joinpath(dir, "logs", "F_LTIME.log.gz"))
    check("a family without a log", run_status(dir, FAMS), 2)

    write(joinpath(dir, "logs", "F_LTIME.log.gz"), "not a gzip file")
    check("a family that cannot be replayed", run_status(dir, FAMS), 2)

    rm(joinpath(dir, "logs", "F_LTIME.log.gz"))
    symlink(joinpath(src, "logs", "F_LTIME.log.gz"), joinpath(dir, "logs", "F_LTIME.log.gz"))
    write(joinpath(dir, "logs", "status.tsv"), "F_LTIME\tfailed\ttimeout\nF_POOL\tpassed\n")
    check("a family that failed to record", run_status(dir, FAMS), 2)

    rm(dir; recursive = true)

    # the oracle keeps no state from one family to the next: the mismatches of F_TANGPT3
    # after F_DWPOOL are those of F_TANGPT3 alone, which has some
    dir = scratch(src, ORACLE_FAMS)
    alone = replay_ids(dir, ["F_TANGPT3"], "F_TANGPT3")
    after = replay_ids(dir, ORACLE_FAMS, "F_TANGPT3")
    same = alone !== nothing && after !== nothing && !isempty(alone) && after == alone
    check("F_TANGPT3 after F_DWPOOL, as alone", same, true)
    rm(dir; recursive = true)

    nfail = count(!last, results)
    println(nfail == 0 ? "the gate rejects every case" : "$nfail case(s) not rejected")
    return nfail == 0 ? 0 : 1
end

exit(main(ARGS))
