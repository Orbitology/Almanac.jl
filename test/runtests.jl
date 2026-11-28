# Almanac.jl verification suite. CSPICE (through SPICE.jl) is the oracle.
#
# Run all groups, with threads (group 7 needs more than one):
#
#     julia --startup-file=no --threads=8 --project=test test/runtests.jl
#
# Run some groups only (names are the file names without `test_` and `.jl`):
#
#     julia --startup-file=no --threads=8 --project=test test/runtests.jl states threads
#
# Run the groups in several processes, which is faster but needs several GB of memory for
# each process that compiles an AD backend:
#
#     julia --startup-file=no --threads=8 test/parallel.jl
#
# ALMANAC_TEST_FULL=1 selects the exhaustive run (large samples, all derivative orders).
# ALMANAC_TEST_EXCLUDE lists groups, separated by spaces, that a run without group names
# leaves out (CI runs the AD and static groups on one or two jobs only).
#
# The kernels are downloaded into `test/data/` on the first run (about 110 MB); the
# environment variables of `test/kernels.jl` change their paths and the download. A
# group whose kernels are missing is skipped with a visible notice, or fails with
# ALMANAC_REQUIRE_KERNELS=1. ALMANAC_TEST_SEED changes the random samples.

include("helpers.jl")
Kernels.fetch_kernels()

const GROUPS = [
    "states", "rotations", "derivatives", "identities", "pool", "management",
    "threads", "quality", "aberration", "spktypes", "time", "sclk", "ck", "dafwrite",
    "windows", "roundtrip", "sgp4", "math_vectors", "math_rotations", "math_coordinates",
    "math_ellipsoid", "math_conics", "math_polynomials", "dynframes", "poolwrite", "frameinfo", "das", "dsk",
    "gf_solver", "gf_quantities", "gf_posc", "gf_threads", "gf_surface", "gf_illum", "gf_occult",
    "gf_fov",
    "geom_subpoints", "geom_illum", "geom_occult", "geom_limbterm", "geom_fov", "geom_constobs",
    "geom_azel", "table",
    "ek_file", "ek_parse", "ek_query", "ek_threads", "ek_write", "jsonkernel", "providers", "queries",
    "memkernels", "emptystring", "alloccheck", "ad_enzyme", "ad_mooncake", "ad_zygote", "ad_tracked", "static",
]

"Wall time of each group of this run, in the order of the run."
const GROUP_SECONDS = Pair{String, Float64}[]

excluded = split(get(ENV, "ALMANAC_TEST_EXCLUDE", ""))
selected = isempty(ARGS) ? setdiff(GROUPS, excluded) : ARGS
unknown = setdiff(union(selected, excluded), GROUPS)
isempty(unknown) || error("unknown test group(s) $(unknown); the groups are $(GROUPS)")
isempty(ARGS) && foreach(g -> push!(SKIPPED, (g, "left out by ALMANAC_TEST_EXCLUDE")), intersect(GROUPS, excluded))

println("Almanac tests: Julia ", VERSION, ", ", Threads.nthreads(), " thread(s), seed ", SEED)
Threads.nthreads() == 1 && printstyled(
    "NOTE  one thread only: group 7 will be skipped; run " *
        "with `julia --threads=8 --project=test test/runtests.jl`\n"; color = :yellow, bold = true
)

try
    @testset "Almanac" begin
        for g in selected
            push!(GROUP_SECONDS, g => @elapsed include(joinpath(@__DIR__, "test_$g.jl")))
        end
    end
finally
    SPICE.kclear()
    kclear()
    println("\nMeasured values (for tightening the thresholds):")
    foreach(println, MEASURED)
    println("\nTime of each group (s), longest first:")
    for (g, t) in sort(GROUP_SECONDS; by = last, rev = true)
        println("  ", rpad(g, 18), lpad(round(t; digits = 1), 7))
    end
    println("  ", rpad("all", 18), lpad(round(sum(last, GROUP_SECONDS; init = 0.0); digits = 1), 7))
    if !isempty(SKIPPED)
        printstyled("\nSkipped:\n"; color = :yellow, bold = true)
        for (name, reason) in SKIPPED
            println("  ", name, ": ", reason)
        end
    end
end
