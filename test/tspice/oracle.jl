# The oracle of the replay: the CSPICE that cspice/build_cspice.sh builds into
# DIR/oracle/lib, which SPICE.jl loads in place of the library of CSPICE_jll through an
# override of the artifact of CSPICE_jll; and the functions of CORE-MATH that the same
# script builds beside it (DIR/oracle/lib/libcrmath.so), which replay.jl calls in place
# of Julia's sin, cos, exp, ... (crmath.jl). run.jl and selftest.jl give replay.jl the
# environments of `oracle_env` and `math_env`.

using TOML: TOML

"The library of the oracle of the recording `dir`."
oracle_lib(dir) = joinpath(abspath(dir), "oracle", "lib", "libcspice.so")

"The library of the functions of CORE-MATH of the recording `dir`."
crmath_lib(dir) = joinpath(abspath(dir), "oracle", "lib", "libcrmath.so")

"""
    cspice_artifact_hash(project)

The git tree hash of the artifact of CSPICE_jll in the environment `project`, read from
the Artifacts.toml of CSPICE_jll by a separate Julia process, so that no override
applies.
"""
function cspice_artifact_hash(project)
    code = """
    import Artifacts, SPICE
    toml = joinpath(pkgdir(SPICE.CSPICE_jll), "Artifacts.toml")
    print(Artifacts.artifact_meta("CSPICE", toml)["git-tree-sha1"])
    """
    hash = readchomp(`$(Base.julia_cmd()) --startup-file=no --project=$project -e $code`)
    @assert occursin(r"^[0-9a-f]{40}$", hash) "not the hash of an artifact: $hash"
    return hash
end

"""
    oracle_env(dir, project)

The environment in which SPICE.jl, in the environment `project`, loads the oracle of the
recording `dir`. Writes DIR/oracle_depot/artifacts/Overrides.toml, which maps the
artifact of CSPICE_jll to DIR/oracle, and returns JULIA_DEPOT_PATH (the depots of this
process, then that depot) and TSPICE_ORACLE_LIB (the real path of the library, which
replay.jl compares with the library that SPICE.jl loads).
"""
function oracle_env(dir, project)
    lib = oracle_lib(dir)
    @assert isfile(lib) "no oracle: $lib"
    depot = joinpath(abspath(dir), "oracle_depot")
    overrides = Dict(cspice_artifact_hash(project) => realpath(dirname(dirname(lib))))
    path = joinpath(mkpath(joinpath(depot, "artifacts")), "Overrides.toml")
    open(io -> TOML.print(io, overrides), path, "w")
    @assert TOML.parsefile(path) == overrides
    # Julia reads Overrides.toml from every depot, and writes its caches to the first, so
    # the depot of the override goes last; the list is given in full, because an empty
    # entry at the end would add only the depots bundled with Julia
    depots = unique([DEPOT_PATH; depot])
    return Dict("JULIA_DEPOT_PATH" => join(depots, ':'), "TSPICE_ORACLE_LIB" => realpath(lib))
end

"""
    math_env(dir, math)

The environment in which replay.jl computes with the elementary functions `math`:
`"core"`, those of CORE-MATH in the library of the recording `dir`, compared at tolerance
zero, or `"julia"`, Julia's own. Returns TSPICE_MATH and, for `"core"`, TSPICE_CRMATH_LIB
(the real path of the library).
"""
function math_env(dir, math)
    @assert math in ("core", "julia") "not a choice of functions: $math"
    math == "julia" && return Dict("TSPICE_MATH" => "julia")
    lib = crmath_lib(dir)
    @assert isfile(lib) "no library of CORE-MATH: $lib"
    return Dict("TSPICE_MATH" => "core", "TSPICE_CRMATH_LIB" => realpath(lib))
end
