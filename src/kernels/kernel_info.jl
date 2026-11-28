# Information calls: the bodies and the coverage of an SPK file (`spkobj`, `spkcov`),
# the load list of a context (`ktotal`, `kdata`, `kinfo`) and the toolkit version.

"CSPICE version whose behavior and file formats this package follows."
const CSPICE_VERSION = "CSPICE_N0067"

"""
    tkvrsn(item="TOOLKIT") -> String

Return the version string of the toolkit, as
`tkvrsn_c`.

For `item = "TOOLKIT"` the result is `"Almanac_v<package version>/CSPICE_N0067"`.
The part after the slash is the CSPICE version whose behavior and file formats the
package follows; `tkvrsn_c` gives that part alone. Letter case is ignored, and the
symbol `:TOOLKIT` is accepted as in SPICE.jl. Any other item gives
`"No version found."`, as in CSPICE, and an empty one throws `SpiceError(:emptystring)`.

# Examples
```jldoctest
julia> endswith(tkvrsn(), "/CSPICE_N0067")
true
```

# References

- [NAIF documentation of `tkvrsn_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/tkvrsn_c.html)
"""
function tkvrsn(item::Union{AbstractString, Symbol} = "TOOLKIT")
    s = String(item)
    _chkfstr("tkvrsn", "string `item`", s)       # CHKFSTR_VAL of tkvrsn_c
    uppercase(strip(s)) == "TOOLKIT" || return "No version found."
    return string("Almanac_v", pkgversion(@__MODULE__), "/", CSPICE_VERSION)
end

# ---------------------------------------------------------------------------------------
# Load list
# ---------------------------------------------------------------------------------------

"""
    kernel_filetype(entry::Loaded) -> String

Return the type of a load-list entry as CSPICE reports it: `"SPK"`, `"CK"`, `"PCK"`,
`"DSK"`, `"EK"`, `"META"` for a text kernel that assigns `KERNELS_TO_LOAD`, or `"TEXT"`.
"""
function kernel_filetype(entry::Loaded)
    content = entry.content
    content isa TextKernel || return String(content.kind)
    return any(a -> a.name == "KERNELS_TO_LOAD", content.assignments) ? "META" : "TEXT"
end

"""
    _kind_filter(kind) -> Union{Nothing,Vector{String}}

Return the kernel types that the string `kind` selects: `nothing` (every type) if one of
its blank-separated words is `"ALL"`, else the words in upper case.
"""
function _kind_filter(kind::AbstractString)
    words = String[uppercase(w) for w in split(kind, ' '; keepempty = false)]
    "ALL" in words && return nothing
    return words
end

function _kernels_of_kind(ks::KernelSet, kind::AbstractString)
    words = _kind_filter(kind)
    words === nothing && return ks.loaded
    return Loaded[e for e in ks.loaded if kernel_filetype(e) in words]
end

"""
    ktotal([ctx,] kind="ALL") -> Int

Return the number of loaded kernels of the types in `kind`, as
`ktotal_c`.

`kind` is `"ALL"` or a blank-separated list of `SPK`, `CK`, `PCK`, `DSK`, `EK`, `TEXT`
and `META`, in any letter case; a list with the word `ALL` selects every kernel. A
meta-kernel counts once as `META`, and each file that it lists counts once. A path that
is loaded twice counts twice. A word that is not a kernel type selects nothing.

Throws `SpiceError(:emptystring)` if `kind` is the empty string.

# Examples
```jldoctest
julia> ktotal(SpiceContext())
0
```

See also [`kdata`](@ref), [`kinfo`](@ref).

# References

- [NAIF documentation of `ktotal_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ktotal_c.html)
"""
function ktotal(ks::KernelSet, kind::AbstractString = "ALL")
    _chkfstr("ktotal", "kind", kind)                # CHKFSTR of ktotal_c
    return length(_kernels_of_kind(ks, kind))
end
ktotal(ctx::SpiceContext, kind::AbstractString = "ALL") = ktotal(kernelset(ctx), kind)
ktotal(kind::AbstractString = "ALL") = ktotal(default_context(), kind)

"""
    kdata([ctx,] which, kind="ALL") -> Union{Nothing,Tuple{String,String,String,Int}}

Return `(file, filtyp, source, handle)` of loaded kernel number `which` among the
kernels of the types in `kind`, as
`kdata_c`, or `nothing` if there is no such kernel.

`which` counts from 1 in load order, as in SPICE.jl; CSPICE counts from 0. `kind` is
as for [`ktotal`](@ref). `file` is the path given to [`furnsh`](@ref), or for a file
of a meta-kernel its `KERNELS_TO_LOAD` entry after the `PATH_SYMBOLS` substitution.
`filtyp` is one of the type words of `ktotal`, and `source` is the meta-kernel that
loaded the file, or `""`. `handle` is always 0, because a loaded kernel has no file
handle.

Throws `SpiceError(:emptystring)` if `kind` is the empty string.

# Examples
```julia
for i in 1:ktotal("SPK")
    file, filtyp, source, _ = kdata(i, "SPK")
    println(file)
end
```

See also [`kinfo`](@ref).

# References

- [NAIF documentation of `kdata_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/kdata_c.html)
"""
function kdata(ks::KernelSet, which::Integer, kind::AbstractString = "ALL")
    _chkfstr("kdata", "kind", kind)                 # CHKFSTR of kdata_c
    list = _kernels_of_kind(ks, kind)
    (1 <= which <= length(list)) || return nothing
    e = list[which]
    return e.path, kernel_filetype(e), e.source, 0
end
kdata(ctx::SpiceContext, which::Integer, kind::AbstractString = "ALL") =
    kdata(kernelset(ctx), which, kind)
kdata(which::Integer, kind::AbstractString = "ALL") = kdata(default_context(), which, kind)

"""
    kinfo([ctx,] file) -> Union{Nothing,Tuple{String,String,Int}}

Return `(filtyp, source, handle)` of the loaded kernel with the path `file`, as
`kinfo_c`, or `nothing` if the file is not loaded.

`file` is compared as a string with the path given to [`furnsh`](@ref), as in CSPICE.
If the path is loaded more than once, the first entry is reported. The elements have
the meaning that they have in [`kdata`](@ref); `handle` is always 0.

Throws `SpiceError(:emptystring)` if `file` is empty.

# References

- [NAIF documentation of `kinfo_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/kinfo_c.html)
"""
function kinfo(ks::KernelSet, file::AbstractString)
    _chkfstr("kinfo", "file name", file)
    i = findfirst(e -> e.path == file, ks.loaded)
    i === nothing && return nothing
    e = ks.loaded[i]
    return kernel_filetype(e), e.source, 0
end
kinfo(ctx::SpiceContext, file::AbstractString) = kinfo(kernelset(ctx), file)
kinfo(file::AbstractString) = kinfo(default_context(), file)

# ---------------------------------------------------------------------------------------
# Bodies and coverage of an SPK file
# ---------------------------------------------------------------------------------------

"""
    _require_binary_kernel(who, path, kind) -> Nothing

Check with `getfat` that `path` is a binary kernel of type `kind` (`"SPK"`, `"CK"`,
`"PCK"`), as `spkobj`, `spkcov`, `ckobj`, `ckcov` and `pckcov` do. Throws
`SpiceError(:invalidformat)` for a transfer file, `SpiceError(:invalidarchtype)` for a
file that is not a DAF and `SpiceError(:invalidfiletype)` for a DAF of another type.
"""
function _require_binary_kernel(who::AbstractString, path::AbstractString, kind::AbstractString)
    arch, typ = getfat(path)
    arch == "XFR" && throw(
        SpiceError(
            :invalidformat, "$who: input file $path has " *
                "architecture $arch; the file must be a binary $kind file"
        )
    )
    arch == "DAF" || throw(
        SpiceError(
            :invalidarchtype, "$who: input file $path has " *
                "architecture $arch; the file must be a binary $kind file (DAF architecture)"
        )
    )
    typ == kind || throw(
        SpiceError(
            :invalidfiletype, "$who: input file $path has file " *
                "type $typ; the file must be a binary $kind file"
        )
    )
    return nothing
end

"""
    _each_summary(f, path) -> Nothing

Call `f(dc, ic, handle)` for the unpacked summary of every array of the DAF `path`, in
file order, on a read handle that is closed afterwards.
"""
function _each_summary(f, path::AbstractString)
    h = dafopr(path)
    try
        dafbfs(h)
        while daffna(h)                  # invariant: f was called for every array before this one
            dc, ic = dafus(dafgs(h), h.nd, h.ni)
            f(dc, ic, h)
        end
    finally
        dafcls(h)
    end
    return nothing
end

"Add the IDs `new` to the set `ids`, which is sorted and without duplicates on return."
function _add_ids!(ids::AbstractVector{<:Integer}, new::AbstractVector{<:Integer})
    append!(ids, new)
    sort!(unique!(ids))
    @assert issorted(ids) && allunique(ids)
    return ids
end

"""
    spkobj!(ids, path) -> ids

Insert the ID codes of the bodies that have data in the SPK file `path` into the
vector `ids`, as
`spkobj_c` does with a SPICE cell.

`ids` is used as a set: on return it is sorted and without duplicates.

!!! warning "Difference from CSPICE"
    The vector comes first, as in SPICE.jl; CSPICE has the cell last.

Throws the errors of [`spkobj`](@ref).

# References

- [NAIF documentation of `spkobj_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkobj_c.html)
"""
function spkobj!(ids::AbstractVector{<:Integer}, path::AbstractString)
    _require_binary_kernel("spkobj", path, "SPK")
    sort!(unique!(ids))
    _each_summary(path) do dc, ic, h
        k = searchsortedfirst(ids, ic[1])
        (k <= length(ids) && ids[k] == ic[1]) || insert!(ids, k, ic[1])
    end
    return ids
end

"""
    spkobj(path) -> Vector{Int}

Return the ID codes of the bodies that have data in the SPK file `path`, sorted and
without duplicates, as
`spkobj_c`.

The file is read from disk and need not be loaded. An open [`DAFHandle`](@ref)
(in memory, or of a file) can take the place of the path: the call then reads the
segments that the handle holds now, without an array that `dafbna` began and `dafena`
did not end, and throws `SpiceError(:dafnosuchhandle)` if the handle is closed.

Throws `SpiceError(:filenotfound)` if the file does not exist, `SpiceError(:invalidformat)`
if it is a transfer file (a binary kernel in text form), `SpiceError(:invalidarchtype)` if
it is not a DAF, and `SpiceError(:invalidfiletype)` if it is a DAF of another kernel type.

# Examples
```julia
for id in spkobj("de440s.bsp")
    w = spkcov("de440s.bsp", id)
    println(bodc2n(id), ": ", collect(w))
end
```

See also [`spkcov`](@ref), [`spkobj!`](@ref).

# References

- [NAIF documentation of `spkobj_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkobj_c.html)
"""
spkobj(path::AbstractString) = spkobj!(Int[], path)

"""
    spkcov!(cover, path, idcode) -> cover

Insert the coverage of the body `idcode` in the SPK file `path` into the window
`cover`, as
`spkcov_c`.

The result is the union of the old content of `cover` with the new intervals. See
[`spkcov`](@ref) for the intervals and the errors.

!!! warning "Difference from CSPICE"
    The window comes first, as in SPICE.jl; CSPICE has the cell last.

# References

- [NAIF documentation of `spkcov_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcov_c.html)
"""
function spkcov!(cover::Window, path::AbstractString, idcode::Integer)
    _require_binary_kernel("spkcov", path, "SPK")
    _each_summary(path) do dc, ic, h
        ic[1] == idcode && wninsd!(cover, dc[1], dc[2])
    end
    return cover
end

"""
    spkcov(path, idcode) -> Window

Return the coverage [`Window`](@ref) of the body `idcode` in the SPK file `path`, as
`spkcov_c`.

The window is the union of the intervals `[start, stop]` of the segments of that body,
in seconds past J2000 TDB; it is empty if the body has no data. `collect(w)` gives the
endpoints as a `Vector{Float64}`. The file need not be loaded, and an open
[`DAFHandle`](@ref) can take the place of the path, as for [`spkobj`](@ref).

Throws the errors of [`spkobj`](@ref).

See also [`spkcov!`](@ref), [`Window`](@ref), [`ckcov`](@ref), [`pckcov`](@ref).

# References

- [NAIF documentation of `spkcov_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcov_c.html)
"""
spkcov(path::AbstractString, idcode::Integer) = spkcov!(Window(), path, idcode)
