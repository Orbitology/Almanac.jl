# FixedString: an ASCII string of bounded length stored inline, and the byte writers of
# the number and time formats. A writer takes a buffer and returns it: a
# `Vector{UInt8}` is changed in place, a `FixedString` is immutable and comes back as a
# new value. So `buf = _putint(buf, n)` builds either one, and a `FixedString` never
# touches the heap.

"""
    FixedString(s::AbstractString)

An ASCII string of at most 64 bytes stored inline, as an `SVector` stores its elements,
so that creating or returning one allocates nothing. [`et2utc`](@ref) and
[`etcal`](@ref) return a `FixedString`.

A `FixedString` is an `AbstractString`: it compares equal to the `String` with the same
characters, and prints, interpolates and indexes as one. `String(s)` converts it, for a
call that requires a `String`.

Throws `SpiceError(:buffertoosmall)` if `s` has more than 64 bytes, and
`SpiceError(:nonprintablechars)` if `s` is not ASCII. This type has no CSPICE
counterpart.

# Examples
```jldoctest
julia> s = Almanac.FixedString("ABC")
"ABC"

julia> s == "ABC", length(s), String(s) isa String
(true, 3, true)
```
"""
struct FixedString <: AbstractString
    words::NTuple{8, UInt64}     # byte `i` is byte `(i - 1) % 8` of word `(i - 1) ÷ 8 + 1`
    len::Int
    function FixedString(words::NTuple{8, UInt64}, len::Int)
        0 <= len <= 64 || _fixedstring_overflow(len)
        return new(words, len)
    end
end

@noinline _fixedstring_overflow(len) = throw(
    SpiceError(:buffertoosmall, "a FixedString holds at most 64 bytes, not $len")
)

FixedString() = FixedString(ntuple(_ -> zero(UInt64), Val(8)), 0)

function FixedString(s::AbstractString)
    isascii(s) || throw(SpiceError(:nonprintablechars, "a FixedString holds ASCII characters only"))
    return _puts(FixedString(), s)
end

@inline _fsbyte(s::FixedString, i::Int) =
    (s.words[((i - 1) >> 3) + 1] >> (((i - 1) & 7) << 3)) % UInt8

Base.ncodeunits(s::FixedString) = s.len
Base.codeunit(::FixedString) = UInt8
Base.@propagate_inbounds function Base.codeunit(s::FixedString, i::Int)
    @boundscheck checkbounds(s, i)
    return _fsbyte(s, i)
end
Base.isvalid(s::FixedString, i::Int) = checkbounds(Bool, s, i)
Base.@propagate_inbounds function Base.iterate(s::FixedString, i::Int = 1)
    i > s.len && return nothing
    return Char(_fsbyte(s, i)), i + 1
end
Base.sizeof(s::FixedString) = s.len
function Base.String(s::FixedString)
    v = Base.StringVector(s.len)
    @inbounds for i in 1:(s.len)
        v[i] = _fsbyte(s, i)
    end
    return String(v)
end

# The buffer operations. Bytes past the length of a `FixedString` are undefined: a
# writer that grows a buffer with `_bresize` writes every new byte.
_blen(b::Vector{UInt8}) = length(b)
_blen(b::FixedString) = b.len
@inline _bget(b::FixedString, i::Int) = _fsbyte(b, i)
@inline _bset(b::Vector{UInt8}, i::Int, c::UInt8) = (b[i] = c; b)
@inline function _bset(b::FixedString, i::Int, c::UInt8)
    k = ((i - 1) >> 3) + 1
    sh = ((i - 1) & 7) << 3
    w = b.words
    words = ntuple(j -> j == k ? (w[j] & ~(UInt64(0xff) << sh)) | (UInt64(c) << sh) : w[j], Val(8))
    return FixedString(words, b.len)
end
@inline _bresize(b::Vector{UInt8}, n::Int) = resize!(b, n)
@inline _bresize(b::FixedString, n::Int) = FixedString(b.words, n)
@inline _bpush(b::Vector{UInt8}, c::UInt8) = push!(b, c)
@inline function _bpush(b::FixedString, c::UInt8)
    n = b.len + 1
    n <= 64 || _fixedstring_overflow(n)
    return _bset(FixedString(b.words, n), n, c)
end

"Copy the `n` bytes at `src` to `dst` in the same buffer; the ranges may overlap."
function _bmove(b, dst::Int, src::Int, n::Int)
    if dst < src
        for k in 0:(n - 1)
            b = _bset(b, dst + k, _bget(b, src + k))
        end
    elseif dst > src
        for k in (n - 1):-1:0
            b = _bset(b, dst + k, _bget(b, src + k))
        end
    end
    return b
end

"Delete the bytes `from:to` of `b` (nothing if `to < from`)."
function _bdelete(b, from::Int, to::Int)
    to < from && return b
    n = _blen(b)
    b = _bmove(b, from, to + 1, n - to)
    return _bresize(b, n - (to - from + 1))
end

"Insert the byte `c` before position `i` of `b`."
function _binsert(b, i::Int, c::UInt8)
    n = _blen(b)
    b = _bresize(b, n + 1)
    b = _bmove(b, i + 1, i, n - i + 1)
    return _bset(b, i, c)
end

"The first position from `i` on that holds no blank; `_blen(b) + 1` if there is none."
function _skip_blanks(b, i::Int)
    while i <= _blen(b) && _bget(b, i) == UInt8(' ')
        i += 1
    end
    return i
end

"Append the byte `c` `n` times to `b` (nothing for `n <= 0`)."
@inline function _putn(b, c::UInt8, n::Int)
    for _ in 1:n
        b = _bpush(b, c)
    end
    return b
end

"Append the bytes of `s`."
@inline _puts(b::Vector{UInt8}, s::AbstractString) = append!(b, codeunits(s))
@inline function _puts(b::FixedString, s::AbstractString)
    for c in codeunits(s)
        b = _bpush(b, c)
    end
    return b
end

"Append the decimal digits of `n >= 0`, zero padded on the left to `width` (`lpad(string(n), width, '0')`)."
function _putint(b, n::Integer, width::Int = 1)
    nd = ndigits(n)
    b = _putn(b, UInt8('0'), width - nd)
    base = _blen(b)
    b = _bresize(b, base + nd)
    # Loop invariant: the `k` last digits of `n` are written.
    for k in 0:(nd - 1)
        n, d = divrem(n, 10)
        b = _bset(b, base + nd - k, UInt8('0') + UInt8(d))
    end
    return b
end
