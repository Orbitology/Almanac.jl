# Time strings to ET: a close port of the CSPICE parser. `zztime.c` gives the token
# list and its editing operations, `tpartv.c` the rules that give a meaning to every
# token, `tcheck.c` the range checks, and `str2et.c`, `utc2et.c`, `tparse.c` and
# `tpictr.c` the public calls. The acceptance rules are those of the C source, step by
# step, so the code keeps its structure (and its fixed sizes: 64 tokens, a picture of
# 320 characters, a 12-character pattern key).

# ---------------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------------

"Recognized substrings of `zztokns`: (text, number of characters that must match, class). The order inside one first letter is the search order."
const _TK_RECOG = (
    (" ", 1, 'b'), ("(", 1, '['), (")", 1, ']'), (",", 1, ','), ("-", 1, '-'), (".", 1, '.'),
    ("//", 2, 'd'), ("/", 1, '/'), ("::", 2, 'd'), (":", 1, ':'),
    ("A.D.", 4, 'e'), ("AD", 2, 'e'), ("A.M.", 4, 'N'), ("AM", 2, 'N'), ("APRIL", 3, 'm'),
    ("AUGUST", 3, 'm'),
    ("B.C.", 4, 'e'), ("BC", 2, 'e'),
    ("CDT", 3, 'Z'), ("CST", 3, 'Z'),
    ("DECEMBER", 3, 'm'), ("D+", 2, 'E'), ("D-", 2, 'E'), ("D", 1, 'E'),
    ("EDT", 3, 'Z'), ("EST", 3, 'Z'), ("E+", 2, 'E'), ("E-", 2, 'E'), ("E", 1, 'E'),
    ("FEBRUARY", 3, 'm'), ("FRIDAY", 3, 'w'),
    ("JANUARY", 3, 'm'), ("JD", 2, 'j'), ("JULY", 3, 'm'), ("JUNE", 3, 'm'),
    ("MARCH", 3, 'm'), ("MAY", 3, 'm'), ("MDT", 3, 'Z'), ("MONDAY", 3, 'w'), ("MST", 3, 'Z'),
    ("NOVEMBER", 3, 'm'),
    ("OCTOBER", 3, 'm'),
    ("P.M.", 4, 'N'), ("PDT", 3, 'Z'), ("PM", 2, 'N'), ("PST", 3, 'Z'),
    ("SATURDAY", 3, 'w'), ("SEPTEMBER", 3, 'm'), ("SUNDAY", 3, 'w'),
    ("TDB", 3, 's'), ("TDT", 3, 's'), ("THURSDAY", 3, 'w'), ("TUESDAY", 3, 'w'), ("T", 1, 't'),
    ("UTC+", 4, 'O'), ("UTC-", 4, 'o'), ("UTC", 3, 's'),
    ("'", 1, '\''),
    ("WEDNESDAY", 3, 'w'),
    ("Z", 1, 'x'),
)

"Names of the token classes for the messages of `zzunpck`."
const _TK_NAMES = Dict{Char, String}(
    '\'' => "\"Year Abbreviation Mark\"", ',' => "comma",
    '-' => "dash", '.' => "period", '/' => "slash", ':' => "colon", 'D' => "Day of Month",
    'H' => "Hour", 'M' => "Minute", 'N' => "AM/PM indicator", 'O' => "UTC-Offset indicator",
    'S' => "Second", 'Y' => "Year", 'Z' => "Time-Zone indicator", '[' => "Left Parenthesis",
    ']' => "Right Parenthesis", 'b' => "White Space", 'd' => "Day-of-Year indicator",
    'e' => "Era", 'i' => "Integer", 'j' => "Julian Date indicator", 'm' => "Month",
    'n' => "Decimal Number", 'o' => "UTC-Offset indicator",
    's' => "Time System specification", 't' => "ISO Time Separator", 'w' => "Weekday",
    'y' => "Day of Year"
)

"Built-in token patterns of `zztpats.c` and their meanings, in the order of the source."
const _TK_PATTERNS = (
    ("Y-i-it", "Y*m*D*"), ("Y-i-iti:i", "Y*m*D*H*M"), ("Y-i-iti:i:i", "Y*m*D*H*M*S"),
    ("Y-i-iti:i:n", "Y*m*D*H*M*S"), ("Y-i-iti:n", "Y*m*D*H*M"), ("Y-i/", "Y*y*"),
    ("Y-i/i:i", "Y*y*H*M"), ("Y-i/i:i:i", "Y*y*H*M*S"), ("Y-i/i:i:n", "Y*y*H*M*S"),
    ("Y-i/i:n", "Y*y*H*M"), ("Y-id", "Y*y*"), ("Y-idi:i", "Y*y*H*M"),
    ("Y-idi:i:i", "Y*y*H*M*S"), ("Y-idi:i:n", "Y*y*H*M*S"), ("Y-idi:n", "Y*y*H*M"),
    ("Y-it", "Y*y*"), ("Y-iti:i", "Y*y*H*M"), ("Y-iti:i:i", "Y*y*H*M*S"),
    ("Y-iti:i:n", "Y*y*H*M*S"), ("Y-iti:n", "Y*y*H*M"), ("Yid", "Yy*"), ("Yidi:i", "Yy*H*M"),
    ("Yidi:i:i", "Yy*H*M*S"), ("Yidi:i:n", "Yy*H*M*S"), ("Yidi:n", "Yy*H*M"), ("Yii", "YmD"),
    ("Yiii", "YmDH"), ("Yiii:i", "YmDH*M"), ("Yiii:i:i", "YmDH*M*S"),
    ("Yiii:i:n", "YmDH*M*S"), ("Yiii:n", "YmDH*M"), ("Yiiii", "YmDHM"), ("Yiiiii", "YmDHMS"),
    ("Yiiiin", "YmDHMS"), ("Yiiin", "YmDHM"), ("Yiin", "YmDH"), ("Yim", "YDm"),
    ("Yimi", "YDmH"), ("Yimi:i", "YDmH*M"), ("Yimi:i:i", "YDmH*M*S"),
    ("Yimi:i:n", "YDmH*M*S"), ("Yimi:n", "YDmH*M"), ("Yimn", "YDmH"), ("Yin", "YmD"),
    ("Ymi", "YmD"), ("Ymii", "YmDH"), ("Ymii:i", "YmDH*M"), ("Ymii:i:i", "YmDH*M*S"),
    ("Ymii:i:n", "YmDH*M*S"), ("Ymii:n", "YmDH*M"), ("Ymin", "YmDH"), ("Ymn", "YmD"),
    ("Ynm", "YDm"), ("i-Y/", "y*Y*"), ("i-Y/i:i", "y*Y*H*M"), ("i-Y/i:i:i", "y*Y*H*M*S"),
    ("i-Y/i:i:n", "y*Y*H*M*S"), ("i-Y/i:n", "y*Y*H*M"), ("i-Yd", "y*Y*"),
    ("i-Ydi:i", "y*Y*H*M"), ("i-Ydi:i:i", "y*Y*H*M*S"), ("i-Ydi:i:n", "y*Y*H*M*S"),
    ("i-Ydi:n", "y*Y*H*M"), ("i-i-Y", "m*D*Y"), ("i-i-Yi:i", "m*D*YH*M"),
    ("i-i-Yi:i:i", "m*D*YH*M*S"), ("i-i-Yi:i:n", "m*D*YH*M*S"), ("i-i-Yi:n", "m*D*YH*M"),
    ("i-i-it", "Y*m*D*"), ("i-i-iti:i", "Y*m*D*H*M"), ("i-i-iti:i:i", "Y*m*D*H*M*S"),
    ("i-i-iti:i:n", "Y*m*D*H*M*S"), ("i-i-iti:n", "Y*m*D*H*M"), ("i-i/i:i", "Y*y*H*M"),
    ("i-i/i:i:i", "Y*y*H*M*S"), ("i-i/i:i:n", "Y*y*H*M*S"), ("i-i/i:n", "Y*y*H*M"),
    ("i-idi:i", "Y*y*H*M"), ("i-idi:i:i", "Y*y*H*M*S"), ("i-idi:i:n", "Y*y*H*M*S"),
    ("i-idi:n", "Y*y*H*M"), ("i-it", "Y*y*"), ("i-iti:i", "Y*y*H*M"),
    ("i-iti:i:i", "Y*y*H*M*S"), ("i-iti:i:n", "Y*y*H*M*S"), ("i-iti:n", "Y*y*H*M"),
    ("i/i/Y", "m*D*Y"), ("i/i/Y/i:n", "m*D*Y*H*M"), ("i/i/Yi:i", "m*D*YH*M"),
    ("i/i/Yi:i:i", "m*D*YH*M*S"), ("i/i/Yi:i:n", "m*D*YH*M*S"), ("i/i/i", "m*D*Y"),
    ("i/i/ii:i", "m*D*YH*M"), ("i/i/ii:i:i", "m*D*YH*M*S"), ("i/i/ii:i:n", "m*D*YH*M*S"),
    ("i/i/ii:n", "m*D*YH*M"), ("i/i/ii:n", "m*D*YH*M"), ("i:i:ii-i-Y", "H*M*Sm*D*Y"),
    ("i:i:ii/i/Y", "H*M*Sm*D*Y"), ("i:i:ii/i/i", "H*M*Sm*D*Y"), ("i:i:iimY", "H*M*SDmY"),
    ("i:i:imiY", "H*M*SmDY"), ("i:i:ni-i-Y", "H*M*Sm*D*Y"), ("i:i:ni/i/Y", "H*M*Sm*D*Y"),
    ("i:i:ni/i/i", "H*M*Sm*D*Y"), ("i:i:nimY", "H*M*SDmY"), ("i:i:nmiY", "H*M*SmDY"),
    ("i:ii-i-Y", "H*Mm*D*Y"), ("i:ii/i/Y", "H*Mm*D*Y"), ("i:ii/i/i", "H*Mm*D*Y"),
    ("i:iimY", "H*MDmY"), ("i:imiY", "H*MmDY"), ("i:ni-i-Y", "H*Mm*D*Y"),
    ("i:ni/i/Y", "H*Mm*D*Y"), ("i:ni/i/i", "H*Mm*D*Y"), ("i:nimY", "H*MDmY"),
    ("i:nmiY", "H*MmDY"), ("iYd", "yY*"), ("iYdi:i", "yY*H*M"), ("iYdi:i:i", "yY*H*M*S"),
    ("iYdi:i:n", "yY*H*M*S"), ("iYdi:n", "yY*H*M"), ("iiY", "mDY"), ("iiYi", "mDYH"),
    ("iiYi:i", "mDYH*M"), ("iiYi:i:i", "mDYH*M*S"), ("iiYi:i:n", "mDYH*M*S"),
    ("iiYi:n", "mDYH*M"), ("iiYn", "mDYH"), ("iid", "Yy*"), ("iidi:i", "Yy*H*M"),
    ("iidi:i:i", "Yy*H*M*S"), ("iidi:i:n", "Yy*H*M*S"), ("iidi:n", "Yy*H*M"), ("iim", "YDm"),
    ("iimi", "YDmH"), ("iimi:i", "YDmH*M"), ("iimi:i:i", "YDmH*M*S"),
    ("iimi:i:n", "YDmH*M*S"), ("iimi:n", "YDmH*M"), ("iimii", "YDmHM"), ("iimiii", "YDmHMS"),
    ("iimiin", "YDmHMS"), ("iimin", "YDmHM"), ("iimn", "YDmH"), ("imY", "DmY"),
    ("imYi", "DmYH"), ("imYi:i", "DmYH*M"), ("imYi:i:i", "DmYH*M*S"),
    ("imYi:i:n", "DmYH*M*S"), ("imYi:n", "DmYH*M"), ("imYn", "DmYH"), ("imi", "YmD"),
    ("imi:i:iY", "DmH*M*SY"), ("imi:i:nY", "DmH*M*SY"), ("imi:iY", "DmH*MY"),
    ("imi:nY", "DmH*MY"), ("imii", "YmDH"), ("imii:i", "YmDH*M"), ("imii:i:i", "YmDH*M*S"),
    ("imii:i:n", "YmDH*M*S"), ("imii:n", "YmDH*M"), ("imiii", "YmDHM"), ("imiiii", "YmDHMS"),
    ("imiiin", "YmDHMS"), ("imiin", "YmDHM"), ("imin", "YmDH"), ("imn", "YmD"),
    ("inY", "mDY"), ("inm", "YDm"), ("miY", "mDY"), ("miYi", "mDYH"), ("miYi:i", "mDYH*M"),
    ("miYi:i:i", "mDYH*M*S"), ("miYi:i:n", "mDYH*M*S"), ("miYi:n", "mDYH*M"),
    ("miYn", "mDYH"), ("mii", "mDY"), ("mii:i:iY", "mDH*M*SY"), ("mii:i:nY", "mDH*M*SY"),
    ("mii:iY", "mDH*MY"), ("mii:nY", "mDH*MY"), ("miii", "mDYH"), ("miii:i", "mDYH*M"),
    ("miii:i:i", "mDYH*M*S"), ("miii:i:n", "mDYH*M*S"), ("miii:n", "mDYH*M"),
    ("miiii", "mDYHM"), ("miiiii", "mDYHMS"), ("miiiin", "mDYHMS"), ("miiin", "mDYHM"),
    ("miin", "mDYH"), ("mnY", "mDY"), ("mni", "mDY"), ("nmY", "DmY"),
    ("Y-i-itx", "Y*m*D**"), ("Y-i-iti:ix", "Y*m*D*H*M*"), ("Y-i-iti:i:ix", "Y*m*D*H*M*S*"),
    ("Y-i-iti:i:nx", "Y*m*D*H*M*S*"), ("Y-i-iti:nx", "Y*m*D*H*M*"), ("Y-itx", "Y*y**"),
    ("Y-iti:ix", "Y*y*H*M*"), ("Y-iti:i:ix", "Y*y*H*M*S*"), ("Y-iti:i:nx", "Y*y*H*M*S*"),
    ("Y-iti:nx", "Y*y*H*M*"), ("i-i-itx", "Y*m*D**"), ("i-i-iti:ix", "Y*m*D*H*M*"),
    ("i-i-iti:i:ix", "Y*m*D*H*M*S*"), ("i-i-iti:i:nx", "Y*m*D*H*M*S*"),
    ("i-i-iti:nx", "Y*m*D*H*M*"), ("i-itx", "Y*y**"), ("i-iti:ix", "Y*y*H*M*"),
    ("i-iti:i:ix", "Y*y*H*M*S*"), ("i-iti:i:nx", "Y*y*H*M*S*"), ("i-iti:nx", "Y*y*H*M*"),
    ("Y-i-iti", "Y*m*D*H"), ("Y-i-itn", "Y*m*D*H"), ("Y-iti", "Y*y*H"), ("Y-itn", "Y*y*H"),
    ("i-i-iti", "Y*m*D*H"), ("i-i-itn", "Y*m*D*H"), ("i-iti", "Y*y*H"), ("i-itn", "Y*y*H"),
    ("Y-i-itix", "Y*m*D*H*"), ("Y-i-itnx", "Y*m*D*H*"), ("Y-itix", "Y*y*H*"),
    ("Y-itnx", "Y*y*H*"), ("i-i-itix", "Y*m*D*H*"), ("i-i-itnx", "Y*m*D*H*"),
    ("i-itix", "Y*y*H*"), ("i-itnx", "Y*y*H*"),
)

const _TK_ZONES = ("EST", "EDT", "CST", "CDT", "MST", "MDT", "PST", "PDT")
const _TK_OFFSETS = ("UTC-5", "UTC-4", "UTC-6", "UTC-5", "UTC-7", "UTC-6", "UTC-8", "UTC-7")

const _MONTH_NAMES = (
    "January", "February", "March", "April", "May", "June", "July", "August",
    "September", "October", "November", "December",
)

# Indices of the five modifiers of a time string.
const _MOD_ERA = 1
const _MOD_WDAY = 2
const _MOD_ZONE = 3
const _MOD_AMPM = 4
const _MOD_SYSTEM = 5

const _BLANK = UInt8(' ')

# ---------------------------------------------------------------------------------
# Token list (ZZTIME)
# ---------------------------------------------------------------------------------

const _MAXTOK = 64
const _MAXPIC = 320

"""
    _Tokens

Token list of one time string, as the saved variables of `zztime.c` hold it: `str` is
the string without trailing blanks, `rep[1:size]` has one class letter for each token
and is blank after `size`, `begs`/`ends` are the positions of the tokens in `str` and
`pbegs`/`pends` in the picture `pictur` (320 characters, blank padded). `meta[1]` is
`size`, and `meta[2]` is 1 when the string has an A.M./P.M. mark.

The buffers are `_StackVec`s, and every function that takes a `_Tokens` is inlined
into `_tpartv`, so that the list stays on the stack.
"""
struct _Tokens{S <: AbstractVector{UInt8}}
    str::S
    rep::_StackVec{_MAXTOK, UInt8}
    begs::_StackVec{_MAXTOK, Int}
    ends::_StackVec{_MAXTOK, Int}
    pbegs::_StackVec{_MAXTOK, Int}
    pends::_StackVec{_MAXTOK, Int}
    pictur::_StackVec{_MAXPIC, UInt8}
    meta::_StackVec{2, Int}
end

"An empty token list of the bytes `str` (`size` 0, `rep` and `pictur` blank)."
@inline _Tokens(str::AbstractVector{UInt8}) = _Tokens(
    str, _StackVec{_MAXTOK, UInt8}(_BLANK), _StackVec{_MAXTOK, Int}(0),
    _StackVec{_MAXTOK, Int}(0), _StackVec{_MAXTOK, Int}(0), _StackVec{_MAXTOK, Int}(0),
    _StackVec{_MAXPIC, UInt8}(_BLANK), _StackVec{2, Int}(0)
)

@inline _ntok(st::_Tokens) = st.meta[1]
@inline _setntok!(st::_Tokens, n::Int) = (st.meta[1] = n; nothing)
@inline _ampm(st::_Tokens) = st.meta[2] == 1

"Return `s[b:e]` as a `String`; empty when the range is empty or outside `s`."
_sub(s::AbstractVector{UInt8}, b::Int, e::Int) =
    (b >= 1 && e <= length(s) && b <= e) ? String(s[b:e]) : ""

"Return the text of `s` with `<` before position `b` and `>` after position `e`, the marking of the CSPICE messages. Requires `1 <= b`."
function _marked(s::AbstractVector{UInt8}, b::Int, e::Int)
    b = clamp(b, 1, length(s) + 1)
    e = clamp(e, b - 1, length(s))
    return string(String(s[1:(b - 1)]), '<', String(s[b:e]), '>', String(s[(e + 1):end]))
end

"Write `text` into `dest[from:to]` of the blank-padded buffer `dest`, cut or blank padded to that length (the Fortran substring assignment)."
@inline function _assign!(dest, from::Int, to::Int, text)
    t = text isa AbstractString ? codeunits(text) : text
    # Loop invariant: `dest[from:i-1]` holds the first characters of `text`, blank padded.
    for i in from:to
        k = i - from + 1
        dest[i] = k <= length(t) ? t[k] : _BLANK
    end
    return dest
end

"""
    _repsub!(buf, left, right, text, ndec=-1) -> buf

Replace `buf[left:right]` of the fixed-length buffer `buf` by `text`, in
place, cut or blank padded to the length of `buf` (`zzrepsub`). With `ndec >= 0`,
`text` is followed by a decimal point and `ndec` characters `#` (the mark of
`zzmkpc`). `text` is a `String` or a tuple of bytes.
"""
@inline function _repsub!(buf, left::Int, right::Int, text, ndec::Int = -1)
    t = text isa AbstractString ? codeunits(text) : text
    n = length(buf)
    head = clamp(left - 1, 0, n)                 # bytes kept before the text
    tail0 = clamp(right, 0, n) + 1               # first byte kept after the text
    tlen = length(t) + (ndec >= 0 ? 1 + ndec : 0)
    dest = head + tlen + 1                       # where the kept tail goes
    ntail = n - tail0 + 1
    # The tail moves first: the text may overlap where it was. The copy runs from the
    # end that it does not overwrite before it reads.
    if dest <= n && ntail > 0
        m = min(ntail, n - dest + 1)
        if dest > tail0
            for k in (m - 1):-1:0
                buf[dest + k] = buf[tail0 + k]
            end
        elseif dest < tail0
            for k in 0:(m - 1)
                buf[dest + k] = buf[tail0 + k]
            end
        end
    end
    for k in 1:min(tlen, n - head)
        buf[head + k] = k <= length(t) ? t[k] : (k == length(t) + 1 ? UInt8('.') : UInt8('#'))
    end
    for k in (head + tlen + ntail + 1):n
        buf[k] = _BLANK
    end
    return buf
end

"Return whether `rep[b:b+length(p)-1]` equals the bytes `p`."
@inline function _matches(rep, b::Int, p)
    for k in 1:length(p)
        rep[b + k - 1] == p[k] || return false
    end
    return true
end

"Return the first position `b` in `1:n` at which `pat` starts inside `rep[1:n]`, or 0 (`pos`)."
@inline function _pos(rep, n::Int, pat)::Int
    p = codeunits(pat)
    for b in 1:(n - length(p) + 1)
        _matches(rep, b, p) && return b
    end
    return 0
end

"Return the last position `b <= start` at which `pat` starts inside `rep[1:n]`, or 0 (`posr`)."
@inline function _posr(rep, n::Int, start::Int, pat)::Int
    p = codeunits(pat)
    for b in min(n - length(p) + 1, start):-1:1
        _matches(rep, b, p) && return b
    end
    return 0
end

"Move the position records of tokens `get:size` down to `put...` (the common loop of `zzcmbt` and `zznote`). Requires `get > put`."
@inline function _shift!(st::_Tokens, put::Int, get::Int)
    # Loop invariant: records `get0:g-1` are at `put0:p-1`.
    for g in get:_ntok(st)
        st.begs[put] = st.begs[g]
        st.ends[put] = st.ends[g]
        st.pbegs[put] = st.pbegs[g]
        st.pends[put] = st.pends[g]
        put += 1
    end
    return nothing
end

"Combine the first (`l2r`) or last occurrence of the token run `pat` into the one token `letter`; return whether it did (`zzcmbt`)."
@inline function _zzcmbt!(st::_Tokens, pat::String, letter::Char, l2r::Bool)::Bool
    n = ncodeunits(pat)
    size = _ntok(st)
    from = l2r ? _pos(st.rep, size, pat) : _posr(st.rep, size, size, pat)
    from > 0 || return false
    to = from + n - 1
    st.ends[from] = st.ends[to]
    st.pends[from] = st.pends[to]
    _repsub!(st.rep, from, to, (UInt8(letter),))
    _shift!(st, from + 1, to + 1)
    _setntok!(st, size - (n - 1))
    return true
end

"Return the 12-character key of the token list (`zzgrep` into the 12-character variable of `tpartv`; longer lists are cut) as an integer, the first character in the high byte."
@inline function _zzgrep_key(st::_Tokens)::UInt128
    key = zero(UInt128)
    for i in 1:12
        key = (key << 8) | st.rep[i]
    end
    return key
end

"The key of `_zzgrep_key` for the text `s` of 12 characters."
function _tk_key(s::String)::UInt128
    @assert ncodeunits(s) == 12
    key = zero(UInt128)
    for c in codeunits(s)
        key = (key << 8) | c
    end
    return key
end

"Pattern key (12 characters, blank padded, as `_tk_key` gives it) to meaning (12 characters, blank padded)."
const _TK_KNOWN = Dict{UInt128, String}(_tk_key(rpad(k, 12)) => rpad(m, 12) for (k, m) in _TK_PATTERNS)

"Return whether the token list has a token of class `letter` (`zzist`)."
@inline function _zzist(st::_Tokens, letter::Char)
    c = UInt8(letter)
    for i in 1:_ntok(st)
        st.rep[i] == c && return true
    end
    return false
end

"Return whether the byte `c` is one of the bytes of `chars`."
@inline function _isin(c::UInt8, chars)
    for k in 1:length(chars)
        chars[k] == c && return true
    end
    return false
end

"Return (found, b, e): whether two tokens in a row both have a class in `chars`, and the range of the pair in the string (`zzispt`)."
@inline function _zzispt(st::_Tokens, chars::String)
    c = codeunits(chars)
    # Loop invariant: no pair starts before `from`.
    from = 1
    while from <= _MAXTOK
        if _isin(st.rep[from], c)
            from < _ntok(st) || return false, 0, 0
            to = from + 1
            _isin(st.rep[to], c) && return true, st.begs[from], st.ends[to]
            from = to
        else
            from += 1
        end
    end
    return false, 0, 0
end

"Remove the first token of class `letter`; return (found, b, e) with its range in the string, or (false, 0, 0) (`zznote`)."
@inline function _zznote!(st::_Tokens, letter::Char)
    c = UInt8(letter)
    put = 0
    for i in 1:_MAXTOK
        if st.rep[i] == c
            put = i
            break
        end
    end
    put == 0 && return false, 0, 0
    b, e = st.begs[put], st.ends[put]
    _shift!(st, put, put + 1)
    size = _ntok(st)
    # Loop invariant: `rep[put0:i-1]` holds the letters that were at `put0+1:i`.
    for i in put:(size - 1)
        st.rep[i] = st.rep[i + 1]
    end
    st.rep[size] = _BLANK
    _setntok!(st, size - 1)
    return true, b, e
end

"Remove every token of class `letter`; return whether one was removed (`zzremt`)."
@inline function _zzremt!(st::_Tokens, letter::Char)::Bool
    put = 0
    did = false
    c = UInt8(letter)
    # Loop invariant: tokens `1:i-1` that are not of class `letter` are at `1:put`.
    for i in 1:_ntok(st)
        if st.rep[i] != c
            put += 1
            st.rep[put] = st.rep[i]
            st.begs[put] = st.begs[i]
            st.ends[put] = st.ends[i]
            st.pbegs[put] = st.pbegs[i]
            st.pends[put] = st.pends[i]
        else
            did = true
        end
    end
    _setntok!(st, put)
    for i in (put + 1):_MAXTOK
        st.rep[i] = _BLANK
    end
    return did
end

"""
    _zzsubt!(st, pat, transl, l2r) -> Bool

Replace the letters of the first (`l2r`) or last occurrence of `pat` by
`transl` and returns whether it did (`zzsubt`). A `pat` that starts with `<` must
match at the start of the list and one that ends with `>` at its end; in these two
forms the replaced range has the length of the pattern without the mark, and a shorter
`transl` is blank padded, as the Fortran assignment does.
"""
@inline function _zzsubt!(st::_Tokens, pat::String, transl::String, l2r::Bool)::Bool
    p = codeunits(pat)
    k = length(p)
    size = _ntok(st)
    if p[1] == UInt8('<') && k > 1
        to = min(k - 1, size)
        to == k - 1 || return false
        for i in 1:to
            st.rep[i] == p[i + 1] || return false
        end
        _assign!(st.rep, 1, to, transl)
        return true
    elseif p[k] == UInt8('>') && k > 1
        from = max(1, size - k + 2)
        to = size
        to - from + 1 == k - 1 || return false
        _matches(st.rep, from, @view p[1:(k - 1)]) || return false
        _assign!(st.rep, from, to, transl)
        return true
    end
    from = l2r ? _pos(st.rep, _MAXTOK, pat) : _posr(st.rep, _MAXTOK, size, pat)
    from > 0 || return false
    _assign!(st.rep, from, min(_MAXTOK, from + ncodeunits(transl) - 1), transl)
    return true
end

"Apply `_zzsubt!`, then remove every `*` token; return whether the substitution was made (`zzrept`)."
@inline function _zzrept!(st::_Tokens, pat::String, transl::String, l2r::Bool)::Bool
    did = _zzsubt!(st, pat, transl, l2r)
    _zzremt!(st, '*')
    return did
end

"""
    _zzknown!(st) -> Bool

Replace the classes of the token list by the meaning of its built-in pattern, if its
12-character key is one, and remove the `*` tokens (`zzrept` with the key of `zzgrep`
and the meaning of `zztpats`); return whether the key is a pattern. The key matches at
the start of the list, where `zzsubt` finds it.
"""
@inline function _zzknown!(st::_Tokens)::Bool
    meaning = get(_TK_KNOWN, _zzgrep_key(st), nothing)
    meaning === nothing && return false
    _assign!(st.rep, 1, 12, meaning)
    _zzremt!(st, '*')
    return true
end

"Give class `letter` to every integer token whose value is in `b..e`; return whether one changed (`zzvalt`)."
@inline function _zzvalt!(st::_Tokens, b::Int, e::Int, letter::Char)::Bool
    did = false
    for i in 1:_ntok(st)
        st.rep[i] == UInt8('i') || continue
        # `nparsi` rejects a value that is not a 32-bit integer.
        x = _try_nparsi(@view st.str[st.begs[i]:st.ends[i]])
        x === nothing && continue
        if b <= x <= e
            st.rep[i] = UInt8(letter)
            did = true
        end
    end
    return did
end

"Return whether characters `from:to` of `s` exist and equal `text[1:to-from+1]` when case is ignored (`samsbi`)."
@inline function _samsbi(s::AbstractVector{UInt8}, from::Int, to::Int, text::String)::Bool
    to <= length(s) || return false
    t = codeunits(text)
    for i in from:to
        _upper(s[i]) == t[i - from + 1] || return false
    end
    return true
end

@noinline _toomany(str, s, b, e) = "The input time string '$str' cannot be processed " *
    "because it contains more than $_MAXTOK recognizable tokens. The token that " *
    "could not be processed was '$(_sub(s, b, e))'."
@noinline _toolong(str, s, b, e) = "The input time string '$str' cannot be processed " *
    "because the internal picture describing it requires more than $_MAXPIC " *
    "characters. The token that could not be processed was '$(_sub(s, b, e))'."

"""
    _zztokns!(st, str) -> String

Fill the empty token list `st` of the bytes `st.str` of `str` and return `""`, or
return a message that says which part of the string is not a recognized token
(`zztokns`). Integers, runs of blanks, tabs, punctuation, month and weekday names
(three letters or more of the name), eras, A.M./P.M., time zones, `UTC+`/`UTC-`,
time systems, `JD`, `T` and `Z` are tokens; the picture of the string is built at the
same time.
"""
@inline function _zztokns!(st::_Tokens, str::AbstractString)
    bytes = st.str
    last = length(bytes)
    next = 1
    pnext = 1
    put = 0
    pto = 0
    # Loop invariant: `bytes[1:next-1]` is tokenized into `put` tokens and its picture
    # is `pictur[1:pnext-1]`.
    while next <= last
        from = next
        pfrom = pnext
        item = bytes[next]
        to = from - 1
        while to < last && _isdigit(bytes[to + 1])
            to += 1
        end
        nchar = to - from + 1
        if nchar > 0
            put += 1
            put > _MAXTOK && return _toomany(str, bytes, from, to)
            st.rep[put] = UInt8('i')
            next = to + 1
            pto = pfrom + nchar - 1
            pto > _MAXPIC && return _toolong(str, bytes, from, to)
            pnext = pto + 1
            _assign!(st.pictur, pfrom, pto, @view bytes[from:to])
        elseif item == _BLANK
            put += 1
            put > _MAXTOK && return _toomany(str, bytes, from, from)
            to = from
            while to < last && bytes[to + 1] == _BLANK
                to += 1
            end
            st.rep[put] = UInt8('b')
            next = to + 1
            pto = pfrom + to - from
            pto > _MAXPIC && return _toolong(str, bytes, from, from)
            pnext = pto + 1
            _assign!(st.pictur, pfrom, pto, "")
        elseif item == 0x09
            put += 1
            put > _MAXTOK && return _toomany(str, bytes, from, from)
            st.rep[put] = UInt8('b')
            to = from
            next += 1
            pto = pfrom
            pto > _MAXPIC && return _toolong(str, bytes, from, from)
            pnext = pto + 1
            st.pictur[pfrom] = _BLANK
        elseif item < 32 || item > 126
            return _nonprinting(bytes, next)
        else
            up = _upper(item)
            found = false
            for (text, w, class) in _TK_RECOG
                codeunit(text, 1) == up || continue
                to = from + w - 1
                _samsbi(bytes, from, to, text) || continue
                if class == 'm' || class == 'w'
                    # A month or a weekday: the token runs as far as the name matches.
                    r = ncodeunits(text)
                    while w < r && to < last && _upper(bytes[to + 1]) == codeunit(text, w + 1)
                        w += 1
                        to += 1
                    end
                    kind = w >= r ? 1 : 2
                    # The letter case: lower, as the name is listed (its first three
                    # letters in upper case) or other.
                    case = up != item ? 1 : (
                            next + 2 <= last && bytes[next + 1] == codeunit(text, 2) &&
                            bytes[next + 2] == codeunit(text, 3) ? 2 : 3
                        )
                    mark = class == 'm' ?
                        (("month", "MONTH", "Month"), ("mon", "MON", "Mon"))[kind][case] :
                        (("weekday", "WEEKDAY", "Weekday"), ("wkd", "WKD", "Wkd"))[kind][case]
                    pto = pfrom + ncodeunits(mark) - 1
                    pto > _MAXPIC && return _toolong(str, bytes, from, to)
                    pnext = pto + 1
                    _assign!(st.pictur, pfrom, pto, mark)
                elseif class == 'e'
                    pto = pfrom + 2
                    pto > _MAXPIC && return _toolong(str, bytes, from, to)
                    pnext = pto + 1
                    _assign!(st.pictur, pfrom, pto, item == up ? "ERA" : "era")
                elseif class == 'N'
                    pto = pfrom + 3
                    pto > _MAXPIC && return _toolong(str, bytes, from, to)
                    pnext = pto + 1
                    _assign!(st.pictur, pfrom, pto, item == up ? "AMPM" : "ampm")
                    st.meta[2] = 1
                elseif class == 'x'
                    # Recognized and ignored: the picture gets nothing.
                else
                    pto = pfrom + to - from
                    pto > _MAXPIC && return _toolong(str, bytes, from, to)
                    pnext = pto + 1
                    _assign!(st.pictur, pfrom, pto, @view bytes[from:to])
                end
                put += 1
                put > _MAXTOK && return _toomany(str, bytes, from, to)
                st.rep[put] = UInt8(class)
                next = to + 1
                found = true
                break
            end
            found || return _unrecognized(bytes, from)
        end
        st.begs[put] = from
        st.ends[put] = to
        st.pbegs[put] = pfrom
        st.pends[put] = pto
    end
    _setntok!(st, put)
    return ""
end

@noinline _nonprinting(bytes, next) = "There is a non-printing, non-tab character (ASCII " *
    "$(Int(bytes[next]))) at position $next of the time string: " *
    String(bytes[1:(next - 1)]) * "<>" * String(bytes[(next + 1):end])
@noinline _unrecognized(bytes, from) = "The input string contains an unrecognizable " *
    "substring beginning at the character marked by <$(Char(bytes[from]))>: \"" *
    _marked(bytes, from, from) * "\""

"Return the number of `#` that `zzmkpc` writes after the picture mark of the token text `token`: one for each decimal (14 at most), or -1 if the token has no decimal point."
function _zzmkpc(token::AbstractVector{UInt8})
    point = findfirst(==(UInt8('.')), token)
    return point === nothing ? -1 : min(14, length(token) - point)
end

"Return the month number of the first three characters of `s[j:k]` (blank padded, any case) as `_MONTH_ABBR` lists them, or `nothing`."
function _month_abbr_index(s::AbstractVector{UInt8}, j::Int, k::Int)
    c1, c2, c3 = ntuple(n -> j + n - 1 <= min(k, j + 2) ? _upper(s[j + n - 1]) : _BLANK, Val(3))
    for (i, m) in enumerate(_MONTH_ABBR)
        codeunit(m, 1) == c1 && codeunit(m, 2) == c2 && codeunit(m, 3) == c3 && return i
    end
    return nothing
end

"Return the value of the number in `s[j:k]` as `nparsd` computes it, or `old` if the text is not a number (`nparsd` leaves its output unchanged then)."
function _tokval(s::AbstractVector{UInt8}, j::Int, k::Int, old::Float64)::Float64
    x = _try_nparsd(@view s[j:k])
    return x === nothing ? old : x
end

"""
    _zzunpck!(pic, st, yabbrv) -> (ok, tvec, ntvec, type, np, error)

Return the time vector `tvec` (eight components, the absent ones 0), the number of
components `ntvec` that the string gives, the type `"YMD"`, `"YD"` or `"JD"`, and the
length `np` of the `timout` picture of the string, which is written into `pic` (the
picture `st.pictur`, which changes); or `ok == false` and a message that says which token has no meaning or
which component is absent or given twice (`zzunpck`). Every token of `st` is one of the
resolved classes `Y m D y H M S J`. The picture is blank (`np == 0`) when a token has no
meaning, and it is set in every other case.
"""
@inline function _zzunpck!(pic, st::_Tokens, yabbrv::Bool)
    s = st.str
    tvec = zero(SVector{8, Float64})
    h = m = sec = 0.0
    nyear = nmon = nday = nhour = nmin = nsec = ndoy = njd = 0
    ntvec = 0
    # Loop invariant: tokens `i+1:size` are counted, their values are stored and their
    # picture marks are written (from the right, so the positions to the left hold).
    for i in _ntok(st):-1:1
        item = Char(st.rep[i])
        j, k = st.begs[i], st.ends[i]
        tok = @view s[j:k]
        if item == 'Y'
            nyear += 1; ntvec += 1
            tvec = Base.setindex(tvec, _tokval(s, j, k, tvec[1]), 1)
            _repsub!(pic, st.pbegs[i], st.pends[i], yabbrv ? "YR" : "YYYY")
        elseif item == 'm'
            nmon += 1; ntvec += 1
            value = _month_abbr_index(s, j, k)
            if value === nothing
                tvec = Base.setindex(tvec, _tokval(s, j, k, tvec[2]), 2)
                _repsub!(pic, st.pbegs[i], st.pends[i], "MM")
            else
                tvec = Base.setindex(tvec, Float64(value), 2)
            end
        elseif item == 'D'
            nday += 1; ntvec += 1
            tvec = Base.setindex(tvec, _tokval(s, j, k, tvec[3]), 3)
            _repsub!(pic, st.pbegs[i], st.pends[i], "DD", _zzmkpc(tok))
        elseif item == 'y'
            ndoy += 1; ntvec += 1
            tvec = Base.setindex(tvec, _tokval(s, j, k, tvec[2]), 2)
            _repsub!(pic, st.pbegs[i], st.pends[i], "DOY", _zzmkpc(tok))
        elseif item == 'H'
            nhour += 1; ntvec += 1
            h = _tokval(s, j, k, h)
            _repsub!(pic, st.pbegs[i], st.pends[i], _ampm(st) ? "AP" : "HR", _zzmkpc(tok))
        elseif item == 'M'
            nmin += 1; ntvec += 1
            m = _tokval(s, j, k, m)
            _repsub!(pic, st.pbegs[i], st.pends[i], "MN", _zzmkpc(tok))
        elseif item == 'S'
            nsec += 1; ntvec += 1
            sec = _tokval(s, j, k, sec)
            _repsub!(pic, st.pbegs[i], st.pends[i], "SC", _zzmkpc(tok))
        elseif item == 'J'
            njd += 1; ntvec += 1
            # A token such as `- jd 2451545.0` is not a number. CSPICE leaves the value
            # of an earlier call in the vector; NaN marks the date as having no value.
            tvec = Base.setindex(tvec, _tokval(s, j, k, NaN), 1)
            _repsub!(pic, st.pbegs[i], st.pends[i], "JULIAND", _zzmkpc(tok))
        else
            return false, tvec, 0, "", 0, _unpck_token_error(s, item, j, k)
        end
    end
    np = length(pic)
    while np > 0 && pic[np] == _BLANK
        np -= 1
    end
    if nyear == 1 && nmon == 1 && nday == 1 && ndoy == 0 && njd == 0 && nhour <= 1 &&
            nmin <= nhour && nsec <= nmin
        tvec = SVector(tvec[1], tvec[2], tvec[3], h, m, sec, tvec[7], tvec[8])
        return true, tvec, ntvec, "YMD", np, ""
    elseif nyear == 1 && nmon == 0 && nday == 0 && njd == 0 && ndoy == 1 && nhour <= 1 &&
            nmin <= nhour && nsec <= nmin
        tvec = SVector(tvec[1], tvec[2], h, m, sec, tvec[6], tvec[7], tvec[8])
        return true, tvec, ntvec, "YD", np, ""
    elseif nyear == 0 && nmon == 0 && nday == 0 && njd == 1 && ndoy == 0 && nhour <= 0 &&
            nmin <= 0 && nsec <= 0
        return true, tvec, ntvec, "JD", np, ""
    end
    error = _unpck_count_error(s, nyear, nmon, nday, nhour, nmin, nsec, ndoy, njd)
    return false, tvec, 0, "", np, error
end

"The message of `zzunpck` for the token `s[j:k]` of class `item`, which has no meaning."
@noinline function _unpck_token_error(s, item::Char, j::Int, k::Int)
    item == 'i' && return "The meaning of the integer <$(_sub(s, j, k))> " *
        "could not be determined: '" * _marked(s, j, k) * "'"
    item == 'n' && return "The meaning of the decimal number " *
        "<$(_sub(s, j, k))> could not be determined: " * _marked(s, j, k)
    return "An unexpected " * get(_TK_NAMES, item, "substring") * " (\"$(_sub(s, j, k))\") was " *
        "encountered in the time string: " * _marked(s, j, k)
end

"The message of `zzunpck` for the counts of the components that do not make a date."
@noinline function _unpck_count_error(s, nyear, nmon, nday, nhour, nmin, nsec, ndoy, njd)
    str = String(s[1:end])
    twice(what) = "Two substrings indicating $what were identified in the input time " *
        "string \"$str\""
    return if nyear == 0 && nday == 0 && njd == 0 && ndoy == 0 && nhour == 0 && nmin == 0 &&
            nsec == 0
        "No numeric components were supplied in the time string."
    elseif njd == 1
        "The string possesses calendar components in addition to Julian Date specifier."
    elseif njd > 1
        "There is more than one Julian Date specified in the epoch string."
    elseif nyear == 0
        "The year associated with the calendar string \"$str\" could not be identified."
    elseif nyear > 1
        twice("a calendar year")
    elseif nmon > 0 && ndoy > 0
        "Both a day of year and month were identified in the input string. \"$str\""
    elseif nmon > 1
        twice("a calendar month")
    elseif ndoy > 1
        twice("a day of year")
    elseif nday > 1
        twice("a day of month")
    elseif nhour > 1
        twice("an hour of the day")
    elseif nmin > 1
        twice("minutes of the hour")
    elseif nsec > 1
        twice("seconds")
    elseif ndoy == 0 && nmon == 0
        "Neither a month nor day of year could be identified in the input time string: \"$str\""
    elseif nmon == 1 && nday == 0
        "A month was identified in the time string \"$str\", but a day of month could " *
            "not be identified."
    elseif nmon == 0 && nday == 1
        "A day of month was identified in the time string \"$str\", but the month it " *
            "belongs to could not be identified."
    elseif nmin > nhour
        "A minutes components of the time was identified in the time string \"$str\", " *
            "but the hours component could not be identified."
    elseif nsec > nmin
        "A seconds components of the time was identified in the time string \"$str\", " *
            "but the minutes component could not be identified."
    else
        "The time string \"$str\" could not be parsed."
    end
end

# ---------------------------------------------------------------------------------
# TPARTV
# ---------------------------------------------------------------------------------

"""
    _TimeParts

Result of `_tpartv`: `ok` and, if it is false, `error`; the time vector `tvec` with
`ntvec` components given and its `type` (`"YMD"`, `"YD"`, `"JD"`); the five modifiers
`modify` (era `A.D.`/`B.C.`, weekday, zone as `UTC+h:mm`, `A.M.`/`P.M.`, system), `""`
when absent, and `mods` (one is present); `yabbrv` (the year has the form `'yy`); and
the `timout` picture `pictur` of the string (`""` if it has none, or if it was not
asked for).
"""
struct _TimeParts
    ok::Bool
    tvec::SVector{8, Float64}
    ntvec::Int
    type::String
    modify::NTuple{5, String63}
    mods::Bool
    yabbrv::Bool
    pictur::String
    error::String
end

"Return `text` after `pictur` without its trailing blanks and `spaces` blanks (`suffix`)."
function _suffix(pictur::String, text::AbstractString, spaces::Int)
    head = rstrip(pictur, ' ')
    buf = Base.StringVector(ncodeunits(head) + spaces + ncodeunits(text))
    copyto!(buf, codeunits(head))
    fill!(view(buf, (ncodeunits(head) + 1):(ncodeunits(head) + spaces)), _BLANK)
    copyto!(buf, ncodeunits(head) + spaces + 1, codeunits(text), 1, ncodeunits(text))
    return String(buf)
end

"""
    _picture_string(pic, np, needle, modify, zone) -> String

Return the picture `pic[1:np]` as a `String`, with the marks that `tpartv` adds after
it: `::RND` when `needle` occurs in it, then the zone (if `zone`) and the system
modifiers.
"""
@inline function _picture_string(pic, np::Int, needle::String, modify::NTuple{5, String63}, zone::Bool)
    return _picture_suffix(_string_of(pic, np), needle, modify, zone)
end

"Return `pictur` with the marks that `tpartv` adds after it: `::RND` when `needle` occurs in it, then the zone (if `zone`) and the system modifiers."
@noinline function _picture_suffix(pictur::String, needle::String, modify::NTuple{5, String63}, zone::Bool)
    occursin(needle, pictur) && (pictur = _suffix(pictur, "::RND", 1))
    if zone && modify[_MOD_ZONE] != ""
        pictur = _suffix(_suffix(pictur, "::", 1), modify[_MOD_ZONE], 0)
    end
    if modify[_MOD_SYSTEM] != ""
        pictur = _suffix(_suffix(pictur, "::", 1), modify[_MOD_SYSTEM], 0)
    end
    return pictur
end

"Return the text of `s[b:e]`, cut to `modlen` characters, in upper case, as a `String63` (`modlen <= 63`)."
@inline function _upcut(s::AbstractVector{UInt8}, b::Int, e::Int, modlen::Int)::String63
    buf = _ByteBuf()
    if b >= 1 && b <= e && e <= length(s)
        for i in b:min(e, b + modlen - 1)
            buf = _bpush(buf, _upper(s[i]))
        end
    end
    return _inline(buf)
end

"The weekday of the token `s[b:e]` (three letters or more of a weekday name, any case) in three upper case letters."
@inline function _weekday_abbr(s::AbstractVector{UInt8}, b::Int)::String63
    c1, c2 = _upper(s[b]), _upper(s[b + 1])
    for name in ("MON", "TUE", "WED", "THU", "FRI", "SAT", "SUN")
        codeunit(name, 1) == c1 && codeunit(name, 2) == c2 && return String63(name)
    end
    return String63("")
end

"The time system of a token of class `s` (`TDB`, `TDT` or `UTC`, any case), in upper case."
@inline function _system_name(s::AbstractVector{UInt8}, b::Int)::String63
    _upper(s[b]) == UInt8('U') && return String63("UTC")
    return _upper(s[b + 2]) == UInt8('B') ? String63("TDB") : String63("TDT")
end

const _NOMODS = ntuple(_ -> String63(""), Val(5))

"The result of `_tpartv` for a string that does not parse; `modify` holds the modifiers found before the error."
@noinline _fail_parts(error::String, modify::NTuple{5, String63} = _NOMODS) =
    _TimeParts(false, zero(SVector{8, Float64}), 0, "", modify, false, false, "", error)

"""
    _tpartv(str, modlen, wantpic = true) -> _TimeParts

Return the components of the time string, its modifiers and its picture, or `ok ==
false` and a message, with the rules of `tpartv` in their order: combine `UTC+h:mm`,
find one decimal number, handle a Julian date, mark integers of 1000 or more as years,
match the ISO patterns if a `T` is present, extract the modifiers, match the built-in
patterns (again after the removal of commas, dashes and slashes, one kind at a time),
and at last resolve hours, minutes, seconds and the integers next to a month name.
A modifier is cut to `modlen` characters (at most 63), the length of the modifier
variables of the CSPICE caller (16 in `str2et`, 8 in the others). With `wantpic`
false the picture is not made (`pictur == ""`), and a string that parses allocates
nothing.
"""
function _tpartv(str::AbstractString, modlen::Int, wantpic::Bool = true)::_TimeParts
    @assert 1 <= modlen <= 63
    cu = codeunits(str)
    last = length(cu)
    # Loop invariant: `cu[last+1:end]` is blank.
    while last > 0 && cu[last] == _BLANK
        last -= 1
    end
    last == 0 && return _fail_parts("The input time string is blank.")
    s = view(cu, 1:last)
    st = _Tokens(s)
    error = _zztokns!(st, str)
    error == "" || return _fail_parts(error)
    l2r, r2l = true, false
    era = wday = zone = ampm = system = String63("")

    for sign in ("Oi", "oi")
        if _zzcmbt!(st, sign, 'z', l2r)
            _zzcmbt!(st, "z:i", 'Z', l2r)
            _zzsubt!(st, "z", "Z", l2r)
        end
    end
    _zzrept!(st, "m.", "m*", l2r)
    _zzrept!(st, "w.", "w*", l2r)
    _zzrept!(st, "w,", "w*", l2r)
    _zzcmbt!(st, "i.i", 'n', r2l) || _zzcmbt!(st, "i.", 'n', r2l)
    _zzremt!(st, 'b')
    # `zzunpck` writes its marks into the picture of the tokens, which is used once.
    pic = st.pictur

    # Julian date.
    if _zzist(st, 'j')
        _zzrept!(st, "[s]", "*s*", l2r)
        mods, b, e = _zznote!(st, 's')
        mods && (system = _system_name(s, b))
        _zzrept!(st, "[j]", "*j*", l2r)
        _zzremt!(st, 'j')
        _zzist(st, 'n') || _zzsubt!(st, "i", "n", l2r)
        _zzcmbt!(st, "-n", 'n', l2r)
        _zzsubt!(st, "n", "J", l2r)
        ok, tvec, ntvec, type, np, error = _zzunpck!(pic, st, false)
        modify = (era, wday, zone, ampm, system)
        pictur = wantpic ? _picture_string(pic, np, "JULIAND.", modify, false) : ""
        return _TimeParts(ok, tvec, ntvec, type, modify, mods, false, pictur, error)
    end

    # Calendar date: an integer of 1000 or more is a year.
    _zzvalt!(st, 1000, Int(typemax(Int32)), 'Y')

    # With the ISO separator only the built-in patterns are tried.
    if _zzist(st, 't')
        _zzknown!(st) || return _fail_parts(
            "The input string uses the ISO \"T\" date/time delimiter but does " *
                "not match any of the accepted ISO formats."
        )
        ok, tvec, ntvec, type, np, error = _zzunpck!(pic, st, false)
        modify = (era, wday, zone, ampm, system)
        pictur = wantpic ? _picture_string(pic, np, ".#", modify, true) : ""
        return _TimeParts(ok, tvec, ntvec, type, modify, false, false, pictur, error)
    end

    # Abbreviated years: a quote and an integer below 100.
    _zzvalt!(st, 100, 1000, 'I')
    yabbrv = _zzrept!(st, "'i", "*Y", l2r)
    while _zzsubt!(st, "I", "i", l2r)
    end

    # Modifiers, with or without parentheses.
    for (pat, transl) in (("[e]", "*e*"), ("[w]", "*w*"), ("[N]", "*N*"), ("[Z]", "*Z*"), ("[s]", "*s*"))
        _zzrept!(st, pat, transl, l2r)
    end
    _zzsubt!(st, "ie", "Ye", l2r)
    havera, be, ee = _zznote!(st, 'e')
    havwdy, bw, ew = _zznote!(st, 'w')
    havzon, bz, ez = _zznote!(st, 'Z')
    havapm, bn, en = _zznote!(st, 'N')
    havsys, bs, es = _zznote!(st, 's')
    mods = havera || havwdy || havzon || havapm || havsys
    havera && (era = _upper(s[be]) == UInt8('A') ? String63("A.D.") : String63("B.C."))
    havwdy && (wday = _weekday_abbr(s, bw))
    havapm && (ampm = _upper(s[bn]) == UInt8('A') ? String63("A.M.") : String63("P.M."))
    havsys && (system = _system_name(s, bs))
    if havzon
        zone = _upcut(s, bz, ez, modlen)
        for k in eachindex(_TK_ZONES)
            if zone == _TK_ZONES[k]
                zone = String63(_TK_OFFSETS[k])
                break
            end
        end
    end
    modify = (era, wday, zone, ampm, system)

    known = _zzknown!(st)
    if !known
        # Two delimiters in a row, or a delimiter at one end, make the string ambiguous.
        pair, from, to = _zzispt(st, ",/-:d.")
        pair && return _fail_parts(
            "There are two successive delimiters <$(_sub(s, from, to))> in " *
                "the input string. This is an ambiguous input. '" *
                _marked(s, from, to) * "'", modify
        )
        # The first and the last of the 12 letters of the key that are not blank.
        r = 1
        for k in 1:12
            st.rep[k] != _BLANK && (r = k)
        end
        edge = _isin(st.rep[1], codeunits(",/-:.")) ? st.rep[1] :
            (_isin(st.rep[r], codeunits(",/-:.")) ? st.rep[r] : _BLANK)
        if edge != _BLANK
            _zzsubt!(st, _delimiter_string(edge), "Q", l2r)
            _, from, to = _zznote!(st, 'Q')
            return _fail_parts(
                "An unexpected delimiter ('$(_sub(s, from, to))') was encountered " *
                    "in the input string. '" * _marked(s, from, to) * "'", modify
            )
        end

        # Remove commas, then dashes, then slashes, and try the patterns after each step.
        for delim in (',', '-', '/')
            _zzremt!(st, delim)
            known = _zzknown!(st)
            known && break
        end
    end

    if !known
        # What is left of the modifier classes is a second modifier of one kind.
        b = e = 0
        for c in ('e', 's', 'Z', 'w', 'N')
            found, b, e = _zznote!(st, c)
            found && break
        end
        b > 0 && return _fail_parts(
            "The substring \"$(_sub(s, b, e))\" is a duplicate modifier of " *
                "the input string: '" * _marked(s, b, e) * "'", modify
        )
        for c in ('[', ']', 'O', 'o', 'z')
            found, b, e = _zznote!(st, c)
            found && break
        end
        b > 0 && return _fail_parts(
            "The substring \"$(_sub(s, b, e))\" could not be resolved in " *
                "the input string: '" * _marked(s, b, e) * "'", modify
        )
        if _zzist(st, 'm') && _zzist(st, 'd')
            _, b1, e1 = _zznote!(st, 'm')
            _, b2, e2 = _zznote!(st, 'd')
            return _fail_parts(
                "Both a month \"$(_sub(s, b1, e1))\" and day of year delimiter " *
                    "\"$(_sub(s, b2, e2))\" appear in the input string: '$str'", modify
            )
        end

        # Hours, minutes and seconds: the rightmost run of numbers with colons.
        for (pat, meaning) in (
                ("i:i:i:n", "D*H*M*S"), ("i:i:i:i", "D*H*M*S"), ("i:i:n", "H*M*S"),
                ("i:i:i", "H*M*S"), ("i:n", "H*M"), ("i:i", "H*M"),
            )
            _zzrept!(st, pat, meaning, r2l) && break
        end
        _zzremt!(st, ':')

        # An integer next to a month name, and the day-of-year forms.
        for (pat, meaning, remove) in (
                ("<miiH", "mDY", false), ("<mi", "mD", false),
                ("Siim>", "SYDm", false), ("im>", "Dm", false), ("miY>", "mDY", false),
                ("Ymi", "YmD", false), ("Smi", "SmD", false), ("Mmi", "MmD", false),
                ("imY", "DmY", false), ("imH", "DmH", false), ("Yid", "Yy*", true),
                ("iYd", "yY*", true), ("Ydi", "Y*y", true),
            )
            (remove ? _zzrept!(st, pat, meaning, l2r) : _zzsubt!(st, pat, meaning, l2r)) && break
        end
    end

    ok, tvec, ntvec, type, np, error = _zzunpck!(pic, st, yabbrv)
    # A string that matches a built-in pattern always gets its picture; another one
    # has none when a token has no meaning.
    pictur = wantpic && (known || np > 0) ? _picture_string(pic, np, ".#", modify, true) : ""
    return _TimeParts(ok, tvec, ntvec, type, modify, mods, yabbrv, pictur, error)
end

"The delimiter `c` (one of `,/-:.`) as a `String` of one character."
_delimiter_string(c::UInt8) = c == UInt8(',') ? "," : c == UInt8('/') ? "/" :
    c == UInt8('-') ? "-" : c == UInt8(':') ? ":" : "."


"Message for a Julian date whose number cannot be read (CSPICE returns the value of an earlier call for such a string)."
const _NOJDVALUE = "The Julian date of the string has no value: the sign and the number " *
    "must be one token, as in 'JD -2451545.5'."

# ---------------------------------------------------------------------------------
# Checks and small rules: TCHECK, TEXPYR, ZZUTCPM
# ---------------------------------------------------------------------------------

"Return a year 0..99 as the year of the 100 years from `lbound` (1969..2068 at start) that it abbreviates, and another year unchanged (`texpyr`; bound)."
function _expand_year(year::Int, lbound::Int)::Int
    (year >= 100 || year < 0) && return year
    y = year + lbound ÷ 100 * 100
    return y < lbound ? y + 100 : y
end

"The year of `utc2et_` and `tparse_`: `1 - year` in the era B.C., else a year below 100 without the era A.D. expanded by `_expand_year`."
function _era_year(year::Int, era::AbstractString, lbound::Int)::Int
    era == "B.C." && return 1 - year
    (era != "A.D." && year < 100) && return _expand_year(year, lbound)
    return year
end

"""
    _zzutcpm(str, start) -> (hoff, moff, last, ok)

Return the signed hour and minute offsets of `str[start:]`, of the form `::UTC+h[:mm]`
or `::UTC-h[:mm]` with one or two digits in each number, `h < 13` and `mm <= 59`,
the position of the last character used and whether the hours were read (`zzutcpm`); minutes that cannot be read count as 0.
`str` is any vector of bytes with `length` and indexing; the call is inlined, so that a
buffer on the stack stays there.
"""
@inline function _zzutcpm(str, start::Int)
    n = length(str)
    at(i) = i <= n ? str[i] : _BLANK
    signat = start + 5
    sign = at(signat) == UInt8('+') ? 1.0 : (at(signat) == UInt8('-') ? -1.0 : 0.0)
    sign == 0.0 && return 0.0, 0.0, start - 1, false
    digits(from) = begin
        to = from - 1
        while _isdigit(at(to + 1))
            to += 1
        end
        to
    end
    unsat = start + 6
    unsto = digits(unsat)
    (1 <= unsto - unsat + 1 <= 2) || return 0.0, 0.0, start - 1, false
    x = _digits_value(str, unsat, unsto)
    x >= 13.0 && return 0.0, 0.0, start - 1, false
    hoff = sign * x
    last = unsto
    at(unsto + 1) == UInt8(':') || return hoff, 0.0, last, true
    unsat = unsto + 2
    unsto = digits(unsat)
    (1 <= unsto - unsat + 1 <= 2) || return hoff, 0.0, last, true
    x = _digits_value(str, unsat, unsto)
    x > 59.0 && return hoff, 0.0, last, true
    return hoff, sign * x, unsto, true
end

"The value of the one or two decimal digits `str[from:to]`, which `nparsd` gives exactly."
@inline function _digits_value(str, from::Int, to::Int)::Float64
    v = 0
    for i in from:to
        v = 10 * v + Int(str[i] - UInt8('0'))
    end
    return Float64(v)
end

"Return (hoff, moff, ok) of `_zzutcpm` for the zone `zone` (`\"\"` or a zone such as `UTC-7:30`) as `str2et` holds it in a 16-character variable after `::`."
function _zone_offsets(zone::AbstractString)
    n = ncodeunits(zone)
    field = SVector(
        ntuple(Val(16)) do i
            i <= 2 ? UInt8(':') : (i - 2 <= n ? codeunit(zone, i - 2) : _BLANK)
        end
    )
    hoff, moff, _, ok = _zzutcpm(field, 1)
    return hoff, moff, ok
end

_num(x::Float64) = isinteger(x) && abs(x) < 1.0e15 ? string(Int(x)) : string(x)

"""
    _tcheck(tvec, type, mods, modify) -> (ok, error)

Return whether the components `tvec` of type `"YMD"` or `"YD"` are in the ranges of a
calendar date (`tcheck` with checking on): a whole year, a month 1..12, a day of that month or of that year, hours 0..23 (1..12 with an
A.M./P.M. modifier), minutes 0..59, seconds below 60, or below 61 in the last minute
of June 30 and December 31; and a fraction only in the last component that is not
zero. `error` says which component is out of range.
"""
function _tcheck(
        tvec::SVector{8, Float64}, type::String, mods::Bool, modify::NTuple{5, String63}
    )
    year = _nint(tvec[1])
    myear = (mods && modify[_MOD_ERA] == "B.C.") ? 1 - year : year
    leapdy = _leap(abs(myear))
    dinmon(m) = m == 2 ? 28.0 + leapdy : Float64((31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[m])
    dinyr = 365.0 + leapdy
    jun30 = 181.0 + leapdy
    ampm = mods && modify[_MOD_AMPM] != ""
    hubnd, hlbnd = ampm ? (13.0, 1.0) : (24.0, 0.0)
    tvec[1] == Float64(year) || return false, "The year value was $(_num(tvec[1])). " *
        "This must be an integral value."
    local day::Int, doy::Float64
    if type == "YD"
        day = 2
        doy = tvec[2]
        (tvec[2] >= dinyr + 1.0 || tvec[2] < 1.0) && return false, "Day $(_num(tvec[2])) " *
            "has been specified for the year $year. The correct range for the day of " *
            "year for this year is from 1 to $(365 + leapdy)."
    else
        month = _nint(tvec[2])
        day = 3
        tvec[2] == Float64(month) || return false, "The month specified, " *
            "$(_num(tvec[2])), was not an integer. The month must be an integer in the " *
            "range from 1 to 12."
        (tvec[2] < 1.0 || tvec[2] > 12.0) && return false, "The month specified was " *
            "$month. The month must be an integer in the range from 1 to 12 (inclusive)."
        (tvec[3] < 1.0 || tvec[3] >= dinmon(month) + 1.0) && return false, "The day of " *
            "the month specified for the month of $(_MONTH_NAMES[month]) was " *
            "$(_num(tvec[3])). For $(_MONTH_NAMES[month]) the day must be at least 1 " *
            "and less than $(_num(dinmon(month) + 1.0))."
        doy = 0.0
        for i in 1:(month - 1)
            doy += dinmon(i)
        end
        doy += tvec[3]
    end
    hour, minute, second = day + 1, day + 2, day + 3
    (tvec[hour] >= hubnd || tvec[hour] < hlbnd) && return false, "The hours component " *
        "of the time specified was $(_num(tvec[hour])). " * (
            ampm ?
            "When either A.M. or P.M. is specified with the time the hours component must " *
            "be at least 1 and less than 13." :
            "The hours component must be greater than or equal to 0 and less than 24."
        )
    (tvec[minute] >= 60.0 || tvec[minute] < 0.0) && return false, "The minutes " *
        "component of the time specified was $(_num(tvec[minute])). This value must " *
        "be greater than or equal to 0 and less than 60."
    if tvec[second] >= 60.0 || tvec[second] < 0.0
        pm = mods && modify[_MOD_AMPM] == "P.M."
        leapsec = tvec[second] < 61.0 && tvec[second] > 0.0 && tvec[minute] == 59.0 &&
            (tvec[hour] == 23.0 || (tvec[hour] == 11.0 && pm)) &&
            (doy == dinyr || doy == jun30)
        leapsec || return false, "The seconds component of time must be at least 0 and " *
            "less than 60 (61 during the last minute of June 30 and December 31). The " *
            "value supplied was $(_num(tvec[second]))."
    end
    names = ("days", "hours", "minutes", "seconds")
    # Loop invariant: no component before `i` has a fraction together with a later
    # component that is not zero.
    for i in day:minute
        tvec[i] == Float64(_nint(tvec[i])) && continue
        for j in (i + 1):second
            tvec[j] == 0.0 || return false, "The '$(names[i - day + 1])' component of the " *
                "date has a fractional component. This is allowed only if all " *
                "components of lesser significance have value 0. However the " *
                "'$(names[j - day + 1])' component has value $(_num(tvec[j]))."
        end
    end
    return true, ""
end

# ---------------------------------------------------------------------------------
# Defaults of the time routines: TIMDEF
# ---------------------------------------------------------------------------------

"""
    _timdef_field(s, n = 16) -> String

`s` as `ljust_` and `ucase_` write it into a variable of `n` bytes in TIMDEF, without
the trailing blanks: the bytes after the first `n` of `s` without leading blanks are
lost, and only ASCII letters change case.
"""
function _timdef_field(s::AbstractString, n::Int = 16)::String
    b = codeunits(lstrip(s, ' '))
    m = min(n, length(b))
    while m > 0 && b[m] == UInt8(' ')     # invariant: b[m+1:min(n, end)] is blank
        m -= 1
    end
    out = String(map(_upper, b[1:m]))
    @assert ncodeunits(out) <= n
    return out
end

"""
    _timdef_set(defs, item, value) -> TimeDefaults

Return `defs` with the default `item` (`"CALENDAR"`, `"SYSTEM"` or `"ZONE"`) set to
`value`, both in any case with leading blanks allowed: a system clears the zone and a
zone clears the system, as in CSPICE. Throws
`SpiceError(:badtimeitem)` for another item and `SpiceError(:baddefaultvalue)` for
another value. Pure.
"""
function _timdef_set(defs::TimeDefaults, item::AbstractString, value::AbstractString)
    it = _timdef_field(item)
    val = _timdef_field(value)
    if it == "SYSTEM"
        val in ("TDB", "TDT", "TT", "UTC") || throw(
            SpiceError(
                :baddefaultvalue,
                "the default time system '$value' is not one of UTC, TDB, TDT, TT"
            )
        )
        return TimeDefaults(defs.calendar, val, "", defs.yearbound, defs.check)
    elseif it == "ZONE"
        zone = findfirst(==(val), _TK_ZONES)
        zone === nothing || (val = _TK_OFFSETS[zone])
        val = _timdef_field(val, 14)
        _, _, ok = _zone_offsets(val)
        ok || throw(
            SpiceError(
                :baddefaultvalue, "the default time zone '$value' is not " *
                    "one of EST, EDT, CST, CDT, MST, MDT, PST, PDT and is not an offset " *
                    "UTC+h[:mm] or UTC-h[:mm] with h below 13"
            )
        )
        return TimeDefaults(defs.calendar, "", String(val), defs.yearbound, defs.check)
    elseif it == "CALENDAR"
        val in ("JULIAN", "GREGORIAN", "MIXED") || throw(
            SpiceError(
                :baddefaultvalue,
                "the default calendar '$value' is not one of GREGORIAN, JULIAN, MIXED"
            )
        )
        return TimeDefaults(val, defs.system, defs.zone, defs.yearbound, defs.check)
    end
    throw(
        SpiceError(
            :badtimeitem,
            "the default item '$item' is not one of CALENDAR, SYSTEM, ZONE"
        )
    )
end

"Return the value of the default `item` (`CALENDAR`, `SYSTEM` or `ZONE`, any case) in `defs`; throws `SpiceError(:badtimeitem)` for another item."
function _timdef_get(defs::TimeDefaults, item::AbstractString)::String
    it = _timdef_field(item)
    it == "CALENDAR" && return defs.calendar
    it == "SYSTEM" && return defs.system
    it == "ZONE" && return defs.zone
    throw(
        SpiceError(
            :badtimeitem,
            "the default item '$item' is not one of CALENDAR, SYSTEM, ZONE"
        )
    )
end

# ---------------------------------------------------------------------------------
# Public calls: STR2ET, UTC2ET, TPARSE, TPICTR
# ---------------------------------------------------------------------------------

# The kinds of CSPICE: for `str2et`, UNPARSEDTIME for a string that does not parse,
# BADTIMESTRING for a component out of range, TIMECONFLICT and TIMEZONEERROR; for
# `utc2et`, INVALIDTIMESTRING for every reason. `tparse_c` and `tpictr_c` return the
# message without an error; here they throw the kinds of `str2et`.
@noinline _badtime(s, why, kind::Symbol) =
    throw(SpiceError(kind, "cannot convert the time string '$s': $why"))

"""
    str2et([ctx,] s) -> Float64

Return the epoch of the time string `s` in seconds past J2000 TDB, as
`str2et_c`.

A leap-second kernel must be loaded. The string forms of CSPICE are accepted:

| Form                             | Examples                                             |
|:---------------------------------|:-----------------------------------------------------|
| ISO 8601                         | `"2025-03-01T12:00:00.5"`, `"2025-060T12:00"`        |
| calendar, month name or number   | `"1 MAR 2025 12:00:00"`, `"March 1, 2025 3:00 PM"`, `"3/1/2025"` |
| day of year                      | `"2025-060 // 12:00"`                                |
| Julian date                      | `"JD 2460736.0"`                                     |

A string can also hold an era (`B.C.`, `A.D.`), `A.M.` or `P.M.`, a weekday, and a time zone
(`PDT`, `UTC-7:30`) or a time system (`UTC`, `TDB`, `TDT`; not `TT`, as in CSPICE), but not
a zone and a system together. A year can be abbreviated as `'97` or given with two digits.

- Without a zone and a system the string is UTC, or what [`timdef`](@ref) set for the
  context; a default system `"TT"` leaves it UTC, as in CSPICE.
- The date is Gregorian, or of the calendar that `timdef` set for the context.
- A year of two digits is in 1969 to 2068, or in the 100 years that
  [`tsetyr`](@ref) set for the context.
- A UTC leap second `23:59:60` is accepted on a day that has one.
- The components are checked against their ranges whatever [`tparch`](@ref) set, so
  `"1996-02-30"` throws.

Throws `SpiceError`, with the reason in the message:
- `:emptystring` if `s` is `""`;
- `:unparsedtime` if `s` is not a time string of these forms;
- `:badtimestring` if a component is out of range;
- `:timeconflict` if `s` names both a time system and a time zone;
- `:timezoneerror` if the time zone is not one of the above;
- `:valueoutofrange` if the day number of the epoch is outside the 32-bit integers;
- `:noleapseconds` if no leap-second kernel is loaded.

!!! warning "Difference from CSPICE"
    CSPICE N0067 keeps day numbers in 32-bit integers and lets them wrap, so
    `"238409829538 JD"` gives a wrong ET. This call throws
    `SpiceError(:valueoutofrange)` for a day number outside the 32-bit integers.

# Examples
```julia
furnsh("naif0012.tls")
et = str2et("2025-03-01T12:00:00")              # UTC
et = str2et("2025 MAR 01 12:00:00 TDB")         # the string names the time system
et = str2et("March 1, 2025 4:00 AM PST")        # the string names a time zone
timout(et, "YYYY-MM-DDTHR:MN:SC.### ::UTC ::RND")    # back to a string, rounded
```

See also [`timout`](@ref), [`et2utc`](@ref), [`utc2et`](@ref), [`tparse`](@ref).

# References

- [NAIF documentation of `str2et_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/str2et_c.html)
"""
function str2et(ks::KernelSet, s::AbstractString, defs::TimeDefaults = TimeDefaults())::Float64
    _chkfstr("str2et", "time string", s)       # CHKFSTR of str2et_c
    return _str2et(ks, s, defs)
end

# `str2et_` without the check of `str2et_c`. The readers of kernel variables call it, as
# SPICELIB calls `str2et_`: there an empty string is a blank one, which does not parse.
function _str2et(ks::KernelSet, s::AbstractString, defs::TimeDefaults = TimeDefaults())::Float64
    l = _lsk(ks)
    dhoff, dmoff, _ = defs.zone == "" ? (0.0, 0.0, true) : _zone_offsets(defs.zone)
    p = _tpartv(s, 16, false)
    p.ok || _badtime(s, p.error, :unparsedtime)
    modify = p.modify
    era, zone, ampm, system = modify[_MOD_ERA], modify[_MOD_ZONE], modify[_MOD_AMPM], modify[_MOD_SYSTEM]
    mods = p.mods
    tvec = p.tvec
    type = p.type
    (zone != "" && system != "") && _badtime(
        s,
        "both a time system '$system' and a time zone '$zone' are specified; specify " *
            "one of them only", :timeconflict
    )
    # Without a zone and a system in the string, those of `defs` apply.
    usedefs = zone == "" && system == ""
    zoned = usedefs ? defs.zone != "" : zone != ""
    hoff = moff = 0.0
    if usedefs
        hoff, moff = dhoff, dmoff
    elseif zone != ""
        hoff, moff, ok = _zone_offsets(zone)
        ok || _badtime(
            s, "the time zone '$zone' is not an offset UTC+h[:mm] " *
                "or UTC-h[:mm] with h below 13", :timezoneerror
        )
    end
    sysname = usedefs ? defs.system : system
    tdb = sysname == "TDB"
    tdt = sysname == "TDT"

    if type == "JD"
        isnan(tvec[1]) && _badtime(s, _NOJDVALUE, :unparsedtime)
        from = tdb ? _TT_JDTDB : (tdt ? _TT_JDTDT : _TT_JDUTC)
        return _ttrans(l, from, _TT_TDB, _tvec(tvec[1]))[1]
    end
    yd = type == "YD"
    yr = 1
    mm = 2
    dy, hr, mn, sc = yd ? (2, 3, 4, 5) : (3, 4, 5, 6)
    forml = yd ? _TT_YDF : _TT_YMDF
    ptype = yd ? _TT_YD : _TT_YMD
    calndr = defs.calendar
    set(v, x, i) = Base.setindex(v, x, i)

    # A leap second in a zone or in the Julian calendar cannot be checked before the
    # date is moved to UTC and to the Gregorian calendar: check with one second less.
    adjust = false
    if (zoned || calndr != "GREGORIAN") && 60.0 <= tvec[sc] < 61.0
        adjust = true
        tvec = set(tvec, tvec[sc] - 1.0, sc)
    end
    if calndr == "JULIAN" || (calndr == "MIXED" && tvec[yr] < 1580.0)
        # The Julian calendar has a leap day every fourth year: accept the date if it
        # is valid in this year or in the Gregorian leap pattern four years later.
        ok1, _ = _tcheck(set(tvec, tvec[1] + 4.0, 1), type, mods, modify)
        ok2, error = _tcheck(tvec, type, mods, modify)
        ok = ok1 || ok2
    else
        ok, error = _tcheck(tvec, type, mods, modify)
    end
    ok || _badtime(s, error, :badtimestring)
    adjust && (tvec = set(tvec, tvec[sc] + 1.0, sc))
    ((tdt || tdb) && tvec[sc] >= 60.0) && _badtime(
        s,
        "the seconds component is $(_num(tvec[sc])), and the time system $sysname has " *
            "no leap seconds: the seconds must be less than 60", :badtimestring
    )
    (era == "B.C." && p.yabbrv) && _badtime(
        s,
        "the year is abbreviated with a quote and the era is B.C.; these do not go together",
        :badtimestring
    )
    era == "B.C." && (tvec = set(tvec, 1.0 - tvec[yr], yr))
    if ampm == "P.M."
        tvec[hr] < 12.0 && (tvec = set(tvec, tvec[hr] + 12.0, hr))
    elseif ampm == "A.M."
        tvec[hr] >= 12.0 && (tvec = set(tvec, tvec[hr] - 12.0, hr))
    end
    year = _nint(tvec[yr])
    if p.yabbrv || (year < 100 && era == "")
        tvec = set(tvec, Float64(_expand_year(year, defs.yearbound)), yr)
    end

    dojul = if calndr == "MIXED"
        yd ? (tvec[yr] < 1582.0 || (tvec[yr] == 1582.0 && tvec[dy] < 279.0)) :
            (
                tvec[yr] < 1582.0 || (tvec[yr] <= 1582.0 && tvec[mm] < 10.0) ||
                (tvec[yr] <= 1582.0 && tvec[mm] <= 10.0 && tvec[dy] < 6.0)
            )
    else
        calndr == "JULIAN"
    end
    orgnyr = _cint(tvec[yr])
    if dojul
        day = _cint(tvec[dy])
        frac = tvec[dy] - Float64(day)
        gy, gm, gd, gdoy = jul2gr(orgnyr, yd ? 1 : _cint(tvec[mm]), day)
        tvec = set(tvec, Float64(gy), yr)
        if yd
            tvec = set(tvec, Float64(gdoy) + frac, dy)
        else
            tvec = set(set(tvec, Float64(gm), mm), Float64(gd) + frac, dy)
        end
    end

    tv(v) = _tvec(v[1], v[2], v[3], v[4], v[5], v[6])
    if tdt
        return _unitim(l, _ttrans(l, forml, _TT_FORMAL, tv(tvec))[1], _US_TDT, _US_TDB)
    elseif tdb
        return _ttrans(l, forml, _TT_FORMAL, tv(tvec))[1]
    end
    if zoned
        tvec = set(set(tvec, tvec[hr] - hoff, hr), tvec[mn] - moff, mn)
        secs = tvec[sc]
        out = _ttrans(l, forml, forml, tv(set(tvec, 0.0, sc)))
        tvec = SVector(out[1], out[2], out[3], out[4], out[5], out[6], tvec[7], tvec[8])
        tvec = set(tvec, secs, sc)
    end
    if adjust
        ok, _ = _tcheck(tvec, type, false, modify)
        ok || _badtime(
            s, "the seconds component is out of range: " * (
                zoned ?
                    "in the time zone $(usedefs ? defs.zone : zone) " : ""
            ) * (
                dojul ? "on the Julian " *
                    "calendar " : ""
            ) * "a leap second can occur in the year $orgnyr only in the " *
                "second that follows 23:59:59 UTC on June 30 and December 31", :badtimestring
        )
    end
    return _ttrans(l, ptype, _TT_TDB, tv(tvec))[1]
end

# The range check of `tcheck_` that `tparch` turns on, for `utc2et` and `tparse`: it
# runs on the components as parsed, before a year of two digits is expanded.
function _tparch_check(s::AbstractString, defs::TimeDefaults, p, kind::Symbol)
    defs.check || return nothing
    ok, why = _tcheck(p.tvec, p.type, p.mods, p.modify)
    ok || _badtime(s, why, kind)
    return nothing
end

"""
    utc2et([ctx,] s) -> Float64

Return the epoch of the UTC time string `s` in seconds past J2000 TDB, as
`utc2et_c`.

[`str2et`](@ref) is the call for new code; this one keeps the behavior of the older
CSPICE routine:

- `s` has a form that `str2et` accepts, without a time zone, without `A.M.` or `P.M.`
  and without a time system other than `UTC`.
- The defaults of [`timdef`](@ref) are not used: the string is UTC and the date is
  Gregorian. The year bound of [`tsetyr`](@ref) is used.
- The components are not checked against their ranges, so `1996-02-30` is March 1,
  unless [`tparch`](@ref) has turned the check on for the context.

Throws `SpiceError(:emptystring)` if `s` is `""`, and
`SpiceError(:invalidtimestring)` if `s` is not a time string of these forms or,
with the check on, if a component is out of range, as CSPICE does for every such
string, `SpiceError(:valueoutofrange)` if the day number of the epoch is outside the
32-bit integers (CSPICE lets it wrap; see [`str2et`](@ref)), and
`SpiceError(:noleapseconds)` if no leap-second kernel is loaded.

# References

- [NAIF documentation of `utc2et_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/utc2et_c.html)
"""
function utc2et(ks::KernelSet, s::AbstractString, defs::TimeDefaults = TimeDefaults())::Float64
    _chkfstr("utc2et", "time string", s)
    l = _lsk(ks)
    p = _tpartv(s, 8)
    p.ok || _badtime(s, p.error, :invalidtimestring)
    modify = p.modify
    (modify[_MOD_SYSTEM] != "" && modify[_MOD_SYSTEM] != "UTC") && _badtime(
        s,
        "utc2et does not support a time system in the string; the time system " *
            "$(modify[_MOD_SYSTEM]) was specified. Use str2et.", :invalidtimestring
    )
    modify[_MOD_ZONE] != "" && _badtime(
        s, "utc2et does not support a time zone in the " *
            "string; the time zone '$(modify[_MOD_ZONE])' was specified. Use str2et.",
        :invalidtimestring
    )
    modify[_MOD_AMPM] != "" && _badtime(
        s,
        "utc2et does not support the A.M./P.M. conventions. Use str2et.", :invalidtimestring
    )
    tvec = p.tvec
    p.type == "JD" && isnan(tvec[1]) && _badtime(s, _NOJDVALUE, :invalidtimestring)
    p.type == "JD" && return _ttrans(l, _TT_JDUTC, _TT_TDB, _tvec(tvec[1]))[1]
    _tparch_check(s, defs, p, :invalidtimestring)
    year = _era_year(_nint(tvec[1]), modify[_MOD_ERA], defs.yearbound)
    tv = _tvec(year, tvec[2], tvec[3], tvec[4], tvec[5], tvec[6])
    return _ttrans(l, p.type == "YD" ? _TT_YD : _TT_YMD, _TT_TDB, tv)[1]
end

"""
    tparse([ctx,] s) -> Float64

Return the seconds past J2000 of the time string `s` on a formal calendar, as
`tparse_c`.

On the formal calendar every day has 86400 seconds: there are no leap seconds, and
the call needs no kernel. The result is the count that [`etcal`](@ref) writes back,
not an epoch of a uniform time scale; for a UTC string it is the `"UTC"` epoch of
[`deltet`](@ref). `s` has a form that [`str2et`](@ref) accepts, without a time zone,
a time system, or `A.M.` and `P.M.`.

A year of two digits follows [`tsetyr`](@ref) of the context: at start the year 69 is
1969 and the year 68 is 2068. The components are not checked against their ranges, so
`"1996-02-30"` is March 1, unless [`tparch`](@ref) has turned the check on for the
context.

Throws `SpiceError(:emptystring)` if `s` is `""`, `SpiceError(:unparsedtime)` if the
string cannot be parsed and, with the check
on, `SpiceError(:badtimestring)` if a component is out of range: the kinds of
[`str2et`](@ref). The message is the one that `tparse_c` returns in `errmsg`;
`tparse_c` signals no error.

# Examples
```jldoctest
julia> tparse("2000-01-01T12:00:00")
0.0

julia> tparse("2 JAN 2000")
43200.0
```

# References

- [NAIF documentation of `tparse_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/tparse_c.html)
"""
function tparse(s::AbstractString, defs::TimeDefaults)::Float64
    _chkfstr("tparse", "time string", s)
    p = _tpartv(s, 8)
    p.ok || _badtime(s, p.error, :unparsedtime)
    modify = p.modify
    modify[_MOD_SYSTEM] != "" && _badtime(
        s, "tparse does not support a time system " *
            "in the string; the time system $(modify[_MOD_SYSTEM]) was specified.", :unparsedtime
    )
    modify[_MOD_ZONE] != "" && _badtime(
        s, "tparse does not support a time zone in " *
            "the string; the time zone '$(modify[_MOD_ZONE])' was specified.", :unparsedtime
    )
    modify[_MOD_AMPM] != "" && _badtime(
        s, "tparse does not support the A.M./P.M. conventions.",
        :unparsedtime
    )
    tvec = p.tvec
    p.type == "JD" && isnan(tvec[1]) && _badtime(s, _NOJDVALUE, :unparsedtime)
    p.type == "JD" && return (tvec[1] - _JD2000) * _SPD
    _tparch_check(s, defs, p, :badtimestring)
    if p.type == "YD"
        tvec = SVector(tvec[1], 1.0, tvec[2], tvec[3], tvec[4], tvec[5], tvec[7], tvec[8])
    end
    year = _era_year(_nint(tvec[1]), modify[_MOD_ERA], defs.yearbound)
    month, day = _nint(tvec[2]), _nint(tvec[3])
    q = 0
    adjust = year < 1
    if adjust
        q, year = fldmod(year, 400)
        year += 400
        q -= 1
    end
    day = year * 367 - div((year + div(month + 9, 12)) * 7, 4) -
        div((div(year + div(month - 9, 7), 100) + 1) * 3, 4) + div(month * 275, 9) + day - 730516
    adjust && (day += q * 146097)
    day = Int(day % Int32)
    return (Float64(day) - 0.5) * _SPD + tvec[4] * 3600.0 + tvec[5] * 60.0 + tvec[6]
end

"""
    tpictr([ctx,] sample) -> String

Return the picture with which [`timout`](@ref) writes epochs in the form of the time
string `sample`, as
`tpictr_c`.

`sample` has a form that [`str2et`](@ref) accepts. The call needs no kernel; `ctx` is
accepted and not used.

Throws `SpiceError(:emptystring)` if `sample` is `""`, and `SpiceError(:unparsedtime)`,
the kind of [`str2et`](@ref), if `sample` is not a
time string; the message is the one that `tpictr_c` returns in `errmsg`, and
`tpictr_c` signals no error.

# Examples
```jldoctest
julia> tpictr("Mar 1, 2025 3:04:05 PM")
"Mon DD, YYYY AP:MN:SC AMPM"

julia> tpictr("2025-03-01T12:00:00.000")
"YYYY-MM-DDTHR:MN:SC.### ::RND"
```

# References

- [NAIF documentation of `tpictr_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/tpictr_c.html)
"""
function tpictr(sample::AbstractString)::String
    _chkfstr("tpictr", "sample", sample)
    p = _tpartv(sample, 8)
    strip(p.pictur, ' ') == "" && _badtime(sample, p.error, :unparsedtime)
    return String(rstrip(p.pictur, ' '))
end
