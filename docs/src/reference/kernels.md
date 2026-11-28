# Kernels, contexts and the kernel pool

## Contexts and errors

```@docs
SpiceContext
default_context
Almanac.kernelset
Almanac.KernelSet
SpiceError
SpiceKindError
SpkInsufficientDataError
Almanac.ERROR_KINDS
```

## Loading and inspecting kernels

```@docs
furnsh
unload
kclear
ktotal
kinfo
kdata
tkvrsn
getfat
spkobj
spkobj!
spkcov
spkcov!
```

## Text kernels in JSON

A file named `*.json` is read as a text kernel written as one JSON object, in the
format of [jsonspice](https://github.com/jacobwilliams/jsonspice), which
[PEPPER](https://github.com/jacobwilliams/pepper) also reads. [`furnsh`](@ref),
[`ldpool`](@ref), [`unload`](@ref), [`kinfo`](@ref) and meta-kernels accept it;
CSPICE does not. The extension decides the format: a JSON file with another extension
is read as a classic text kernel, and loads nothing, because it has no `\begindata`.

```json
// earth.json
{
    "BODY399_RADII": [6378.1366, 6378.1366, 6356.7519],
    "+NAIF_BODY_NAME": "SPOCK",
    "+NAIF_BODY_CODE": -1000,
    "EPOCH": "@2026-NOV-11/12:00",
}
```

- Each member is one assignment. A key `"+NAME"` appends to the variable, like `+=`.
- Numbers, `true` (1) and `false` (0) are numeric values. Strings are string values,
  and a string that starts with `@` is a time.
- An array holds numbers and times, or strings.
- `//` line comments and trailing commas are allowed; any other JSON5 syntax is an
  error.

JSON numbers are read correctly rounded, while the numbers of a text kernel are read
as CSPICE reads them, which can be a few units in the last place off, so one value
written in the two formats can give pools that differ in the last bits.

## Reading the kernel pool

```@docs
gdpool
gcpool
gipool
bodvrd
bodvcd
bodfnd
expool
dtpool
gnpool
szpool
stpool
badkpv
```

## Changing the kernel pool

```@docs
pdpool
pcpool
pipool
dvpool
lmpool
ldpool
clpool
swpool
cvpool
dwpool
poolversion
```

## Body and frame names

```@docs
bodn2c
bodc2n
bodc2s
bods2c
scn2id
scid2n
boddef
namfrm
frmnam
```
