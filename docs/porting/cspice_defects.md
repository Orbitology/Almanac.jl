# CSPICE defects

Where CSPICE N0067 has a defect, Almanac.jl either copies it, so that its results equal
those of CSPICE, or fixes it, so that its results differ for the inputs that the defect
concerns. This page lists the known defects of the routines that Almanac.jl ports, and
what Almanac.jl does about each.

The defects come from the tests of Almanac.jl, from [NAIF's test suite](@ref "NAIF's test suite")
replayed against Almanac.jl, and from the
[list of SPICELIB bugs](https://github.com/jacobwilliams/pepper/blob/master/spicelib_bugs.md)
of [PEPPER](https://github.com/jacobwilliams/pepper), a Fortran port of SPICELIB. Not
listed: defects of routines that Almanac.jl does not port, and defects of the Fortran
source that cannot occur in Julia, such as a string passed as both input and output.
CSPICE is SPICELIB translated from Fortran to C by f2c, and some defects come from the
run-time library of f2c. In the tables:

- **Where** names the routine of SPICELIB in lower case, `_c` for a routine written in
  C. Parentheses name the public calls of Almanac.jl through which the defect shows,
  where the routine is not one of them.
- **Test** names the test group (`test/test_<group>.jl`) that covers the behavior, or
  the family of NAIF's test suite (`tspice F_NAME`) whose difference from CSPICE is in
  the baseline of that suite, the list of known differences that a replay accepts.

## Defects that Almanac.jl copies

| Where | Defect | In Almanac.jl | Test |
|---|---|---|---|
| `nparsd` (text kernels, SCLK kernels, TLE fields) | Converts a decimal string to a double that is not always the nearest one: the GM of the Sun in `gm_de431.tpc` comes out 5 units in the last place low | The same doubles. The EK query scanner fixes it (below) | `pool` |
| C casts of f2c (`(integer) x`, `i_dnnt`; `str2et`, `ekfind`) | A double outside the 32-bit range, or NaN, is undefined behavior in C: x86 gives the smallest 32-bit integer, ARM saturates | The x86 result on every machine, so it differs from CSPICE built for ARM | `time`, `ek_query` |
| `hx2dp` | The exponent -2³¹ overflows when negated, and the routine divides by 16 once instead of returning zero | The same | `math_vectors` |
| `gfsntc_c` | Checks the string `abcorr` twice and `obsrvr` never, so an empty observer name is not `SPICE(EMPTYSTRING)`: it goes on to the lookup of the body | The same checks | `emptystring` |

## Defects that Almanac.jl fixes

Most fixes replace a crash, a read of undefined memory or an endless loop with an
error or the correct result.

### Time, text kernels, frames

| Where | Defect | In Almanac.jl | Test |
|---|---|---|---|
| `tisbod` (`sxform`, states in a body-fixed frame) | For cubic angles (`BODY<sys>_MAX_PHASE_DEGREE = 3`, `<sys>` the ID of the system barycenter), divides the rate of the cubic term by the square of the seconds per century in place of the seconds per century, so the derivative block misses that term | The exact derivative; see [Derivatives of rotations](@ref) | `frameinfo` |
| `str2et`, `tparch` | Converting a second 60 in a time zone, or a Julian-calendar date, turns on the range check of `tparch` for the rest of the process, so later conversions depend on earlier ones | Conversions do not change the setting of `tparch`, which belongs to the context | `time` |
| `ttrans` (`str2et`) | Day numbers are 32-bit integers and wrap: the epoch `"238409829538 JD"` gives a wrong ET | `:valueoutofrange` | `time` |
| `ttrans` (`str2et`) | A month outside 1 to 12 in a year-month-weekday input aborts with a subscript error | `:badsubscript` | `time` |
| `zzrefch0`, `zzfrmch0` (`pxform`, `sxform`, every query in a frame) | Frame definitions that form a cycle (frame A relative to B, B relative to A) make the search of the chain loop forever | `:framechain` after 1000 links | `rotations` |
| `et2lst` | The short message is cut to `SPICE(CANTGETROTATIONTYPE`, without its closing parenthesis | `:cantgetrotationtype` | `frameinfo` |
| `ccifrm` | Finds the `FRAME_<ID>_NAME` variable by replacing the first `_CLASS` of the key `FRAME_<name>_CLASS`, so a frame name that contains `_CLASS` gives a wrong variable | Replaces the trailing `_CLASS` | `frameinfo` |

### SPK, CK and binary PCK

| Where | Defect | In Almanac.jl | Test |
|---|---|---|---|
| `spkpvn` (types 1, 2, 3, 5, 10, 20, 21) | No check of the time bounds: the reader extrapolates, reads outside the segment, or stops the process (type 1) | `:timeoutofbounds`. Types 15 and 17 hold one set of elements and are evaluated at any epoch, as in CSPICE | `dafwrite` |
| `spke01` (`spkezr` and the other state calls, on SPK type 1) | No check of the step sizes of a difference line: a zero step divides by zero (`spke21` signals `SPICE(ZEROSTEP)` for type 21) | `:baddaf` | `spktypes` |
| `spks19` (`spksub` of type 19) | Reads `data[0]` and aborts when the last group of 100 epochs in the epoch directory starts after the `stop` epoch given to `spksub` | `:spkstructureerror` | `dafwrite` |
| `spkgeo` | An SPK chain with a cycle never returns | `SpkInsufficientDataError`, kind `:spkinsuffdata` | `management` |
| `ckw04a` | With more than one record in one call, packs the array in place and copies the midpoint and the radius of the second and later records after their coefficients have overwritten them | The values given | `dafwrite` |
| `ckr06` (`ckgp`, `ckgpav`, on CK type 6) | On moving to the next mini-segment, rounds its clock rate to a whole number and keeps it, so later derivatives are wrong | The stored rate | `ck` |
| `ckbsr` (`ckgp`, `ckgpav`) | Keeps the last segment found: `ckgp` right after `ckgpav` reuses it, although a segment of higher priority covers the request | No state between calls | `ck` |
| `ckr05` (`ckgp`, `ckgpav`, on CK type 5) | After a request in a gap between interpolation intervals, keeps that interval, so a later request finds no pointing where a new process finds some | No state between calls | `ck` |
| `zzsgp4` (`evsgp4`, SPK type 10) | An eccentricity above 1 gives no error, and the constants of the previous element set are kept | `:invalidelements` | `sgp4` |
| f2c `open` (`dafopw`, `dasopw`) | A file that the process may not write is opened for read without an error, so the open for write succeeds and the first write fails | `:fileopenfail` at the open | `dafwrite` |

### Geometry, events, vectors

| Where | Defect | In Almanac.jl | Test |
|---|---|---|---|
| `limbpt` | Does not count the room left, so with a DSK method it writes past `maxn` points | `:outofroom` | `geom_limbterm` |
| `termpt` | The room left in the output arrays (`maxn`) goes down by n² for each half plane with n points, not by n, so with a DSK method `OUTOFROOM` can come before `maxn` points are found | The room goes down by n | `geom_limbterm` |
| `zzgfocin` (`gfoclt`, `gfocce`, `occult`) | The surface names given for the back target are translated to IDs as surfaces of the front target. Surfaces given by number are not affected | The surfaces of the back target | `geom_occult`, `gf_occult` |
| `dskxsi` | Returns the plate of the last segment searched that had an intercept, which need not be the segment of the returned point | The plate of the returned point | `dsk` |
| `dskxsi`, `dskxv` | A call whose frame is not centered at the target (`INVALIDFRAME`) keeps the ID of that frame; until the pool changes, the next call with the same target and the frame of the last call that passed uses the failed frame, so the ray is in the wrong frame, or the call fails if that frame has no orientation data | Each call uses the frame it names | `dsk`, `tspice F_DSKXSI` |
| `dski02` | For item 1 (the number of vertices) with `start > 1`, gives `INDEXOUTOFRANGE`, `DASNOSUCHADDRESS` or no error, depending on earlier calls | `:indexoutofrange` | `dsk` |
| `dskmi2` | Does not check that the vertex indices are in 1 to `nv` | `:badvertexindex` | `dsk` |
| `zzmkspin`, `dskw02` (`dskmi2`, `dskw02`) | The fine voxel count is the 32-bit product of the three grid extents, so a grid of more than 2³¹ voxels can wrap to a count that passes the check against 100000000: 1000 × 1000 × 4295 voxels pass as 32704 | The product in 64 bits: `:valueoutofrange` | `dsk` |
| `zzdsksba` (`dskxsi`, `dskxv`, `srfnrm`, `latsrf` and the geometry calls with a DSK shape) | The segment table of a body looks up the SPK offset of every segment whose frame is not centered on the body, also of segments that no search uses (other surfaces, other times): missing ephemeris data for one of them makes every call for the body fail | The offset of a segment is looked up when a search uses it | `dsk` |
| `zzdsksph`, `zzdskbux`, `zzdskbun` (the same calls) | A segment whose frame is that of the segment examined before it reuses the offset of that segment, which is taken at the middle of the coverage of that segment, not of its own. In `zzdskbux` this affects the test of which segments the ray meets, so a ray can miss a segment that it crosses; in `zzdskbun` a point on a segment can be found in no segment (`POINTNOTINSEGMENT`) | The offset of each segment | `dsk` |
| `zzgflong` (`gfposc`, `gfsubc`, `gfsntc`) | For ABSMIN or ABSMAX with `adjust` on a longitude or right ascension, reads the first interval of an empty result window | Returns the empty window | `gf_posc` |
| `gfevnt` | A blank quantity name matches one of the two blank entries of the table of quantities, whose search does nothing: the call returns with no error and the result is left as it was | `:notrecognized` | `gf_quantities` |
| GF searches | A step below the spacing of doubles at the epoch makes the search sample the same epoch forever | `:invalidstep` | `gf_solver` |
| `zzgfrelx`, `zzgflong` (`gfuds`, `gfdist`, `gfposc` and the other GF searches with `adjust`) | Take a NaN or infinite `adjust` for ABSMIN and ABSMAX, since only `adjust < 0` and, for the other relations, `adjust != 0` are tested | `:invalidvalue` for every relation, the kind that CSPICE gives for the other relations | `gf_solver` |
| `zzgfrelx`, `zzgflong` (`gfuds`, `gfdist`, `gfposc` and the other GF searches with `refval`) | Take a NaN reference value: every comparison with it is false, so `<`, `>` and `=` find no event | `:invalidvalue` for every relation | `gf_solver`, `gf_posc` |
| `wnsumd` | Where rounding makes the variance of the interval lengths negative (intervals of equal length), the standard deviation is NaN | 0 | `windows` |
| `vperp_c` | The header says that a zero `b` gives `a`; the code gives the zero vector | `a`, as the header says | `math_vectors` |
| `vsepg_c` | Loses up to half the digits near 0 and π | The method of `vsep` | `math_vectors` |

### EK queries

The `zz` routines and `nparsd` in this table are steps of [`ekfind`](@ref), and the
parse errors show in [`ekpsel`](@ref) too. The first five rows change the rows that
`ekfind` returns. The tests reproduce each of the five with an internal switch of the
search, so that every difference from CSPICE is traced to one of them.

| Where | Defect | In Almanac.jl | Test |
|---|---|---|---|
| `zzekweed` | A row that two OR branches give can appear twice, and two different rows of two segments can be merged into one | Each row once | `ek_query` |
| `zzekillt`, `zzekille` | A double compared with an indexed integer column is rounded: `ID < 17.4` loses `ID = 17` | The exact comparison | `ek_query` |
| `zzekjsqz` | After `zzekweed` empties a segment, the rows of the next segment are lost and later rows are read from the wrong segment | The correct rows | `ek_query` |
| `nparsd` (query scanner) | `2.756` in a query gives the next double above the one nearest to 2.756, so it is not equal to a stored 2.756 | The nearest double | `ek_query` |
| `zzekqcon`, `zzekscmp` | An integer value outside the 32-bit range becomes the smallest 32-bit integer (`i_dnnt`): `DCOL > 10000000000` holds for every row that is not null | The value compared as a double | `ek_query` |
| `zzeksemc` | The check of a query stops at the first `IS NULL` or `IS NOT NULL`: the constraints after it and the ORDER BY columns are not checked, so a type mismatch passes, and a non-scalar column in a later constraint makes the search signal `SPICE(BUG)` | Checks every constraint and every ORDER BY column | `ek_parse` |
| `zzekjsrt` | ORDER BY or a join on a fixed-count string column longer than 32 characters fails with `STRINGTRUNCATED` | The sorted rows | `ek_query` |
| `zzekscan` | Reads past its buffer for a query over 2000 characters | A parse error, reported in `errmsg` | `ek_parse` |
| `zzeknrml` | A WHERE clause that ends with NOT reads one token past the end of the query | Reads no token past the end | `ek_parse` |
| `zzeknrml`, `zzekpars` | A WHERE clause with more than 1000 conjunctions, or that normalizes to slightly more than 1000 relations, overruns arrays of 1000 elements and stops the process | `:toomanyconjunctions` for more than 1000 conjunctions; the rows otherwise | `ek_parse` |
| `zzeksemc` | The message for a non-scalar column prints the query from an offset that is never set | Prints the column name | `ek_parse` |
| `zzekpars` | A qualified ORDER BY column (`T.COL`) directly before FROM, WHERE or SELECT is encoded one integer short, so the SELECT columns shift and a later step rejects the query | `:bug`, before the query runs | `ek_parse` |
| `ekgc`, `ekgd`, `ekgi`, `eknelt` | The message for a SELECT column index out of range prints the number of tables as the range | The number of SELECT columns | `ek_query` |

### EK files

| Where | Defect | In Almanac.jl | Test |
|---|---|---|---|
| `ekrcei` | Does not check the segment number, and searches pages that were never written | `:invalidindex` | `ek_file` |
| `ekucec`, `ekuced`, `ekucei` | Delete the old entry before checking the new one, so a failing update loses the entry | The entry stays | `ek_write` |
| `ekacec`, `ekaced`, `ekacei` | A second add to an entry is refused only for scalar numeric columns; for other columns the pages of the first entry are never freed | `:nonemptyentry` | `ek_write` |
| `ekacei`, `ekucei` | A bad segment number reads the column descriptor through a segment descriptor that is not set | The segment is checked first | `ek_write` |
| `ekaclc`, `ekacld`, `ekacli` | A column stored twice in one fast load gives an unreadable segment | `:columnstoredtwice` | `ek_write` |
| `ekffld` | Finishes a fast load in which a column was not stored; the segment cannot be read | Refuses it | `ek_write` |
| `ekacld`, `ekacli` (`zzekac04`, `zzekac05`) | An entry size below 1 for a non-null row aborts with a subscript error | `:badsubscript` | `ek_write` |
| `ekaclc`, `ekacld`, `ekacli` | An entry size below 0 for a null row reads outside the values | `:invalidcount` | `ek_write` |
| `ekdelr` on a fast-loaded segment | The record pages have a link count of 0, so the first deleted row frees a page and new rows overwrite the others | Correct link counts | `ek_write` |
| `zzekac04` to `zzekac06`, `zzekde04` to `zzekde06` (`ekaclc`, `ekacld`, `ekacli`, then `ekdelr`, `ekucec`, `ekuced`, `ekucei`) | A data page on which an array entry continues has no link, so an update frees it and the next entry overwrites another | Correct links | `ek_write` |
| `zzekfrx` (`ekdelr`, `ekucec`, `ekuced`, `ekucei`) | Looks for a key at one position of the index only. The index can hold equal keys in another order (the fast writers sort unstably, and strings longer than the column are cut), so the call gives `ITEMNOTFOUND` | Looks at all equal keys | `ek_write` |
| `zzektr1s` (`ekffld`) | For a fast load of 5167 or 5168 rows (the smallest trees of depth 3), sizes the subtrees of the middle nodes by dividing by MNKEYC in place of MNKIDC, so their leaves overlap: the record tree and the index trees hold more keys than they state, some of them twice | A valid tree | `ek_write` |

### Fixes that change only the bytes of a file

The values read back are the same; a byte-by-byte comparison with a file that CSPICE
writes differs.

| Where | Defect | In Almanac.jl | Test |
|---|---|---|---|
| `dafps` (binary PCK writers) | The 4 padding bytes of a summary with 5 integers keep the sixth integer of the last summary packed, often an SPK or CK address | Zero | `dafwrite` |
| `pck03b` | Leaves the unused half of the last integer word of the segment descriptor uninitialized, so the file holds what the stack held | Zero | `dafwrite` |
| `dafarr` (`dafac`) | Extends the file with a buffer that is not set | Zeros | `dafwrite` |
| DAS writers (DSK, EK) | The padding of a record is whatever the buffer held before | Blanks in a character record, zeros in a double precision or integer record | `dsk` |
| `zzekpage` (the EK writers) | Unused bytes of new character pages are NUL until an `ekopn` in the same process | Blanks | `ek_write` |

### Fixes that change no computed value

- `scs2e`, `wninsd`, `wncomd`, `lstcld`, `lxqstr`, `zzrvar`: a test that does not
  short-circuit reads one element past the end.
- Slow: `namfrm` never hits its cache; `ccifrm` examines each name up to 100 times.
- `stmp03`: SPICELIB built with Intel Fortran on Windows sums 10 terms; CSPICE and
  Almanac.jl sum 11.
- `zzekjtst`: the element indices of a non-key join constraint are never set. Only a
  non-scalar column reaches them, which the check of `zzeksemc` lets through after
  `IS NULL`; Almanac.jl refuses such a query.
- `zzekad06`: a message prints the record count as the range of column indices; the
  check cannot fail through a public call.

## Limits that Almanac.jl does not have

Where CSPICE signals an error past the size of a fixed table, Almanac.jl makes the call.

| Where | Limit of CSPICE | In Almanac.jl | Test |
|---|---|---|---|
| `zzddhopn`, `rdker`, `rdtext` (`furnsh`, `ldpool`, `dafopr`, `dasopr`) | A path is cut to 255 characters, so two paths with the same start are confused | The whole path | `dafwrite` |
| `ttrans` (`str2et`) | Only the first 140 leap seconds of a leap-seconds kernel are used; later ones are ignored with no error | All of them | `time` |
| `framex` (`namfrm`, `frmnam` and every query in a frame) | `HASHISFULL` past 5209 frames defined in kernels | No limit | `frameinfo` |
| `zzboddef` (`boddef`) | `TOOMANYPAIRS` past 161 pairs, where the header says 150 | No limit | `poolwrite` |
| `zzrvbf` (`lmpool`) | A string value over 80 characters is written past its slot | Keeps the first 80 characters and writes nothing past them | `poolwrite` |
| `keeper` (`furnsh`) | `FILENAMETOOLONG`, `PATHTOOLONG` past 255 characters, after the path symbols of a meta-kernel are replaced | The whole name | `tspice F_KEEPER` |
| `zzddhman`, `dafah` (`spkopn`, `ckopn`, `pckopn`, `furnsh`, the DSK loader) | `FTFULL`, `DAFFTFULL` past 5000 open files | No limit | `tspice F_DDHFTSIZE`, `F_ZZDSKBSR`, `F_PCKBSR` |
| `dafana` (`dafbna`) | `STFULL` past 20 arrays written at once | No limit | `tspice F_DAFANA` |
| `pool` (`swpool`) | `KERVARSETOVERFLOW`, `TOOMANYWATCHES`, `AGENTLISTOVERFLOW` past the sizes of the tables of watched names and agents | No limit | `tspice F_DWPOOL` |
| `zzbodker` (`bodc2n`, `bodn2c`) | `KERVARTOOBIG` past 14983 values of `NAIF_BODY_NAME` | No limit | `tspice F_BODCOD` |
| `zzscad01` (`scdecd`, `sce2c` and the other SCLK calls) | `TOOMANYCOEFFS` past 100000 coefficients, `TOOMANYPARTITIONS` past 9999 partitions of a type 1 SCLK kernel | No limit | `tspice F_ZZSC01` |
| `zzswffet` (`sxform`, `pxform` of switch frames) | `TOOMANYBASEFRAMES` past 15000 base frames of switch frames | No limit | `tspice F_ZZSWFXFM` |
