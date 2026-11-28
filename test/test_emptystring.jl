# Group: empty string arguments. Each C function of CSPICE N0067 checks its string
# arguments first, with the macro CHKFSTR, and signals SPICE(EMPTYSTRING) for an empty
# one, whatever its other arguments are. Each entry of `EMPTYSTRING_CALLS` calls a port
# with the strings that its C function checks as the arguments of a closure, in the order
# of the Almanac signature; the test empties each in turn and gives "X" to the others.
# The other arguments are invalid where they can be (NaN, -1, `true`, empty arrays,
# closed handles), so that a check made before the strings gives another kind. The lists
# of checked strings come from the `_c.c` files of CSPICE N0067.

FULL || @eval Base.Experimental.@optlevel 0

const EMPTYSTRING_DIR = mktempdir()

"A handle of an SPK file that was created in `EMPTYSTRING_DIR` and closed."
function emptystring_daf()
    h = spkopn(joinpath(EMPTYSTRING_DIR, "f$(time_ns()).bsp"), "test", 0)
    dafcls(h)
    return h
end

"A handle of an EK file that was created in `EMPTYSTRING_DIR` and closed."
function emptystring_das()
    h = ekopn(joinpath(EMPTYSTRING_DIR, "f$(time_ns()).bes"), "test", 0)
    ekcls(h)
    return h
end

"The kind of the error of `f(args...)`, its type if it is not a `SpiceError`, or `nothing`."
function emptystring_kind(f, args)
    try
        f(args...)
    catch err
        return err isa SpiceError ? err.kind : typeof(err)
    end
    return nothing
end

"Fast-load pointers of an EK file that was created in `EMPTYSTRING_DIR` and closed."
function emptystring_fastload()
    h = ekopn(joinpath(EMPTYSTRING_DIR, "f$(time_ns()).bes"), "test", 0)
    _, rcptrs = ekifld(h, "T", 1, ["S"], ["DATATYPE = CHARACTER*(4)"])
    ekcls(h)
    return rcptrs
end

const EMPTYSTRING_CALLS = let KS = Almanac.kernelset(SpiceContext()), CTX = SpiceContext(),
        DAF = emptystring_daf(), DAS = emptystring_das(), FASTLOAD = emptystring_fastload(),
        NOVEC = SVector{3, Float64}[]
    [
        "azlcpo" => (method, target, abcorr, obsctr, obsref) -> azlcpo(KS, method, target, NaN, abcorr, true, true, Float64[], obsctr, obsref),
        "badkpv" => (caller, name, comp) -> badkpv(KS, caller, name, comp, -1, -1, ' '),
        "boddef" => (name) -> boddef(CTX, name, -1),
        "bodfnd" => (item) -> bodfnd(KS, -1, item),
        "bodn2c" => (name) -> bodn2c(KS, name),
        "bods2c" => (name) -> bods2c(KS, name),
        "bodvcd" => (item) -> bodvcd(KS, -1, item),
        "bodvrd" => (bodynm, item) -> bodvrd(KS, bodynm, item),
        "ckcov" => (ck, level, timsys) -> ckcov(KS, ck, -1, true, level, NaN, timsys),
        "ckgp" => (ref) -> ckgp(KS, -1, NaN, -1.0, ref),
        "ckgpav" => (ref) -> ckgpav(KS, -1, NaN, -1.0, ref),
        "ckmeta" => (meta) -> ckmeta(KS, -1, meta),
        "ckobj" => (ck) -> ckobj(ck),
        "ckopn" => (fname, ifname) -> ckopn(fname, ifname),
        "ckw01" => (ref, segid) -> ckw01(KS, DAF, NaN, NaN, -1, ref, segid, NaN, NaN, NaN),
        "ckw02" => (ref, segid) -> ckw02(KS, DAF, NaN, NaN, -1, ref, segid, NaN, NaN, NaN, NaN, NaN),
        "ckw03" => (ref, segid) -> ckw03(KS, DAF, NaN, NaN, -1, ref, segid, NaN, NaN, NaN, NaN),
        "ckw05" => (ref, segid) -> ckw05(KS, DAF, -1, -1, NaN, NaN, -1, ref, true, segid, NaN, NaN, NaN, -1, NaN),
        "cnmfrm" => (cname) -> cnmfrm(KS, cname),
        "convrt" => (in, out) -> convrt(NaN, in, out),
        "cvpool" => (agent) -> cvpool(CTX, agent),
        "dafopr" => (fname) -> dafopr(fname),
        "dafopw" => (fname) -> dafopw(fname),
        "dasonw" => (fname, ftype, ifname) -> dasonw(fname, ftype, ifname, -1),
        "dasopr" => (fname) -> dasopr(fname),
        "dasopw" => (fname) -> dasopw(fname),
        "deltet" => (eptype) -> deltet(KS, NaN, eptype),
        "dlaopn" => (fname, ftype, ifname) -> dlaopn(fname, ftype, ifname, -1),
        "dpgrdr" => (body) -> dpgrdr(KS, body, NaN, NaN, NaN, NaN, NaN),
        "drdpgr" => (body) -> drdpgr(KS, body, NaN, NaN, NaN, NaN, NaN),
        "dskobj" => (dskfnm) -> dskobj(dskfnm),
        "dskopn" => (fname, ifname) -> dskopn(fname, ifname, -1),
        "dsksrf" => (dskfnm) -> dsksrf(dskfnm, -1),
        "dskw02" => (frame) -> dskw02(CTX, DAS, -1, -1, -1, frame, -1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN),
        "dskxsi" => (target, fixref) -> dskxsi(KS, true, target, NaN, NaN, fixref, NaN, NaN),
        "dskxv" => (target, fixref) -> dskxv(KS, true, target, NaN, NaN, fixref, zeros(2, 1), zeros(2, 1)),
        "dskxv!" => (target, fixref) -> dskxv!(SVector{3, Float64}[], Bool[], KS, true, target, NaN, NaN, fixref, NOVEC, NOVEC),
        "dtpool" => (name) -> dtpool(KS, name),
        "dvpool" => (name) -> dvpool(CTX, name),
        "edterm" => (trmtyp, source, target, fixref, abcorr, obsrvr) -> edterm(KS, trmtyp, source, target, NaN, fixref, abcorr, obsrvr, -1),
        "ekacec" => (column) -> ekacec(DAS, -1, -1, column, Float64[]),
        "ekaced" => (column) -> ekaced(DAS, -1, -1, column, Float64[]),
        "ekacei" => (column) -> ekacei(DAS, -1, -1, column, Float64[]),
        "ekaclc" => (column) -> ekaclc(DAS, -1, column, Float64[], Float64[], Float64[], FASTLOAD),
        "ekacld" => (column) -> ekacld(DAS, -1, column, Float64[], Float64[], Float64[], FASTLOAD),
        "ekacli" => (column) -> ekacli(DAS, -1, column, Float64[], Float64[], Float64[], FASTLOAD),
        "ekbseg" => (tabnam) -> ekbseg(DAS, tabnam, String[], String[]),
        "ekccnt" => (table) -> ekccnt(KS, table),
        "ekfind" => (query) -> ekfind(CTX, query),
        "ekifld" => (tabnam) -> ekifld(DAS, tabnam, -1, String[], String[]),
        "eklef" => (fname) -> eklef(CTX, fname),
        "ekopn" => (fname, ifname) -> ekopn(fname, ifname),
        "ekopr" => (fname) -> ekopr(fname),
        "ekopw" => (fname) -> ekopw(fname),
        "ekpsel" => (query) -> ekpsel(CTX, query),
        "ekrcec" => (column) -> ekrcec(DAS, -1, -1, column),
        "ekrced" => (column) -> ekrced(DAS, -1, -1, column),
        "ekrcei" => (column) -> ekrcei(DAS, -1, -1, column),
        "ekucec" => (column) -> ekucec(DAS, -1, -1, column, Float64[]),
        "ekuced" => (column) -> ekuced(DAS, -1, -1, column, Float64[]),
        "ekucei" => (column) -> ekucei(DAS, -1, -1, column, Float64[]),
        "et2lst" => (type) -> et2lst(KS, NaN, -1, NaN, type),
        "et2utc" => (format) -> et2utc(KS, NaN, format, -1),
        "expool" => (name) -> expool(KS, name),
        "fovray" => (inst, abcorr, obsrvr) -> fovray(KS, inst, Float64[], "X", abcorr, obsrvr, NaN),
        "fovtrg" => (inst, target, tshape, abcorr, obsrvr) -> fovtrg(KS, inst, target, tshape, "X", abcorr, obsrvr, NaN),
        "furnsh" => (file) -> furnsh(CTX, file),
        "gcpool" => (name) -> gcpool(KS, name),
        "gdpool" => (name) -> gdpool(KS, name),
        "getfat" => (file) -> getfat(file),
        "getfvn" => (inst) -> getfvn(KS, inst),
        "gfdist" => (target, abcorr, obsrvr, relate) -> gfdist(KS, target, abcorr, obsrvr, relate, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gfevnt" => (gquant, op) -> gfevnt(KS, NaN, NaN, gquant, String[], String[], Float64[], Int[], Bool[], op, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gffove" => (inst, tshape, abcorr) -> gffove(KS, inst, tshape, Float64[], "X", "X", abcorr, "X", NaN, NaN, NaN, Window([0.0, 1.0])),
        "gfilum" => (method, angtyp, target, illmn, fixref, abcorr, obsrvr, relate) -> gfilum(KS, method, angtyp, target, illmn, fixref, abcorr, obsrvr, Float64[], relate, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gfocce" => (occtyp, front, fshape, back, bshape, abcorr, obsrvr) -> gfocce(KS, occtyp, front, fshape, "X", back, bshape, "X", abcorr, obsrvr, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gfocce, context form" => (occtyp, front, fshape, back, bshape, abcorr, obsrvr) -> gfocce(CTX, occtyp, front, fshape, "X", back, bshape, "X", abcorr, obsrvr, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gfoclt" => (occtyp, front, fshape, back, bshape, abcorr, obsrvr) -> gfoclt(KS, occtyp, front, fshape, "X", back, bshape, "X", abcorr, obsrvr, NaN, Window([0.0, 1.0])),
        "gfoclt, context form" => (occtyp, front, fshape, back, bshape, abcorr, obsrvr) -> gfoclt(CTX, occtyp, front, fshape, "X", back, bshape, "X", abcorr, obsrvr, NaN, Window([0.0, 1.0])),
        "gfpa" => (target, illmn, abcorr, obsrvr, relate) -> gfpa(KS, target, illmn, abcorr, obsrvr, relate, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gfposc" => (target, frame, abcorr, obsrvr, crdsys, coord, relate) -> gfposc(KS, target, frame, abcorr, obsrvr, crdsys, coord, relate, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gfrfov" => (inst, rframe, abcorr) -> gfrfov(KS, inst, Float64[], rframe, abcorr, "X", NaN, Window([0.0, 1.0])),
        "gfrr" => (target, abcorr, obsrvr, relate) -> gfrr(KS, target, abcorr, obsrvr, relate, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gfsep" => (targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, relate) -> gfsep(KS, targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, relate, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gfsntc" => (target, fixref, method, abcorr, dref, crdsys, coord, relate) -> gfsntc(KS, target, fixref, method, abcorr, "X", dref, Float64[], crdsys, coord, relate, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gfsubc" => (target, fixref, method, abcorr, obsrvr, crdsys, coord, relate) -> gfsubc(KS, target, fixref, method, abcorr, obsrvr, crdsys, coord, relate, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gftfov" => (inst, target, tshape, abcorr, obsrvr) -> gftfov(KS, inst, target, tshape, "X", abcorr, obsrvr, NaN, Window([0.0, 1.0])),
        "gfuds" => (relate) -> gfuds(KS, NaN, NaN, relate, NaN, NaN, NaN, Window([0.0, 1.0])),
        "gipool" => (name) -> gipool(KS, name),
        "gnpool" => (name) -> gnpool(KS, name),
        "hx2dp" => (string) -> hx2dp(string),
        "illumf" => (method, target, ilusrc, fixref, abcorr, obsrvr) -> illumf(KS, method, target, ilusrc, NaN, fixref, abcorr, obsrvr, Float64[]),
        "illumg" => (method, target, ilusrc, fixref, abcorr, obsrvr) -> illumg(KS, method, target, ilusrc, NaN, fixref, abcorr, obsrvr, Float64[]),
        "ilumin" => (method, target, fixref, abcorr, obsrvr) -> ilumin(KS, method, target, NaN, fixref, abcorr, obsrvr, Float64[]),
        "kdata" => (kind) -> kdata(KS, -1, kind),
        "kinfo" => (file) -> kinfo(KS, file),
        "ktotal" => (kind) -> ktotal(KS, kind),
        "latsrf" => (method, target, fixref) -> latsrf(KS, method, target, NaN, fixref, NaN),
        "ldpool" => (fname) -> ldpool(CTX, fname),
        "limbpt" => (method, target, fixref, abcorr, corloc, obsrvr) -> limbpt(KS, method, target, NaN, fixref, abcorr, corloc, obsrvr, Float64[], NaN, -1, NaN, NaN, -1),
        "lspcn" => (body) -> lspcn(KS, body, NaN, "X"),
        "ltime" => (dir) -> ltime(KS, NaN, -1, dir, -1),
        "namfrm" => (frname) -> namfrm(KS, frname),
        "occult" => (targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr) -> occult(KS, targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, NaN),
        "pckcov" => (pckfnm) -> pckcov(pckfnm, -1),
        "pckfrm" => (pckfnm) -> pckfrm(pckfnm),
        "pckopn" => (name, ifname) -> pckopn(name, ifname),
        "pckw02" => (frame, segid) -> pckw02(DAF, -1, frame, NaN, NaN, segid, NaN, NaN, NaN),
        "pcpool" => (name) -> pcpool(CTX, name, String[]),
        "pdpool" => (name) -> pdpool(CTX, name, Float64[]),
        "pgrrec" => (body) -> pgrrec(KS, body, NaN, NaN, NaN, NaN, NaN),
        "phaseq" => (target, illmn, obsrvr, abcorr) -> phaseq(KS, NaN, target, illmn, obsrvr, abcorr),
        "pipool" => (name) -> pipool(CTX, name, Int[]),
        "pxform" => (from, to) -> pxform(KS, from, to, NaN),
        "pxfrm2" => (from, to) -> pxfrm2(KS, from, to, NaN, NaN),
        "recpgr" => (body) -> recpgr(KS, body, Float64[], NaN, NaN),
        "scencd" => (sclkch) -> scencd(KS, -1, sclkch),
        "scs2e" => (sclkch) -> scs2e(KS, -1, sclkch),
        "sctiks" => (clkstr) -> sctiks(KS, -1, clkstr),
        "sincpt" => (method, target, fixref, abcorr, obsrvr, dref) -> sincpt(KS, method, target, NaN, fixref, abcorr, obsrvr, dref, Float64[]),
        "sincpt, context form" => (method, target, fixref, abcorr, obsrvr, dref) -> sincpt(CTX, method, target, NaN, fixref, abcorr, obsrvr, dref, Float64[]),
        "spk14b" => (segid, frame) -> spk14b(KS, DAF, segid, -1, -1, frame, NaN, NaN, -1),
        "spkacs" => (ref, abcorr) -> spkacs(KS, -1, NaN, ref, abcorr, -1),
        "spkapo" => (ref, abcorr) -> spkapo(KS, -1, NaN, ref, Float64[], abcorr),
        "spkaps" => (ref, abcorr) -> spkaps(KS, -1, NaN, ref, abcorr, Float64[], Float64[]),
        "spkcov" => (spkfnm) -> spkcov(spkfnm, -1),
        "spkcpo" => (target, outref, refloc, abcorr, obsctr, obsref) -> spkcpo(KS, target, NaN, outref, refloc, abcorr, Float64[], obsctr, obsref),
        "spkcpt" => (trgctr, trgref, outref, refloc, abcorr, obsrvr) -> spkcpt(KS, Float64[], trgctr, trgref, NaN, outref, refloc, abcorr, obsrvr),
        "spkcvo" => (target, outref, refloc, abcorr, obsctr, obsref) -> spkcvo(KS, target, NaN, outref, refloc, abcorr, Float64[], NaN, obsctr, obsref),
        "spkcvt" => (trgctr, trgref, outref, refloc, abcorr, obsrvr) -> spkcvt(KS, Float64[], NaN, trgctr, trgref, NaN, outref, refloc, abcorr, obsrvr),
        "spkez" => (ref, abcorr) -> spkez(KS, -1, NaN, ref, abcorr, -1),
        "spkezp" => (ref, abcorr) -> spkezp(KS, -1, NaN, ref, abcorr, -1),
        "spkezr" => (targ, ref, abcorr, obs) -> spkezr(KS, targ, NaN, ref, abcorr, obs),
        "spkgeo" => (ref) -> spkgeo(KS, -1, NaN, ref, -1),
        "spkgps" => (ref) -> spkgps(KS, -1, NaN, ref, -1),
        "spkltc" => (ref, abcorr) -> spkltc(KS, -1, NaN, ref, abcorr, Float64[]),
        "spkobj" => (spkfnm) -> spkobj(spkfnm),
        "spkopa" => (file) -> spkopa(file),
        "spkopn" => (fname, ifname) -> spkopn(fname, ifname),
        "spkpds" => (frame) -> spkpds(KS, -1, -1, frame, -1, NaN, NaN),
        "spkpos" => (targ, ref, abcorr, obs) -> spkpos(KS, targ, NaN, ref, abcorr, obs),
        "spkssb" => (ref) -> spkssb(KS, -1, NaN, ref),
        "spksub" => (ident) -> spksub(DAF, Float64[], ident, NaN, NaN, DAF),
        "spkw02" => (frame, segid) -> spkw02(KS, DAF, -1, -1, frame, NaN, NaN, segid, NaN, NaN, NaN),
        "spkw03" => (frame, segid) -> spkw03(KS, DAF, -1, -1, frame, NaN, NaN, segid, NaN, NaN, NaN),
        "spkw05" => (frame, segid) -> spkw05(KS, DAF, -1, -1, frame, NaN, NaN, segid, NaN, NaN, NaN),
        "spkw08" => (frame, segid) -> spkw08(KS, DAF, -1, -1, frame, NaN, NaN, segid, -1, NaN, NaN, NaN),
        "spkw09" => (frame, segid) -> spkw09(KS, DAF, -1, -1, frame, NaN, NaN, segid, -1, NaN, NaN),
        "spkw10" => (frame, segid) -> spkw10(KS, DAF, -1, -1, frame, NaN, NaN, segid, NaN, NaN, NaN),
        "spkw12" => (frame, segid) -> spkw12(KS, DAF, -1, -1, frame, NaN, NaN, segid, -1, NaN, NaN, NaN),
        "spkw13" => (frame, segid) -> spkw13(KS, DAF, -1, -1, frame, NaN, NaN, segid, -1, NaN, NaN),
        "spkw15" => (frame, segid) -> spkw15(KS, DAF, -1, -1, frame, NaN, NaN, segid, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN),
        "spkw17" => (frame, segid) -> spkw17(KS, DAF, -1, -1, frame, NaN, NaN, segid, NaN, NaN, NaN, NaN),
        "spkw18" => (frame, segid) -> spkw18(KS, DAF, -1, -1, -1, frame, NaN, NaN, segid, -1, NaN, NaN),
        "spkw20" => (frame, segid) -> spkw20(KS, DAF, -1, -1, frame, NaN, NaN, segid, NaN, NaN, NaN, NaN, NaN, NaN),
        "srfcss" => (bodstr) -> srfcss(KS, -1, bodstr),
        "srfnrm" => (method, target, fixref) -> srfnrm(KS, method, target, NaN, fixref, NaN),
        "srfs2c" => (srfstr, bodstr) -> srfs2c(KS, srfstr, bodstr),
        "srfscc" => (srfstr) -> srfscc(KS, srfstr, -1),
        "stpool" => (item, contin) -> stpool(KS, item, -1, contin),
        "str2et" => (timstr) -> str2et(KS, timstr),
        "subpnt" => (method, target, fixref, abcorr, obsrvr) -> subpnt(KS, method, target, NaN, fixref, abcorr, obsrvr),
        "subpnt, context form" => (method, target, fixref, abcorr, obsrvr) -> subpnt(CTX, method, target, NaN, fixref, abcorr, obsrvr),
        "subslr" => (method, target, fixref, abcorr, obsrvr) -> subslr(KS, method, target, NaN, fixref, abcorr, obsrvr),
        "subslr, context form" => (method, target, fixref, abcorr, obsrvr) -> subslr(CTX, method, target, NaN, fixref, abcorr, obsrvr),
        "swpool" => (agent) -> swpool(CTX, agent, String[]),
        "sxform" => (from, to) -> sxform(KS, from, to, NaN),
        "szpool" => (name) -> szpool(name),
        "tangpt" => (method, target, fixref, abcorr, corloc, obsrvr, dref) -> tangpt(KS, method, target, NaN, fixref, abcorr, corloc, obsrvr, dref, Float64[]),
        "termpt" => (method, ilusrc, target, fixref, abcorr, corloc, obsrvr) -> termpt(KS, method, ilusrc, target, NaN, fixref, abcorr, corloc, obsrvr, Float64[], NaN, -1, NaN, NaN, -1),
        "timdef" => (action, item) -> timdef(CTX, action, item),
        "timdef" => (value) -> timdef(CTX, "SET", "CALENDAR", value),   # `value` for SET only
        "timout" => (pictur) -> timout(KS, NaN, pictur),
        "tipbod" => (ref) -> tipbod(KS, ref, -1, NaN),
        "tisbod" => (ref) -> tisbod(KS, ref, -1, NaN),
        "tkvrsn" => (item) -> tkvrsn(item),
        "tparch" => (type) -> tparch(CTX, type),
        "tparse" => (string) -> tparse(KS, string),
        "tpictr" => (sample) -> tpictr(sample),
        "trgsep" => (targ1, shape1, frame1, targ2, shape2, frame2, obsrvr, abcorr) -> trgsep(KS, NaN, targ1, shape1, frame1, targ2, shape2, frame2, obsrvr, abcorr),
        "unitim" => (insys, outsys) -> unitim(KS, NaN, insys, outsys),
        "unload" => (file) -> unload(CTX, file),
        "utc2et" => (utcstr) -> utc2et(KS, utcstr),
        "wnreld" => (op) -> wnreld(Window([0.0, 1.0]), op, Window([0.0, 1.0])),
        "xfmsta" => (icosys, ocosys, body) -> xfmsta(KS, Float64[], icosys, ocosys, body),
    ]
end

group("empty string arguments") do
    cd(EMPTYSTRING_DIR) do       # a call that does not check would create its file here
        @testset "each checked string, emptied" begin
            for (name, f) in EMPTYSTRING_CALLS
                n = only(methods(f)).nargs - 1
                for i in 1:n
                    args = ntuple(j -> j == i ? "" : "X", n)
                    @test (name, i, emptystring_kind(f, args)) == (name, i, :emptystring)
                end
            end
        end
        # The defect of `gfsntc_c` that Almanac keeps: it checks `abcorr` twice and never
        # `obsrvr`, so an empty observer goes on to the lookup of the body.
        @testset "gfsntc: an empty observer is a body that is not found" begin
            ks = Almanac.kernelset(SpiceContext())
            f = () -> gfsntc(
                ks, "EARTH", "IAU_EARTH", "ELLIPSOID", "NONE", "", "J2000", [1.0, 0.0, 0.0],
                "LATITUDINAL", "LATITUDE", "=", 0.0, 0.0, 3600.0, Window([0.0, 1.0])
            )
            @test emptystring_kind(f, ()) === :idcodenotfound
        end
        # `gfevnt_c` checks the number of parameters before its strings.
        @testset "gfevnt: too many parameters before an empty quantity" begin
            ks = Almanac.kernelset(SpiceContext())
            f = () -> gfevnt(
                ks, NaN, NaN, "", fill("X", 11), String[], Float64[], Int[], Bool[], "X",
                NaN, NaN, NaN, Window([0.0, 1.0])
            )
            @test emptystring_kind(f, ()) === :invalidcount
        end
        # Regressions: `timdef_c` compares the action with `eqstr_c`, which ignores
        # every blank; a blank shape is a shape that is not recognized, not an empty
        # string; an empty unit in a kernel is a unit that is not recognized.
        @testset "blank and kernel strings are not empty arguments" begin
            ctx = SpiceContext()
            @test emptystring_kind(() -> timdef(ctx, "S ET", "CALENDAR", ""), ()) === :emptystring
            ks = Almanac.kernelset(ctx)
            f = () -> occult(ks, "MOON", "   ", "IAU_MOON", "SUN", "ELLIPSOID", "IAU_SUN", "LT", "EARTH", 1.0e8)
            @test emptystring_kind(f, ()) === :badmethodsyntax
            for (nm, v) in (
                    ("INS-999_FOV_SHAPE", "CIRCLE"), ("INS-999_FOV_FRAME", "J2000"),
                    ("INS-999_FOV_CLASS_SPEC", "ANGLES"), ("INS-999_FOV_ANGLE_UNITS", ""),
                )
                pcpool(ctx, nm, [v])
            end
            for (nm, v) in (
                    ("INS-999_BORESIGHT", [0.0, 0.0, 1.0]),
                    ("INS-999_FOV_REF_VECTOR", [1.0, 0.0, 0.0]), ("INS-999_FOV_REF_ANGLE", [1.0]),
                )
                pdpool(ctx, nm, v)
            end
            @test emptystring_kind(() -> getfov(ctx, -999), ()) === :unitsnotrec
        end
    end
end
