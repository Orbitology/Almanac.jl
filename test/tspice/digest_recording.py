"""Write the digest of a recording of record.sh to stdout: one line per file of DIR/logs,
sorted, with its path relative to DIR/logs and its SHA-256, separated by a tab. Two
recordings are the same when their digests are equal:

    python3 -I test/tspice/digest_recording.py DIR > digest.tsv
    diff <(python3 -I test/tspice/digest_recording.py A) <(python3 -I test/tspice/digest_recording.py B)

The digest of a file is that of its content, with these exceptions:

- a log (.log.gz) is digested as the replay reads it, decompressed and without the
  output arguments of the calls that failed (`! T`, at any depth), which CSPICE leaves
  partly written; a binary file of values (`B` line) that only such an output names is
  left out too;
- another .gz file is digested decompressed, because two versions of gzip can compress
  the same data differently;
- in the `S` and `B` lines of a log, the path of a snapshot up to its directory
  NAME.log.files is replaced by `<logs>/NAME.log.files`, as replay.jl finds the snapshots
  of a recording that moved, so that the digest does not depend on the directory in
  which the recording was made;
- a symbolic link gives `-> TARGET`;
- status.tsv is sorted, because the families end in any order.

The files that record the run, not its calls, are left out: output.txt and cspice*.log
in the snapshot directories (*.log.files/), which hold the output and the dates of the
family, and the files .last.* of the logger in *.log.files/latest/, which hold the
absolute paths of copies. A FIFO in DIR/logs is the log of a family that did not end,
and is an error. JOBS sets the number of parallel processes (default: the CPU threads).
"""

import fnmatch
import gzip
import hashlib
import multiprocessing
import os
import re
import stat
import sys

SNAPSHOT_PATH = re.compile(rb"^([SB]) .*/([^/]+\.log\.files/)")
BIN_LINE = b"B <logs>/"
LEFT_OUT = ("*.log.files/*/output.txt", "*.log.files/*/cspice*.log", "*.log.files/latest/.last.*")


def left_out(rel):
    """True if the file at the path `rel` records the run, not its calls (LEFT_OUT)."""
    return any(fnmatch.fnmatchcase(rel, p) for p in LEFT_OUT)


class Lines:
    """The lines of a log without their newline, with one line of look-ahead."""

    def __init__(self, f):
        self.it = iter(f)
        self.pending = None

    def next(self):
        if self.pending is not None:
            l, self.pending = self.pending, None
            return l
        l = next(self.it, None)
        return l[:-1] if l is not None and l.endswith(b"\n") else l

    def back(self, l):
        assert self.pending is None and l is not None
        self.pending = l


def integer(b):
    try:
        return int(b)
    except ValueError:
        return None


def read_call(lines, head, bins):
    """The lines of the call that begins with `head`, in the order of the log, without its
    outputs if it failed; a deeper call that begins before it ends is read first, as
    foreach_call of replay.jl reads it. `bins` collects the binary files that the kept and
    the left-out `B` lines name."""
    assert head.startswith(b"@ ") and integer(head.split()[2]) is not None
    depth = int(head.split()[2])
    items = [head]  # a line, or (is_output, lines of an argument)
    binline = None  # a `B` line: the values of the next argument of n < 0
    end = None
    while True:
        l = lines.next()
        if l is None:
            break
        if l.startswith(b"! "):
            end = l
            break
        if l.startswith(b"@ "):
            q = l.split()
            d = integer(q[2]) if len(q) == 4 else None
            if d is not None and d > depth:
                items.extend(read_call(lines, l, bins))
                continue
            lines.back(l)  # a call that begins before this one ended: this one ends here
            break
        if l[:1] in (b"<", b">"):
            p = l.split()
            n = integer(p[2]) if len(p) == 3 else None
            if n is not None:
                arg = [l] + [v for v in (lines.next() for _ in range(max(n, 0))) if v is not None]
                if n < 0 and binline is not None:
                    arg.insert(0, binline)
                    binline = None
                items.append((l[:1] == b">", arg))
                continue
        if l.startswith(b"S ") or l.startswith(b"B "):
            l = SNAPSHOT_PATH.sub(rb"\1 <logs>/\2", l)
        if l.startswith(b"B "):
            if binline is not None:
                items.append(binline)
            binline = l
            continue
        items.append(l)
    if binline is not None:
        items.append(binline)
    failed = end is not None and end.split()[1:2] == [b"T"]
    out = []
    for it in items:
        if isinstance(it, tuple):
            is_output, arg = it
            kept = not (is_output and failed)
            if arg[0].startswith(BIN_LINE):
                bins[kept].add(arg[0][len(BIN_LINE) :].decode("utf-8", "surrogateescape"))
            if kept:
                out.extend(arg)
        else:
            out.append(it)
            if it.startswith(BIN_LINE):
                bins[True].add(it[len(BIN_LINE) :].decode("utf-8", "surrogateescape"))
    if end is not None:
        out.append(end)
    assert out[0] == head and (end is None or out[-1] == end)
    return out


def digest_log(path):
    """The SHA-256 of the log at `path` as the replay reads it, and the binary files that
    its kept and its left-out `B` lines name."""
    assert path.endswith(".log.gz")
    h = hashlib.sha256()
    bins = {True: set(), False: set()}
    with gzip.open(path, "rb") as f:
        lines = Lines(f)
        while True:
            l = lines.next()
            if l is None:
                break
            if l.startswith(b"@ ") and len(l.split()) == 4 and integer(l.split()[2]) is not None:
                h.update(b"\n".join(read_call(lines, l, bins)) + b"\n")
            else:
                h.update(SNAPSHOT_PATH.sub(rb"\1 <logs>/\2", l) + b"\n")
    return h.hexdigest(), bins


def digest_file(rel):
    """The SHA-256 of the file at the path `rel`, by the rules of the module docstring."""
    assert not rel.endswith(".log.gz") and not left_out(rel)
    h = hashlib.sha256()
    if rel == "status.tsv":
        with open(rel, "rb") as f:
            h.update(b"".join(l + b"\n" for l in sorted(f.read().splitlines())))
    else:
        with gzip.open(rel, "rb") if rel.endswith(".gz") else open(rel, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
    return h.hexdigest()


def task(job):
    kind, rel = job
    try:
        return rel, digest_log(rel) if kind == "log" else (digest_file(rel), None)
    except (OSError, EOFError, ValueError) as e:
        raise RuntimeError(f"{rel}: {e}") from None


def walk(top):
    """The paths below `top`, relative to it: the files and the symbolic links."""
    files, links = [], []
    stack = [""]
    while stack:
        d = stack.pop()
        with os.scandir(os.path.join(top, d) if d else top) as it:
            for e in it:
                rel = d + "/" + e.name if d else e.name
                mode = e.stat(follow_symlinks=False).st_mode
                if stat.S_ISLNK(mode):
                    links.append(rel)
                elif stat.S_ISDIR(mode):
                    stack.append(rel)
                elif stat.S_ISREG(mode):
                    files.append(rel)
                elif stat.S_ISFIFO(mode):
                    sys.exit(f"{os.path.join(top, rel)} is a FIFO: its family did not end")
    assert not set(files) & set(links)
    return files, links


def main(args):
    if len(args) != 1:
        sys.exit("usage: digest_recording.py DIR")
    logs = os.path.join(args[0], "logs")
    if not os.path.isdir(logs):
        sys.exit(f"{args[0]} is not a recording of record.sh: it has no directory logs")
    files, links = walk(logs)
    os.chdir(logs)
    result = {rel: "-> " + os.readlink(rel) for rel in links}
    jobs = int(os.environ.get("JOBS") or os.cpu_count() or 1)
    with multiprocessing.get_context("fork").Pool(jobs) as pool:
        try:
            # the logs first: they name the binary files to leave out
            logfiles = [rel for rel in files if rel.endswith(".log.gz")]
            kept, dropped = set(), set()
            for rel, (h, bins) in pool.imap_unordered(task, [("log", rel) for rel in logfiles]):
                result[rel] = h
                kept |= bins[True]
                dropped |= bins[False]
            others = [
                rel
                for rel in files
                if not rel.endswith(".log.gz") and not left_out(rel) and not (rel in dropped and rel not in kept)
            ]
            for rel, (h, _) in pool.imap_unordered(task, [("file", rel) for rel in others], chunksize=64):
                result[rel] = h
        except RuntimeError as e:
            sys.exit(str(e))
    out = sys.stdout.buffer
    for rel in sorted(result, key=os.fsencode):
        out.write(os.fsencode(rel) + b"\t" + result[rel].encode() + b"\n")


if __name__ == "__main__":
    main(sys.argv[1:])
