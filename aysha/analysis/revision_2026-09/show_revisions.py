#!/usr/bin/env python
"""Report prose revisions to a markdown manuscript as a reviewable changeset.

The manuscript is one paragraph per line, so a line diff marks whole paragraphs
as changed and is useless for review. This runs a word-level diff and prints
only the changed spans with surrounding context, grouped under the section each
change falls in, and flags the two classes of edit that need re-verification:
a changed numeric value, and a changed bracketed citation marker.

Run from the directory holding the manuscript.

Usage:
    python show_revisions.py                     # HEAD vs working tree
    python show_revisions.py --old <ref>         # <ref> vs working tree
    python show_revisions.py --old A --new B     # commit A vs commit B
    python show_revisions.py --file other.md
"""
from __future__ import annotations
import argparse
import re
import subprocess
import sys

CTX = 100
CITE = re.compile(r"\[\d+(?:\s*,\s*\d+)*\]")
NUMTOK = re.compile(r"-?\d+(?:\.\d+)*")   # multi-dot versions match as one token
# tokens that are navigation or software versions, not measured quantities
SKIP_BEFORE = re.compile(
    r"(?:section|sections|figure|figures|table|tables|appendix|panel|ref|refs|"
    r"python|numpy|scipy|pandas|matplotlib|version|v)\W*$", re.I)
CODE_SPAN = re.compile(r"`[^`]*`")


def data_numbers(text):
    """Numeric tokens that plausibly denote measured or derived quantities.

    Drops anything inside an inline code span (filenames, keys) and anything
    immediately preceded by a navigation or package word, so a reworded
    cross-reference does not masquerade as a changed result.
    """
    clean = CODE_SPAN.sub(" ", text)
    out = set()
    for m in NUMTOK.finditer(clean):
        if SKIP_BEFORE.search(clean[max(0, m.start() - 24):m.start()]):
            continue
        out.add(m.group(0))
    return out


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True)


def word_diff(path, old, new):
    cmd = ["git", "diff", "--word-diff=plain", "--unified=0", old]
    if new:
        cmd.append(new)
    cmd += ["--", path]
    r = run(cmd)
    if r.returncode != 0:
        sys.exit("git diff failed:\n" + r.stderr)
    return r.stdout


def section_map(path, old):
    """Line number -> nearest preceding heading, taken from the OLD revision."""
    r = run(["git", "show", old + ":./" + path])
    text = r.stdout if r.returncode == 0 else open(path, encoding="utf-8").read()
    out, cur = {}, "(front matter)"
    for i, line in enumerate(text.split("\n"), start=1):
        if line.startswith("#"):
            cur = line.lstrip("#").strip()
        out[i] = cur
    return out


def merge(spans, n):
    """Merge overlapping context windows so adjacent edits print once."""
    out = []
    for a, b in spans:
        a, b = max(0, a - CTX), min(n, b + CTX)
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="manuscript_revised_v3.md")
    ap.add_argument("--old", default="HEAD")
    ap.add_argument("--new", default=None)
    args = ap.parse_args()

    raw = word_diff(args.file, args.old, args.new)
    target = args.new or "working tree"
    if not raw.strip() and args.new is None and args.old == "HEAD":
        # Nothing uncommitted: the revisions were committed, so fall back to the
        # last two commits that touched this file (2026-09-08).
        revs = run(["git", "log", "--format=%H", "--", args.file]).stdout.split()
        if len(revs) >= 2:
            args.old, args.new = revs[1], revs[0]
            target = args.new
            raw = word_diff(args.file, args.old, args.new)
            print("Working tree is clean; comparing the last two commits of this file.\n")
    if not raw.strip():
        print("No changes to %s (%s -> %s)." % (args.file, args.old, target))
        return
    sects = section_map(args.file, args.old)

    heads = re.findall(r"^@@ (.*?) @@", raw, flags=re.M)
    bodies = re.split(r"^@@ .*?@@", raw, flags=re.M)[1:]
    n_num = n_cite = n_para = 0
    print("%s: %d changed paragraph(s), %s -> %s\n"
          % (args.file, len(bodies), args.old, target))

    for head, body in zip(heads, bodies):
        m = re.match(r"-(\d+)", head.strip())
        ln = int(m.group(1)) if m else 0
        marks = [(x.start(), x.end()) for x in
                 re.finditer(r"(?:\[-.*?-\]|\{\+.*?\+\})+", body, re.S)]
        if not marks:
            continue
        n_para += 1
        removed = " ".join(re.findall(r"\[-(.*?)-\]", body, re.S))
        added = " ".join(re.findall(r"\{\+(.*?)\+\}", body, re.S))
        num_hit = sorted(data_numbers(removed) ^ data_numbers(added))
        cite_hit = sorted(set(CITE.findall(removed)) ^ set(CITE.findall(added)))
        n_num += len(num_hit)
        n_cite += len(cite_hit)

        print("--- line %d | %s" % (ln, sects.get(ln, "?")))
        if num_hit:
            print("    NUMBERS: %s   <- re-verify against results.json"
                  % ", ".join(num_hit))
        if cite_hit:
            print("    CITATIONS: %s   <- re-check the reference list"
                  % ", ".join(cite_hit))
        for a, b in merge(marks, len(body)):
            print("    ...%s..." % re.sub(r"\s+", " ", body[a:b]).strip())
        print()

    print("Totals: %d changed paragraph(s), %d numeric token(s), %d citation marker(s)."
          % (n_para, n_num, n_cite))

    notes = [(i, ln.strip()) for i, ln in
             enumerate(open(args.file, encoding="utf-8").read().split("\n"), 1)
             if "<!--" in ln]
    if notes:
        print("\nInline notes to the reviewer (%d):" % len(notes))
        for i, ln in notes:
            print("  line %d: %s" % (i, ln[:200]))


if __name__ == "__main__":
    main()
