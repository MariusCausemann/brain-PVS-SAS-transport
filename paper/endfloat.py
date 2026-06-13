#!/usr/bin/env python3
r"""
Move all float environments (figure, table, and their starred forms) to the end
of a LaTeX document, just before \end{document}. Leaves a placeholder comment
where each float was, so references and surrounding text are untouched.

Usage:
    python3 endfloat_source.py master.tex submission.tex
    python3 endfloat_source.py master.tex submission.tex --envs figure table algorithm
    python3 endfloat_source.py master.tex submission.tex --clearpage --markers

Maintain ONE source (master.tex). Regenerate submission.tex on demand.
"""

import argparse
import re
import sys

DEFAULT_ENVS = ["figure", "table"]


def line_is_comment_start(text, pos):
    """True if an unescaped % appears before `pos` on the same line."""
    line_start = text.rfind("\n", 0, pos) + 1
    i = line_start
    while i < pos:
        if text[i] == "%" and (i == 0 or text[i - 1] != "\\"):
            return True
        i += 1
    return False


def find_floats(text, envs):
    """Return list of (start, end, env, body) for top-level float environments,
    skipping any that begin inside a comment. Handles nesting of the same env."""
    names = "|".join(re.escape(e) for e in envs)
    # match \begin{env} or \begin{env*}
    token = re.compile(r"\\(begin|end)\{(" + names + r")(\*?)\}")
    floats = []
    stack = []
    for m in token.finditer(text):
        if line_is_comment_start(text, m.start()):
            continue
        kind, base, star = m.group(1), m.group(2), m.group(3)
        if kind == "begin":
            stack.append((m.start(), base + star))
        else:
            if stack and stack[-1][1] == base + star:
                start, env = stack.pop()
                if not stack:  # only record top-level floats
                    floats.append((start, m.end(), env, text[start:m.end()]))
    return floats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--envs", nargs="+", default=DEFAULT_ENVS,
                    help="environment names to relocate (default: figure table)")
    ap.add_argument("--clearpage", action="store_true",
                    help="insert \\clearpage before each relocated float")
    ap.add_argument("--markers", action="store_true",
                    help="leave a visible [Figure/Table about here] marker, not just a comment")
    ap.add_argument("--figures-first", action="store_true",
                    help="group all figures before all tables instead of source order")
    args = ap.parse_args()

    with open(args.infile, encoding="utf-8") as f:
        text = f.read()

    floats = find_floats(text, args.envs)
    if not floats:
        print("No floats found; copying unchanged.", file=sys.stderr)

    # Build the relocated block and placeholders first, then cut the floats
    # from the body in REVERSE order so earlier offsets stay valid.
    collected = []
    edits = []  # (start, end, placeholder)
    for n, (start, end, env, body) in enumerate(floats, 1):
        base = env.rstrip("*")
        if args.markers:
            placeholder = "\n%% [{} about here]\n[{} {} about here]\n".format(
                base.capitalize(), base.capitalize(), n)
        else:
            placeholder = "\n%% {} moved to end of document\n".format(env)
        collected.append((base, body))
        edits.append((start, end, placeholder))

    for start, end, placeholder in reversed(edits):
        text = text[:start] + placeholder + text[end:]

    if args.figures_first:
        collected.sort(key=lambda x: 0 if x[0] == "figure" else 1)

    block_parts = ["\n\n%% ===== Floats relocated to end of document =====\n"]
    for base, body in collected:
        if args.clearpage:
            block_parts.append("\\clearpage\n")
        block_parts.append(body + "\n\n")
    block = "".join(block_parts)

    # Insert before \end{document}; if absent, append at the very end.
    m = re.search(r"\\end\{document\}", text)
    if m:
        text = text[:m.start()] + block + text[m.start():]
    else:
        text = text + block

    with open(args.outfile, "w", encoding="utf-8") as f:
        f.write(text)

    print("Relocated {} float(s) -> {}".format(len(collected), args.outfile),
          file=sys.stderr)


if __name__ == "__main__":
    main()