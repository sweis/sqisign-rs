#!/usr/bin/env python3
"""Rewrite the `let s = X; fp2_op(&mut X, &s, ...)` aliasing-temp pattern
into in-place `fp2_op_ip(&mut X, ...)` calls."""

import re
import sys

# Match `let s = EXPR;` + `OP(&mut EXPR, ...);` either on the same line
# or on consecutive lines.
PAT = re.compile(
    r"^(?P<ind>[ \t]*)let s = (?P<expr>[^;]+);"
    r"(?P<sep>\s*)"
    r"(?P<op>fp2_(?:mul|sqr|add|sub|neg|half)|hadamard|pointwise_square|to_squared_theta)"
    r"\((?P<args>[^;]+)\);\n",
    re.MULTILINE,
)


def rewrite_one(m: re.Match) -> str:
    ind = m["ind"]
    expr = m["expr"]
    op = m["op"]
    args = [a.strip() for a in split_args(m["args"])]
    out = args[0]
    rest = args[1:]
    # Output may be `&mut X` (X==expr) or just `X` when expr is `*X`
    # (e.g. `let s = *p; op(p, &s)` where p: &mut T).
    if out == f"&mut {expr}":
        out_expr = out
    elif expr.startswith("*") and out == expr[1:]:
        out_expr = out  # already `p` which is &mut
    else:
        return m.group(0)  # output isn't the captured var — leave alone

    s = "&s"
    def emit(ip_op: str, *extra: str) -> str:
        if extra:
            return f"{ind}{ip_op}({out_expr}, {', '.join(extra)});\n"
        return f"{ind}{ip_op}({out_expr});\n"
    if op in ("hadamard", "pointwise_square", "to_squared_theta"):
        if rest == [s]:
            return emit(f"{op}_ip")
        return m.group(0)
    if op == "fp2_sqr":
        if rest == [s]:
            return emit("fp2_sqr_ip")
    elif op == "fp2_neg":
        if rest == [s]:
            return emit("fp2_neg_ip")
    elif op == "fp2_half":
        return m.group(0)
    elif op == "fp2_mul":
        if rest == [s, s]:
            return emit("fp2_sqr_ip")
        if rest[0] == s:
            return emit("fp2_mul_ip", rest[1])
        if rest[1] == s:
            return emit("fp2_mul_ip", rest[0])
    elif op == "fp2_add":
        if rest == [s, s]:
            return emit("fp2_dbl_ip")
        if rest[0] == s:
            return emit("fp2_add_ip", rest[1])
        if rest[1] == s:
            return emit("fp2_add_ip", rest[0])
    elif op == "fp2_sub":
        if rest == [s, s]:
            return m.group(0)
        if rest[0] == s:
            return emit("fp2_sub_ip", rest[1])
        if rest[1] == s:
            return emit("fp2_rsub_ip", rest[0])
    return m.group(0)


def split_args(s: str) -> list[str]:
    """Split top-level commas (no nested parens in our case, but be safe)."""
    out, depth, cur = [], 0, []
    for c in s:
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        if c == "," and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(c)
    out.append("".join(cur))
    return out


def process(path: str) -> int:
    with open(path) as f:
        src = f.read()
    new = src
    # Iterate until fixpoint (a rewrite can expose another adjacent pattern).
    while True:
        nxt = PAT.sub(rewrite_one, new)
        if nxt == new:
            break
        new = nxt
    if new != src:
        with open(path, "w") as f:
            f.write(new)
    # Count rewrites by line shrinkage (each rewrite removes one `let s =` line).
    return src.count("let s = ") - new.count("let s = ")


if __name__ == "__main__":
    total = 0
    for p in sys.argv[1:]:
        n = process(p)
        print(f"{p}: {n} rewrites")
        total += n
    print(f"total: {total}")
