#!/usr/bin/env python3
"""Generate fuzz_verify seed corpus from KAT vectors (pk || sm)."""
import os, re

KAT = os.path.join(os.path.dirname(__file__), "..", "KAT", "PQCsignKAT_353_SQIsign_lvl1.rsp")
OUT = os.path.join(os.path.dirname(__file__), "corpus", "fuzz_verify")
os.makedirs(OUT, exist_ok=True)

with open(KAT) as f:
    text = f.read()

blocks = re.split(r"\n\n", text)
n = 0
for blk in blocks:
    pk = re.search(r"^pk = ([0-9A-Fa-f]+)$", blk, re.M)
    sm = re.search(r"^sm = ([0-9A-Fa-f]+)$", blk, re.M)
    if pk and sm:
        data = bytes.fromhex(pk.group(1)) + bytes.fromhex(sm.group(1))
        with open(os.path.join(OUT, f"kat_{n:03d}"), "wb") as out:
            out.write(data)
        n += 1
        if n >= 5:
            break
print(f"wrote {n} seeds to {OUT}")
