#!/usr/bin/env python3
"""Simulate a SAM of evenly spread reads over one reference, for examples.

    python examples/make_reads.py REFNAME LENGTH OUT.sam [--depth 20] [--circular]

Coverage is flat with mild noise; with --circular, reads starting near the
end run on past it (wrap-aware coordinates, as after padding + wrapping).
SEQ/QUAL are '*'. Deterministic.
"""
import argparse
import random

ap = argparse.ArgumentParser()
ap.add_argument("ref"); ap.add_argument("length", type=int); ap.add_argument("out")
ap.add_argument("--depth", type=float, default=20); ap.add_argument("--read", type=int, default=150)
ap.add_argument("--circular", action="store_true"); ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
random.seed(a.seed)
n = round(a.length * a.depth / a.read)
hi = a.length if a.circular else max(1, a.length - a.read + 1)
recs = [f"r{i:06d}\t{random.choice([0, 16])}\t{a.ref}\t{random.randint(1, hi)}\t60\t{a.read}M\t*\t0\t0\t*\t*"
        for i in range(n)]
with open(a.out, "w") as fp:
    fp.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{a.ref}\tLN:{a.length}\n"
             "@PG\tID:simulate\tPN:make_reads.py\n" + "\n".join(recs) + "\n")
print(f"{n} reads -> {a.out}")
