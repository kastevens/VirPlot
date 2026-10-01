#!/usr/bin/env python3
"""Simulate a SAM of evenly spread reads over every molecule in a GFF3.

    python examples/make_example_sam.py GFF OUT.sam [--depth 20] [--read 150]

One @SQ line per GFF `region` row (name = seqid, length = region end), reads
spread with mild noise over each molecule in proportion to its length. For a
region with Is_circular=true, reads may start anywhere and run on past the end
(wrap-aware coordinates, as after padding the reference and wrapping back).
SEQ/QUAL are '*'. Deterministic. The single-reference make_reads.py does the
same for one molecule given on the command line.
"""
import argparse
import random

ap = argparse.ArgumentParser()
ap.add_argument("gff"); ap.add_argument("out")
ap.add_argument("--depth", type=float, default=20); ap.add_argument("--read", type=int, default=150)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
random.seed(a.seed)

regions = []
with open(a.gff) as fp:
    for line in fp:
        parts = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(parts) != 9 or parts[2] != "region":
            continue
        attrs = dict(kv.split("=", 1) for kv in parts[8].split(";") if "=" in kv)
        regions.append((parts[0], int(parts[4]), attrs.get("Is_circular", "").lower() == "true"))

head = ["@HD\tVN:1.6\tSO:unsorted"] + [f"@SQ\tSN:{s}\tLN:{n}" for s, n, _ in regions]
head.append("@PG\tID:simulate\tPN:make_example_sam.py")
recs, k = [], 0
for seqid, length, circular in regions:
    n = round(length * a.depth / a.read)
    hi = length if circular else max(1, length - a.read + 1)
    for _ in range(n):
        recs.append(f"r{k:06d}\t{random.choice([0, 16])}\t{seqid}\t{random.randint(1, hi)}"
                    f"\t60\t{a.read}M\t*\t0\t0\t*\t*")
        k += 1
with open(a.out, "w") as fp:
    fp.write("\n".join(head + recs) + "\n")
print(f"{k} reads over {len(regions)} molecule(s) -> {a.out}")
