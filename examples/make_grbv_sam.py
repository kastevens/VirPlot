#!/usr/bin/env python3
"""Generate examples/grbv.sam: simulated reads over the GRBV circular genome.

Grapevine red blotch virus (RefSeq NC_022002.1) is a 3206 nt circular ssDNA
geminivirus. Coverage here is simulated, not real: it is deliberately even
(as rolling-circle replication tends to give) and, importantly, *continuous
across the origin* — reads starting near position 3206 run on into position 1.

Those origin-spanning reads are written the way a wrap-aware pipeline emits
them: POS stays within 1..3206 and the CIGAR runs past the end of the
reference. That is what you get by padding the reference by a read length,
aligning with bowtie2, and wrapping the coordinates back (bowtie2 itself has
no notion of a circular reference, so without that step these reads are
soft-clipped or lost and the depth plot shows a false dip at the origin).

VirPlot reads SEQ/QUAL not at all, so both are '*' rather than inventing
sequence we do not have.

Deterministic.  Usage:  python examples/make_grbv_sam.py examples/grbv.sam
"""
import math
import random
import sys

REF = "NC_022002.1"
L = 3206
READ = 150
TARGET_DEPTH = 30

random.seed(20260929)

out = sys.argv[1] if len(sys.argv) > 1 else "examples/grbv.sam"
n_reads = round(L * TARGET_DEPTH / READ)

records = []


def emit(name, flag, pos, mapq, cigar):
    records.append(f"{name}\t{flag}\t{REF}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t*\t*")


# Even coverage with a gentle undulation, sampled uniformly around the whole
# circle: start positions run 1..L, so reads near the end overrun the origin.
i = 0
for _ in range(n_reads):
    while True:
        pos = random.randint(1, L)
        # mild sinusoidal variation (~±25%) so the trace is not a flat line
        if random.random() < 0.75 + 0.25 * math.sin(2 * math.pi * pos / L):
            break
    i += 1
    flag = random.choice([0, 16])
    r = random.random()
    if r < 0.88:
        cigar = f"{READ}M"
    elif r < 0.94:
        cigar = f"8S{READ - 8}M"          # soft clip: 142 ref bases
    elif r < 0.98:
        cigar = f"70M2I{READ - 72}M"      # small insertion
    else:
        cigar = f"60M3D{READ - 60}M"      # small deletion
    emit(f"grbv{i:05d}", flag, pos, 60, cigar)

# a few records that must not be counted (filter check)
emit("grbv_sec", 256, 500, 60, f"{READ}M")
emit("grbv_dup", 1024, 900, 60, f"{READ}M")
emit("grbv_unmapped", 4, 0, 0, "*")

random.shuffle(records)

with open(out, "w") as fp:
    fp.write("@HD\tVN:1.6\tSO:unsorted\n")
    fp.write(f"@SQ\tSN:{REF}\tLN:{L}\n")
    fp.write("@PG\tID:simulate\tPN:make_grbv_sam.py\tDS:simulated reads, wrap-aware coordinates\n")
    fp.write("\n".join(records) + "\n")

spanning = sum(1 for r in records if int(r.split("\t")[3]) > L - READ)
print(f"{len(records)} records, {spanning} crossing the origin -> {out}")
