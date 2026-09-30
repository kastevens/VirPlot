"""Generate examples/sample.sam: synthetic reads over a two-segment genome.

SyntheticVirus1 (8000 bp) gets most reads, shaped so the plot is interesting
(a peak and a gap); SyntheticVirus2 (3000 bp, see sample_multi.gff3) gets a
sparser, flatter set so the two-RNA example has something to compare. A
handful of records exercise filters: secondary, duplicate, unmapped, low MAPQ.

Deterministic; run as  python examples/make_sample_sam.py examples/sample.sam
"""
import random, sys

random.seed(20260929)
L = 8000
out = sys.argv[1]
reads = []

def seq(n): return "".join(random.choice("ACGT") for _ in range(n))

def rec(name, flag, ref, pos, mapq, cigar, qlen):
    return f"{name}\t{flag}\t{ref}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq(qlen)}\t{'I'*qlen}"

i = 0  # noqa
def next_name():
    global i; i += 1; return f"read{i:05d}"

# background coverage, weighted towards the middle, with a hole at 5200-5400
for _ in range(600):
    while True:
        pos = int(random.triangular(1, L - 150, 3000))
        if not (5000 < pos < 5400):
            break
    r = random.random()
    if r < 0.75:
        cigar, qlen = "150M", 150
    elif r < 0.85:
        cigar, qlen = "10S140M", 150            # soft clip: 140 ref bases
    elif r < 0.92:
        cigar, qlen = "70M3I77M", 150           # insertion: 147 ref bases
    elif r < 0.97:
        cigar, qlen = "60M5D90M", 150           # deletion: 155 ref span, 150 counted
    else:
        cigar, qlen = "50M400N100M", 150        # spliced: two blocks
    reads.append(rec(next_name(), random.choice([0, 16]), "SyntheticVirus1", pos, 60, cigar, qlen))

# records that must NOT count
reads.append(rec(next_name(), 256, "SyntheticVirus1", 100, 60, "150M", 150))    # secondary
reads.append(rec(next_name(), 1024, "SyntheticVirus1", 200, 60, "150M", 150))   # duplicate
reads.append(rec(next_name(), 512, "SyntheticVirus1", 300, 60, "150M", 150))    # QC fail
reads.append(rec(next_name(), 4, "*", 0, 0, "*", 150))                          # unmapped
# second segment: ~1/4 the depth of segment 1, one shallow dip around 1500
for _ in range(90):
    while True:
        pos = random.randint(1, 3000 - 150)
        if not (1450 < pos < 1550) or random.random() < 0.3:
            break
    reads.append(rec(next_name(), random.choice([0, 16]), "SyntheticVirus2", pos, 60, "150M", 150))
# low MAPQ reads: counted by default, dropped with --min-mapq 20
for _ in range(40):
    reads.append(rec(next_name(), 0, "SyntheticVirus1", random.randint(7000, 7800), 3, "150M", 150))
# a read hanging off the 3' end (clipped to the genome length)
reads.append(rec(next_name(), 0, "SyntheticVirus1", 7950, 60, "150M", 150))

random.shuffle(reads)
with open(out, "w") as fp:
    fp.write("@HD\tVN:1.6\tSO:unsorted\n")
    fp.write("@SQ\tSN:SyntheticVirus1\tLN:8000\n")
    fp.write("@SQ\tSN:SyntheticVirus2\tLN:3000\n")
    fp.write("@PG\tID:synth\tPN:make_sample_sam.py\n")
    fp.write("\n".join(reads) + "\n")
print(len(reads), "records")
