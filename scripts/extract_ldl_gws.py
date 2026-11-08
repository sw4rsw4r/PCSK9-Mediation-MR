"""Stream GLGC LDL-C and keep variants with p < 5e-8 that are in the 1000G EUR panel."""
from pathlib import Path

from config import BFILE, GLGC as LDL, LDL_GWS_DIR as OUT

BIM = Path(str(BFILE) + ".bim")
OUT.mkdir(parents=True, exist_ok=True)
P_CUT = 5e-8

print("Loading 1000G EUR ids", flush=True)
ref = set()
with BIM.open() as handle:
    for line in handle:
        ref.add(line.split()[1])
print(f"EUR variants: {len(ref)}", flush=True)

kept = {}
n_rows = 0
n_gws = 0
n_not_in_ref = 0
n_bad = 0
with LDL.open(encoding="utf-8", errors="replace") as handle:
    header = handle.readline().rstrip("\n").split("\t")
    idx = {name: i for i, name in enumerate(header)}
    need = [
        "variant_id", "chromosome", "base_pair_location", "other_allele",
        "effect_allele", "n", "effect_allele_frequency", "beta",
        "standard_error", "p_value",
    ]
    col = [idx[name] for name in need]
    for line in handle:
        n_rows += 1
        if n_rows % 5000000 == 0:
            print(f"  read {n_rows:,}  gws {n_gws:,}  kept {len(kept):,}", flush=True)
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= col[-1]:
            n_bad += 1
            continue
        p_raw = parts[col[9]].strip()
        try:
            p_val = float(p_raw)
        except ValueError:
            n_bad += 1
            continue
        if p_val >= P_CUT:
            continue
        n_gws += 1
        rsid = parts[col[0]].strip()
        if rsid not in ref:
            n_not_in_ref += 1
            continue
        prev = kept.get(rsid)
        if prev is not None and float(prev[9]) <= p_val:
            continue
        chrom = parts[col[1]].strip()
        if chrom in {"X", "Y", "MT", "M"}:
            continue
        kept[rsid] = [
            rsid,
            chrom,
            parts[col[2]].strip(),
            parts[col[3]].strip(),
            parts[col[4]].strip(),
            parts[col[5]].strip(),
            parts[col[6]].strip(),
            parts[col[7]].strip(),
            parts[col[8]].strip(),
            p_raw,
        ]

stats = OUT / "ldl_gws.tsv"
clump = OUT / "ldl_gws_for_clump.txt"
with stats.open("w", encoding="utf-8", newline="\n") as out, clump.open("w", encoding="utf-8", newline="\n") as cl:
    out.write("rsid\tchromosome\tposition\tother_allele\teffect_allele\tn\teaf\tbeta\tse\tp\n")
    cl.write("SNP\tP\n")
    for row in kept.values():
        out.write("\t".join(row) + "\n")
        cl.write(f"{row[0]}\t{row[9]}\n")

print(
    f"rows={n_rows:,} gws={n_gws:,} not_in_EUR={n_not_in_ref:,} "
    f"bad={n_bad:,} kept={len(kept):,}",
    flush=True,
)
