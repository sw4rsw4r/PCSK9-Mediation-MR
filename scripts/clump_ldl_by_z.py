"""Genome-wide LDL-C clump for the PCSK9 mediation analysis, ranked by |z|.

Input is results/ldl_gws/ldl_gws.tsv from extract_ldl_gws.py. The thresholds
are p < 5e-8, r2 < 0.001 and 10 Mb in the 1000 Genomes European reference.
GLGC stores P values that become 0 in double precision for 827 variants with
|z| above about 38.5, so ranking by P leaves the index variant at the strongest
loci to file order. PLINK receives 10^-(|z| + 3), which orders by |z| and keeps
every variant below the p1/p2 threshold.

Variants that cannot be aligned to CAD are removed before clumping, so a locus
is not left without an instrument when its strongest variant is unusable.
"""
import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import extract_pcsk9_setting1 as base
from config import BFILE, CLUMP_DIR as OUT, LDL_GWS as LDL, PLINK

OUT.mkdir(parents=True, exist_ok=True)


def main():
    with LDL.open(encoding="utf-8") as src:
        rows = list(csv.DictReader(src, delimiter="\t"))
    print(f"LDL-C p < 5e-8 {len(rows):,}", flush=True)
    print("streaming CAD", flush=True)
    cad = base.stream_cad({row["rsid"] for row in rows})

    usable, dropped = [], []
    for row in rows:
        if base.harmonize_to_ldl(row, cad.get(row["rsid"], [])) is None:
            dropped.append(row)
        else:
            usable.append(row)
    print(f"alignable with CAD {len(usable):,} not alignable {len(dropped):,}", flush=True)
    with (OUT / "not_alignable_with_cad.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("rsid\tchromosome\tposition\teffect_allele\tother_allele\teaf\tabs_z\tin_cad\n")
        for row in dropped:
            z = abs(float(row["beta"]) / float(row["se"]))
            handle.write("\t".join([
                row["rsid"], row["chromosome"], row["position"], row["effect_allele"],
                row["other_allele"], row["eaf"], f"{z:.2f}", str(bool(cad.get(row["rsid"]))),
            ]) + "\n")

    clump_in = OUT / "ldl_gws_for_clump.txt"
    with clump_in.open("w", encoding="utf-8", newline="\n") as dst:
        dst.write("SNP\tP\n")
        for row in usable:
            z = abs(float(row["beta"]) / float(row["se"]))
            dst.write(f"{row['rsid']}\t{10 ** -(z + 3):.6g}\n")
    subprocess.run([
        str(PLINK), "--bfile", str(BFILE), "--clump", str(clump_in),
        "--clump-kb", "10000", "--clump-p1", "5e-8", "--clump-p2", "5e-8",
        "--clump-r2", "0.001", "--out", str(OUT / "ldl_clump"),
    ], check=True, capture_output=True)
    print("wrote", OUT / "ldl_clump.clumped", flush=True)


if __name__ == "__main__":
    main()
