"""Setting 2 LDL-C instruments, including the PCSK9 region.

Same clump as Setting 1 (p < 5e-8, r2 < 0.001, 10 Mb, 1000 Genomes EUR, ranked by |z|).
Setting 1 drops the PCSK9 window and LD with the cis PCSK9 instruments. This set keeps
those variants. Allele mismatches are still dropped.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import extract_pcsk9_setting1 as base
from config import SETTING2 as OUT

OUT.mkdir(parents=True, exist_ok=True)


def main():
    clumped = base.read_clumped()
    print(f"clumped {len(clumped)}", flush=True)
    ldl = base.load_ldl({row["snp"] for row in clumped})
    print(f"LDL stats found {len(ldl)}", flush=True)
    print("streaming CAD", flush=True)
    cad = base.stream_cad(set(ldl))
    print("streaming deCODE", flush=True)
    pq = base.stream_decode(set(ldl))

    kept_rows = []
    removed = []
    for row in clumped:
        rs = row["snp"]
        src = ldl.get(rs)
        if src is None:
            removed.append((rs, row["chr"], row["bp"], "missing LDL stats"))
            continue
        cad_hit = base.harmonize_to_ldl(src, cad.get(rs, []))
        pq_hit = base.harmonize_to_ldl(src, pq.get(rs, []))
        if cad_hit is None:
            removed.append((rs, row["chr"], row["bp"], "CAD allele mismatch"))
            continue
        kept_rows.append((src, cad_hit, pq_hit, row))

    print(f"kept {len(kept_rows)} removed {len(removed)}", flush=True)
    with (OUT / "removed.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("rsid\tchromosome\tposition\treason\n")
        for rec in removed:
            handle.write("\t".join(map(str, rec)) + "\n")
    with (OUT / "ldl_ivs_setting2_new.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "rsid\tchromosome\tposition\teffect_allele\tother_allele\t"
            "beta_ldl\tse_ldl\teaf_ldl\tp_ldl\tn_ldl\t"
            "beta_cad\tse_cad\teaf_cad\tp_cad\t"
            "beta_pqtl\tse_pqtl\tp_pqtl\tn_pqtl\n"
        )
        for src, cad_hit, pq_hit, row in kept_rows:
            pq_beta = "" if pq_hit is None else f"{pq_hit['beta']:.8g}"
            pq_se = "" if pq_hit is None else f"{pq_hit['se']:.8g}"
            pq_p = "" if pq_hit is None else str(pq_hit["p"])
            pq_n = "" if pq_hit is None else str(pq_hit["n"])
            handle.write("\t".join([
                src["rsid"], src["chromosome"], src["position"],
                src["effect_allele"], src["other_allele"],
                src["beta"], src["se"], src["eaf"], src["p"], src["n"],
                f"{cad_hit['beta']:.8g}", f"{cad_hit['se']:.8g}",
                base.fmt(cad_hit["eaf"], ".6g"), str(cad_hit["p"]),
                pq_beta, pq_se, pq_p, pq_n,
            ]) + "\n")
    print("wrote", OUT / "ldl_ivs_setting2_new.tsv", flush=True)


if __name__ == "__main__":
    main()
