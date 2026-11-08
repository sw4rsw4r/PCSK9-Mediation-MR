"""Extract the PCSK9 gene body +/- 50 kb from GLGC LDL-C, CAD, and deCODE pQTL."""
from config import CAD, DECODE as PQ, GLGC as LDL, WINDOW_DIR as OUT

OUT.mkdir(parents=True, exist_ok=True)

# GRCh37 ENSG00000169174 55,505,221-55,530,525. GRCh38 55,039,445-55,064,852.
B37_LO, B37_HI = 55_455_221, 55_580_525
B38_LO, B38_HI = 54_989_445, 55_114_852


def in_window(chrom, pos, lo, hi):
    return chrom in {"1", "chr1"} and lo <= pos <= hi


def stream_table(path, lo, hi, opener, chrom_name, pos_name, keep_names):
    rows = {}
    n = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        ic, ip = idx[chrom_name], idx[pos_name]
        cols = [idx[name] for name in keep_names]
        for line in handle:
            n += 1
            if n % 2_000_000 == 0:
                print(f"  {path.name} {n:,}", flush=True)
            parts = line.rstrip("\n").split("\t")
            chrom = parts[ic].replace("chr", "")
            if chrom not in {"1", ""} and chrom.isdigit() and int(chrom) > 1:
                break
            if chrom != "1":
                continue
            try:
                pos = int(float(parts[ip]))
            except ValueError:
                continue
            if pos < lo:
                continue
            if pos > hi:
                break
            rows[n] = [parts[c] if c < len(parts) else "" for c in cols]
    print(f"{path.name}: scanned {n:,}, kept {len(rows):,}", flush=True)
    return rows


def main():
    import gzip

    print("LDL-C", flush=True)
    ldl_keep = [
        "variant_id", "base_pair_location", "effect_allele", "other_allele",
        "n", "effect_allele_frequency", "beta", "standard_error", "p_value",
    ]
    ldl = stream_table(LDL, B37_LO, B37_HI, open, "chromosome", "base_pair_location", ldl_keep)

    # GCST90132314.h is on GRCh38. rs11591147 sits at 55,039,974.
    print("CAD", flush=True)
    cad_keep = [
        "rsid", "base_pair_location", "effect_allele", "other_allele",
        "beta", "standard_error", "effect_allele_frequency", "p_value", "n",
    ]
    cad = stream_table(
        CAD, B38_LO, B38_HI, gzip.open, "chromosome", "base_pair_location", cad_keep
    )

    print("deCODE", flush=True)
    pq_keep = [
        "rsids", "Pos", "effectAllele", "otherAllele", "Beta", "SE", "Pval", "N", "ImpMAF",
    ]
    pq = stream_table(PQ, B38_LO, B38_HI, gzip.open, "Chrom", "Pos", pq_keep)

    def dump(name, rows, id_index):
        path = OUT / f"_ids_{name}.tsv"
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            for rec in rows.values():
                handle.write("\t".join(rec) + "\n")
        ids = []
        for rec in rows.values():
            raw = rec[id_index].strip()
            ids.append(raw)
        print(name, "examples", ids[:5], flush=True)
        return path

    dump("ldl", ldl, 0)
    dump("cad", cad, 0)
    dump("pqtl", pq, 0)

    ldl_by = {}
    for rec in ldl.values():
        rsid = rec[0].strip()
        if rsid and rsid not in {".", "NA"}:
            ldl_by[rsid] = rec
    cad_by = {}
    for rec in cad.values():
        rsid = rec[0].strip()
        if rsid and rsid not in {".", "NA"}:
            cad_by.setdefault(rsid, rec)
    pq_by = {}
    for rec in pq.values():
        raw = rec[0].strip()
        if not raw or raw in {".", "NA"}:
            continue
        for rsid in raw.split(","):
            rsid = rsid.strip()
            if rsid:
                pq_by.setdefault(rsid, rec)

    print(
        f"rsids LDL {len(ldl_by):,} CAD {len(cad_by):,} pQTL {len(pq_by):,} "
        f"LDL&CAD {len(set(ldl_by)&set(cad_by)):,} "
        f"LDL&pQTL {len(set(ldl_by)&set(pq_by)):,} "
        f"CAD&pQTL {len(set(cad_by)&set(pq_by)):,}",
        flush=True,
    )
    shared = sorted(set(ldl_by) & set(cad_by) & set(pq_by))
    print(f"overlap {len(shared):,}", flush=True)

    out = OUT / "pcsk9_window_overlap.tsv"
    header = [
        "rsid", "position_b37", "effect_allele_ldl", "other_allele_ldl",
        "beta_LDL_C", "se_LDL_C", "pval_LDL_C", "maf_LDL_C", "n_LDL_C",
        "beta_CAD", "se_CAD", "pval_CAD", "maf_CAD", "n_CAD",
        "effect_allele_cad", "other_allele_cad", "position_cad_b38",
        "position_b38", "effect_allele_pqtl", "other_allele_pqtl",
        "beta_PCSK9_pQTL", "se_PCSK9_pQTL", "pval_PCSK9_pQTL", "maf_PCSK9_pQTL", "n_PCSK9_pQTL",
    ]
    with out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(header) + "\n")
        for rsid in shared:
            a, c, q = ldl_by[rsid], cad_by[rsid], pq_by[rsid]
            row = [
                rsid, a[1], a[2], a[3], a[6], a[7], a[8], a[5], a[4],
                c[4], c[5], c[7], c[6], c[8], c[2], c[3], c[1],
                q[1], q[2], q[3], q[4], q[5], q[6], q[8], q[7],
            ]
            handle.write("\t".join(row) + "\n")
    (OUT / "overlap_rsids.txt").write_text("\n".join(shared) + "\n", encoding="utf-8")
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()
