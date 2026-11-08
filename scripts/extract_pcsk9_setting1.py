"""Setting 1 LDL-C instruments for PCSK9.

Start from the genome-wide LDL-C clump ranked by |z| (clump_ldl_by_z.py;
p < 5e-8, r2 < 0.001, 10 Mb, 1000 Genomes EUR). Drop a variant if it is inside
the PCSK9 gene body +/- 50 kb or if r2 with a cis PCSK9 instrument is at least
0.001. A genome-wide association with plasma PCSK9 outside that region is
not a reason to drop the variant. CAD and deCODE betas are harmonized to
the LDL-C effect allele. Measured pQTL betas are kept. Nothing is filled
with 0 or 1.

ANCHORS is the list of 16 cis PCSK9 instruments written by select_cis_instruments.R.
"""
import gzip
import subprocess
from pathlib import Path

from config import BFILE, CAD, CLUMPED, DECODE, LDL_GWS as LDL, PLINK, SETTING1 as OUT

OUT.mkdir(parents=True, exist_ok=True)

# GRCh37 PCSK9 gene body 55,505,221-55,530,525, plus 50 kb.
WIN_LO, WIN_HI = 55_455_221, 55_580_525
GENE = 55_505_221
R2_CUT = 0.001
ANCHORS = [
    "rs11591147", "rs472495", "rs505151", "rs7525503", "rs2479409",
    "rs6681159", "rs41297885", "rs4927191", "rs11583974", "rs45613943",
    "rs56349475", "rs150119739", "rs11810371", "rs143646411",
    "rs75050571", "rs17192725",
]
COMP = str.maketrans("ACGT", "TGCA")


def read_clumped():
    snps = []
    for line in CLUMPED.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 5 or not parts[0].isdigit() or parts[2] == "SNP":
            continue
        snps.append({"chr": int(parts[0]), "snp": parts[2], "bp": int(parts[3]), "p": parts[4]})
    return snps


def load_ldl(wanted):
    rows = {}
    with LDL.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            if rec["rsid"] in wanted:
                rows[rec["rsid"]] = rec
    return rows


def r2_against_anchors(candidates):
    use = sorted({row["snp"] for row in candidates} | set(ANCHORS))
    (OUT / "ld_extract.txt").write_text("\n".join(use) + "\n", encoding="utf-8")
    (OUT / "pcsk9_anchors.txt").write_text("\n".join(ANCHORS) + "\n", encoding="utf-8")
    subprocess.check_call([
        str(PLINK), "--bfile", str(BFILE),
        "--extract", str(OUT / "ld_extract.txt"),
        "--r2",
        "--ld-snp-list", str(OUT / "pcsk9_anchors.txt"),
        "--ld-window-kb", "10000",
        "--ld-window", "999999",
        "--ld-window-r2", str(R2_CUT),
        "--out", str(OUT / "ld_vs_pcsk9"),
    ])
    tagged = {}
    ld_path = OUT / "ld_vs_pcsk9.ld"
    anchors = set(ANCHORS)
    if not ld_path.exists():
        return tagged
    for line in ld_path.read_text(encoding="utf-8", errors="replace").splitlines()[1:]:
        parts = line.split()
        if len(parts) < 7:
            continue
        a, b, r2 = parts[2], parts[5], float(parts[6])
        if a in anchors and b not in anchors:
            other, anchor = b, a
        elif b in anchors and a not in anchors:
            other, anchor = a, b
        else:
            continue
        prev = tagged.get(other)
        if prev is None or r2 > prev[0]:
            tagged[other] = (r2, anchor)
    return tagged


def complement(allele):
    allele = allele.upper()
    if any(base not in "ACGT" for base in allele):
        return None
    return allele.translate(COMP)


def palindromic(a, b):
    return {a, b} in ({"A", "T"}, {"C", "G"})


def harmonize_to_ldl(ldl, rows, beta_key="beta", se_key="se", ea_key="ea", oa_key="oa", eaf_key="eaf", p_key="p"):
    ea, oa = ldl["effect_allele"].upper(), ldl["other_allele"].upper()
    try:
        eaf = float(ldl["eaf"])
    except ValueError:
        eaf = None
    for rec in rows:
        cea, coa = rec[ea_key].upper(), rec[oa_key].upper()
        flip = None
        if cea == ea and coa == oa:
            flip = False
        elif cea == oa and coa == ea:
            flip = True
        else:
            cea_c, coa_c = complement(cea), complement(coa)
            if cea_c is None:
                continue
            if palindromic(ea, oa):
                continue
            if cea_c == ea and coa_c == oa:
                flip = False
            elif cea_c == oa and coa_c == ea:
                flip = True
        if flip is None:
            continue
        if palindromic(ea, oa) and eaf is not None and 0.42 <= eaf <= 0.58:
            return None
        beta = float(rec[beta_key])
        se = float(rec[se_key])
        raw_eaf = rec.get(eaf_key, "")
        ceaf = float(raw_eaf) if raw_eaf not in {"", "NA", None} else None
        if flip:
            beta = -beta
            if ceaf is not None:
                ceaf = 1 - ceaf
        if palindromic(ea, oa) and ceaf is not None and eaf is not None:
            if abs(eaf - ceaf) > abs(eaf - (1 - ceaf)):
                beta = -beta
                ceaf = 1 - ceaf
        return {
            "beta": beta,
            "se": se,
            "eaf": ceaf,
            "p": rec.get(p_key, ""),
            "n": rec.get("n", ""),
        }
    return None


def stream_cad(wanted):
    found = {rs: [] for rs in wanted}
    with gzip.open(CAD, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        n = 0
        left = set(wanted)
        for line in handle:
            n += 1
            if n % 5000000 == 0:
                print(f"  CAD rows {n:,} still {len(left):,}", flush=True)
            parts = line.rstrip("\n").split("\t")
            rs = parts[idx["rsid"]].strip()
            if rs not in found:
                continue
            found[rs].append({
                "ea": parts[idx["effect_allele"]].strip(),
                "oa": parts[idx["other_allele"]].strip(),
                "beta": parts[idx["beta"]].strip(),
                "se": parts[idx["standard_error"]].strip(),
                "eaf": parts[idx["effect_allele_frequency"]].strip(),
                "p": parts[idx["p_value"]].strip(),
                "n": parts[idx["n"]].strip() if "n" in idx else "",
            })
    return found


def stream_decode(wanted):
    found = {rs: [] for rs in wanted}
    with gzip.open(DECODE, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        n = 0
        for line in handle:
            n += 1
            if n % 2000000 == 0:
                print(f"  deCODE rows {n:,}", flush=True)
            parts = line.rstrip("\n").split("\t")
            raw = parts[idx["rsids"]].strip()
            if raw in {"", ".", "NA"}:
                continue
            ids = [x for x in raw.replace(";", ",").split(",") if x]
            hit = [x for x in ids if x in found]
            if not hit:
                continue
            rec = {
                "ea": parts[idx["effectAllele"]].strip(),
                "oa": parts[idx["otherAllele"]].strip(),
                "beta": parts[idx["Beta"]].strip(),
                "se": parts[idx["SE"]].strip(),
                "eaf": "",
                "p": parts[idx["Pval"]].strip(),
                "n": parts[idx["N"]].strip(),
            }
            for rs in hit:
                found[rs].append(rec)
    return found


def fmt(value, spec):
    if value is None or value == "":
        return ""
    if isinstance(value, float):
        return format(value, spec)
    return str(value)


def main():
    clumped = read_clumped()
    print(f"clumped {len(clumped)}", flush=True)
    near = [
        row for row in clumped
        if row["chr"] == 1 and abs(row["bp"] - GENE) <= 10_000_000
    ]
    print(f"chr1 within 10 Mb of PCSK9: {len(near)}", flush=True)
    tagged = r2_against_anchors(near)
    print(f"tagged by r2>={R2_CUT}: {len(tagged)}", flush=True)

    ldl = load_ldl({row["snp"] for row in clumped})
    print(f"LDL stats found {len(ldl)}", flush=True)
    print("streaming CAD", flush=True)
    cad = stream_cad(set(ldl))
    print("streaming deCODE", flush=True)
    pq = stream_decode(set(ldl))

    removed = []
    kept_rows = []
    n_no_cad = 0
    n_no_pq = 0
    for row in clumped:
        rs = row["snp"]
        src = ldl.get(rs)
        if src is None:
            removed.append((rs, row["chr"], row["bp"], "missing LDL stats"))
            continue
        reason = None
        if rs in set(ANCHORS) or (row["chr"] == 1 and WIN_LO <= row["bp"] <= WIN_HI):
            reason = "inside PCSK9 +/-50 kb or is a cis pQTL instrument"
        elif rs in tagged:
            r2, anchor = tagged[rs]
            reason = f"r2={r2:.4f} with {anchor}"
        cad_hit = harmonize_to_ldl(src, cad.get(rs, []))
        pq_hit = harmonize_to_ldl(src, pq.get(rs, []))
        if reason:
            removed.append((rs, row["chr"], row["bp"], reason))
            continue
        if cad_hit is None:
            n_no_cad += 1
            removed.append((rs, row["chr"], row["bp"], "CAD allele mismatch"))
            continue
        if pq_hit is None:
            n_no_pq += 1
        kept_rows.append((src, cad_hit, pq_hit, row))

    print(
        f"kept {len(kept_rows)} removed {len(removed)} "
        f"no_cad {n_no_cad} no_pqtl {n_no_pq}",
        flush=True,
    )
    with (OUT / "removed_pcsk9_ld.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("rsid\tchromosome\tposition\treason\n")
        for rec in removed:
            handle.write("\t".join(map(str, rec)) + "\n")
    with (OUT / "ldl_ivs_setting1.tsv").open("w", encoding="utf-8", newline="\n") as handle:
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
                fmt(cad_hit["eaf"], ".6g"), str(cad_hit["p"]),
                pq_beta, pq_se, pq_p, pq_n,
            ]) + "\n")


if __name__ == "__main__":
    main()
