"""Total, indirect and direct effects of PCSK9 on CAD through LDL-C, from summary_results.

Inputs (summary_results/)
  pcsk9_instruments.tsv           the 16 PCSK9 instruments (PCSK9 gene +/- 50 kb, r2 < 0.1)
  ldl_setting1.tsv, ldl_setting2.tsv   LDL-C instruments of Settings 1 and 2
  region_corr.ld, region_snps.bim      1000 Genomes EUR correlations (r, PLINK --r square) of the
                                       16 PCSK9 instruments and the Setting 2 LDL-C instruments
                                       within 10 Mb of PCSK9 (24 variants)
Outputs (results/estimates/)
  estimates.tsv    one row per analysis x setting x LDL-C instrument set
  region_ivw.tsv   LDL-C -> CAD by IVW within the APOE and LPA regions and outside them

Analyses
  main         16 PCSK9 instruments; LD within 10 Mb of PCSK9 modelled; instruments elsewhere
               independent (LDL-C instruments were clumped at r2 < 0.001)
  independent  the PCSK9 instruments re-clumped at r2 < 0.001 (3 remain); LD not modelled;
               Setting 2 MVMR drops LDL-C instruments with r2 >= 0.001 to them
LDL-C instrument sets: "APOE+LPA excluded" (main analysis of the paper) and "all".

Estimators
  IVW and multivariable IVW in generalised least squares form, Omega = diag(se_Y) R diag(se_Y),
  with multiplicative random effects (residual standard error floored at 1). With R = I they are
  the MendelianRandomization estimators; with R given they match mr_ivw / mr_mvivw(correl = TRUE).
  Two-step: indirect = theta_XM * theta_MY (delta method); direct = total - indirect with
  Var = Var(total) + Var(indirect) - 2 theta_XM Cov(total, theta_MY). The covariance comes from
  CAD associations shared by both estimates (rs11591147 in Setting 2) or linked by LD; it is 0 in
  Setting 1. MVMR: direct from MVMR; indirect = total - direct, with Cov(total, direct).
  Covariances use the two weight vectors and Omega, each scaled by its random-effects factor.
  Conditional F (Sanderson et al.) with the LD matrix, minimised over delta; with two exposures it is
  the same for both.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import norm

from config import ESTIMATES_DIR as OUT, SUMMARY

COMP = str.maketrans("ACGT", "TGCA")
APOE = (19, 44_409_011, 46_412_650)   # gene +/- 1 Mb, GRCh37
LPA = (6, 159_952_515, 162_087_407)


def in_region(d, reg):
    c, lo, hi = reg
    return (d.chromosome == c) & (d.position >= lo) & (d.position <= hi)


EXCLUSIONS = {
    "APOE+LPA excluded": lambda d: (in_region(d, APOE) | in_region(d, LPA)).values,
    "all": lambda d: np.zeros(len(d), dtype=bool),
}


def allele_sign(ea, a1, a2):
    """+1 if the LDL-C effect allele is PLINK's A1, -1 if A2 (either strand), else nan."""
    ea, a1, a2 = ea.upper(), a1.upper(), a2.upper()
    if ea in (a1, a1.translate(COMP)):
        return 1.0
    if ea in (a2, a2.translate(COMP)):
        return -1.0
    return np.nan


# ---------- data ----------
def load_inputs():
    pc = pd.read_csv(SUMMARY / "pcsk9_instruments.tsv", sep="\t")
    cis = pc[["rsid", "beta_pqtl", "se_pqtl", "beta_ldl", "se_ldl", "beta_cad", "se_cad"]].copy()
    s1 = pd.read_csv(SUMMARY / "ldl_setting1.tsv", sep="\t")
    s2 = pd.read_csv(SUMMARY / "ldl_setting2.tsv", sep="\t")
    ea = dict(zip(pc.rsid, pc.effect_allele))
    for d in (s1, s2):
        ea.update(dict(zip(d.rsid, d.effect_allele)))
    return cis, s1, s2, ea


def ld_panel(ea):
    bim = pd.read_csv(SUMMARY / "region_snps.bim", sep="\t", header=None)
    m = np.loadtxt(SUMMARY / "region_corr.ld")
    m = np.tril(m) + np.tril(m, -1).T
    s = np.array([allele_sign(ea[r], a1, a2) for r, a1, a2 in zip(bim[1], bim[4], bim[5])])
    if np.isnan(s).any():
        raise ValueError("allele not aligned in the LD panel")
    return list(bim[1]), m * np.outer(s, s)


def corr_from(panel, ids):
    c = np.eye(len(ids))
    if panel is None:
        return c
    pos = {r: i for i, r in enumerate(panel[0])}
    for a, ra in enumerate(ids):
        for b, rb in enumerate(ids):
            if a != b and ra in pos and rb in pos:
                c[a, b] = panel[1][pos[ra], pos[rb]]
    return c


def cross_corr(panel, a_ids, b_ids):
    """Correlations between two instrument lists: 1 for the same variant, else LD from the panel."""
    pos = {} if panel is None else {r: i for i, r in enumerate(panel[0])}
    c = np.zeros((len(a_ids), len(b_ids)))
    for i, a in enumerate(a_ids):
        for j, b in enumerate(b_ids):
            if a == b:
                c[i, j] = 1.0
            elif a in pos and b in pos:
                c[i, j] = panel[1][pos[a], pos[b]]
    return c


def max_r2_to(panel, rsid, targets):
    pos = {r: i for i, r in enumerate(panel[0])}
    if rsid not in pos:
        return 0.0
    return max((panel[1][pos[rsid], pos[t]] ** 2 for t in targets if t in pos), default=0.0)


# ---------- estimators ----------
def gls(X, y, sy, C):
    """Returns estimates, random-effects SEs, weight matrix (theta = A @ y), RE factor, Omega."""
    om = np.outer(sy, sy) * C
    w = np.linalg.inv(om)
    v = np.linalg.inv(X.T @ w @ X)
    a = v @ X.T @ w
    th = a @ y
    r = y - X @ th
    k = max(1.0, float(np.sqrt(r @ w @ r / (len(y) - X.shape[1]))))
    return th, np.sqrt(np.diag(v)) * k, a, k, om


def cond_f(x, sx, z, sz, C):
    def q(de):
        s = np.outer(sx, sx) * C + de ** 2 * np.outer(sz, sz) * C
        e = x - de * z
        return e @ np.linalg.solve(s, e)
    r = minimize_scalar(q, bounds=(-50, 50), method="bounded")
    return q(r.x) / (len(x) - 1)


def run(label, cis, panel, est_panel, s1, s2, mv_r2_cut, exclusions):
    """panel chooses instruments (LD filters); est_panel is the LD used in estimation (None = none)."""
    cis_ids = list(cis.rsid)
    cc = corr_from(est_panel, cis_ids)
    p1 = cis[["beta_pqtl"]].values
    tot = gls(p1, cis.beta_cad.values, cis.se_cad.values, cc)
    xm = gls(p1, cis.beta_ldl.values, cis.se_ldl.values, cc)
    T, sT, A, sA = tot[0][0], tot[1][0], xm[0][0], xm[1][0]
    rows = []
    for setting, d0 in (("Setting 1", s1), ("Setting 2", s2)):
        if setting == "Setting 1":   # no LDL-C instrument with r2 >= 0.001 to a PCSK9 instrument
            d0 = d0[[max_r2_to(panel, r, cis_ids) < 0.001 for r in d0.rsid]]
        for ex_name in exclusions:
            d = d0[~EXCLUSIONS[ex_name](d0)]
            u = gls(d[["beta_ldl"]].values, d.beta_cad.values, d.se_cad.values, corr_from(est_panel, list(d.rsid)))
            B, sB, kU = u[0][0], u[1][0], u[3]
            ind = A * B
            s_ind = np.sqrt(B ** 2 * sA ** 2 + A ** 2 * sB ** 2)
            om_x = cross_corr(est_panel, cis_ids, list(d.rsid)) * np.outer(cis.se_cad.values, d.se_cad.values)
            cov_tb = tot[3] * kU * tot[2][0] @ om_x @ u[2][0]
            dir2, s_dir2 = T - ind, np.sqrt(max(sT ** 2 + s_ind ** 2 - 2 * A * cov_tb, 0.0))
            g = d[d.beta_pqtl.notna() & ~d.rsid.isin(cis_ids)]
            dropped = [r for r in g.rsid if max_r2_to(panel, r, cis_ids) >= mv_r2_cut]
            g = g[~g.rsid.isin(dropped)]
            mv = pd.concat([g[["rsid", "beta_ldl", "se_ldl", "beta_cad", "se_cad", "beta_pqtl", "se_pqtl"]], cis])
            mid = list(mv.rsid)
            cm = corr_from(est_panel, mid)
            m = gls(mv[["beta_ldl", "beta_pqtl"]].values, mv.beta_cad.values, mv.se_cad.values, cm)
            Bm, sBm, Dm, sDm = m[0][0], m[1][0], m[0][1], m[1][1]
            wt = np.zeros(len(mid))
            wt[[mid.index(r) for r in cis_ids]] = tot[2][0]
            cov = wt @ m[4] @ m[2][1] * tot[3] * m[3]
            indm, s_indm = T - Dm, np.sqrt(max(sT ** 2 + sDm ** 2 - 2 * cov, 0.0))
            f = cond_f(mv.beta_pqtl.values, mv.se_pqtl.values, mv.beta_ldl.values, mv.se_ldl.values, cm)
            rec = dict(analysis=label, setting=setting, exclusion=ex_name, n_cis=len(cis_ids),
                       n_uni=len(d), n_mvmr=len(mid), mvmr_dropped_ld=",".join(dropped), cond_f=f, rse_ldl_cad=kU,
                       corr_total_ldl_cad=cov_tb / (sT * sB))
            for q, e, s in (("total", T, sT), ("pcsk9_ldl", A, sA), ("ldl_cad_uni", B, sB),
                            ("ldl_cad_mvmr", Bm, sBm), ("indirect_twostep", ind, s_ind),
                            ("direct_twostep", dir2, s_dir2), ("indirect_mvmr", indm, s_indm),
                            ("direct_mvmr", Dm, sDm)):
                rec[q], rec[q + "_se"] = e, s
                rec[q + "_p"] = 2 * norm.sf(abs(e / s))
            rec["pct_mediated_twostep"] = 100 * ind / T
            rec["pct_mediated_mvmr"] = 100 * indm / T
            rows.append(rec)
    return rows


def clump_independent(cis, panel, cut):
    order = cis.assign(z=np.abs(cis.beta_pqtl / cis.se_pqtl)).sort_values("z", ascending=False).rsid
    keep = []
    for r in order:
        if max_r2_to(panel, r, keep) < cut:
            keep.append(r)
    return cis[cis.rsid.isin(keep)].reset_index(drop=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cis, s1, s2, ea = load_inputs()
    panel = ld_panel(ea)
    cis3 = clump_independent(cis, panel, 0.001)
    print("independent PCSK9 instruments:", list(cis3.rsid))
    rows = run("main", cis, panel, panel, s1, s2, 0.1, ["APOE+LPA excluded", "all"])
    rows += run("independent", cis3, panel, None, s1, s2, 0.001, ["APOE+LPA excluded"])
    est = pd.DataFrame(rows)
    est.to_csv(OUT / "estimates.tsv", sep="\t", index=False, float_format="%.6g")

    reg = []
    for name, mask in (("APOE region", in_region(s1, APOE)), ("LPA region", in_region(s1, LPA)),
                       ("outside APOE and LPA", ~(in_region(s1, APOE) | in_region(s1, LPA)))):
        d = s1[mask]
        u = gls(d[["beta_ldl"]].values, d.beta_cad.values, d.se_cad.values, np.eye(len(d)))
        reg.append(dict(setting="Setting 1", instruments=name, n=len(d), ldl_cad_ivw=u[0][0], se=u[1][0]))
    pd.DataFrame(reg).to_csv(OUT / "region_ivw.tsv", sep="\t", index=False, float_format="%.6g")

    show = est[["analysis", "setting", "exclusion", "n_uni", "n_mvmr", "total", "ldl_cad_uni",
                "direct_twostep", "direct_twostep_p", "direct_mvmr", "direct_mvmr_p", "cond_f"]]
    print(show.round(4).to_string(index=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
