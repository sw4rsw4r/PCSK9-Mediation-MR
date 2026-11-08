# Local data

The variant-level tables used to draw the paper figures are already in `summary_results`. The files below are needed only to rebuild those tables from the original studies. They are not included in this repository.

Place the downloads in `data/raw`, or set `PCSK9_RAW` to another directory.

- `GCST90239658_buildGRCh37.tsv`, the GLGC LDL-C GWAS (GCST90239658).
- `GCST90132314.h.tsv.gz`, the CAD GWAS (GCST90132314).
- `Proteomics_SMP_PC0_5231_79_PCSK9_PCSK9_10032022.txt.gz`, the deCODE plasma PCSK9 SomaScan summary for aptamer 5231_79.
- 1000 Genomes Phase 3 European PLINK files. The prefix defaults to `data/raw/1000G_EUR`. Set `PCSK9_BFILE` if it is somewhere else.
- PLINK 1.9 on the path, or set `PLINK` to the executable.

Scripts write intermediate files under `results`. That directory is not part of the repository.

```bash
python scripts/extract_ldl_gws.py
python scripts/clump_ldl_by_z.py
python scripts/extract_pcsk9_window.py
python scripts/extract_pcsk9_setting1.py
python scripts/extract_pcsk9_setting2.py
Rscript scripts/select_cis_instruments.R
```

Before the cis-instrument script, build the linkage-disequilibrium matrix for the variants in `results/window/pcsk9_window_overlap.tsv`.

```bash
plink --bfile data/raw/1000G_EUR --extract results/window/overlap_rsids.txt --make-bed --out results/window/eur_window
plink --bfile results/window/eur_window --r square --out results/window/eur_window
```

`scripts/extract_pcsk9_window.py` writes `overlap_rsids.txt` next to the overlap table.
