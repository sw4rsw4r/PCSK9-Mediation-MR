"""Local paths for files that are not redistributed with this repository.

Set PCSK9_RAW to the directory that holds the GWAS downloads, PCSK9_BFILE to the
1000 Genomes European PLINK prefix, and PLINK to the PLINK 1.9 executable.
PCSK9_RESULTS is where the scripts write intermediate files.
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("PCSK9_RAW", ROOT / "data" / "raw"))
RESULTS = Path(os.environ.get("PCSK9_RESULTS", ROOT / "results"))

GLGC = RAW / "GCST90239658_buildGRCh37.tsv"
CAD = RAW / "GCST90132314.h.tsv.gz"
DECODE = RAW / "Proteomics_SMP_PC0_5231_79_PCSK9_PCSK9_10032022.txt.gz"
BFILE = Path(os.environ.get("PCSK9_BFILE", RAW / "1000G_EUR"))
PLINK = Path(os.environ.get("PLINK", "plink"))

LDL_GWS_DIR = RESULTS / "ldl_gws"
LDL_GWS = LDL_GWS_DIR / "ldl_gws.tsv"
CLUMP_DIR = RESULTS / "ldl_clump"
CLUMPED = CLUMP_DIR / "ldl_clump.clumped"
WINDOW_DIR = RESULTS / "window"
SETTING1 = RESULTS / "setting1"
SETTING2 = RESULTS / "setting2"
