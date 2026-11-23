# Finish the conda environment on Windows (not needed on Linux or macOS).
#
# conda-forge's Windows builds of several R packages with compiled code (for example robustbase,
# quantreg, KernSmooth and RcppEigen) fail to load at random with "Mingw-w64 runtime failure: 32 bit
# pseudo relocation ... out of range", depending on where Windows places the DLLs; the more packages
# are loaded, the more often it happens. The CRAN Windows binaries load reliably. This script installs
# the R packages used here, with their dependencies, as CRAN Windows binaries from the Posit Package
# Manager snapshot of 2026-11-20 into a separate library inside the environment (conda's own files
# are not changed), and sets R_LIBS so that R uses them first. R itself still comes from conda.
#
#   conda env create -f environment.yml
#   powershell -ExecutionPolicy Bypass -File scripts/setup_env_windows.ps1
param([string]$EnvName = 'pcsk9-mr', [string]$Conda = 'conda')
$ErrorActionPreference = 'Stop'
$line = (& $Conda env list) | Where-Object { $_ -match "^$([regex]::Escape($EnvName))\s" }
if (-not $line) { throw "conda environment '$EnvName' not found; create it with: conda env create -f environment.yml" }
$prefix = ($line -split '\s+')[-1]
$lib = (Join-Path $prefix 'Lib\R\library-cran') -replace '\\', '/'
if (Test-Path $lib) { Remove-Item -Recurse -Force $lib }
New-Item -ItemType Directory -Force $lib | Out-Null
$r = @"
options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2026-11-20'))
pins <- c(MendelianRandomization = '0.10.0', ggplot2 = '4.0.4', rmarkdown = '2.32', knitr = '1.52', jsonlite = '2.0.0')
# all recursive dependencies, except the base packages that come with R
db <- available.packages(type = 'binary')
deps <- unique(c(names(pins), unlist(tools::package_dependencies(names(pins), db = db, recursive = TRUE,
                                                                  which = c('Depends', 'Imports', 'LinkingTo')))))
deps <- setdiff(deps, rownames(installed.packages(priority = 'base')))
install.packages(deps, lib = '$lib', type = 'binary', dependencies = FALSE, quiet = TRUE)
got <- sapply(names(pins), function(p) as.character(packageVersion(p, lib.loc = '$lib')))
print(rbind(pinned = pins, installed = got))
if (!all(got == pins)) stop('installed versions differ from the pinned versions')
cat(length(deps), 'packages installed from the 2026-11-20 snapshot\n')
"@
$tmp = Join-Path ([IO.Path]::GetTempPath()) 'pcsk9_setup_cran.R'
[IO.File]::WriteAllText($tmp, $r)
& $Conda run -n $EnvName Rscript $tmp
& $Conda env config vars set -n $EnvName "R_LIBS=$lib" | Out-Null
& $Conda run -n $EnvName Rscript -e "for (p in c('MendelianRandomization', 'ggplot2', 'rmarkdown')) suppressPackageStartupMessages(library(p, character.only = TRUE)); cat('R packages load from', dirname(find.package('MendelianRandomization')), '\n')"
"$EnvName ready (R_LIBS=$lib)"
