# Max-ID spatial regression tutorial
> Ongoing research by Nathan and Likun

A worked R tutorial for a spatial regression model built on a max-infinitely-divisible
(max-id) latent field: an exponentially-tilted positive-stable (expPS) spatial process
observed on the log scale with independent log-Laplace noise, fit by Monte Carlo EM
(MCEM). Start at **`maxid_regression_tutorial.html`** — it covers the model definition,
spatial design, simulation, the expPS/log-Laplace distribution theory, the E-step and
M-step derivations, bias correction, uncertainty quantification, and a real worked fit
at T=400 time replicates.

## Contents

| File | Purpose |
|---|---|
| `maxid_regression_tutorial.Rmd` | Source document for the tutorial |
| `maxid_regression_tutorial.html` | Rendered tutorial (open this) |
| `tutorial_figures/` | Figures embedded in the tutorial |
| `maxid_regression_mcem.R` | Core simulator: `simulate_model()`, expPS sampling (`rexpPS_half`), log-Laplace noise |
| `maxid_regression_fast_fit.R` | Production MCEM fitter: `fit_mcem_fast()`, `plot_mcem_traces()`, the Louis-information uncertainty routine, and the leave-one-replicate-out jackknife bias correction for theta |
| `maxid_tutorial_examples.R` | Spatial design helpers used by the tutorial: `make_spatial_design()`, `simulate_spatial()` |
| `check_mcem_performance.R` | Independent posterior re-checks (`check_final_mcem()`) and a batch-means Monte Carlo SE helper |
| `maxid_theta_jackknife_correction.R` | Standalone jackknife-correction re-check for a fit from before the automatic correction was added, or for an independent-seed verification |
| `elevation_grid.csv` | Saved AWS zoom-6 elevation samples, so the tutorial's spatial design doesn't need a network download to reproduce |
| `mcem_results_tutorial_T400.rds` | Saved simulation + fit object from the tutorial's own T=400 worked example |

## Reproducing the tutorial

```r
source("maxid_regression_mcem.R")
source("maxid_regression_fast_fit.R")
source("maxid_tutorial_examples.R")
```

reproduces every function the tutorial calls. The `.Rmd` itself only evaluates its
heavy MCEM fit when knit with `params$run_fit = TRUE`; by default it renders quickly
and shows that code without re-running it, since the worked numbers and most figures
are already baked into the prose.

To re-render the HTML:

```r
rmarkdown::render("maxid_regression_tutorial.Rmd", output_options = list(mathjax = NULL))
```

The tutorial's math is static, pre-rendered KaTeX (not live MathJax), produced by
post-processing pandoc's output with the `katex` R package:

```r
# install.packages(c("V8", "katex"))
katex::render_math_in_html("maxid_regression_tutorial.html")
```

This requires a pandoc binary that `rmarkdown` can find (set `RSTUDIO_PANDOC` to its
directory if `rmarkdown::pandoc_available()` returns `FALSE`).
