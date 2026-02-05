# evalinger

**E-Values for Adaptive Clinical Trial Monitoring**

An R package implementing e-value and e-process methodology for interim
monitoring of adaptive clinical trials. Companion software to:

> Sokolova, A. and Sokolov, V. (2026). "E-values for Adaptive Clinical Trials:
> A Practical Methodology for Anytime-Valid Monitoring and Adaptation."

## Overview

E-values and e-processes provide anytime-valid inference under optional stopping
and optional continuation. Unlike group sequential boundaries that require a
prespecified monitoring schedule, e-value monitoring uses a single time-uniform
threshold that controls Type I error regardless of how many times the data are
inspected or whether the trial continues beyond its planned horizon.

The `evalinger` package provides a complete toolkit:

| Module | Functions | Purpose |
|---|---|---|
| **E-process construction** | `eprocess_binary()`, `eprocess_logrank()` | Betting martingale for binary endpoints; safe logrank for survival |
| **Design calibration** | `grow_lambda()`, `expected_growth_rate()`, `edesign_binary()` | GROW-optimal betting fraction, power, expected sample size |
| **Real-time monitoring** | `emonitor()`, `update()` | Stateful streaming monitor with batch updates |
| **Confidence sequences** | `confseq_binary()` | Always-valid confidence intervals for treatment effects |
| **Futility** | `futility_cs()`, `futility_eprocess()` | Futility via confidence sequences or reciprocal e-process |
| **Platform trials** | `platform_monitor()`, `ebh()` | Multi-arm monitoring with e-Benjamini-Hochberg control |
| **Hybrid monitoring** | `hybrid_monitor()` | Joint e-process + group sequential monitoring |
| **Method comparison** | `simulate_comparison()` | Five-method Monte Carlo comparison (e-value, group sequential, naive p, naive Bayes, calibrated Bayes) |

## Installation

```r
# Install from GitHub
devtools::install_github("vsokolov/evalinger")

# Or using pak
pak::pak("vsokolov/evalinger")
```

## Quick start

```r
library(evalinger)

# 1. Design: find the optimal betting fraction
design <- edesign_binary(p_C = 0.30, delta = 0.15, alpha = 0.025, power = 0.80)
print(design)

# 2. Construct an e-process from trial data
x_T <- rbinom(150, 1, 0.45)  # treatment outcomes
x_C <- rbinom(150, 1, 0.30)  # control outcomes
ep <- eprocess_binary(x_T, x_C, lambda = design$lambda, alpha = 0.025)
plot(ep)

# 3. Real-time monitoring (batch-by-batch)
mon <- emonitor(alpha = 0.025, lambda = design$lambda)
mon <- update(mon, x_T = x_T[1:50],  x_C = x_C[1:50])   # interim 1
mon <- update(mon, x_T = x_T[51:100], x_C = x_C[51:100]) # interim 2
print(mon)

# 4. Always-valid confidence sequence
cs <- confseq_binary(x_T, x_C, alpha = 0.05)
plot(cs)

# 5. Compare methods head-to-head
cmp <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 200,
                           n_looks = 20, alpha = 0.025, nrep = 5000)
print(cmp)
plot(cmp)
```

## Interactive web application

The package includes a Shiny web application with three dashboards:

- **Design Calculator** -- calibrate the GROW-optimal betting fraction and
  visualize expected stopping time as a function of lambda
- **Monitoring Dashboard** -- simulate batch-by-batch monitoring with live
  e-process trajectory and confidence sequence visualization
- **Method Comparison** -- run head-to-head comparisons of e-value, group
  sequential, and Bayesian monitoring under user-specified scenarios

### Run locally

```r
library(evalinger)
shiny::runApp(system.file("shiny", package = "evalinger"))
```

### Deploy to shinyapps.io

1. Install required packages:

```r
install.packages(c("shiny", "bslib", "rsconnect"))
devtools::install_github("vsokolov/evalinger")
```

2. Configure your shinyapps.io account:

```r
rsconnect::setAccountInfo(
  name   = "your-account",
  token  = "your-token",
  secret = "your-secret"
)
```

3. Deploy:

```r
rsconnect::deployApp(
  appDir  = system.file("shiny", package = "evalinger"),
  appName = "evalinger",
  appTitle = "evalinger: E-Values for Clinical Trials"
)
```

### Deploy to Posit Connect / Shiny Server

Copy the app directory to your server's content area:

```bash
cp -r "$(Rscript -e 'cat(system.file("shiny", package="evalinger"))')" /srv/shiny-server/evalinger
```

Ensure the `evalinger` package is installed in the R library used by Shiny
Server. The app will source the package functions automatically if the package
is not installed, but installing the package is recommended for production.

### Deploy via Docker

```dockerfile
FROM rocker/shiny:4.4.0

RUN R -e "install.packages('devtools'); \
          devtools::install_github('vsokolov/evalinger'); \
          install.packages(c('bslib'))"

COPY inst/shiny /srv/shiny-server/evalinger

EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
```

Build and run:

```bash
docker build -t evalinger-app .
docker run -p 3838:3838 evalinger-app
```

The app will be available at `http://localhost:3838/evalinger`.

## Regulatory context

This package implements methodology aligned with the regulatory framework
established by:

- **FDA Draft Guidance (January 2026)**: "Use of Bayesian Methodology in
  Clinical Trials of Drug and Biological Products" (CDER/CBER)
- **FDA Guidance (2019)**: "Adaptive Designs for Clinical Trials of Drugs and
  Biologics"
- **FDA Complex Innovative Design (CID) Program**

The interactive web application is designed to generate the operating
characteristics documentation that the 2026 FDA draft Bayesian guidance
requires for adaptive designs: Type I error rates, power curves, expected
sample sizes under design priors, and sensitivity to calibration parameters.

## Testing

```r
devtools::test()
```

## Citation

```bibtex
@article{SokolovaSokolov2026,
  author  = {Sokolova, Alexandra and Sokolov, Vadim},
  title   = {E-values for Adaptive Clinical Trials: A Practical Methodology
             for Anytime-Valid Monitoring and Adaptation},
  year    = {2026},
  note    = {Companion R package: https://github.com/vsokolov/evalinger}
}
```

## License

MIT
