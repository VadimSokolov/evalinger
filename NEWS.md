# evalinger 0.1.0

* Initial CRAN release.
* Core e-process construction for two-arm binary trials (`eprocess_binary()`).
* Safe logrank e-process for survival endpoints (`eprocess_logrank()`).
* GROW-optimal betting fraction computation (`grow_lambda()`,
  `expected_growth_rate()`, `expected_stopping_time()`).
* Design calibration with simulated power and Type I error
  (`edesign_binary()`).
* Stateful incremental monitor for real-time DSMB use (`emonitor()`,
  `update.emonitor()`).
* Always-valid confidence sequences for treatment effects
  (`confseq_binary()`).
* CS-based and reciprocal e-process futility monitoring (`futility_cs()`,
  `futility_eprocess()`).
* Simulation tools comparing e-value, group sequential, naive, and
  Bayesian monitoring (`simulate_comparison()`).
* Hybrid e-process + group sequential monitoring (`hybrid_monitor()`).
* Multi-arm platform trial monitoring with e-BH multiplicity control
  (`platform_monitor()`, `ebh()`).
* Bridge to gsDesign (`as_eprocess_design()`).
* Shiny interactive demo (`inst/shiny/app.R`).
* Two vignettes: "Getting Started" and "Reproducing the Paper's Numerical
  Study".
