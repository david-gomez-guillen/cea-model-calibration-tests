# cea-model-calibration-tests

A THALASSA test model built to show that a calibration needs constraints, not
only a good fit. It is a Markov model of screening for a hypothetical cancer,
with the same interface as `cea-model-test`, calibrated to a single target,
clinical incidence in six age groups, that can be met exactly in two ways:

- few preclinical lesions, nearly all of them progressive: the sensible way;
- many more lesions, a large share of them indolent (never symptomatic, but
  found and treated by screening): a way screening trials rule out.

Both fit the target with zero error, and the cost-effectiveness analysis
concludes differently on each:

- uncalibrated base values: **annual screening** is optimal;
- calibrated with the constraints (at most 20% of the lesions indolent, onset
  rising with age): **screening every 2 years** is optimal;
- calibrated without them: the fit is just as exact, but overdiagnosis makes
  every screening strategy a bad buy, and the **new treatment** is optimal.

The DSA of each analysis points at different parameters too.

In the Calibration tab, run both schemes with simulated annealing and 1,000
iterations (a few seconds each); the default Bayesian optimization cannot search
the twelve calibrated values, see the overview.

- [overview.md](overview.md): the model, its parameters, the calibration and the
  results of the three analyses. It is also the Overview tab of the app.
- [model.R](model.R): the model.
- [thalassa_interface.R](thalassa_interface.R): what THALASSA loads, including the
  two calibration schemes, `unconstrained` and `constrained`.
- [analysis.R](analysis.R): the whole demonstration outside the app.
  `Rscript analysis.R` writes it to `results/`.
