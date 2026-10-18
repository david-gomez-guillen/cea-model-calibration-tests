# Constrained calibration test model

A cost-effectiveness model of screening for a hypothetical cancer, built to show
why a calibration has to respect what is known about the parameters, and not only
fit its targets.

The model is calibrated to one target, clinical cancer incidence by age, and
that target can be met exactly in more than one way:

- a **sensible** way: few preclinical lesions, nearly all of which become cancers
  if left alone, and more of them the older people are;
- a **clinically inconsistent** way: many more lesions, a third or more of them
  indolent, never causing symptoms, which screening trials rule out, and often a
  risk that falls and rises again with age.

Both fit the registry equally well, error zero, and they lead to opposite
decisions. The uncalibrated values lead to a third one:

| Parameters | Error | Indolent lesions, by age group | Optimal at 25,000 €/QALY |
|---|---|---|---|
| Uncalibrated (base values) | 5.9 | 10% everywhere | **Screening every year** |
| Calibrated **without** constraints | 0 | 36% to 58% | **New treatment, no screening** |
| Calibrated **with** constraints | 0 | 1% to 19% | **Screening every 2 years** |

Nothing in the fit tells the last two apart. Only the constraints do.

## Model structure

The model is a seven-state Markov cohort model with annual cycles:

- **Healthy**: alive and cancer-free. The whole cohort starts here.
- **Progressive lesion**: a preclinical cancer that causes no symptoms yet, but
  will, and that a screening test can find.
- **Indolent lesion**: a preclinical lesion that a screening test finds just as
  well, but that never causes symptoms and never kills.
- **Clinical**: a cancer diagnosed because of its symptoms.
- **Screen-detected**: a progressive lesion found by screening and treated early.
  Much less deadly than a clinical cancer.
- **Overdiagnosed**: an indolent lesion found by screening. Nobody can tell it from
  a screen-detected cancer, so it is treated and followed up the same way, at the
  same cost and loss of quality of life, with no benefit.
- **Dead**: absorbing state. Everyone dies of other causes at a background rate
  that rises with age (a Gompertz curve, from 0.08% a year at 40 to 2.8% at 79);
  clinical and screen-detected cancers can also kill.

Each year a share `p.onset` of the healthy develops a lesion, a share
`p.progressive` of which are progressive and the rest indolent. Progressive
lesions become clinical with probability `p.progression` a year, a mean sojourn
time of 3 years. The cohort is followed from age 40 to age 79 (40 cycles), with
no new entries. Screening, when a strategy has it, happens at the start of a
cycle; costs and utilities of the cycle are then accrued, and the transitions
applied at its end.

## Strategies

- `no_intervention`: every progressive lesion is diagnosed when it becomes
  symptomatic and receives standard care; indolent lesions are never found. The
  reference.
- `screening_5y`, `screening_2y`, `screening_1y`: everyone alive and undiagnosed is
  tested every 5 years, every 2 years or every year, from 50 to 74. The test finds
  a share `screening.sensitivity` of the lesions, progressive and indolent alike;
  a share `1 - screening.specificity` of the healthy test positive and need a
  work-up.
- `new_treatment`: no screening, and clinical cancers receive a new treatment that
  multiplies their annual probability of cancer death by `hr.new.treatment`, at
  `cost.new.treatment` a year.

The strategies are alternatives: the analysis picks one of them.

## Parameters

| Parameter | Base value | DSA range | Class | Meaning |
|---|---|---|---|---|
| `p.onset` | 0.0025 to 0.005 | ± 20% | Natural history | Annual probability of developing a preclinical lesion while healthy, by age group: 0.0025 at 40-49, rising by 0.0005 per age group to 0.005 at 70-79. Calibrated. |
| `p.progressive` | 0.9 | ± 20% | Natural history | Share of new lesions that are progressive; the rest are indolent. Calibrated, by age group. |
| `p.progression` | 0.33 | 0.25 – 0.5 | Natural history | Annual probability of a progressive lesion becoming symptomatic. |
| `p.clinical.death` | 0.12 | 0.08 – 0.16 | Natural history | Annual probability of cancer death of a clinical cancer under standard care. |
| `p.screen.detected.death` | 0.02 | 0.01 – 0.04 | Natural history | Annual probability of cancer death of a screen-detected cancer. |
| `screening.sensitivity` | 0.75 | 0.6 – 0.9 | Screening | Share of the lesions a screening test finds. |
| `screening.specificity` | 0.95 | 0.9 – 0.98 | Screening | Share of the healthy who test negative. |
| `cost.screening` | 70 | 50 – 100 | Screening | Cost of a screening test. |
| `cost.false.positive` | 1000 | 500 – 2000 | Screening | Cost of the work-up of a false positive. |
| `cost.early.treatment` | 15000 | 10000 – 20000 | Screening | One-off cost of treating a screen-detected or overdiagnosed cancer. |
| `cost.screen.detected` | 1000 | 500 – 1500 | Screening | Annual follow-up cost after a screen-detected or overdiagnosed cancer. |
| `utility.screen.detected` | 0.9 | 0.85 – 0.95 | Screening | Utility of a year lived after a screen-detected or overdiagnosed cancer. |
| `hr.new.treatment` | 0.6 | 0.45 – 0.75 | Treatment | Hazard ratio of cancer death of the new treatment. |
| `cost.new.treatment` | 1500 | 1000 – 2500 | Treatment | Annual cost of the new treatment per clinical cancer. |
| `cost.clinical` | 6000 | 4000 – 8000 | Treatment | Annual cost of standard care per clinical cancer. |
| `utility.clinical` | 0.65 | 0.55 – 0.75 | Treatment | Utility of a year lived with a clinical cancer. |
| `discount` | 0.03 | | General | Discount rate of costs and effects. |

The base values of `p.onset` and `p.progressive` are the uncalibrated starting
point: a risk that rises with age and 10% of indolent lesions, the kind of guess
a first look at the literature gives, with far too many lesions in the younger
age groups. They are also the initial guess of both calibration schemes, which
start strictly inside the feasible region of the constrained one.

## Strata

Results are reported, and the calibrated parameters estimated, by age group:
40-49, five-year groups through the screening ages (50-54, 55-59, 60-64, 65-69)
and 70-79.

## Outputs

- `summary`: one row per strategy with the discounted cost (`C`, €) and the
  discounted quality-adjusted life years (`E`) per person over the whole horizon.
- `incidence`: clinical cancers diagnosed per 100,000 people alive and per year,
  averaged over the years of each age group.
- `mortality`: cancer deaths per 100,000 people alive and per year.
- `outcomes`: screen-detected cancers, overdiagnosed cancers, clinical cancers and
  cancer deaths per 1,000 people of the cohort over the horizon, undiscounted.

## Calibration

Both schemes calibrate `p.onset` and `p.progressive` in every age group, twelve
values, against clinical incidence from a cancer registry, measured in the
`no_intervention` strategy, i.e. before screening existed:

| Age group | 40-49 | 50-54 | 55-59 | 60-64 | 65-69 | 70-79 |
|---|---|---|---|---|---|---|
| Clinical incidence (per 100,000) | 56 | 107 | 158 | 212 | 269 | 348 |

The error is the sum of the squared relative differences between the simulated
and the target incidence. Parameter sets that make the simulation fail get an
infinite error.

### The target can be met in more than one way

Only progressive lesions ever become clinical, so incidence measures how many
progressive lesions there are, `p.onset * p.progressive`, and nothing else. Any
share of progressive lesions fits the registry exactly, provided the onset rises
to make up for it, in every age group independently of the others. This table
shows some of the parameter sets that fit the target exactly, with the same
share in every age group:

| Indolent lesions | Onset per 1,000, 40-49 to 70-79 | Error | Optimal | Overdiagnosed per 1,000 screened every 2 years |
|---|---|---|---|---|
| 0% | 0.80 – 3.80 | 0 | Screening every 2 years | 0 |
| 10% | 0.89 – 4.25 | 0 | Screening every 2 years | 6.5 |
| 20% | 1.00 – 4.83 | 0 | Screening every 2 years | 14.6 |
| 30% | 1.14 – 5.59 | 0 | Screening every 2 years | 25.0 |
| 40% | 1.33 – 6.63 | 0 | New treatment | 38.9 |

What decides between them is overdiagnosis. Every indolent lesion screening
finds is a person treated as a cancer patient for nothing. Once a third or so of
the lesions are indolent, screening costs more than the cancers it catches early
are worth, and the new treatment, which only reaches cancers that would have
become clinical anyway, is the better buy.

The targets are synthetic: they are the incidence the model produces with 5% of
the lesions indolent and an onset rising steadily with age, the "truth" of this
example.

### What the constraints add

Two things are known about the parameters besides the registry:

- Screening trials measure overdiagnosis directly, as the excess of cancers in the
  screened arm that is still there years after screening stopped, and put the
  share of indolent lesions below 20%: `p.progressive` is at least 0.8 in every
  age group.
- The risk of developing a solid cancer only grows with age: `p.onset` does not
  decrease from one age group to the next.

That is what the `constrained` scheme is told, in `calibration.constraints()`,
each of whose eleven values must be zero or below. The constraints accept a
partial vector, so they also apply step by step to a stepwise run. The latent
space of the constrained scheme is trained on feasible samples only, drawn with
`p.progressive` in the plausible range and `p.onset` sorted by age.

### Schemes

- `unconstrained`: the twelve values within the search box, and nothing else.
- `constrained`: the same, plus the constraints above.

Both start from the base values, and with the default scope of 100% search
between 0 and twice them, capped at 1. Within that box, `p.progressive` can go
as low as about 0.15 in the youngest age group and 0.4 in the oldest while still
fitting the target, so most of the surface of exact fits lies outside the
plausible region, and an optimizer that searches it globally, with nothing to
choose between its points, lands there.

Two evaluation plots come with the schemes: the share of indolent lesions each
evaluation implies, against the 20% that trials allow, and the onset of all
lesions next to that of the progressive ones, which is the part the target pins
down.

## What the three analyses conclude

At 25,000 €/QALY, with the base values of every parameter that is not calibrated:

**Uncalibrated.** The base onset puts far too many cancers in the younger age
groups (157 per 100,000 a year at 40-49 against 56, 241 at 50-54 against 107),
so screening has a lot to find. The efficient frontier runs from no intervention
to screening every 5, 2 and 1 years, the last at 18,000 €/QALY: **screening
every year** is optimal. The decision turns on the screening parameters:
`p.screen.detected.death`, `screening.specificity` and `cost.false.positive`
each switch it to screening every 2 years within their range.

**Unconstrained calibration.** The run `analysis.R` reports fits the target
exactly with 36% to 58% of the lesions indolent, and an onset that falls from
65-69 to 70-79. Screening every 2 years overdiagnoses 61 people per 1,000, more
than the 39 cancers it finds, at 64,000 €/QALY over no intervention. Every
screening strategy is dominated or extendedly dominated, and the **new
treatment**, at 19,000 €/QALY, is optimal. Only the parameters of the treatment,
`hr.new.treatment` and `cost.new.treatment`, change the decision, to no
intervention: screening is not even in contention.

**Constrained calibration.** The run `analysis.R` reports fits the target just as
exactly, with 1% to 19% of the lesions indolent and an onset rising with age.
Overdiagnosis falls to 6 per 1,000, the frontier runs through screening every 5
years (10,000 €/QALY), every 2 years (14,000 €/QALY) and every year (31,000
€/QALY), and **screening every 2 years** is optimal, with the new treatment
dominated. The decision turns on the sojourn time (`p.progression`) and on the
cost of false positives (`screening.specificity`, `cost.false.positive`), which
move it to screening every year; and on `p.progressive`, 20% less of which,
breaking the constraints, makes the new treatment optimal again.

Every one of the ten unconstrained runs of `analysis.R`, from random starting
points, fits the target to within 1e-13, ends outside the plausible region and
recommends the new treatment; every one of the ten constrained runs fits it just
as well and recommends screening every 2 years.

## In THALASSA

To see it in the Calibration tab, run each scheme with **simulated annealing and
1,000 iterations**, the other settings at their defaults; a run takes a few
seconds. Over ten seeds, every `unconstrained` run fits the target to an error
below 2e-4 with 57% to 74% of the lesions indolent in its worst age group, and
its parameters make the new treatment optimal; every `constrained` run fits it as
well with at most 20% indolent in every age group, and makes screening every 2
years optimal. Pin the best evaluation and use the calibrated parameters in the
Base case and DSA tabs to see the conclusions above.

Twelve values are more than the defaults of the tab can search:

- The default, Bayesian optimization with 60 iterations, leaves the constrained
  scheme where it started: almost nothing it draws meets all eleven constraints
  at once, and it never improves on the initial guess. On the unconstrained
  scheme it reaches errors around 1, with parameters as implausible as those
  above. TuRBO, SCBO and constraint modeling do not fix the constrained scheme
  at that budget either.
- Nelder-Mead with 1,000 iterations fits both schemes, but searches locally: most
  of its unconstrained runs end close to where they started, still implausible
  in some age group but not enough to change the decision. That is a lesson too:
  without constraints, which of the equally good fits a calibration returns, and
  so what it recommends, depends on the optimizer.

## Reproducing it outside THALASSA

`analysis.R` runs the whole demonstration with base R: the exact fits with the
same share of indolent lesions in every age group, both schemes calibrated from
ten starting points each with `stats::optim()` and `stats::constrOptim()` within
the same search box, and the cost-effectiveness analysis and DSA of the three
parameter sets. Run `Rscript analysis.R` from this directory; it prints the
conclusions above and writes them, with the calibrated parameters, the CE planes,
the tornado plots and the exact fits, to `results/`.
