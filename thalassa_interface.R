library(ggplot2)

source('model.R')

get.overview <- function() {
  # Markdown shown in the Overview tab (see get.overview.markdown() in shiny-cea).
  # The text lives in overview.md. Everything described there is derived from
  # simulate() in model.R and from the rest of this interface, so it must be kept
  # in sync with them.
  return(paste(readLines('overview.md'), collapse='\n'))
}

get.model.settings <- function() {
  return(list(cost.unit='€', effectiveness.unit='QALY'))
}

get.strategies <- function() {
  # Hardcoded strategies for the model. In a real application, these could be loaded from a file or database.
  # The descriptions and attributes describe what each strategy changes in
  # simulate() in model.R, and must be kept in sync with it.
  screening <- function(interval) {
    every <- if (interval == 1) 'year' else sprintf('%d years', interval)
    list(
      name=sprintf('screening_%dy', interval),
      display.name=sprintf('Screening every %s', every),
      description=sprintf('Everyone alive and undiagnosed is tested every %s from age 50 to 74. A share screening.sensitivity of the preclinical lesions is found and treated at once: the progressive ones become screen-detected cancers, the indolent ones overdiagnosed cancers. A share 1 - screening.specificity of the healthy need a work-up.', every),
      attributes=list(intervention='Screening', interval=sprintf('%d year%s', interval, if (interval == 1) '' else 's'))
    )
  }
  return(list(
    list(
      name='no_intervention',
      display.name='No intervention',
      description='Natural history of the cancer: every progressive lesion is diagnosed when it becomes symptomatic and receives standard care, and the indolent ones are never found. The reference the other strategies are compared against.',
      attributes=list(intervention='None')
    ),
    screening(5),
    screening(2),
    screening(1),
    list(
      name='new_treatment',
      display.name='New treatment',
      description='No screening. Clinical cancers receive a new treatment on top of standard care, which multiplies their annual probability of cancer death by hr.new.treatment, at cost.new.treatment per year of treatment.',
      attributes=list(intervention='Treatment')
    )
  ))
}

get.strategy.attributes <- function() {
  # What describes the strategies, shown as columns of the Strategies tab. The
  # intervention is drawn as the shape of the point on the base case plot and
  # the screening interval as its fill, which is why every shape is one with an
  # inside to fill.
  return(list(
    intervention=list(
      label='Intervention',
      plot='shape',
      values=c(None='triangle', Screening='circle', Treatment='square')
    ),
    interval=list(
      label='Screening interval',
      plot='fill',
      ordered=TRUE,
      levels=c('5 years', '2 years', '1 year')
    )
  ))
}

# Constraints: each returns TRUE when the value is acceptable, or the message
# shown for it when it is not.
screen.detected.death.below.clinical.death <- function(par.name, params) {
  if (any(unlist(params[['p.screen.detected.death']]) >= unlist(params[['p.clinical.death']])))
    'The cancer death probability of a screen-detected cancer must be below that of a clinical one'
  else TRUE
}

# The uncalibrated onset by age group: a risk that rises with age, the kind of
# guess a first look at the literature gives, and too high in the younger age
# groups.
BASE.ONSET <- c(`40-49`=0.0025, `50-54`=0.003, `55-59`=0.0035, `60-64`=0.004, `65-69`=0.0045, `70-79`=0.005)

get.parameters <- function() {
  # Hardcoded parameters for the model. In a real application, these could be loaded from a file or database.
  #
  # p.onset and p.progressive are the parameters the calibration estimates. Their
  # base values are the uncalibrated starting point, and they have no dsa.range,
  # so that the DSA moves each age group of a calibrated profile by the same
  # proportion rather than flattening it onto two values. p.onset is stratified,
  # one entry per age group.
  onset <- lapply(AGE.GROUPS, function(group) list(
    name='p.onset',
    display.name='Annual probability of developing a preclinical lesion while healthy',
    stratum=group,
    base.value=BASE.ONSET[[group]],
    distribution='beta',
    class='Natural history',
    min.value=0,
    max.value=1
  ))
  return(c(onset, list(
    list(
      name='p.progressive',
      display.name='Share of new preclinical lesions that are progressive (the rest are indolent)',
      base.value=0.9,
      distribution='beta',
      class='Natural history',
      min.value=0,
      max.value=1
    ),
    list(
      name='p.progression',
      display.name='Annual probability of a progressive lesion becoming symptomatic (1 / mean sojourn time)',
      base.value=0.33,
      distribution='beta',
      class='Natural history',
      min.value=0,
      max.value=1,
      dsa.range=c(0.25, 0.5)
    ),
    list(
      name='p.clinical.death',
      display.name='Annual probability of cancer death of a clinical cancer under standard care',
      base.value=0.12,
      distribution='beta',
      class='Natural history',
      min.value=0,
      max.value=1,
      dsa.range=c(0.08, 0.16),
      constraints=list(screen.detected.death.below.clinical.death)
    ),
    list(
      name='p.screen.detected.death',
      display.name='Annual probability of cancer death of a screen-detected cancer',
      base.value=0.02,
      distribution='beta',
      class='Natural history',
      min.value=0,
      max.value=1,
      dsa.range=c(0.01, 0.04),
      constraints=list(screen.detected.death.below.clinical.death)
    ),
    list(
      name='screening.sensitivity',
      display.name='Sensitivity of the screening test for a preclinical lesion',
      base.value=0.75,
      distribution='beta',
      class='Screening',
      min.value=0,
      max.value=1,
      dsa.range=c(0.6, 0.9)
    ),
    list(
      name='screening.specificity',
      display.name='Specificity of the screening test',
      base.value=0.95,
      distribution='beta',
      class='Screening',
      min.value=0,
      max.value=1,
      dsa.range=c(0.9, 0.98)
    ),
    list(
      name='cost.screening',
      display.name='Cost of a screening test',
      base.value=70,
      distribution='gamma',
      class='Screening',
      min.value=0,
      dsa.range=c(50, 100)
    ),
    list(
      name='cost.false.positive',
      display.name='Cost of the work-up of a false positive test',
      base.value=1000,
      distribution='gamma',
      class='Screening',
      min.value=0,
      dsa.range=c(500, 2000)
    ),
    list(
      name='cost.early.treatment',
      display.name='One-off cost of treating a screen-detected cancer',
      base.value=15000,
      distribution='gamma',
      class='Screening',
      min.value=0,
      dsa.range=c(10000, 20000)
    ),
    list(
      name='cost.screen.detected',
      display.name='Annual follow-up cost of a screen-detected cancer',
      base.value=1000,
      distribution='gamma',
      class='Screening',
      min.value=0,
      dsa.range=c(500, 1500)
    ),
    list(
      name='utility.screen.detected',
      display.name='Utility of a year lived after a screen-detected cancer',
      base.value=0.9,
      distribution='beta',
      class='Screening',
      min.value=0,
      max.value=1,
      dsa.range=c(0.85, 0.95)
    ),
    list(
      name='hr.new.treatment',
      display.name='Hazard ratio of cancer death of the new treatment',
      base.value=0.6,
      distribution='lognormal',
      class='Treatment',
      min.value=0,
      max.value=1,
      dsa.range=c(0.45, 0.75)
    ),
    list(
      name='cost.new.treatment',
      display.name='Annual cost of the new treatment per person with a clinical cancer',
      base.value=1500,
      distribution='gamma',
      class='Treatment',
      min.value=0,
      dsa.range=c(1000, 2500)
    ),
    list(
      name='cost.clinical',
      display.name='Annual cost of standard care per person with a clinical cancer',
      base.value=6000,
      distribution='gamma',
      class='Treatment',
      min.value=0,
      dsa.range=c(4000, 8000)
    ),
    list(
      name='utility.clinical',
      display.name='Utility of a year lived with a clinical cancer',
      base.value=0.65,
      distribution='beta',
      class='Treatment',
      min.value=0,
      max.value=1,
      dsa.range=c(0.55, 0.75)
    ),
    list(
      name='discount',
      display.name='Discount rate',
      base.value=0.03,
      class='General',
      min.value=0,
      max.value=1
    )
  )))
}

get.strata <- function() {
  # Hardcoded strata for the model. In a real application, these could be loaded from a file or database.
  return(AGE.GROUPS)
}

get.model.states <- function() {
  # State diagram shown in the Overview tab, next to the text of overview.md. It
  # describes the transition matrix of simulate() in model.R and must be kept in
  # sync with it.
  return(list(
    title='Model states',
    description='The seven states of the cohort and the annual transitions between them. Hover a state or a transition for what drives it.',
    nodes=data.frame(
      id=c('healthy', 'progressive', 'indolent', 'clinical', 'screen.detected', 'overdiagnosed', 'dead'),
      label=c('Healthy', 'Progressive\nlesion', 'Indolent\nlesion', 'Clinical', 'Screen-\ndetected', 'Over-\ndiagnosed', 'Dead'),
      description=c(
        'Alive and cancer-free. The whole cohort starts here, at age 40.',
        'Alive with a preclinical lesion that will become symptomatic unless found first.',
        'Alive with a preclinical lesion that will never become symptomatic.',
        'Alive with a cancer diagnosed because of its symptoms.',
        'Alive with a progressive lesion found by screening and treated early.',
        'Alive with an indolent lesion found by screening and treated as a cancer, with no benefit.',
        'Absorbing state, reached from every other state.'
      ),
      group=c('Undiagnosed', 'Undiagnosed', 'Undiagnosed', 'Diagnosed', 'Diagnosed', 'Diagnosed', 'Dead')
    ),
    edges=data.frame(
      from=c('healthy', 'healthy', 'progressive', 'progressive', 'indolent',
             'healthy', 'progressive', 'indolent', 'clinical', 'screen.detected', 'overdiagnosed'),
      to=c('progressive', 'indolent', 'clinical', 'screen.detected', 'overdiagnosed',
           'dead', 'dead', 'dead', 'dead', 'dead', 'dead'),
      label=c('Onset', 'Onset', 'Symptoms', 'Detection', 'Detection',
              'Other causes', 'Other causes', 'Other causes', 'Cancer or other causes', 'Cancer or other causes', 'Other causes'),
      description=c(
        'p.onset * p.progressive, by age group once calibrated.',
        'p.onset * (1 - p.progressive), by age group once calibrated.',
        'p.progression. Its inverse is the mean sojourn time.',
        'At each screening round, a share screening.sensitivity of the progressive lesions.',
        'At each screening round, a share screening.sensitivity of the indolent lesions.',
        'Background mortality by age.',
        'Background mortality by age.',
        'Background mortality by age.',
        'Background mortality plus p.clinical.death, multiplied by hr.new.treatment under new_treatment.',
        'Background mortality plus p.screen.detected.death.',
        'Background mortality by age: an indolent lesion never kills.'
      )
    )
  ))
}

get.code.sample <- function() {
  # Code shown in the Code panel of the Overview tab, one tab per entry, for
  # whoever wants to see how the model works. The whole of model.R is written to
  # be read, so it is shown as it is, and the constraints of the calibration
  # next to it, since they are what this model is about: the lines of this file
  # from the largest plausible share of indolent lesions to the end of
  # calibration.constraints().
  interface <- readLines('thalassa_interface.R')
  first <- grep('^# The largest share of indolent lesions', interface)
  constraints.start <- grep('^calibration.constraints <- function', interface)
  last <- constraints.start + which(interface[constraints.start:length(interface)] == '}')[1] - 1
  return(list(
    `model.R`=paste(readLines('model.R'), collapse='\n'),
    `Calibration constraints`=paste(interface[first:last], collapse='\n')
  ))
}

run.simulation <- function(strategies, pars) {
  # The pars vector should be transformed to the format expected by the simulate function.
  # This is a simple mapping based on the parameter names.
  results <- simulate(strategies,
                      p.onset=pars[['p.onset']],
                      p.progressive=pars[['p.progressive']],
                      p.progression=pars[['p.progression']],
                      p.clinical.death=pars[['p.clinical.death']],
                      p.screen.detected.death=pars[['p.screen.detected.death']],
                      screening.sensitivity=pars[['screening.sensitivity']],
                      screening.specificity=pars[['screening.specificity']],
                      hr.new.treatment=pars[['hr.new.treatment']],
                      cost.screening=pars[['cost.screening']],
                      cost.false.positive=pars[['cost.false.positive']],
                      cost.early.treatment=pars[['cost.early.treatment']],
                      cost.screen.detected=pars[['cost.screen.detected']],
                      cost.clinical=pars[['cost.clinical']],
                      cost.new.treatment=pars[['cost.new.treatment']],
                      utility.clinical=pars[['utility.clinical']],
                      utility.screen.detected=pars[['utility.screen.detected']],
                      discount=pars[['discount']])
  return(results)
}

# ---- Calibration ---------------------------------------------------------------

# What the calibration is fitted to: clinical incidence by age group, from a
# cancer registry, before any screening existed. Only progressive lesions ever
# become clinical, so this target says how many progressive lesions there are,
# p.onset * p.progressive, and nothing about the indolent ones: twice the onset
# with half of the lesions indolent fits it just as well.
CALIBRATION.TARGETS <- list(
  `Clinical incidence (per 100,000)`=list(`40-49`=56, `50-54`=107, `55-59`=158, `60-64`=212, `65-69`=269, `70-79`=348)
)

calibration.error <- function(pars, target) {
  calibration.strategy <- 'no_intervention'
  result <- tryCatch({
    results <- run.simulation(calibration.strategy, pars)
    simulated <- list(`Clinical incidence (per 100,000)`=results$incidence[[calibration.strategy]])
    # The target is a named list with one entry per stratum, so it is flattened
    # into a named numeric vector to match the simulated values by stratum name.
    # The error is the sum of the squared relative differences.
    error <- 0
    for (name in names(target)) {
      observed <- unlist(target[[name]])
      error <- error + sum((simulated[[name]][names(observed)] / observed - 1)^2)
    }
    if (!is.finite(error)) stop('The simulation produced no usable output')
    list(error=error, output=simulated)
  }, error=function(e) {
    list(error=Inf,
         output=list(`Clinical incidence (per 100,000)`=setNames(rep(NA_real_, length(get.strata())), get.strata())))
  })
  return(result)
}

# The largest share of indolent lesions that external evidence allows. Screening
# trials measure it as the excess of cancers diagnosed in the screened arm that
# is still there years after screening stopped, and put it below 20%.
INDOLENT.MAX <- 0.2

# The feasible region of the constrained scheme. x is the calibration vector,
# p.onset of every age group followed by p.progressive of every age group, or,
# under a stepwise run (partial = TRUE), the same layout over the first age
# groups only. Every value must be zero or below:
#   - at most INDOLENT.MAX of the lesions are indolent, i.e. p.progressive is
#     at least 1 - INDOLENT.MAX, in every age group;
#   - p.onset does not decrease with age, since the risk of developing a
#     solid cancer only grows as people age.
# The onset differences are expressed per 1,000 people, so that both kinds of
# constraint are numbers of a similar size.
calibration.constraints <- function(x, partial = FALSE) {
  n <- length(x) %/% 2
  onset <- x[seq_len(n)]
  progressive <- x[n + seq_len(n)]
  c((1 - INDOLENT.MAX) - progressive,
    if (n > 1) 1000 * (onset[-n] - onset[-1]))
}

calibration.constraints.description <- function() {
  sprintf(paste(
    'The calibration vector is p.onset for the age groups %s, followed by p.progressive for the same age groups.',
    'In every age group p.progressive must be at least %.2f: screening trials show that at most %d%% of the lesions are indolent, i.e. never become symptomatic.',
    'p.onset must not decrease from one age group to the next, since the risk of developing the cancer grows with age.'),
    paste(AGE.GROUPS, collapse=', '), 1 - INDOLENT.MAX, round(100 * INDOLENT.MAX))
}

generate.training.dataset <- function(initial_guess, n, ...) {
  # Calibration vectors to train a latent space on, drawn around the initial
  # guess within the search box, initial_guess * (1 ± variation), and below 1.
  f.pars <- list(...)
  variation <- f.pars$variation

  n_params <- length(initial_guess)

  dataset <- matrix(NA, nrow=n, ncol=n_params)

  for(i in 1:n) {
    factors <- runif(n_params, min=1-variation, max=1+variation)
    dataset[i,] <- pmin(1, initial_guess * factors)
  }

  dataset <- dataset[sample(nrow(dataset)),]

  return(dataset)
}

generate.feasible.training.dataset <- function(initial_guess, n, ...) {
  # The same box, but only vectors that meet calibration.constraints(), so that
  # a latent space trained on them carries the constraints: p.progressive
  # within the plausible range, and p.onset sorted so that it rises with age.
  f.pars <- list(...)
  variation <- f.pars$variation

  n.groups <- length(initial_guess) %/% 2
  onset.guess <- initial_guess[seq_len(n.groups)]
  progressive.guess <- initial_guess[n.groups + seq_len(n.groups)]
  lower <- function(guess) pmax(0, guess * (1 - variation))
  upper <- function(guess) pmin(1, guess * (1 + variation))

  # Where the box and the plausible range do not overlap, the end of the box
  # closest to the range is used, and the vector is as feasible as the box allows.
  progressive.lower <- pmin(pmax(lower(progressive.guess), 1 - INDOLENT.MAX), upper(progressive.guess))
  progressive.upper <- upper(progressive.guess)

  dataset <- matrix(NA, nrow=n, ncol=length(initial_guess))
  for (i in 1:n) {
    onset <- sort(runif(n.groups, min=lower(onset.guess), max=upper(onset.guess)))
    onset <- cummax(pmin(pmax(onset, lower(onset.guess)), upper(onset.guess)))
    progressive <- runif(n.groups, min=progressive.lower, max=progressive.upper)
    dataset[i,] <- c(onset, progressive)
  }

  return(dataset)
}

# Additional plots of the calibration tab, drawn for the evaluation hovered or
# pinned on the error plot. Each receives the values of that evaluation that
# belong to an age group, named '<parameter>.<age group>' for the calibrated
# parameters and 'sim_<target>.<age group>' for the simulated outputs.
#
# The tab draws one evaluation and moves the plot to another by replacing the y
# values of its traces, one by one, and fitting the y axis to the values it
# moves in. So both plots draw the same traces whatever the evaluation, and
# carry an invisible one spanning the whole scale, which keeps the axis, and
# its ticks, where they are as the plot moves.
#
# The ggplot is built here, in the process of the model, but drawn by the app,
# in a process of its own where none of this file exists. ggplot evaluates the
# aesthetics only when it draws, so everything they read has to be a column of
# the data of the plot, never a variable of the model.
invisible.range <- function(x, low, high)
  geom_point(data=data.frame(x=x, y=c(low, high)), aes(x=x, y=y), inherit.aes=FALSE, alpha=0)

# The share of indolent lesions the evaluation implies in every age group,
# against the largest share screening trials allow. It is what tells a
# plausible fit apart from another that fits the target just as well. Each age
# group is a bar in both layers, one of them of height 0, so that the colour of
# a bar can change without the traces changing.
indolent.share.plot <- function(data) {
  if (is.null(data)) return(NULL)
  progressive <- data[paste0('p.progressive.', AGE.GROUPS)]
  if (all(is.na(progressive))) return(NULL)
  indolent <- 100 * (1 - as.numeric(progressive))
  plausible <- indolent <= 100 * INDOLENT.MAX
  plt.df <- data.frame(stratum=factor(AGE.GROUPS, levels=AGE.GROUPS),
                       plausible=ifelse(plausible, indolent, 0),
                       implausible=ifelse(plausible, 0, indolent),
                       text=sprintf('%s: %.0f%% of the lesions indolent (%s)', AGE.GROUPS, indolent,
                                    ifelse(plausible, 'plausible', 'implausible')))
  ggplot(plt.df, aes(x=stratum, text=text)) +
    geom_col(aes(y=plausible, fill='Plausible'), width=0.5) +
    geom_col(aes(y=implausible, fill='Implausible'), width=0.5) +
    geom_hline(data=data.frame(limit=100 * INDOLENT.MAX, label='Largest share trials allow'),
               aes(yintercept=limit, linetype=label), color='grey30', inherit.aes=FALSE) +
    invisible.range(AGE.GROUPS[1], 0, 100) +
    scale_fill_manual(values=c(Plausible='forestgreen', Implausible='firebrick'), breaks=c('Plausible', 'Implausible')) +
    scale_linetype_manual(values=c(`Largest share trials allow`='dashed')) +
    scale_y_continuous(breaks=seq(0, 100, 20)) +
    labs(x='Age group', y='Indolent lesions (%)', fill='', linetype='') +
    theme_minimal()
}

# The onset of all lesions, and of the progressive ones alone, in every age
# group. The second is what the target pins down; the first is what the
# calibration is free to choose. The scale reaches the top of the search box at
# the default scope, twice the base onset of the oldest age group.
onset.plot <- function(data) {
  if (is.null(data)) return(NULL)
  onset <- as.numeric(data[paste0('p.onset.', AGE.GROUPS)])
  progressive <- as.numeric(data[paste0('p.progressive.', AGE.GROUPS)])
  if (all(is.na(onset)) || all(is.na(progressive))) return(NULL)
  plt.df <- data.frame(stratum=factor(rep(AGE.GROUPS, 2), levels=AGE.GROUPS),
                       kind=rep(c('All lesions', 'Progressive lesions'), each=length(AGE.GROUPS)),
                       onset=1000 * c(onset, onset * progressive))
  plt.df$text <- sprintf('%s, %s: %.2f per 1,000 a year', plt.df$kind, plt.df$stratum, plt.df$onset)
  top <- 2000 * max(BASE.ONSET)
  ggplot(plt.df, aes(x=stratum, y=onset, color=kind, group=kind, text=text)) +
    geom_line() +
    geom_point(size=3) +
    invisible.range(AGE.GROUPS[1], 0, top) +
    scale_y_continuous(breaks=pretty(c(0, top))) +
    labs(x='Age group', y='Onset (per 1,000 a year)', color='') +
    theme_minimal()
}

get.calibration.schemes <- function() {
  # Both schemes search the same twelve values against the same target with the
  # same error, and the target is met exactly all along a surface of them. The
  # only difference is that the second one is told what external evidence says
  # about the parameters, which rules most of that surface out (see overview.md).
  common <- list(
    parameters=c('p.onset', 'p.progressive'),
    target=CALIBRATION.TARGETS,
    strata=get.strata(),
    initial_guess=c(unname(BASE.ONSET[AGE.GROUPS]), rep(0.9, length(AGE.GROUPS))),
    error_function=calibration.error,
    evaluation.plots=list(
      `Indolent lesions`=indolent.share.plot,
      `Onset`=onset.plot
    )
  )
  return(list(
    unconstrained=c(list(
      description='Onset and share of progressive lesions by age group against incidence, unconstrained (12 values, 6 targets)',
      latent_space_training_set=generate.training.dataset
    ), common),
    constrained=c(list(
      description='Onset and share of progressive lesions by age group against incidence, with plausibility constraints (12 values, 6 targets)',
      constraints=calibration.constraints,
      constraints_description=calibration.constraints.description,
      latent_space_training_set=generate.feasible.training.dataset
    ), common)
  ))
}
