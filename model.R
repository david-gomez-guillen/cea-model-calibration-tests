# A cost-effectiveness model of screening for a hypothetical cancer.
#
# It is a Markov cohort model. Rather than following people one by one, it keeps
# track of which share of a group of people, the cohort, is in each health state,
# and moves those shares from state to state once a year.
#
#   States:      Healthy, Progressive, Indolent, Clinical, Screen-detected,
#                Overdiagnosed, Dead
#   Cycle:       one year
#   Horizon:     from age 40 to age 79
#   Strategies:  no_intervention, screening_5y, screening_2y, screening_1y,
#                new_treatment
#
# The cancer starts silently, as a preclinical lesion that causes no symptoms
# but that a screening test can find. Lesions are of two kinds:
#
#   - progressive ones become symptomatic after some years, which is when the
#     cancer is diagnosed if nobody looked for it earlier (Clinical);
#   - indolent ones never do. Left alone they harm nobody, but a screening test
#     finds them as well as the progressive ones, and whoever has one found is
#     treated as a cancer patient with no benefit from it: overdiagnosis.
#
# A cancer registry only sees the progressive lesions, once they are clinical,
# so it cannot tell how many indolent ones there are. That is what this model
# is about (see overview.md).
#
# For each strategy the model adds up what the cohort costs and how much health
# it enjoys. Comparing those two totals between strategies is what a
# cost-effectiveness analysis does.

START.AGE <- 40
END.AGE <- 79

# Results by age are reported per age group, and the calibrated parameters take
# one value per age group: five-year groups through the screening ages, wider
# ones before and after them. A group runs from its start to the start of the
# next one.
AGE.GROUPS <- c('40-49', '50-54', '55-59', '60-64', '65-69', '70-79')
AGE.GROUP.STARTS <- c(40, 50, 55, 60, 65, 70)
age.group <- function(age) findInterval(age, AGE.GROUP.STARTS)

# Screening is offered from age 50 to age 74, every 5, 2 or 1 years.
SCREENING.START.AGE <- 50
SCREENING.END.AGE <- 74

STRATEGIES <- c('no_intervention', 'screening_5y', 'screening_2y', 'screening_1y', 'new_treatment')

STATES <- c('healthy', 'progressive', 'indolent', 'clinical', 'screen.detected', 'overdiagnosed', 'dead')
ALIVE <- setdiff(STATES, 'dead')

# Annual probability of dying of anything other than this cancer, by age. A
# Gompertz curve, close to the mortality of a European population.
background.mortality <- function(age) 0.0008 * exp(0.095 * (age - START.AGE))

# A parameter is either one value for every age or a list with one value per
# age group, in which case the one of the given age is used. A calibrated
# parameter arrives as such a list, keyed by age group.
value.at.age <- function(parameter, age) {
  if (!is.list(parameter)) return(parameter)
  group <- age.group(age)
  if (!is.null(names(parameter))) parameter[[AGE.GROUPS[group]]] else parameter[[group]]
}

simulate <- function(strategies,
                     p.onset,
                     p.progressive,
                     p.progression,
                     p.clinical.death,
                     p.screen.detected.death,
                     screening.sensitivity,
                     screening.specificity,
                     hr.new.treatment,
                     cost.screening,
                     cost.false.positive,
                     cost.early.treatment,
                     cost.screen.detected,
                     cost.clinical,
                     cost.new.treatment,
                     utility.clinical,
                     utility.screen.detected,
                     discount) {

  ages <- START.AGE:END.AGE

  summary <- data.frame()
  incidence <- list()
  mortality <- list()
  outcomes <- data.frame()
  cohort.info <- list()

  for (strategy in strategies) {

    # ---- What the strategy changes ----------------------------------------
    # Screening tests everyone alive and undiagnosed at the screening ages, and
    # treats the lesions it finds at once: the progressive ones become
    # Screen-detected cancers, much less deadly than clinical ones, and the
    # indolent ones Overdiagnosed cancers, treated just the same. The new
    # treatment lowers the mortality of clinical cancers, at a cost per year of
    # treatment.
    if (!strategy %in% STRATEGIES)
      stop('Unknown strategy: ', strategy)

    screening.interval <- switch(strategy, screening_5y = 5, screening_2y = 2, screening_1y = 1, NA)
    new.treatment <- strategy == 'new_treatment'

    p.death.clinical <- if (new.treatment) p.clinical.death * hr.new.treatment else p.clinical.death

    # What a person in each state costs in a year, and the quality of life of
    # that year: 1 is a year in full health and 0 is being dead. An overdiagnosed
    # cancer costs and weighs what a screen-detected one does: nobody can tell
    # them apart, so both are treated and followed up alike.
    state.costs <- c(healthy = 0, progressive = 0, indolent = 0,
                     clinical = cost.clinical + if (new.treatment) cost.new.treatment else 0,
                     screen.detected = cost.screen.detected,
                     overdiagnosed = cost.screen.detected,
                     dead = 0)
    state.utilities <- c(healthy = 1, progressive = 1, indolent = 1,
                         clinical = utility.clinical,
                         screen.detected = utility.screen.detected,
                         overdiagnosed = utility.screen.detected,
                         dead = 0)

    # ---- The cohort ---------------------------------------------------------
    # The share of the cohort in each state. Everyone starts healthy, and the
    # shares always add up to 1.
    cohort <- setNames(c(1, rep(0, length(STATES) - 1)), STATES)

    # One value per year, filled in as the cohort ages.
    costs <- c()
    utilities <- c()
    yearly.incidence <- c()
    yearly.mortality <- c()
    total <- c(screen.detected = 0, overdiagnosed = 0, clinical = 0, cancer.deaths = 0)
    cohort.trace <- list(cohort)

    for (age in ages) {
      years.elapsed <- age - START.AGE
      discount.factor <- 1 / (1 + discount)^years.elapsed
      screening.cost <- 0

      # ---- 1. Screening round ---------------------------------------------
      # Everyone not yet diagnosed is tested. A share screening.sensitivity of
      # the lesions, progressive or indolent, is found and treated at once; a
      # share 1 - screening.specificity of the healthy test positive and need a
      # work-up to rule the cancer out.
      if (!is.na(screening.interval) &&
          age >= SCREENING.START.AGE && age <= SCREENING.END.AGE &&
          (age - SCREENING.START.AGE) %% screening.interval == 0) {
        found.progressive <- cohort[['progressive']] * screening.sensitivity
        found.indolent <- cohort[['indolent']] * screening.sensitivity
        false.positives <- cohort[['healthy']] * (1 - screening.specificity)
        screening.cost <- sum(cohort[c('healthy', 'progressive', 'indolent')]) * cost.screening +
          false.positives * cost.false.positive +
          (found.progressive + found.indolent) * cost.early.treatment
        cohort[['progressive']] <- cohort[['progressive']] - found.progressive
        cohort[['indolent']] <- cohort[['indolent']] - found.indolent
        cohort[['screen.detected']] <- cohort[['screen.detected']] + found.progressive
        cohort[['overdiagnosed']] <- cohort[['overdiagnosed']] + found.indolent
        total[['screen.detected']] <- total[['screen.detected']] + found.progressive
        total[['overdiagnosed']] <- total[['overdiagnosed']] + found.indolent
      }

      # ---- 2. Transition matrix -------------------------------------------
      # Row: the state a person is in now. Column: the state a year later.
      # Each row adds up to 1, since everyone has to end up somewhere, so the
      # probability of staying is whatever the other transitions leave.
      # Dead is an absorbing state: once there, nobody leaves.
      p.other <- background.mortality(age)
      p.on <- value.at.age(p.onset, age)
      share.progressive <- value.at.age(p.progressive, age)
      p.prog <- value.at.age(p.progression, age)
      p.die.clinical <- min(1, p.other + p.death.clinical)
      p.die.screen <- min(1, p.other + p.screen.detected.death)

      transitions <- matrix(0, nrow = length(STATES), ncol = length(STATES),
                            dimnames = list(STATES, STATES))
      transitions['healthy', 'progressive'] <- p.on * share.progressive
      transitions['healthy', 'indolent'] <- p.on * (1 - share.progressive)
      transitions['progressive', 'clinical'] <- p.prog
      transitions[c('healthy', 'progressive', 'indolent', 'overdiagnosed'), 'dead'] <- p.other
      transitions['clinical', 'dead'] <- p.die.clinical
      transitions['screen.detected', 'dead'] <- p.die.screen
      transitions['dead', 'dead'] <- 1
      for (state in ALIVE) transitions[state, state] <- 1 - sum(transitions[state, ])

      if (any(transitions < 0) || any(transitions > 1))
        stop('Transition probabilities out of [0, 1] at age ', age)

      # ---- 3. Costs and health of this year ---------------------------------
      # Each state contributes its cost and its utility in proportion to the
      # share of the cohort in it. Money and health count for less the further
      # in the future they are, which is what discounting expresses: at a rate
      # of 3%, what happens a year from now is worth 1 / 1.03 of the same today.
      costs <- c(costs, (sum(cohort * state.costs) + screening.cost) * discount.factor)
      utilities <- c(utilities, sum(cohort * state.utilities) * discount.factor)

      # ---- 4. Incidence and mortality of this year ----------------------------
      # New clinical cases and cancer deaths among the people alive, as a
      # cancer registry reports them: per 100,000 people a year.
      alive <- sum(cohort[ALIVE])
      new.clinical <- cohort[['progressive']] * p.prog
      cancer.deaths <- cohort[['clinical']] * p.death.clinical +
        cohort[['screen.detected']] * p.screen.detected.death
      yearly.incidence <- c(yearly.incidence, 1e5 * new.clinical / alive)
      yearly.mortality <- c(yearly.mortality, 1e5 * cancer.deaths / alive)
      total[['clinical']] <- total[['clinical']] + new.clinical
      total[['cancer.deaths']] <- total[['cancer.deaths']] + cancer.deaths

      # ---- 5. Move the cohort one year forward ------------------------------
      # Multiplying the shares by the matrix sends each share to where its row
      # says, which gives the shares at the start of next year.
      cohort <- drop(cohort %*% transitions)
      cohort.trace <- c(cohort.trace, list(cohort))
    }

    # ---- Results of the strategy ----------------------------------------------
    # C is the discounted cost per person and E the discounted quality-adjusted
    # life years (QALYs) per person, both accumulated over the whole horizon.
    summary <- rbind(summary, data.frame(strategy = strategy,
                                         C = sum(costs),
                                         E = sum(utilities)))

    # Incidence and mortality per age group: the average over the years of the
    # group.
    by.group <- function(x) setNames(as.numeric(tapply(x, age.group(ages), mean)), AGE.GROUPS)
    incidence[[strategy]] <- by.group(yearly.incidence)
    mortality[[strategy]] <- by.group(yearly.mortality)

    # Undiscounted counts per 1,000 people of the cohort.
    outcomes <- rbind(outcomes, data.frame(strategy = strategy, t(1000 * total)))

    cohort.info[[strategy]] <- cohort.trace
  }

  list(summary = summary,
       incidence = incidence,
       mortality = mortality,
       outcomes = outcomes,
       cohort.info = cohort.info)
}
