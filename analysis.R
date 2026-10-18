# The whole demonstration of this model, outside THALASSA.
#
#   Rscript analysis.R            # from the directory of the model
#
# 1. The ridge: parameter sets that fit the calibration target exactly, from
#    one where every lesion is progressive to one where half are indolent, and
#    what the cost-effectiveness analysis concludes at each.
# 2. The calibrations: both schemes run from several starting points, as an
#    optimizer would, and where on the ridge each run ends.
# 3. The cost-effectiveness analysis on the uncalibrated parameters, on the
#    unconstrained calibration and on the constrained one, each with its
#    deterministic sensitivity analysis.
#
# It prints the conclusions and writes them to results/. Only base R is needed;
# ggplot2, when installed, adds the cost-effectiveness planes and the tornados.
#
# The calibrations are run with stats::optim() on the problem the Calibration
# tab of THALASSA builds with its default search scope of 100%: the error
# function and the target of the schemes, every value searched between 0 and
# twice its initial guess, and at most 1. The constraints of the constrained
# scheme are linear in the calibration vector, so stats::constrOptim() takes
# them as they are.

source('thalassa_interface.R')

WTP <- 25000                 # €/QALY, the default of THALASSA
DSA.VARIATION <- 0.2         # for the parameters with no dsa.range, as in the DSA tab
OUTPUT.DIR <- 'results'
N.STARTS <- 10               # starting points of each calibration

set.seed(12345)
dir.create(OUTPUT.DIR, showWarnings=FALSE)

# ---- Parameters ----------------------------------------------------------------

# One entry per parameter: a stratified one, declared once per age group, is
# merged into a list of its values keyed by age group, as the app passes it.
entries <- get.parameters()
parameter.names <- unique(sapply(entries, function(p) p$name))
param.info <- lapply(parameter.names, function(name) Filter(function(p) p$name == name, entries)[[1]])
base.parameters <- setNames(lapply(parameter.names, function(name) {
  declared <- Filter(function(p) p$name == name, entries)
  if (is.null(declared[[1]]$stratum)) declared[[1]]$base.value
  else setNames(lapply(declared, function(p) p$base.value), sapply(declared, function(p) p$stratum))
}), parameter.names)
schemes <- get.calibration.schemes()
scheme <- schemes$constrained
strata <- get.strata()
n.groups <- length(strata)
onset.slots <- seq_len(n.groups)
progressive.slots <- n.groups + seq_len(n.groups)

# The calibration vector, laid out as the app lays it out, made into the
# parameter list the model runs with.
with.calibration <- function(parameters, x) {
  for (i in seq_along(scheme$parameters))
    parameters[[scheme$parameters[i]]] <- setNames(as.list(x[(i - 1) * n.groups + seq_len(n.groups)]), strata)
  parameters
}

calibration.vector <- function(parameters)
  unlist(lapply(scheme$parameters, function(p) rep_len(unlist(parameters[[p]]), n.groups)))

label <- function(strategy) {
  names <- c(no_intervention='No intervention', screening_5y='Screening every 5 years',
             screening_2y='Screening every 2 years', screening_1y='Screening every year',
             new_treatment='New treatment')
  unname(names[strategy])
}

# ---- Cost-effectiveness ------------------------------------------------------------

net.monetary.benefit <- function(summary) summary$E * WTP - summary$C

# The efficient frontier: strategies sorted by effect, leaving out the dominated
# ones (more costly and less effective than another) and the extendedly
# dominated ones (whose ICER is above that of a more effective strategy).
efficient.frontier <- function(summary) {
  s <- summary[order(summary$E, summary$C), c('strategy', 'C', 'E')]
  dominated <- sapply(seq_len(nrow(s)), function(i)
    any(s$C <= s$C[i] & s$E >= s$E[i] & (s$C < s$C[i] | s$E > s$E[i])))
  s <- s[!dominated, ]
  repeat {
    icer <- diff(s$C) / diff(s$E)
    extended <- which(diff(icer) < 0)
    if (length(extended) == 0) break
    s <- s[-(extended[1] + 1), ]
  }
  s$ICER <- c(NA, diff(s$C) / diff(s$E))
  s
}

cea <- function(parameters) {
  results <- run.simulation(STRATEGIES, parameters)
  s <- results$summary
  s$NMB <- net.monetary.benefit(s)
  s$IC <- s$C - s$C[1]
  s$IE <- s$E - s$E[1]
  s$ICER.vs.no_intervention <- ifelse(seq_len(nrow(s)) == 1, NA, s$IC / s$IE)
  s <- merge(s, results$outcomes, by='strategy', sort=FALSE)
  ranked <- s$strategy[order(-s$NMB)]
  list(summary=s, frontier=efficient.frontier(results$summary),
       optimal=ranked[1], margin=max(s$NMB) - max(s$NMB[s$strategy != ranked[1]]))
}

optimal.strategy <- function(parameters) cea(parameters)$optimal

# ---- Calibration -----------------------------------------------------------------

x0 <- scheme$initial_guess
search.lower <- pmax(0, x0 * 0)
search.upper <- pmin(1, x0 * 2)

# A vector the model cannot simulate (a probability of leaving a state above 1)
# has an infinite error, which the optimizers of stats::optim() need finite.
error.of <- function(x) {
  error <- scheme$error_function(with.calibration(base.parameters, x), scheme$target)$error
  if (is.finite(error)) error else 1e6
}

# Values are optimized on the unit scale of the search box, where every
# coordinate counts the same.
to.box <- function(u) search.lower + u * (search.upper - search.lower)
from.box <- function(x) (x - search.lower) / (search.upper - search.lower)

# Central differences, for the methods that need a gradient.
numerical.gradient <- function(u, h=1e-7) {
  sapply(seq_along(u), function(j) {
    step <- replace(rep(0, length(u)), j, h)
    (error.of(to.box(u + step)) - error.of(to.box(u - step))) / (2 * h)
  })
}

# Unconstrained: from a point drawn anywhere in the search box, as the initial
# design of an optimizer draws them.
calibrate.unconstrained <- function(u0) {
  fit <- optim(u0, function(u) error.of(to.box(u)), gr=numerical.gradient, method='L-BFGS-B',
               lower=0, upper=1, control=list(maxit=1000, factr=1e2))
  list(x=to.box(fit$par), error=fit$value)
}

# Constrained: the constraints, c(x) <= 0, are linear, c(x) = A x + b, and
# constrOptim() takes them as ui u - ci >= 0 on the unit scale, together with
# the search box. It starts strictly inside the feasible region, from points
# the training set generator of the scheme draws.
constraint.matrix <- local({
  b <- scheme$constraints(rep(0, length(x0)))
  A <- sapply(seq_along(x0), function(j) scheme$constraints(replace(rep(0, length(x0)), j, 1)) - b)
  D <- diag(search.upper - search.lower)
  list(ui=rbind(-A %*% D, diag(length(x0)), -diag(length(x0))),
       ci=c(b + A %*% search.lower, rep(0, length(x0)), rep(-1, length(x0))))
})

calibrate.constrained <- function(x.start) {
  u0 <- from.box(x.start)
  # Strictly inside: away from the bounds and with the onset strictly rising.
  u0 <- pmin(pmax(u0, 0.02), 0.98)
  u0[onset.slots] <- sort(u0[onset.slots]) + seq(0, 1e-3, length.out=n.groups)
  u0[progressive.slots] <- pmax(u0[progressive.slots], from.box(rep(1 - INDOLENT.MAX + 0.01, length(x0)))[progressive.slots])
  stopifnot(all(constraint.matrix$ui %*% u0 - constraint.matrix$ci > 0))
  fit <- constrOptim(u0, function(u) error.of(to.box(u)), grad=numerical.gradient,
                     ui=constraint.matrix$ui, ci=constraint.matrix$ci,
                     control=list(maxit=1000, reltol=1e-14), outer.iterations=100, outer.eps=1e-10)
  list(x=to.box(fit$par), error=fit$value)
}

describe.run <- function(x, error) {
  data.frame(error=signif(error, 3),
             onset.per.1000=paste(sprintf('%.2f', 1000 * x[onset.slots]), collapse=' / '),
             indolent.percent=paste(sprintf('%.0f', 100 * (1 - x[progressive.slots])), collapse=' / '),
             plausible=all(scheme$constraints(x) <= 1e-9),
             optimal=label(optimal.strategy(with.calibration(base.parameters, x))))
}

# ---- 1. The ridge --------------------------------------------------------------------
# For each share of progressive lesions, the same in both age groups, the onset
# that fits the target exactly.

cat('\nExact fits of the target, from all lesions progressive to half of them indolent\n\n')
ridge <- do.call(rbind, lapply(c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5), function(share) {
  error.at <- function(onset) error.of(c(onset, rep(share, n.groups)))
  fit <- optim(rep(0.5, n.groups), function(u) error.at(to.box(c(u, rep(0, n.groups)))[onset.slots]),
               method='L-BFGS-B', lower=0, upper=1, control=list(factr=1e2))
  x <- c(to.box(c(fit$par, rep(0, n.groups)))[onset.slots], rep(share, n.groups))
  a <- cea(with.calibration(base.parameters, x))
  data.frame(progressive=share, describe.run(x, fit$value)[, -5],
             optimal=label(a$optimal), margin=round(a$margin),
             overdiagnosed.per.1000=round(a$summary$overdiagnosed[a$summary$strategy == 'screening_2y'], 1))
}))
print(ridge, row.names=FALSE)
write.csv(ridge, file.path(OUTPUT.DIR, 'ridge.csv'), row.names=FALSE)

# ---- 2. The calibrations -------------------------------------------------------------

cat('\nCalibrating, ', N.STARTS, ' starting points per scheme...\n', sep='')
unconstrained.runs <- lapply(seq_len(N.STARTS), function(i) calibrate.unconstrained(runif(length(x0))))
feasible.starts <- generate.feasible.training.dataset(x0, N.STARTS, variation=1)
constrained.runs <- lapply(seq_len(N.STARTS), function(i) calibrate.constrained(feasible.starts[i, ]))

runs <- rbind(
  data.frame(scheme='unconstrained', run=seq_len(N.STARTS),
             do.call(rbind, lapply(unconstrained.runs, function(r) describe.run(r$x, r$error)))),
  data.frame(scheme='constrained', run=seq_len(N.STARTS),
             do.call(rbind, lapply(constrained.runs, function(r) describe.run(r$x, r$error)))))
cat('\nEvery run, and the strategy its parameters make optimal\n\n')
print(runs, row.names=FALSE)
write.csv(runs, file.path(OUTPUT.DIR, 'calibration_runs.csv'), row.names=FALSE)
cat('\nOptimal strategy across the runs of each scheme:\n')
print(table(runs$scheme, runs$optimal))

# The result of a calibration is its best run, as an optimizer reports it.
best.of <- function(runs) runs[[which.min(sapply(runs, function(r) r$error))]]$x
x.unconstrained <- best.of(unconstrained.runs)
x.constrained <- best.of(constrained.runs)

scenarios <- list(
  `Uncalibrated`=base.parameters,
  `Unconstrained calibration`=with.calibration(base.parameters, x.unconstrained),
  `Constrained calibration`=with.calibration(base.parameters, x.constrained)
)

# ---- 3. Deterministic sensitivity analysis ---------------------------------------
# Every parameter is moved to both ends of its dsa.range, or by ± DSA.VARIATION
# when it declares none (the calibrated ones are moved in every age group
# together, which keeps the shape of their profile). Reported per parameter, at
# both ends: by how much the NMB of the optimal strategy beats that of the best
# of the others, which is below zero where the decision changes, and which
# strategy is optimal there.

dsa <- function(parameters, base.cea) {
  optimal <- base.cea$optimal
  rows <- list()
  for (p in param.info) {
    if (p$name == 'discount') next
    value <- parameters[[p$name]]
    ends <- if (!is.null(p$dsa.range)) as.list(p$dsa.range) else
      list(lapply(value, function(v) v * (1 - DSA.VARIATION)), lapply(value, function(v) v * (1 + DSA.VARIATION)))
    if (!is.list(value)) ends <- lapply(ends, unlist)
    nmb <- c(); best <- c()
    for (end in ends) {
      varied <- parameters
      # A proportional move is kept within the range of the parameter.
      cap <- function(v) if (is.null(p$max.value)) v else pmin(p$max.value, v)
      varied[[p$name]] <- if (is.list(value)) lapply(end, cap) else cap(end)
      s <- run.simulation(STRATEGIES, varied)$summary
      s$NMB <- net.monetary.benefit(s)
      nmb <- c(nmb, s$NMB[s$strategy == optimal] - max(s$NMB[s$strategy != optimal]))
      best <- c(best, s$strategy[which.max(s$NMB)])
    }
    rows[[p$name]] <- data.frame(parameter=p$name, low.margin=nmb[1], high.margin=nmb[2],
                                 swing=abs(nmb[2] - nmb[1]), optimal.low=best[1], optimal.high=best[2],
                                 decision.changes=any(best != optimal))
  }
  d <- do.call(rbind, rows)
  d <- d[order(-d$swing), ]
  rownames(d) <- NULL
  attr(d, 'optimal') <- optimal
  d
}

# ---- Run and report ----------------------------------------------------------------

overview <- list()
for (name in names(scenarios)) {
  parameters <- scenarios[[name]]
  analysis <- cea(parameters)
  sensitivity <- dsa(parameters, analysis)
  x <- calibration.vector(parameters)
  calibration <- calibration.error(parameters, CALIBRATION.TARGETS)
  feasible <- all(scheme$constraints(x) <= 1e-9)

  cat('\n', strrep('=', 78), '\n', name, '\n', strrep('=', 78), '\n', sep='')
  cat(sprintf('Calibration error %.2g; onset by age group %s per 1,000; indolent lesions %s%% (%s)\n',
              calibration$error, paste(sprintf('%.2f', 1000 * x[onset.slots]), collapse=' / '),
              paste(sprintf('%.0f', 100 * (1 - x[progressive.slots])), collapse=' / '),
              if (feasible) 'meets the constraints' else 'breaks the constraints'))
  cat('Incidence (per 100,000):', paste(sprintf('%.0f', calibration$output[[1]]), collapse=' / '),
      '  target', paste(unlist(CALIBRATION.TARGETS[[1]]), collapse=' / '), '\n\n')
  s <- analysis$summary
  print(data.frame(strategy=s$strategy, cost=round(s$C), QALYs=round(s$E, 3), IC=round(s$IC), IE=round(s$IE, 4),
                   `ICER vs none`=round(s$ICER.vs.no_intervention), NMB=round(s$NMB),
                   `screen-detected /1000`=round(s$screen.detected, 1), `overdiagnosed /1000`=round(s$overdiagnosed, 1),
                   `cancer deaths /1000`=round(s$cancer.deaths, 1), check.names=FALSE), row.names=FALSE)
  f <- analysis$frontier
  cat('\nEfficient frontier:', paste(sprintf('%s%s', f$strategy, ifelse(is.na(f$ICER), '', sprintf(' (%s €/QALY)', format(round(f$ICER), big.mark=',')))), collapse=' -> '), '\n')
  cat(sprintf('Optimal at %s €/QALY: %s, by %.0f € of NMB over the next best\n',
              format(WTP, big.mark=','), label(analysis$optimal), analysis$margin))
  cat(sprintf('\nDSA: NMB of %s minus that of the best other strategy, most influential first\n', analysis$optimal))
  print(head(sensitivity[, c('parameter', 'low.margin', 'high.margin', 'optimal.low', 'optimal.high')], 6), row.names=FALSE, digits=4)
  switches <- sensitivity$parameter[sensitivity$decision.changes]
  cat('Parameters whose range changes the decision:', if (length(switches)) paste(switches, collapse=', ') else 'none', '\n')

  overview[[name]] <- data.frame(
    analysis=name,
    calibration.error=signif(calibration$error, 2),
    indolent.percent=paste(sprintf('%.0f', 100 * (1 - x[progressive.slots])), collapse=' / '),
    meets.constraints=feasible,
    optimal.strategy=label(analysis$optimal),
    frontier=paste(f$strategy, collapse=' -> '),
    most.influential=paste(head(sensitivity$parameter, 3), collapse=', '),
    decision.changed.by=if (length(switches)) paste(switches, collapse=', ') else 'nothing')

  slug <- gsub(' ', '_', tolower(name))
  write.csv(s, file.path(OUTPUT.DIR, sprintf('cea_%s.csv', slug)), row.names=FALSE)
  write.csv(sensitivity, file.path(OUTPUT.DIR, sprintf('dsa_%s.csv', slug)), row.names=FALSE)
  attr(scenarios[[name]], 'analysis') <- analysis
  attr(scenarios[[name]], 'dsa') <- sensitivity
}

overview <- do.call(rbind, overview)
rownames(overview) <- NULL
cat('\n', strrep('=', 78), '\nSummary\n', strrep('=', 78), '\n', sep='')
for (i in seq_len(nrow(overview))) {
  o <- overview[i, ]
  cat(sprintf('\n%s\n  error %.2g, indolent lesions %s%% (%s)\n  optimal: %s\n  most influential: %s\n  decision changed by: %s\n',
              o$analysis, o$calibration.error, o$indolent.percent,
              if (o$meets.constraints) 'plausible' else 'implausible',
              o$optimal.strategy, o$most.influential, o$decision.changed.by))
}
write.csv(overview, file.path(OUTPUT.DIR, 'summary.csv'), row.names=FALSE)

calibrated <- data.frame(parameter=rep(scheme$parameters, each=n.groups), age.group=rep(strata, 2),
                         initial.guess=x0, unconstrained=x.unconstrained, constrained=x.constrained)
write.csv(calibrated, file.path(OUTPUT.DIR, 'calibrated_parameters.csv'), row.names=FALSE)

# ---- Plots -------------------------------------------------------------------------

if (requireNamespace('ggplot2', quietly=TRUE)) {
  planes <- do.call(rbind, lapply(names(scenarios), function(name) {
    s <- attr(scenarios[[name]], 'analysis')$summary
    data.frame(analysis=name, strategy=label(s$strategy), IC=s$IC, IE=s$IE)
  }))
  planes$analysis <- factor(planes$analysis, levels=names(scenarios))
  plane <- ggplot(planes, aes(x=IE, y=IC, color=strategy, shape=strategy)) +
    geom_abline(slope=WTP, intercept=0, linetype='dashed', color='grey50') +
    geom_point(size=3) +
    facet_wrap(~analysis, scales='free') +
    labs(x='Incremental QALYs vs no intervention', y='Incremental cost vs no intervention (€)',
         color='', shape='', caption=sprintf('Dashed line: %s €/QALY', format(WTP, big.mark=','))) +
    theme_minimal() + theme(legend.position='bottom')
  ggsave(file.path(OUTPUT.DIR, 'ce_planes.png'), plane, width=12, height=4.5, dpi=120)

  tornados <- do.call(rbind, lapply(names(scenarios), function(name) {
    d <- head(attr(scenarios[[name]], 'dsa'), 8)
    a <- attr(scenarios[[name]], 'analysis')
    data.frame(analysis=sprintf('%s\noptimal: %s', name, label(a$optimal)), parameter=d$parameter,
               low=d$low.margin, high=d$high.margin, base=a$margin, order=rev(seq_len(nrow(d))))
  }))
  # One facet per analysis, each with its own order of the parameters: the
  # bars are keyed by parameter and analysis, and labelled by the parameter.
  tornados$key <- reorder(paste(tornados$parameter, tornados$analysis, sep='@'), tornados$order)
  tornados$analysis <- factor(tornados$analysis, levels=unique(tornados$analysis))
  tornado <- ggplot(tornados) +
    geom_segment(aes(x=base, xend=low, y=key, yend=key, color='Low end of the range'), linewidth=4) +
    geom_segment(aes(x=base, xend=high, y=key, yend=key, color='High end of the range'), linewidth=4) +
    geom_vline(xintercept=0, color='firebrick') +
    scale_y_discrete(labels=function(l) sub('@.*$', '', l)) +
    facet_wrap(~analysis, scales='free') +
    labs(x='NMB of the optimal strategy minus that of the best other one (€); below 0 the decision changes',
         y='', color='') +
    theme_minimal() + theme(legend.position='bottom')
  ggsave(file.path(OUTPUT.DIR, 'tornados.png'), tornado, width=13, height=5, dpi=120)

  ridge.plot <- ggplot(ridge, aes(x=100 * (1 - progressive), y=margin, fill=optimal)) +
    annotate('rect', xmin=-Inf, xmax=100 * INDOLENT.MAX, ymin=-Inf, ymax=Inf, alpha=0.12, fill='forestgreen') +
    geom_col(width=3.5) +
    labs(x='Indolent lesions (%), every parameter set fitting the target exactly',
         y='NMB margin of the optimal strategy (€)', fill='Optimal strategy',
         caption='Shaded: what screening trials allow') +
    theme_minimal() + theme(legend.position='bottom')
  ggsave(file.path(OUTPUT.DIR, 'ridge.png'), ridge.plot, width=8, height=4.5, dpi=120)
}

cat(sprintf('\nResults written to %s/\n', OUTPUT.DIR))
