# cots_shared.R
# Shared helpers for the COTS calibration models (models/cots/*.R): prediction curves
# (NegBin, Binomial), priors, calibration plots, and brms fit-one-pair wrappers that
# save each posterior to output/posteriors/<label>.RDS. Sourced once by the driver
# script before the per-model files.

# Seed for every COTS brm() call.
cots_seed <- 6173

predict_negbin_curve <- function(post, x_seq, log_offset, eps = 0.01, n_draws = 2000) {
  alpha <- post[, "alpha"]; beta <- post[, "beta"]; r <- post[, "r"]
  n <- length(alpha); idx <- if (n_draws < n) sample.int(n, n_draws) else seq_len(n)
  x_lin <- log(x_seq + eps)
  mu <- sapply(x_lin, function(xi) exp(alpha[idx] + beta[idx] * xi + log_offset))
  y_pred <- matrix(rnbinom(length(mu), size = r[idx], mu = mu), nrow = nrow(mu))
  lo95 <- apply(mu, 2, quantile, 0.025); hi95 <- apply(mu, 2, quantile, 0.975)
  pi_lo95 <- apply(y_pred, 2, quantile, 0.025); pi_hi95 <- apply(y_pred, 2, quantile, 0.975)
  # Union the predictive quantiles with the CI so the PI is never narrower than the CI
  # (overdispersed NegBin draws can otherwise give a PI inside the mean-curve CI).
  pi_lo95 <- pmin(pi_lo95, lo95); pi_hi95 <- pmax(pi_hi95, hi95)
  tibble(
    x = x_seq,
    mean    = apply(mu, 2, mean, trim = 0.01),
    lo95    = lo95, hi95 = hi95,
    pi_lo95 = pi_lo95, pi_hi95 = pi_hi95
  )
}

# Prediction interval for a new Binomial observation: simulate Binomial(size, mu) at
# each posterior draw (on the proportion scale, y_pred/size). `size` is the number of
# trials behind one new observation (1 for tow-level presence; a typical trial count
# for cull-site/reef models) and is stored as attr(, "pi_size").
# Reef-level Binomial models carry an observation-level reef effect; its SD column name
# is stored in attr(post, "obs_re_sd"), and a new-reef effect is added to the PI draws.
# mean/lo95/hi95 stay at the average reef (effect = 0).
obs_re_sd_draws <- function(post) {
  col <- attr(post, "obs_re_sd")
  if (is.null(col)) NULL else as.numeric(post[, col])
}

predict_binomial_curve <- function(post, x_seq, offset = NULL, eps = 0.01, link = "logit", n_draws = 2000, size = 1,
                                   obs_sd = obs_re_sd_draws(post)) {
  alpha <- post[, "alpha"]; beta <- post[, "beta"]
  n <- length(alpha); idx <- if (n_draws < n) sample.int(n, n_draws) else seq_len(n)
  x_lin <- log(x_seq + eps)
  linv <- function(lp) if (link == "cloglog") 1 - exp(-exp(lp)) else plogis(lp)
  lp <- sapply(x_lin, function(xi) {
    l <- alpha[idx] + beta[idx] * xi
    if (!is.null(offset)) l <- l + offset
    l
  })
  mu <- linv(lp)
  # One new-reef effect per posterior draw, shared across x
  lp_new <- if (is.null(obs_sd)) lp else lp + rnorm(length(idx), 0, obs_sd[idx])
  y_pred <- matrix(rbinom(length(mu), size = size, prob = linv(lp_new)), nrow = nrow(mu)) / size
  lo95 <- apply(mu, 2, quantile, 0.025); hi95 <- apply(mu, 2, quantile, 0.975)
  pi_lo95 <- apply(y_pred, 2, quantile, 0.025); pi_hi95 <- apply(y_pred, 2, quantile, 0.975)
  # Union with the CI so the PI is never narrower than the CI (as in predict_negbin_curve()).
  pi_lo95 <- pmin(pi_lo95, lo95); pi_hi95 <- pmax(pi_hi95, hi95)
  out <- tibble(
    x = x_seq,
    mean    = apply(mu, 2, mean, trim = 0.01),
    lo95    = lo95, hi95 = hi95,
    pi_lo95 = pi_lo95, pi_hi95 = pi_hi95
  )
  attr(out, "pi_size") <- size
  out
}

# D-02 priors for the grouped eDNA-concentration NegBin models, translated to brms:
# alpha, beta ~ N(0, 10); shape ~ lognormal(log(10), 0.75); reef SD ~ lognormal(log(0.5), 0.5).
d02_priors_eDNAconc_negbin <- c(
  brms::prior(normal(0, 10), class = "Intercept"),
  brms::prior(normal(0, 10), class = "b"),
  brms::prior(lognormal(log(10), 0.75), class = "shape"),
  brms::prior(lognormal(log(0.5), 0.5), class = "sd")
)

# Reef-grain variant without the group-level `sd` prior: one row per reef, so no reef
# random effect is fitted. Also used by the surface-from-eDNA reef models, where the
# informative shape prior stops the NegBin shape collapsing at small n.
d02_priors_eDNAconc_negbin_reef <- c(
  brms::prior(normal(0, 10), class = "Intercept"),
  brms::prior(normal(0, 10), class = "b"),
  brms::prior(lognormal(log(10), 0.75), class = "shape")
)

# D-02 priors for the Binomial/cloglog manta-presence models: alpha, beta ~ N(0, 10);
# one lognormal(log(0.5), 0.5) prior on the group-level SD (D-02 uses it for reef and site).
d02_priors_manta_binary <- c(
  brms::prior(normal(0, 10), class = "Intercept"),
  brms::prior(normal(0, 10), class = "b"),
  brms::prior(lognormal(log(0.5), 0.5), class = "sd")
)

# Weakly regularizing priors for Binomial models without a group term.
d02_priors_binary_reef <- c(
  brms::prior(normal(0, 10), class = "Intercept"),
  brms::prior(normal(0, 10), class = "b")
)

plot_calib_base <- function(pred, obs_x, obs_y, xlab, ylab, main, has_pi = TRUE) {
  cols <- if (has_pi) cbind(pred$mean, pred$lo95, pred$hi95, pred$pi_lo95, pred$pi_hi95)
          else cbind(pred$mean, pred$lo95, pred$hi95)
  lty  <- if (has_pi) c(1, 2, 2, 3, 3) else c(1, 2, 2)
  col  <- if (has_pi) c("black", "blue", "blue", "grey45", "grey45") else c("black", "blue", "blue")
  # Widen the plot limits to cover the observations, which can lie far from a narrow fitted curve.
  ylim <- range(c(cols, obs_y), na.rm = TRUE)
  xlim <- range(c(pred$x, obs_x), na.rm = TRUE)
  matplot(pred$x, cols, type = "l", lty = lty, lwd = c(2, rep(1, length(lty) - 1)), col = col,
          xlab = xlab, ylab = ylab, main = main, ylim = ylim, xlim = xlim)
  points(obs_x, obs_y, pch = 1)
  legend("topleft", bty = "n", cex = 0.7,
         legend = c("Mean", "95% CrI", if (has_pi) "95% PI", "Observations"),
         lty = c(1, 2, if (has_pi) 3, NA), pch = c(NA, NA, if (has_pi) NA, 1),
         col = c("black", "blue", if (has_pi) "grey45", "black"), lwd = 2)
}

# Regularizing priors for the cross-method NegBin models: D-02's N(0, 10) intercept/slope,
# Gamma(2, 1) on the NegBin shape and N(0, 1) on the group-level SD (aids convergence on
# zero-inflated data).
negbin_priors <- function(grouped) {
  p <- list(brms::prior(normal(0, 10), class = "Intercept"),
            brms::prior(normal(0, 10), class = "b"),
            brms::prior(gamma(2, 1), class = "shape"))
  if (grouped) p <- c(p, list(brms::prior(normal(0, 1), class = "sd")))
  do.call(c, p)
}

# Plot one tow-level cross-method direction: x = other method's density, y = own raw count
# (different units, so no 1:1 line). Does not rescale y: tow-level effort can be a per-row
# covariate (e.g. tow distance). Use plot_crossmethod_pair_density() for cull-site/reef level.
plot_crossmethod_pair <- function(post, y_own, effort_own, x_other, effort_other,
                                   xlab, ylab, main, eps = 0.01) {
  x_seq <- seq(0, max(x_other / effort_other), length.out = 100)
  plot_calib_base(predict_negbin_curve(post, x_seq, log_offset = mean(log(effort_own)), eps = eps),
                   x_other / effort_other, y_own, xlab, ylab, main)
}

# NegBin curve on a per-unit-effort scale: predict at a representative effort, then divide
# by it (D-02 convention; predicting at effort = 1 would misstate the dispersion).
predict_negbin_curve_density <- function(post, x_seq, effort_rep, eps = 0.01, n_draws = 2000) {
  predict_negbin_curve(post, x_seq, log_offset = log(effort_rep), eps = eps, n_draws = n_draws) |>
    dplyr::mutate(dplyr::across(-x, ~ .x / effort_rep))
}

# Cull-site/reef counterpart of plot_crossmethod_pair(): effort_own is the number of tows
# aggregated into a row, so both the curve and the observations are shown per tow.
plot_crossmethod_pair_density <- function(post, y_own, effort_own, x_other, effort_other,
                                           xlab, ylab, main, eps = 0.01) {
  x_seq <- seq(0, max(x_other / effort_other), length.out = 100)
  effort_rep <- mean(effort_own)
  plot_calib_base(predict_negbin_curve_density(post, x_seq, effort_rep, eps = eps),
                   x_other / effort_other, y_own / effort_own, xlab, ylab, main)
}

# Fit one NegBin direction: y ~ log(other density + eps) + offset(log(effort_own)), log link,
# optional reef random effect. Saves the posterior (alpha, beta, r = shape) and drops the
# brmsfit. `prior_override` replaces the default negbin_priors(grouped).
fit_crossmethod_pair <- function(y, effort_own, x_other, effort_other, group_id, eps, label, level, prior_override = NULL) {
  # x_dens_eps is its own column because brms needs formula variables in `data`.
  dat <- tibble(y = y, x_dens_eps = x_other / effort_other + eps, .offset = log(effort_own))
  grouped <- !is.null(group_id)
  if (grouped) dat$.group <- factor(group_id)
  form <- if (grouped) bf(y ~ log(x_dens_eps) + offset(.offset) + (1 | .group))
          else          bf(y ~ log(x_dens_eps) + offset(.offset))
  cat("\nFitting:", label, " | level =", level, " | n =", nrow(dat), "\n")
  use_prior <- if (!is.null(prior_override)) prior_override else negbin_priors(grouped)
  fit <- brm(form, family = negbinomial(link = "log"), data = dat, prior = use_prior,
             chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12))
  post <- as.matrix(as_draws_matrix(fit))
  colnames(post)[colnames(post) == "b_Intercept"] <- "alpha"
  colnames(post)[startsWith(colnames(post), "b_") & colnames(post) != "b_Intercept"]  <- "beta"
  colnames(post)[colnames(post) == "shape"]        <- "r"
  saveRDS(post, file.path("output/posteriors", paste0(label, ".RDS")))
  rm(fit)  # free memory
  post
}

# Fit the D-02 tow-level Binomial/cloglog presence model:
# has_scar | trials(1) ~ log(density + eps) + offset(log(tow_dist)) + (1 | group).
# Priors: d02_priors_manta_binary. Saves the posterior and drops the brmsfit.
fit_binary_crossmethod_pair <- function(has_scar, tow_dist, x_dens, eps, group_id, label) {
  # Centre the offset: log(tow_dist) is far from 0, which makes Stan's default inits fail.
  offset_center <- mean(log(tow_dist))
  dat <- tibble(has_scar = has_scar, .offset = log(tow_dist) - offset_center,
                x_dens_eps = x_dens + eps, .group = factor(group_id))
  fit <- brm(bf(has_scar | trials(1) ~ log(x_dens_eps) + offset(.offset) + (1 | .group)),
             family = binomial(link = "cloglog"), data = dat,
             prior = d02_priors_manta_binary,
             chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.97, max_treedepth = 13))
  post <- as.matrix(as_draws_matrix(fit))
  colnames(post)[colnames(post) == "b_Intercept"] <- "alpha"
  colnames(post)[startsWith(colnames(post), "b_") & colnames(post) != "b_Intercept"] <- "beta"
  post[, "alpha"] <- post[, "alpha"] - offset_center
  attr(post, "link") <- "cloglog"
  saveRDS(post, file.path("output/posteriors", paste0(label, ".RDS")))
  rm(fit)
  post
}

# Cull-site/reef counterpart: n_success | trials(n_trials) ~ log(density + eps), logit link
# (trials() carries the exposure, so no offset). `group_id` is an optional reef random
# effect; if every group has one row it is an observation-level effect, flagged via
# attr(post, "obs_re_sd") so predict_binomial_curve() adds it to the PI.
fit_binary_crossmethod_cull_reef_pair <- function(n_success, n_trials, x_dens, eps, group_id, label, level) {
  dat <- tibble(n_success = n_success, n_trials = n_trials, x_dens_eps = x_dens + eps)
  grouped <- !is.null(group_id)
  if (grouped) dat$.group <- factor(group_id)
  obs_level_re <- grouped && !anyDuplicated(group_id)
  form <- if (grouped) bf(n_success | trials(n_trials) ~ log(x_dens_eps) + (1 | .group))
          else          bf(n_success | trials(n_trials) ~ log(x_dens_eps))
  cat("\nFitting:", label, " | level =", level, " | n =", nrow(dat), "\n")
  fit <- brm(form, family = binomial(link = "logit"), data = dat,
             prior = if (grouped) d02_priors_manta_binary else d02_priors_binary_reef,
             chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12))
  post <- as.matrix(as_draws_matrix(fit))
  colnames(post)[colnames(post) == "b_Intercept"] <- "alpha"
  colnames(post)[startsWith(colnames(post), "b_") & colnames(post) != "b_Intercept"] <- "beta"
  attr(post, "link") <- "logit"
  if (obs_level_re) attr(post, "obs_re_sd") <- "sd_.group__Intercept"
  saveRDS(post, file.path("output/posteriors", paste0(label, ".RDS")))
  rm(fit)  # free memory
  post
}
