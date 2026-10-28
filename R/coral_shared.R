# coral_shared.R
# Shared helpers for the coral-cover Beta-regression models (models/coral/*.R):
# prediction curve, base-graphics calibration plot, dispersion prior, and a JAGS
# fitting wrapper. Sourced once by the driver script before the per-model files.

# Beta likelihood needs a response strictly inside (0, 1); cover can be exactly 0.
clip01 <- function(y, eps = 1e-5) pmin(pmax(y, eps), 1 - eps)

# Seed for every coral fit (brms `seed =`, JAGS per-chain .RNG.seed).
coral_seed <- 2958

# Population-level calibration curve (mean, 95% CrI, 95% PI) for Beta models.
# Link is read from attr(post, "link"), default "logit"; loglog is also supported.
# alpha/beta are stored on the direct link scale, so only the inverse link is needed.
predict_beta_curve <- function(post, x_seq, n_draws = 2000) {
  # post: matrix with columns alpha, beta, phi; reef effect ignored (population curve).
  alpha <- post[, "alpha"]; beta <- post[, "beta"]; phi <- post[, "phi"]
  link <- attr(post, "link"); if (is.null(link)) link <- "logit"
  inv_link <- switch(link,
    logit  = plogis,
    loglog = function(eta) exp(-exp(-eta)),
    stop("predict_beta_curve: unsupported link '", link, "'")
  )
  n <- length(alpha); idx <- if (n_draws < n) sample.int(n, n_draws) else seq_len(n)
  mu <- sapply(x_seq, function(xi) inv_link(alpha[idx] + beta[idx] * xi))
  mu_adj <- pmin(pmax(mu, 1e-5), 1 - 1e-5)
  y_pred <- matrix(rbeta(length(mu_adj), mu_adj * phi[idx], (1 - mu_adj) * phi[idx]), nrow = nrow(mu_adj))
  tibble(
    x = x_seq,
    mean    = apply(mu, 2, mean, trim = 0.01),
    lo95    = apply(mu, 2, quantile, 0.025), hi95 = apply(mu, 2, quantile, 0.975),
    pi_lo95 = apply(y_pred, 2, quantile, 0.025), pi_hi95 = apply(y_pred, 2, quantile, 0.975)
  )
}

plot_calib_base <- function(pred, obs_x, obs_y, xlab, ylab, main) {
  matplot(pred$x * 100, cbind(pred$mean, pred$lo95, pred$hi95, pred$pi_lo95, pred$pi_hi95) * 100,
          type = "l", lty = c(1, 2, 2, 3, 3), lwd = c(2, 1, 1, 1, 1),
          col = c("black", "blue", "blue", "grey45", "grey45"),
          xlab = xlab, ylab = ylab, main = main)
  points(obs_x * 100, obs_y * 100, pch = 1)
  legend("topleft", bty = "n", cex = 0.7,
         legend = c("Mean", "95% CrI", "95% PI", "Observations"),
         lty = c(1, 2, 3, NA), pch = c(NA, NA, NA, 1),
         col = c("black", "blue", "grey45", "black"), lwd = 2)
}

# D-02 prior for the Beta dispersion phi. Set explicitly in brms fits, whose default
# phi prior (gamma(0.01, 0.01)) is near-improper; JAGS models below use the same prior.
d02_prior_coral_phi <- brms::prior(gamma(2, 0.01), class = "phi")

# JAGS Beta-regression model (logit link), shared by the JAGS pair fits.
jags_beta_model_string <- "model{
  for (i in 1:N) {
    logit(p[i]) <- alpha + beta * x[i]
    y[i]    ~ dbeta(a1[i], a2[i])
    a1[i]   <- p_ad[i] * phi
    a2[i]   <- (1 - p_ad[i]) * phi
    p_ad[i] <- max(1e-5, min(1 - 1e-5, p[i]))
  }
  alpha ~ dnorm(0, 0.01)
  beta  ~ dnorm(0, 0.01)
  phi   ~ dgamma(2, 0.01)
}"

# loglog-link variant: inverse link written directly, mu = exp(-exp(-eta)).
jags_beta_model_string_loglog <- "model{
  for (i in 1:N) {
    p[i]    <- exp(-exp(-(alpha + beta * x[i])))
    y[i]    ~ dbeta(a1[i], a2[i])
    a1[i]   <- p_ad[i] * phi
    a2[i]   <- (1 - p_ad[i]) * phi
    p_ad[i] <- max(1e-5, min(1 - 1e-5, p[i]))
  }
  alpha ~ dnorm(0, 0.01)
  beta  ~ dnorm(0, 0.01)
  phi   ~ dgamma(2, 0.01)
}"

# Fit one direction with JAGS, save the posterior to output/posteriors/<label>.RDS and
# return it. Fit objects are removed on exit. link = "loglog" tags the posterior
# with attr(post, "link") so predict_beta_curve() uses the right inverse link.
fit_jags_beta_pair <- function(y, x, label, link = c("logit", "loglog")) {
  link <- match.arg(link)
  model_string <- if (link == "logit") jags_beta_model_string else jags_beta_model_string_loglog
  jagsDat <- list(y = pmin(pmax(y, 1e-5), 1 - 1e-5), x = x, N = length(y))
  inits <- lapply(1:5, function(i) list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = coral_seed + i))
  fm <- jags.model(textConnection(model_string), data = jagsDat, inits = inits, n.chains = 5, quiet = TRUE)
  update(fm, n.iter = 10000, progress.bar = "none")
  fm.post <- coda.samples(fm, variable.names = c("alpha", "beta", "phi"), n.iter = 200000, thin = 200, progress.bar = "none")
  post <- as.matrix(fm.post)
  if (link == "loglog") attr(post, "link") <- "loglog"
  saveRDS(post, file.path("output/posteriors", paste0(label, ".RDS")))
  rm(fm, fm.post)
  post
}
