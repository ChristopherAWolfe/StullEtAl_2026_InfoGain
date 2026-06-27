# =============================================================================
# kl_prior_posterior_overview.R
# =============================================================================
#
# This script illustrates information gain in Bayesian inference. It compares
# two probability densities on the interval [0, 2*mu]: a uniform prior Q(x) and
# a truncated Gaussian posterior P(x) with mean mu.
#
# The script computes the Kullback-Leibler divergence D_KL(P || Q) in bits:
#
#   D_KL(P || Q) = integral_0^(2*mu) P(x) log2(P(x) / Q(x)) dx
#
# This divergence measures how much the posterior narrows relative to the prior.
# A smaller value of sigma makes P more concentrated and increases D_KL(P || Q).
# The script evaluates KL divergence from sampled densities, inverts that
# relationship to find sigma for a target bit count, plots the resulting
# density curves, and runs simple checks that stop execution if anything fails.
#
# Run the script with Rscript kl_prior_posterior_overview.R or source it from
# an interactive R session. The symbol mu denotes the center of the posterior
# and half the width of the support interval. The symbol sigma controls the
# spread of P before truncation. The symbol dx is the spacing between
# consecutive grid points from 0 to 2*mu.
# =============================================================================

# --- parameters ---

# mu is the center of the posterior and the upper support bound is 2*mu.
mu <- 10

# dx is the spacing between grid points used for trapezoidal integration. The
# ratio 2*mu/dx must be an integer so that the grid runs from 0 to 2*mu in
# steps of dx, inclusive of both endpoints.
dx <- 0.001

# target_bits lists the information levels, in bits, used for inversion and
# plotting.
target_bits <- c(0.5, 1, 2)

# plot_file is the path of the PDF written by the plotting code.
plot_file <- "kl_density_curves.pdf"

# --- helpers ---

# assert stops immediately with a readable message when condition is FALSE.
assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
  return(invisible(NULL))
}

# integrate_on_grid approximates the integral of a function sampled on the grid
# x using the trapezoidal rule. The arguments values and x must have the same
# length, and the grid spacing is inferred from x and mu.
integrate_on_grid <- function(values, x, mu) {
  if (length(values) != length(x)) {
    stop("values and x must have the same length.")
  }
  if (length(values) < 2) {
    stop("At least two grid points are required for trapezoidal integration.")
  }

  dx <- compute_grid_spacing(x, mu)
  n <- length(values)
  return(dx * (0.5 * values[1] + sum(values[2:(n - 1)]) + 0.5 * values[n]))
}

# compute_kl_divergence_discrete approximates D_KL(P || Q) in bits from density
# samples p and q on the grid x. Both densities must integrate to 1 within
# tol, and Q must be positive wherever P is positive.
compute_kl_divergence_discrete <- function(p, q, x, mu, tol = 1e-6) {
  if (length(p) != length(q)) {
    stop("p and q must have the same length.")
  }
  if (length(p) != length(x)) {
    stop("p and x must have the same length.")
  }
  if (any(p < 0, na.rm = TRUE) || any(q < 0, na.rm = TRUE)) {
    stop("Density samples must be non-negative.")
  }

  p_mass <- integrate_on_grid(p, x, mu)
  q_mass <- integrate_on_grid(q, x, mu)
  if (abs(p_mass - 1) > tol) {
    stop(sprintf("P does not integrate to 1 (integral = %.10f).", p_mass))
  }
  if (abs(q_mass - 1) > tol) {
    stop(sprintf("Q does not integrate to 1 (integral = %.10f).", q_mass))
  }

  # By the usual convention for KL divergence, terms where P is zero contribute
  # nothing, so those grid points are omitted from the sum.
  support <- p > 0
  if (!any(support)) {
    stop("P has no positive mass on the grid.")
  }
  if (any(q[support] <= 0, na.rm = TRUE)) {
    stop("Q must be positive wherever P is positive.")
  }

  integrand <- rep(0, length(p))
  integrand[support] <- p[support] * log2(p[support] / q[support])
  return(integrate_on_grid(integrand, x, mu))
}

# make_support_grid returns the sample points 0, dx, 2*dx, ..., 2*mu, including
# both endpoints. The result has n_bins + 1 points, where n_bins = 2*mu/dx.
make_support_grid <- function(mu, dx) {
  if (mu <= 0) {
    stop("mu must be positive.")
  }
  if (dx <= 0) {
    stop("dx must be positive.")
  }

  n_bins <- 2 * mu / dx
  if (abs(n_bins - round(n_bins)) > 1e-10) {
    stop("2 * mu / dx must be an integer for the requested grid spacing.")
  }

  x <- seq(0, 2 * mu, by = dx)
  if (abs(tail(x, 1) - 2 * mu) > 1e-10) {
    stop("Grid construction failed: last point is not 2 * mu.")
  }
  return(x)
}

# validate_support_grid checks that x is an evenly spaced grid from 0 to 2*mu.
# It raises an error if x[1] is not 0, if the final value is not 2*mu, or if
# the spacing between consecutive points is not uniform.
validate_support_grid <- function(x, mu, tol = 1e-8) {
  if (mu <= 0) {
    stop("mu must be positive.")
  }
  if (length(x) < 2) {
    stop("x must contain at least two grid points.")
  }
  if (any(is.na(x)) || is.na(mu)) {
    stop("x and mu must be finite.")
  }
  if (abs(x[1]) > tol) {
    stop("x[1] must be 0.")
  }
  if (abs(x[length(x)] - 2 * mu) > tol) {
    stop(sprintf("Final x value must be 2 * mu (expected %.10f, got %.10f).", 2 * mu, x[length(x)]))
  }

  steps <- diff(x)
  if (any(abs(steps - steps[1]) > tol)) {
    stop("x must be evenly spaced.")
  }

  expected_spacing <- (2 * mu) / (length(x) - 1)
  if (abs(steps[1] - expected_spacing) > tol) {
    stop(sprintf(
      "Grid spacing must equal 2 * mu / (length(x) - 1) (expected %.10f, got %.10f).",
      expected_spacing,
      steps[1]
    ))
  }

  return(invisible(NULL))
}

# compute_grid_spacing validates x and returns the spacing implied by the grid.
compute_grid_spacing <- function(x, mu, tol = 1e-8) {
  validate_support_grid(x, mu, tol)
  return((2 * mu) / (length(x) - 1))
}

# make_uniform_density evaluates the uniform prior Q(x) = 1 / (2*mu) at each
# point of the validated grid x.
make_uniform_density <- function(x, mu) {
  validate_support_grid(x, mu)
  return(rep(1 / (2 * mu), length(x)))
}

# make_truncated_gaussian_density evaluates the posterior density: a Gaussian
# N(mu, sigma^2) renormalized to integrate to 1 on [0, 2*mu]. Because that
# interval is symmetric about mu, the mean of the truncated distribution is
# also mu.
make_truncated_gaussian_density <- function(x, mu, sigma) {
  validate_support_grid(x, mu)
  if (sigma <= 0) {
    stop("sigma must be positive.")
  }

  lower <- 0
  upper <- 2 * mu
  norm_const <- pnorm((upper - mu) / sigma) - pnorm((lower - mu) / sigma)
  if (norm_const <= 0) {
    stop("Truncation interval has zero probability mass for the given parameters.")
  }

  return(dnorm((x - mu) / sigma) / (sigma * norm_const))
}

# compute_kl_bits_truncated_gaussian computes D_KL(P || Q) in bits for the
# truncated-Gaussian posterior P and uniform prior Q defined by mu, sigma, and dx.
compute_kl_bits_truncated_gaussian <- function(mu, sigma, dx, tol = 1e-6) {
  x <- make_support_grid(mu, dx)
  p <- make_truncated_gaussian_density(x, mu, sigma)
  q <- make_uniform_density(x, mu)
  return(compute_kl_divergence_discrete(p, q, x, mu, tol = tol))
}

# solve_sigma_for_bits finds the value of sigma such that D_KL(P || Q) equals
# target_bits. The divergence decreases as sigma increases, so there is at most
# one solution. The function uses uniroot after bracketing sigma automatically.
# A target of 0 bits corresponds to P equal to Q and has no finite sigma solution,
# so target_bits must be strictly positive.
solve_sigma_for_bits <- function(mu, target_bits, dx, tol = 1e-6, sigma_bounds = NULL) {
  if (mu <= 0) {
    stop("mu must be positive.")
  }
  if (target_bits <= 0) {
    stop("target_bits must be strictly positive.")
  }
  if (dx <= 0) {
    stop("dx must be positive.")
  }

  objective <- function(sigma) {
    return(compute_kl_bits_truncated_gaussian(mu, sigma, dx, tol = tol) - target_bits)
  }

  if (is.null(sigma_bounds)) {
    # A very large sigma makes P nearly flat, which gives a KL divergence below
    # the target, so sigma_hi is increased until the objective is non-positive.
    sigma_hi <- mu
    while (objective(sigma_hi) > 0 && sigma_hi < 1e6 * mu) {
      sigma_hi <- sigma_hi * 2
    }
    if (objective(sigma_hi) > 0) {
      stop("Could not bracket sigma: target_bits may be below the minimum achievable divergence.")
    }

    # A very small sigma makes P sharply peaked, which gives a KL divergence
    # above the target, so sigma_lo is increased until the objective is non-negative.
    # sigma_lo starts at dx because narrower peaks are not resolved reliably below that.
    sigma_lo <- dx
    while (objective(sigma_lo) < 0 && sigma_lo < sigma_hi / 2) {
      sigma_lo <- sigma_lo * 2
    }
    if (objective(sigma_lo) < 0) {
      stop(
        "Could not bracket sigma: target_bits may exceed the maximum achievable divergence ",
        "for the chosen dx."
      )
    }

    sigma_bounds <- c(sigma_lo, sigma_hi)
  }

  f_lo <- objective(sigma_bounds[1])
  f_hi <- objective(sigma_bounds[2])
  if (f_lo * f_hi > 0) {
    stop("Supplied sigma_bounds do not bracket the target_bits solution.")
  }

  root <- uniroot(objective, interval = sigma_bounds, tol = .Machine$double.eps^0.5)$root
  return(root)
}

# make_p_density_for_bits builds the posterior P that achieves the requested KL
# divergence from the prior in bits. It returns the grid x, the density p, the
# fitted sigma, and the achieved bit count.
make_p_density_for_bits <- function(mu, target_bits, dx) {
  x <- make_support_grid(mu, dx)
  sigma <- solve_sigma_for_bits(mu, target_bits, dx)
  return(list(
    x = x,
    p = make_truncated_gaussian_density(x, mu, sigma),
    sigma = sigma,
    bits = compute_kl_bits_truncated_gaussian(mu, sigma, dx)
  ))
}

# --- plots ---

# plot_kl_density_curves draws the uniform prior Q and one truncated Gaussian
# posterior P for each target bit level. The prior uses gray40; each posterior
# uses a distinct Dark 3 color. Bit labels appear on the curves, and the legend
# lists the prior plus one Posterior P entry per curve color. If output_file is NULL,
# the plot is drawn on the active graphics device; otherwise it is written to a PDF.
plot_kl_density_curves <- function(mu, dx, target_bits, output_file = NULL) {
  q_level <- 1 / (2 * mu)
  curves <- lapply(target_bits, function(bits) make_p_density_for_bits(mu, bits, dx))
  y_max <- max(c(q_level, unlist(lapply(curves, function(curve) curve$p))))
  prior_lwd <- 2.5
  posterior_lwd <- 3
  prior_col <- "gray40"
  posterior_cols <- hcl.colors(length(curves), palette = "Dark 3")
  bit_labels <- paste0(target_bits, " bit", ifelse(target_bits == 1, "", "s"))

  y_upper <- y_max * 1.05
  y_ticks <- pretty(c(0, y_upper), n = 4)

  draw <- function() {
    plot(
      NA,
      xlim = c(0, 2 * mu),
      ylim = c(0, y_upper),
      xlab = expression(x),
      ylab = "Density",
      main = "KL divergence of posterior P from prior Q",
      yaxt = "n"
    )

    axis(2, at = y_ticks, las = 1)

    abline(h = q_level, lty = 1, col = prior_col, lwd = prior_lwd)

    for (i in seq_along(curves)) {
      curve <- curves[[i]]
      curve_col <- posterior_cols[i]
      lines(curve$x, curve$p, lty = 1, lwd = posterior_lwd, col = curve_col)

      peak_idx <- which.min(abs(curve$x - mu))
      label_gap <- 0.015 * y_max
      text(
        mu,
        curve$p[peak_idx] + label_gap,
        labels = bit_labels[i],
        adj = c(0.5, 0),
        cex = 0.9,
        col = curve_col
      )
    }

    legend(
      "topright",
      legend = c("Prior Q", rep("Posterior P", length(curves))),
      col = c(prior_col, posterior_cols),
      lty = 1,
      lwd = c(prior_lwd, rep(posterior_lwd, length(curves))),
      bty = "n"
    )

    return(invisible(NULL))
  }

  if (is.null(output_file)) {
    draw()
  } else {
    pdf(output_file, width = 8, height = 5)
    on.exit(dev.off(), add = TRUE)
    draw()
  }

  return(invisible(curves))
}

# --- tests ---
#
# These checks run every time the script is sourced. Any failure stops execution.

x <- make_support_grid(mu, dx)
q <- make_uniform_density(x, mu)
assert(abs(integrate_on_grid(q, x, mu) - 1) < 1e-6, "prior Q should integrate to 1")

bad_x <- x
bad_x[length(bad_x)] <- 2 * mu + compute_grid_spacing(x, mu)
err_grid <- tryCatch(
  make_uniform_density(bad_x, mu),
  error = function(e) conditionMessage(e)
)
assert(grepl("Final x value must be 2 * mu", err_grid, fixed = TRUE), "make_uniform_density should reject invalid grids")

for (bits in target_bits) {
  sigma <- solve_sigma_for_bits(mu, bits, dx)
  achieved <- compute_kl_bits_truncated_gaussian(mu, sigma, dx)
  assert(
    abs(achieved - bits) < 1e-4,
    sprintf(
      "solve_sigma_for_bits failed for %.1f bits (achieved %.6f, sigma %.6f)",
      bits, achieved, sigma
    )
  )
}

# The next check deliberately breaks normalization and confirms that
# compute_kl_divergence_discrete rejects the invalid density.
p_bad <- make_truncated_gaussian_density(x, mu, 2) * 1.01
err <- tryCatch(
  compute_kl_divergence_discrete(p_bad, q, x, mu),
  error = function(e) conditionMessage(e)
)
assert(
  grepl("does not integrate to 1", err),
  "compute_kl_divergence_discrete should reject non-normalized P"
)

# --- run ---

curves <- plot_kl_density_curves(mu, dx, target_bits, output_file = plot_file)

cat("\nPosterior standard deviations (sigma) for each target KL divergence:\n")
cat(sprintf("  mu = %.4f\n", mu))
for (i in seq_along(curves)) {
  cat(sprintf(
    "  %.1f bits: sigma = %.6f (achieved %.6f bits)\n",
    target_bits[i],
    curves[[i]]$sigma,
    curves[[i]]$bits
  ))
}

cat("\nAll checks passed. Wrote plot to", plot_file, "\n")
