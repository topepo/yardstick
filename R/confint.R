
###################################################################
## Bootstrap sampling of tables

# Creates a series of bootstrap sampled tables
boot_tab <- function(x, B = 10, zero_offset = 0) {
  n <- sum(x)
  probs <- as.vector(x) / n
  if (zero_offset > 0) {
    offset <- min(probs[probs > 0]) * zero_offset
    probs[probs == 0] <- offset
    probs <- probs / sum(probs)
  }
  resample_table(probs, n, dimnames(x), times = B)
}

#' @importFrom purrr map 
resample_table <- function(probs, n, dnames, times = 1) {
  dims <- unlist(lapply(dnames, length))
  p <- prod(dims)
  samps <-
    purrr::map(
      rep(p, times),
      sample.int,
      prob = probs,
      size = n,
      replace = TRUE
    )
  res <- purrr::map(samps, samp2table, dims = dims)
  res <- purrr::map(
    res,
    function(x, dn) {
      dimnames(x) <- dn
      x
    },
    dn = dnames
  )
  res
}

# convert sampled integers to tables
samp2table <- function(x, dims) {
  res <- rep(0, prod(dims))
  xtab <- table(x)
  ind <- as.integer(names(xtab))
  res[ind] <- xtab
  res <- matrix(res, ncol = dims[2], nrow = dims[1])
  as.table(res)
}

#' importFrom stats quantile
get_pctl <- function(x, alt, level) {
  if (alt == "less") {
    probs <- level
  } else {
    if (alt == "greater") {
      probs <- 1 - level
    } else {
      conf <- (1 - level) / 2
      probs = c(conf, 1 - conf)
    }
  }
  res <- unname(quantile(x, probs = probs))
  if (alt == "less")
    res <- c(NA, res)
  if (alt == "greater")
    res <- c(res, NA)
  names(res) <- c("lower", "upper")
  res
}

# use with_seed?
# should be able to pass a set of parameters too

#' @importFrom purrr map_dbl
confint_boot_table <- function(
  xtab,
  statistic,
  level = 0.95,
  alternative = c("two.sided", "less", "greater"),
  B = 5000,
  eps = 0
) {
  sim_tabs <- boot_tab(xtab, eps = eps, B = B)
  sim_stats <- map_dbl(sim_tabs, statistic)
  get_pctl(sim_stats, alt = match.arg(alternative), level = level)
}

###################################################################
## Exact assymetric intervals for any binomial proportions

#' @importFrom stats binom.test
confint_binom <- function(
  numer,
  denom,
  level = 0.95,
  alternative = c("two.sided", "less", "greater")
  ) {
  alternative <- match.arg(alternative)
  res <- binom.test(
    x = numer,
    n = denom,
    alternative = alternative,
    conf.level = level
  )
  res <- as.vector(res$conf.int)
  if (alternative == "less")
    res[1] <- NA
  if (alternative == "greater")
    res[2] <- NA
  names(res) <- c("lower", "upper")
  res
}

###################################################################
## Boostrap intervals for non-table data structures (rmse, pr curve)

# not working yet; need to figure out ... with NSE (again)

#' @importFrom rsample bootstraps
confint_boot_df <- function(
  df,
  statistic,
  level = 0.95,
  alternative = c("two.sided", "less", "greater"),
  B = 50,
  ...
) {
  
  samp_df <- rsample::bootstraps(df, times = B)
  sim_stats <- map_dbl(samp_df, stat_wrap, stat = statistic, ...)
  get_pctl(sim_stats, alt = match.arg(alternative), level = level)
}

#' @importFrom rsample analysis
stat_wrap <- function(dat, stat, ...) {
  stat(analysis(dat), ...)
}
