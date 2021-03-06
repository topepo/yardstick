#' Symmetric mean absolute percentage error
#'
#' Calculate the symmetric mean absolute percentage error. This metric is in
#' _relative units_.
#'
#' This implementation of `smape()` is the "usual definition" where the
#' denominator is divided by two.
#'
#' @family numeric metrics
#' @family accuracy metrics
#' @templateVar metric_fn smape
#' @template return
#'
#' @inheritParams rmse
#'
#' @author Max Kuhn, Riaz Hedayati
#'
#' @template examples-numeric
#'
#' @export
#'
smape <- function(data, ...) {
  UseMethod("smape")
}
smape <- new_numeric_metric(
  smape,
  direction = "minimize"
)

#' @rdname smape
#' @export
smape.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {

  metric_summarizer(
    metric_nm = "smape",
    metric_fn = smape_vec,
    data = data,
    truth = !!enquo(truth),
    estimate = !!enquo(estimate),
    na_rm = na_rm
  )

}

#' @export
#' @rdname smape
smape_vec <- function(truth, estimate, na_rm = TRUE, ...) {

  smape_impl <- function(truth, estimate) {

    percent_scale <- 100
    numer <- abs(estimate - truth)
    denom <- (abs(truth) + abs(estimate)) / 2
    mean(numer / denom) * percent_scale

  }

  metric_vec_template(
    metric_impl = smape_impl,
    truth = truth,
    estimate = estimate,
    na_rm = na_rm,
    cls = "numeric"
  )

}
