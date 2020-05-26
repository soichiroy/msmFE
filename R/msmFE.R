#' Main function  for estimating treatment effects
#' @export
msmFE <- function(
  formula_outcome = NULL,
  formula_treatment,
  data, index,
  estimator = "HT",
  fixed_effects = TRUE,
  options = list()
) {


  ## data format
  ## remove no-variation observations if no-variation in the treatment
  pscore_fit <- estimate_pscore()
}



#' Function to estimate propensity scores
#' @export
estimate_pscore <- function(formula, data) {

  ## estimate
  fit <- bife::bife(formula, data = data, model = "logit")

  ## propensity score
  fitted <- predict(fit, type = "response")
  return(list(fit = fit, ps = fitted))
}


estimate_outcome <- function() {

  ## some local change from groupFE repo 
}
