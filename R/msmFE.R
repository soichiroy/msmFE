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

  id_time_vec <- data %>% pull(index[2]) %>% unique() %>% sort()

  ## data format
  ## remove no-variation observations if no-variation in the treatment
  pscore_fit <- estimate_pscore(formula_treatment, data)

  if (!is.null(formula_outcome)) {

  }

  ## estimate treatment effect


}



#' Function to estimate propensity scores
#' @import tidyverse
#' @export
estimate_pscore <- function(formula, data, id_time_vec = NULL) {

  if (is.null(id_time_vec) & ("panel_data" %in% class(data))) {
    id_time_name <- panelr::get_wave(data)
    id_time_vec  <- sort(unique(data %>% pull(id_time_name)))
  }

  ## estimate
  fit <- bife::bife(formula, data = data, model = "logit")

  ## propensity score
  fitted <- predict(fit, type = "response")

  ## format data
  pscore_dat <- enframe(fitted) %>%
    tibble::as_tibble() %>%
    rename(id_unit = name) %>%
    mutate(id_unit = as.numeric(id_unit)) %>% ## this step need to be fixed
    group_by(id_unit) %>%
    mutate(id_time = id_time_vec) %>%
    ungroup()

  return(list(fit = fit, ps = fitted, dat = pscore_dat))
}


estimate_outcome <- function() {

  ## some local change from groupFE repo
}

#' Weighting estimator
#' @import tidyverse
#' @export
estimate_effect <- function(Y, D, ps, trim = TRUE, na_omit = TRUE) {
  dat_complete <- tibble::tibble(Yc = Y, Dc = D, psc = ps)

  if (isTRUE(na_omit)) {
    dat_complete <- na.omit(dat_complete)
  } else {
    dat_complete <- dat_complete %>%
      mutate(psc = if_else(is.na(psc),
                            if_else(Dc == 1, 0.9999, 0.0001),
                            psc
                          ))
  }

  tau_ht <- with(dat_complete, estimator_ht(Yc, Dc, psc, trim = TRUE))
  tau_hj <- with(dat_complete, estimator_hajek(Yc, Dc, psc, trim = TRUE))
  tau_vec <- c(tau_ht, tau_hj)
  names(tau_vec) <- c("HT", "Hajek")

  return(tau_vec)

}


#' HT estimator
#' @export
estimator_ht <- function(Y, D, ps, trim = TRUE) {
  if(isTRUE(trim)) {
    ps <- pmin(ps, 0.9999)
    ps <- pmax(ps, 0.0001)
  }
  tmp1 <- (D * Y) / ps - ((1 - D) * Y / (1 - ps))
  return(sum(tmp1) / length(D))
}


#' HT estimator
#' @export
estimator_hajek <- function(Y, D, ps, trim = TRUE) {
  if(isTRUE(trim)) {
    ps <- pmin(ps, 0.9999)
    ps <- pmax(ps, 0.0001)
  }
  tmp_up1   <- sum((D * Y) / ps)
  tmp_down1 <- sum(D / ps)
  tmp_up0   <- sum(((1-D) * Y) /(1 - ps))
  tmp_down0 <- sum((1 - D) / (1 - ps))

  return(tmp_up1/tmp_down1 - tmp_up0/tmp_down0)
}
