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
#' @import Formula
#' @export
estimate_pscore <- function(
    formula, data, id_time_vec = NULL, bias_correct = FALSE,
is_FE = TRUE, impute_FE = FALSE, model = "logit"
) {

  if (is.null(id_time_vec) & ("panel_data" %in% class(data))) {
    id_time_name <- panelr::get_wave(data)
    id_time_vec  <- sort(unique(data %>% pull(id_time_name)))
  }

  lin_pred <- NULL

  ## estimate
  if (isTRUE(is_FE)) {
    fit <- bife::bife(formula, data = data, model = model)
    ## bias correction
    if (isTRUE(bias_correct)) {
      fit <- bife::bias_corr(fit)
    }

    ## predicted prob
    if (isTRUE(impute_FE)) {
      ## fixed effect estimates (observations are dropped)
      FEest        <- fit$alpha

      ## compute the linear predictor X'β
      coefs        <- as.vector(fit$coef)
      formula_main <- formula(Formula(formula), rhs = 1, lhs = 0)
      formula_main <- update(formula_main, ~ . - 1)
      Xmat         <- model.matrix(formula_main, data = data)
      Xb           <- as.vector(Xmat %*% coefs)

      ## obtain treatment (i.e., outcome variable)
      treatment <- all.vars(formula)[1]
      lin_pred  <- tibble::tibble(
        treatment = pull(data, !!sym(treatment)),
        lin_pred = Xb,
        id_unit  = pull(data, !!sym(panelr::get_id(data))),
        id_time  = pull(data, !!sym(panelr::get_wave(data)))
      ) %>%
      left_join(
        enframe(FEest, name = "id_unit", value = 'FEest'), by = 'id_unit'
      ) %>%
      mutate(id_unit = as.numeric(id_unit))
    }
    ## ps
    fitted <- predict(fit, type = "response")
  } else {
    fit <- glm(formula, data = data, family = binomial(link = model))
    fitted <- predict(fit, type = "response")
    names(fitted) <- pull(data, panelr::get_id(data))
  }

  ## propensity score


  ## format data
  pscore_dat <- enframe(fitted) %>%
    tibble::as_tibble() %>%
    rename(id_unit = name) %>%
    mutate(id_unit = as.numeric(id_unit)) %>% ## this step need to be fixed
    group_by(id_unit) %>%
    mutate(id_time = id_time_vec) %>%
    ungroup()

  return(list(fit = fit, ps = fitted, dat = pscore_dat, lin_pred = lin_pred))
}

#' MSM treatment estimation
#' @export
estimate_outcome <- function(formula, data, weights) {
  ## estimate effect by the weighted ls
  fit <- lm(formula, data = data, weights = weights )
}

#' Construct weigths
#' @param treatment A matrix of treatments where rows correspond to units.
#' @param ps A matrix of propensity scores.
#' @export
construct_weights <- function(treatment, ps) {
  weights <- rep(1, nrow(treatment))
  for(i in seq_len(nrow(treatment))) {
    weights <- weights * (treatment[,i] / ps[,i] + (1 - treatment[,i]) / (1 - ps[,i]))
  }
  return(weights)
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


#' Variance estmator for IPW estimator
#' @export
estimate_variance <- function(Y, D, ps, na_omit = TRUE) {
  ## remove NA
  dat_complete <- tibble::tibble(Yc = Y, Dc = D, psc = ps)
  if (isTRUE(na_omit)) {
    dat_complete <- na.omit(dat_complete)
  } else {
    ## otherwise replace NA pscore with threshold values
    dat_complete <- dat_complete %>%
      mutate(psc = if_else(is.na(psc),
                            if_else(Dc == 1, 0.9999, 0.0001),
                            psc
                          ))
  }


  ##
  ## Compute super-population variance
  ##
  ##    Var(tau) = Var(Y(1)) + Var(Y(0))
  ##
  ## Hajek
  ##    Var(Y(1)) = E[U^2] with U = D * (Y - E[Y(1)]) / ps
  ##    Var(Y(0)) = E[U^2] with U = (1 - D) * (Y - E[Y(0)]) / (1 - ps)
  ##
  ## HT
  ##    Var(Y(1)) = E[U^2] with U = D * Y / ps - E[Y(1)]
  ##    Var(Y(0)) = E[U^2] with U = (1 - D) * Y / (1 - ps) - E[Y(0)]
  ##
  ## where we replace E[Y(d)] with the estimate
  ##

  var_ht <- with(dat_complete, estimator_var_ht(Yc, Dc, psc, trim = TRUE))
  var_hj <- with(dat_complete, estimator_var_hajek(Yc, Dc, psc, trim = TRUE))
  var_est <- c(var_ht, var_hj)
  names(var_est) <- c("HT", "Hajek")
  return(var_est)
}


estimator_var_hajek <- function(Y, D, ps, trim = TRUE) {
  if(isTRUE(trim)) {
    ps <- pmin(ps, 0.9999)
    ps <- pmax(ps, 0.0001)
  }
  tmp_up1   <- sum((D * Y) / ps)
  tmp_down1 <- sum(D / ps)
  tmp_up0   <- sum(((1-D) * Y) /(1 - ps))
  tmp_down0 <- sum((1 - D) / (1 - ps))

  ## mean potential outcome
  EY1 <- tmp_up1 / tmp_down1
  EY0 <- tmp_up0 / tmp_down0

  ## score
  U1 <- D * (Y - EY1) / ps
  U0 <- (1 - D) * (Y - EY0) / (1 - ps)

  var_est  <- mean(U1^2) + mean(U0^2)
  return(var_est / length(D))
}


estimator_var_ht <- function(Y, D, ps, trim = TRUE) {
  if(isTRUE(trim)) {
    ps <- pmin(ps, 0.9999)
    ps <- pmax(ps, 0.0001)
  }

  ## mean potential outcome
  EY1 <- sum(D * Y / ps) / length(D)
  EY0 <- sum((1 - D) * Y / (1 - ps)) / length(D)

  ## score
  U1 <- D * Y / ps - EY1
  U0 <- (1 - D) * Y / (1 - ps) - EY0

  ## variance
  var_est <- mean(U1^2) + mean(U0^2)
  return(var_est / length(D))
}
