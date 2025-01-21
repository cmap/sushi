#' Computing area under a 4 parameter log-logistic dose response curve
#'
#' @param LL : Lower asymptote
#' @param UL : Upper asymptote
#' @param Inflection : Inflection point (EC50)
#' @param Slope : Hill slope ( > 0 for decreasing curves)
#' @param md : Minimum dose
#' @param MD : Maximum dose
#'
#' @return AUC value of the dose-response function (x: log2(dose), y: response) scaled to the dose range.
#' @export
#'
#' @examples
compute_auc <- function(LL, UL, Inflection, Slope, md, MD) {
  f1 = function(x) pmax(pmin((UL + (LL - UL)/(1 + (2^x/Inflection)^Slope)), 1, na.rm = T), 0, na.rm = T)
  return(tryCatch(integrate(f1, log2(md), log2(MD))$value/(log2(MD/md)),
                  error = function(e) {print(e); NA}))
}


#' Computing IC50 for a 4 parameter log-logistic dose response curve
#'
#' @param LL : Lower asymptote
#' @param UL : Upper asymptote
#' @param Inflection : Inflection point (EC50)
#' @param Slope : Hill slope ( > 0 for decreasing curves)
#' @param md : Minimum dose
#' @param MD : Maximum dose
#'
#' @return IC50 : The dose value where the curve intersects with y = 0.5, NA returned if they don't intersect in the given dose range.
#' @export
#'
#' @examples
compute_log_ic50 <- function(LL, UL, Inflection, Slope, md, MD) {
  print("Computing IC50 parameters ...")
  print(list(LL, UL, Inflection, Slope, md, MD))
  if((LL >= 0.5) | (UL <= 0.5)) {
    return(NA)
  } else {
    f1 = function(x) (UL + (LL - UL)/(1 + (2^x/Inflection)^Slope)- 0.5)
    return(tryCatch(uniroot(f1, c(log2(md), log2(MD)))$root,  error = function(x) NA))
  }
}

#' Computing mean square and median absolute errors for for a 4 parameter log-logistic dose response curve and the corresponding (dose,viability) pairs.
#'
#' @param FC : Measured viability vector
#' @param dose : Dose vector corresponding to FC
#' @param UL : Upper asymptote
#' @param LL : Lower asymptote
#' @param Slope : Hill slope
#' @param Inflection : Inflection point (EC50)
#'
#' @return List of mse and mad values.
#' @export
#'
#' @examples
compute_MSE_MAD <- function(FC, dose,  UL, LL,  Slope, Inflection) {
  FC.pred = UL  + (LL -UL )/(1 + (dose/Inflection)^Slope)
  residuals = FC - FC.pred
  return(list(mse = mean(residuals^2), mad = median(abs(residuals))))
}


#' Fitting a dose response curve to the given dose and viability (FC) values.
#' This function fits 5 different dose-response functiosn to the given dose, viability pairs using dr4pl and drc packages and returns
#' the best one (lowest MSE) among them.
#'
#' @param FC : Measured viability vector
#' @param dose : Dose vector corresponding to FC
#' @param UL_low : Lower limit for the upper asymptote
#' @param UL_up : Upper limit for the upper asympotote
#' @param slope_decreasing: Should the curve to be constrained to be decreasing or not.
#'
#' @return Returns a single row data-frame with following columns:
#'          fit_name : Name of the fitting method with the highest explained variance (lowest MSE)
#'          Lower_Limit : Lower asmpytote for the selected fit
#'          Upper_Limit : Upper asymptote for the selected fit
#'          Slope : The Hill slope for the selected fit
#'          Inflection : Inflection point of the selected fit (EC50)
#'          MSE : MSE of the selected fit
#'          MAD : MAD of the selected fit
#'          frac_var_explained : Explained variance for the best fit
#'          successful_fit: If any of the fitting methods has positive frac_var_explained
#'          AUC_Riemann : The average measured viability value across doses
#'          md : Minimum measured dose
#'          MD : Maximum measured dose
#'          AUC : AUC of the log-dose vs viability curve normalized to the measured dose-range (takes values between 0 and 1)
#'          log2.IC50 : Log2 of IC50 for the fitted curve
#' @export
#'
#' @examples
get_best_fit <- function(FC, dose, UL_low=0.8, UL_up=1.01, slope_decreasing=TRUE) {
  require(dr4pl)
  require(drc)
  require(tidyverse)
  require(magrittr)


  # Fits a number of alternate models  to the DRC and chooses the best fit.

  # UL low is the lowerbound of UL we pass to the optimizer and UL_up is the upper bound of UL that we pass to the optimizer
  # fomat of output will be:-
  # results.df <- data.frame("fit_name"=character(),
  #                          "Lower_Limit"=double(),
  #                          "Upper_Limit"=double(),
  #                          "Slope"=double(),
  #                          "Inflection"=double(),
  #                          "MSE"=double(), "MAD" =double(),
  #                          "frac_var_explained"=double())

  dose = as.numeric(dose)
  FC = as.numeric(FC)
  slope_bound <- ifelse(slope_decreasing, 1e-5, Inf)  # bound the slopes by default unless passed another option
  riemann_AUC <- mean(pmin(1,FC)) ## mean fold-change after rounding FC to 1.
  var_data = var(FC)

  md = min(dose); MD = max(dose)

  results.df <- list(); ix = 1


  # FIT 1 ---
  drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                  fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-slope_bound,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)),
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package


  if (drc_model$fit$convergence & all(is.finite(unlist(drc_model$coefficients)))){
    print("fit 1 parameters ...")
    print(drc_model$coefficients)

    mse_mad <- compute_MSE_MAD(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                               -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
    # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.

    results.df[[ix]] <- tibble( fit_name = "drc_drm_constrained",
                                Lower_Limit = as.numeric(drc_model$coefficients[[2]]),
                                Upper_Limit = as.numeric(drc_model$coefficients[[3]]),
                                Slope = -as.numeric(drc_model$coefficients[[1]]),
                                Inflection = as.numeric(drc_model$coefficients[[4]]),
                                MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }



  # FIT 2 ---
  drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                  fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                  lowerl = c(-Inf,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)),
                         error = function(e)
                         {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
  # "slope" in drc package is -ve of slope in dr4pl package


  if (drc_model$fit$convergence & all(is.finite(unlist(drc_model$coefficients)))){
    print("fit 2 parameters ...")
    print(drc_model$coefficients)

    if((!slope_decreasing) | (as.numeric(drc_model$coefficients[[1]]) > 0)){
      mse_mad <- compute_MSE_MAD(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                                 -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
      # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.

      results.df[[ix]] <- tibble( fit_name = "drc_drm_unconstrained",
                                  Lower_Limit = as.numeric(drc_model$coefficients[[2]]),
                                  Upper_Limit = as.numeric(drc_model$coefficients[[3]]),
                                  Slope = -as.numeric(drc_model$coefficients[[1]]),
                                  Inflection = as.numeric(drc_model$coefficients[[4]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
  }


  # FIT 3 ---
  dr4pl_initMan_optNM <- tryCatch(dr4pl(dose, FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_2 = 8*min(dose),
                                                                       theta_3= -3, theta_4 = 0.01),
                                        lowerl = c(UL_low, -Inf, -Inf, 0),
                                        upperl = c(UL_up, Inf, slope_bound, 1.01),
                                        method.optim="Nelder-Mead"),
                                  error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )

  if (!dr4pl_initMan_optNM$convergence){
    if (!is.null(dr4pl_initMan_optNM$dr4pl.robust)) {
      dr4pl_initMan_optNM <- dr4pl_initMan_optNM$dr4pl.robust
    }
  }

  if (dr4pl_initMan_optNM$convergence & all(is.finite(unlist(dr4pl_initMan_optNM$parameters)))){
    print("fit 3 parameters ...")
    print(dr4pl_initMan_optNM$parameters)

    mse_mad <- compute_MSE_MAD(FC, dose, dr4pl_initMan_optNM$parameters[[1]], dr4pl_initMan_optNM$parameters[[4]],
                               dr4pl_initMan_optNM$parameters[[3]], dr4pl_initMan_optNM$parameters[[2]])

    results.df[[ix]] <- tibble( fit_name = "dr4pl_initMan_constrained_optNM",
                                Lower_Limit = as.numeric(dr4pl_initMan_optNM$parameters[[4]]),
                                Upper_Limit = as.numeric(dr4pl_initMan_optNM$parameters[[1]]),
                                Slope = as.numeric(dr4pl_initMan_optNM$parameters[[3]]),
                                Inflection = as.numeric(dr4pl_initMan_optNM$parameters[[2]]),
                                MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }

  # FIT 4 ---
  dr4pl_unconstrained <- tryCatch(dr4pl(dose, FC,
                                        init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.3),
                                        method.init = "logistic",
                                        lowerl = c(0.99, -Inf, -Inf, 0),
                                        upperl = c(1.01, Inf, Inf, 1.01)),
                                  error = function(e) {print(e); return(NA)})

  if (!all(is.na(dr4pl_unconstrained))) {
    if (!dr4pl_unconstrained$convergence) {
      dr4pl_unconstrained <- dr4pl_unconstrained$dr4pl.robust
    }
  }


  param <- tryCatch(dr4pl_unconstrained$parameters, error = function(e) return(NA))
  if (!all(is.na(param))){
    print("fit 4 parameters ...")
    print(dr4pl_unconstrained$parameters)

    if(as.numeric(dr4pl_unconstrained$parameters[[3]])<slope_bound){ ### while slope bound is not passed to this last optimizer, we do not accept a solution not within the bound
      mse_mad <- compute_MSE_MAD(FC, dose, dr4pl_unconstrained$parameters[[1]], dr4pl_unconstrained$parameters[[4]],
                                 dr4pl_unconstrained$parameters[[3]], dr4pl_unconstrained$parameters[[2]])
      results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_unconstrained",
                                  Lower_Limit = as.numeric(dr4pl_unconstrained$parameters[[4]]),
                                  Upper_Limit = as.numeric(dr4pl_unconstrained$parameters[[1]]),
                                  Slope = as.numeric(dr4pl_unconstrained$parameters[[3]]),
                                  Inflection = as.numeric(dr4pl_unconstrained$parameters[[2]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
  }

  # FIT 5 ---
  dr4pl_initL <- tryCatch(dr4pl(dose, FC,
                                init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.005),
                                method.init = "logistic",
                                lowerl = c(UL_low, -Inf, -Inf, 0),
                                upperl = c(UL_up, Inf, slope_bound, 1.01)),
                          error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )

  if (!dr4pl_initL$convergence){
    if (!is.null(dr4pl_initL$dr4pl.robust)) {
      dr4pl_initL <- dr4pl_initL$dr4pl.robust
    }
  }

  if (dr4pl_initL$convergence & all(is.finite(unlist(dr4pl_initL$parameters)))){
    print("fit 5 parameters ...")
    print(dr4pl_initL$parameters)

    mse_mad <- compute_MSE_MAD(FC,dose, dr4pl_initL$parameters[[1]], dr4pl_initL$parameters[[4]],
                               dr4pl_initL$parameters[[3]], dr4pl_initL$parameters[[2]])

    results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_constrained",
                                Lower_Limit = as.numeric(dr4pl_initL$parameters[[4]]),
                                U1pper_Limit = as.numeric(dr4pl_initL$parameters[[1]]),
                                Slope = as.numeric(dr4pl_initL$parameters[[3]]),
                                Inflection = as.numeric(dr4pl_initL$parameters[[2]]),
                                MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
    ix = ix + 1
  }


  # Choose the best fit among the successful fits---
  results.df <- dplyr::bind_rows(results.df)

  if (nrow(results.df)>0){
    results.df <- results.df %>%
      dplyr::arrange(desc(frac_var_explained)) %>%
      head(1) %>%
      dplyr::mutate(successful_fit = TRUE,
                    AUC_Riemann = as.numeric(riemann_AUC),
                    md = md, MD = MD) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(AUC = compute_auc(Lower_Limit, Upper_Limit, Inflection, Slope, md, MD),
                    log2.IC50 = compute_log_ic50(Lower_Limit, Upper_Limit, Inflection, Slope, md, MD))

  }else{
    results.df  <- data.frame(successful_fit=FALSE,
                              AUC_Riemann = as.numeric(riemann_AUC),
                              md = md, MD = MD, AUC = NA, log2.IC50 = NA)
  }

  return (results.df)
}




#' Creating the dose-response table for a single treatment viability screens.
#'
#' @param l2fc : The path to look for the log fold change input file.
#' @param out_path : Where to write the output file, make sure it ends with a "/". If left NULL it will write the output file to in_path.
#' @param cell_line_cols : List of columns to identify each cell line (model). Default: c("depmap_id", "cell_set", "pool_id").
#' @param treatment_cols : List of columns to idenfity each treatment. Default:  c("pert_id", "x_project_id", "pert_name", "pert_plate")
#' @param dose_column : Name of the column that represents the dose information. Default: "pert_dose"
#' @param l2fc_column : Name of the column that holds log2-viability values. Default: "l2fc"
#' @param type_column : Name of the column that contains the perturbation type. Only the rows that has value "trt_cp" will be kept in the calculations. Default: "pert_type"
#' @param cap_for_viability : The upper threshold for the viability values for before fitting the curves. Default is 1.5, any viability value above this value will be made equal to 1.5.
#'
#' @return It writes the calculated dose-response parameters to a file named as output_file_name in the out_path folder.
#' @export
#'
#' @examples create_drc_table(in_path = "data/sushi_io/testing_MTS-SEQ002-KF/test/")
create_drc_table <- function(LFC = l2fc,
                              cell_line_cols =  c("depmap_id", "cell_set", "pool_id"),
                              treatment_cols = c("pert_id", "x_project_id", "pert_name", "pert_plate"),
                              dose_column = "pert_dose", l2fc_column = "l2fc", type_column = "pert_type",
                              cap_for_viability = 1.5) {
  require(data.table)
  require(tidyverse)
  require(rlang)

  # check if the input table populated
  if(nrow(LFC) == 0){
    stop(paste0(LFC, " is empty!"))
  }

  # check for missing critical columns
  necessary_columns <- unique(c(dose_column, l2fc_column, type_column, treatment_cols, cell_line_cols))
  if(any(!necessary_columns %in% colnames(LFC))){
    stop(paste0(paste0(setdiff(necessary_columns, colnames(LFC)), collapse = ", "),
                " are missing in ", LFC))
  }

  # Check pert_type contains trt_cp
  LFC <- dplyr::filter(LFC, .data[[type_column]] == "trt_cp")

  if(nrow(LFC) == 0){
    stop(paste0(type_column, " doesn't contain any trt_cp!"))
  }


  # Fit the curves for single compounds ----
  DRC_SINGLE <- LFC %>%
    dplyr::filter(.data[[type_column]] == "trt_cp") %>%
    dplyr::mutate(dose_ = as.numeric(.data[[dose_column]]),
                  l2fc_ = as.numeric(.data[[l2fc_column]])) %>%
    dplyr::filter(is.finite(dose_), is.finite(l2fc_)) %>%
    dplyr::group_by(across(all_of(c(cell_line_cols, treatment_cols)))) %>%
    dplyr::filter(length(unique(dose_)) > 4) %>%
    dplyr::summarise(get_best_fit(FC = pmin(2^l2fc_, cap_for_viability),
                                  dose = dose_)) %>%
    dplyr::ungroup()
}