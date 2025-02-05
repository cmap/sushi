#' Fits a simple linear model by regressing y over each column of X.
#'
#' @param X : matrix of n by m, it can have NA's in it, columns should be named.
#' @param y : vector of m, elements are in the same order of with the rows of X, there cannot be any NAs.
#' @param v.th : Minimum variance for the columns of X, columns with smaller variances will be dropped.
#' @param n.min : Minimum number of finite pairs between column of X and y, columns not satisfying this condition will be dropped.
#'
#' @return A data frame with: x (corresponding column of X), rho (Pearson correlation), beta (regression coefficient),
#'         p.val.rob / q.val.rob (robust p-value / q-value), p.val / q.val (homoskedastic / traditional p-value / q-value),
#'         n (number of non-na samples), ns (a stability score to represent how many points of y needs to be replaced with
#'         the mean value to flip the sign of regression coefficient).
#' @export
#'
#' @examples
#'
robust_linear_model <- function(X, y, v.th = 0.0025, n.min = 25) {
  X <- scale(X, center = TRUE, scale = FALSE)
  n = colSums(!is.na(X))
  X <- X[, (colMeans(X^2, na.rm = T) >= v.th) & (n >= n.min), drop = FALSE]
  n = colSums(!is.na(X))
  Y = matrix(y, dim(X)[1], dim(X)[2])
  Y <- scale(Y, center = TRUE, scale = FALSE)
  S = X * Y
  muS = colMeans(S, na.rm = T)
  sigmaS = sqrt(rowMeans((t(S) - muS)^2, na.rm = T))
  varX = colMeans(X^2, na.rm = T)
  varY = colMeans(Y^2, na.rm = T)

  # Estimation and inference
  beta = muS / varX
  rho = muS / sqrt(varX * varY)
  p.val <- 2*pt(abs(sqrt(n-2) * muS/sigmaS),
                df = n-2, lower.tail = FALSE)
  p.val.homoskedastic <- 2*pt(sqrt( (n-2) /(1/rho^2-1)),
                              df = n-2, lower.tail = FALSE)
  p.val.homoskedastic[near(abs(rho), 1)] <- 0

  # Experimental robustness score
  ns = apply(S,2, function(s) sum(cumsum(sort(s/sum(s, na.rm = T), decreasing = T)) < 1, na.rm = T)) + 1
  ns = pmin(ns, n - apply(X, 2, function(x) sort(table(x), decreasing = TRUE)[1]))

  res <- data.frame(x = names(beta), rho = rho, beta = beta,
                    p.val.rob = p.val,
                    p.val = p.val.homoskedastic,
                    q.val.rob = p.adjust(p.val, method = "BH"),
                    q.val = p.adjust(p.val.homoskedastic, method = "BH"),
                    n = n,
                    ns = ns)

  return(res)
}



#' Fits a random forest from X to y using k-fold cross validation. Returns
#'
#' @param X : Feature matrix.
#' @param y : Response vector.
#' @param k : Number of cross-validation folds
#' @param vc : Variance threshold, columns of X with lower variance than this threshold are dropped.
#' @param lm : Parameter for gausscov package for feature selection within each fold.
#' @param p0 : Parameter for gausscov package for feature selection within each fold.
#' @param folds : Fold structure for the cross-validation. If NULL, stratified folds are created
#' @param X.test : If provided the the predictions are provided as part of the outputs.
#'
#' @return model_table: Model performance parameters and feature importances.
#'         predictions: Predictions.
#'         folds = Cross-validation folds.
#'         y.test.hat: Test predictions if X.test is provided.
#'
#' @export
#'
#' @examples
random_forest <- function (X, y, k = 5, vc = 0.01, lm = 25, p0 = 0.01, folds = NULL, X.test = NULL) {
  y.clean <- y[is.finite(y)]
  X.clean <- X[, apply(X, 2, function(x) all(is.finite(x)))]
  cl <- sample(dplyr::intersect(rownames(X.clean), names(y.clean)))
  if(is.null(folds)){
    N = floor(length(cl)/k)
    folds <- sapply(1:k, function(kx) cl[((kx - 1) * N + 1):(kx * N)])
  }

  X.clean <- X.clean[cl, ]
  y.clean <- y.clean[cl]
  colnames(X.clean) %<>% make.names()
  yhat_rf <- rep(NA, length(y.clean))
  names(yhat_rf) <- names(y.clean)
  SS = tibble()

  if(!is.null(X.test)){
    y.test.hat = list()
    px = 1
  }

  for (kx in 1:k) {
    test <- intersect(folds[,kx], cl)
    train <- dplyr::setdiff(cl, test)
    X_train <- X.clean[train, , drop = F]
    X_test <- X.clean[test, , drop = F]
    X_train <- X_train[, apply(X_train, 2, var) > vc, drop = F]
    b <- tryCatch(gausscov::f2st(y = y.clean[train], x = X_train, lm = lm, p0 = p0, sub = FALSE), error = function(e) NULL)
    if(is.null(b)){
      cc = cor(X_train, y.clean[train])
      features <- union(rownames(cc[rank(-cc) <= 100, ,drop = FALSE]), rownames(cc[rank(cc) <= 100, ,drop = FALSE]))
    }else{
      features <- colnames(X_train)[b[[1]][, 2]]
    }

    if (length(features) > 0) {
      X_train <- X_train[, features, drop = F]
      rf <- ranger::ranger(y = y.clean[train], x = X_train,
                           importance = "impurity")
      yhat_rf[test] <- predict(rf, data = as.data.frame(X_test[, colnames(X_train), drop = F]))$predictions

      if(!is.null(X.test)){
        X.test.temp <- X.test[, colnames(X_train), drop = F]
        X.test.temp <- X.test.temp[apply(is.finite(X.test.temp), 1, all), , drop = FALSE]
        y.test.hat[[px]] = tibble(y.hat = predict(rf, data = as.data.frame(X.test.temp))$predictions,
                                  sample_id = rownames(X.test))
        px = px + 1
      }

      ss <- tibble::tibble(feature = names(rf$variable.importance),
                           RF.imp = rf$variable.importance/sum(rf$variable.importance),
                           fold = kx)
      SS %<>% dplyr::bind_rows(ss)
    }
  }
  if (nrow(SS) > 0) {
    RF.importances <- SS %>% dplyr::distinct(feature, RF.imp,
                                             fold) %>% reshape2::acast(feature ~ fold, value.var = "RF.imp",
                                                                       fill = 0)
    RF.table <- tibble::tibble(feature = rownames(RF.importances),
                               RF.imp.mean = RF.importances %>% apply(1, mean),
                               RF.imp.sd = RF.importances %>% apply(1, function(x) sd(x,
                                                                                      na.rm = T)), RF.imp.stability = RF.importances %>%
                                 apply(1, function(x) sum(x > 0)/k)) %>% dplyr::filter(feature !=
                                                                                         "(Intercept)") %>% dplyr::arrange(desc(RF.imp.mean)) %>%
      dplyr::mutate(rank = 1:n())
    mse <- mean((yhat_rf - y.clean)^2, na.rm = T)
    mse.se <- sqrt(var((yhat_rf - y.clean)^2, na.rm = T))/sqrt(length(y.clean))
    r2 <- 1 - (mse/var(y.clean, na.rm = T))
    ps <- cor(y.clean, yhat_rf, use = "pairwise.complete.obs")
    RF.table %<>% dplyr::mutate(MSE = mse, MSE.se = mse.se,
                                R2 = r2, PearsonScore = ps)

    if(is.null(X.test)){
      return(list(model_table = RF.table, predictions = yhat_rf, folds = folds))
    }else{
      y.test.hat <- dplyr::bind_rows(y.test.hat)
      return(list(model_table = RF.table, predictions = yhat_rf, folds = folds, y.test.hat = y.test.hat))
    }
  }
  else {
    if(is.null(X.test)){
      return(list(model_table = tibble(), predictions = tibble(), folds = folds))
    }
    else{
      return(list(model_table = tibble(), predictions = tibble(), folds = folds, y.test.hat = tibble()))
    }
  }
}



#' Calculates and returns univariate analysis results by correlating each column of Y with features sets from depmap.org.
#'
#' @param Y : Matrix n x m, make sure rownames are depmap_id's and columns are named. There can be NA's.
#' @param file : Please point out to the downloaded depmap_datasets.h5 file.
#' @param features : You can give a subset of available feature sets, but if left as NULL, it computes for everything.
#' @param homoskedastic : Should homoskedastic or robust q-values be used for ranking and filtering. Default is homoskedastic.
#' @param n.X.min : Results with less thana given sample size are dropped, default is 100
#' @param ns.min : Results that can be flipped with less than ns.min data points are dropped, the default is 3
#' @param q.val.max : Results with q-values less than q.val.max are dropped, default is 0.2
#' @param parallel : Should multiple cores to be used to decrease the compute time. Default is FALSE.
#'
#' @return Returns a data-table with each row corresponds to a particular (feature, feature set, y) triplet. See robust_linear_model for the other columns.
#' @export
#'
#' @examples
#'
univariate_biomarker_table <- function(Y, file = '/data/biomarker/current/depmap_datasets_public.h5',
                                       features = NULL,
                                       homoskedastic = TRUE, n.X.min = 100,
                                       ns.min = 3, q.val.max = .2,
                                       parallel = FALSE){
  require(tidyverse)
  require(magrittr)
  require(rhdf5)


  if(!is.matrix(Y)){
    Y <- as.matrix(Y)
  }

  if(is.null(features)){
    features <-  substr(setdiff(h5ls(file)$group, "/"),2,100)
  }else{
    features <- intersect(features, substr(setdiff(h5ls(file)$group, "/"),2,100))
  }
  print(features)

  RESULTS <- list()
  for(feat in features){

    X <- read_dataset(file , feat)
    cl = intersect(rownames(X), rownames(Y))


    X <- X[cl, ]
    X <- X[ , apply(X, 2, function(x) sum(is.finite(x))) >= n.X.min, drop = FALSE]

    print(paste0(feat, " - ", dim(X)[1] ,'x', dim(X)[2]))

    if((dim(X)[1] >= n.X.min) & dim(X)[2] > 0){
      f <- function(ix){
        y <- Y[cl, ix]; y <- y[is.finite(y)]

        res <- robust_linear_model(X = X[names(y), ], y = y, v.th = 0.0025, n.min = n.X.min) %>%
          as.tibble()%>%
          dplyr::rename(feature = x) %>%
          dplyr::mutate(feature.set = feat,
                        y = colnames(Y)[ix]) %>%
          dplyr::filter(ns >= ns.min)

        if(homoskedastic){
          res <- res %>%
            dplyr::filter(q.val <= q.val.max)
        } else{
          res <- res %>%
            dplyr::filter(q.val.rob <= q.val.max)
        }

        if(nrow(res) > 0){
          res <- res %>%
            dplyr::arrange(rho) %>%
            dplyr::mutate(rank = ifelse(sign(rho) < 0, 1:n(), NA)) %>%
            dplyr::arrange(-rho) %>%
            dplyr::mutate(rank = pmin(rank, ifelse(sign(rho > 0), 1:n(), NA), na.rm = T))
        }
        return(res)
      }

      if(parallel){
        require(parallel)
        require(parallelly)
        RES <- mclapply(1:dim(Y)[2], f, mc.cores = parallelly::availableCores() )
      }else{
        RES <- lapply(1:dim(Y)[2], f)
      }

      RESULTS[[feat]] <- dplyr::bind_rows(RES)
    }
  }

  RESULTS <- dplyr::bind_rows(RESULTS)
  return(RESULTS)
}



#' Exports individual datasets from depmap_datasets.h5 file. You can specify the row names either as ccle_names or depmap_ids.
#'
#' @param file : Please point out to the downloaded depmap_datasets.h5 file.
#' @param dataset : The dataset you want to export. You can check the available datasets with substr(setdiff(rhdf5::h5ls(file)$group, "/"),2,100)
#' @param rownames_depmap_ids : Default TRUE, you can get the rownames as ccle_names by switching to FALSE.
#'
#' @return The requested data matrix.
#' @export
#'
#' @examples
read_dataset <- function(file = '/data/biomarker/current/depmap_datasets_public.h5', dataset,
                         rownames_depmap_ids = TRUE) {
  require(rhdf5)
  if(word(file, sep = fixed("://")) %in% c("s3", "http", "https")){
    s3 = TRUE
    print(paste0("Reading ", file, " from S3"))
  } else{
    s3 = FALSE
    print(paste0("Reading ", file, " from local"))
  }
  X <- h5read(file, name = paste0(dataset, "/mat"), s3 = s3)
  row_meta <- h5read(file, name = paste0(dataset, "/row_meta"), s3 = s3)
  column_meta <- h5read(file, name = paste0(dataset, "/column_meta"), s3 = s3)
  colnames(X) <- column_meta$column_name
  if(rownames_depmap_ids){
    rownames(X) <- row_meta$ModelID
  }else{
    rownames(X) <- row_meta$CCLEName
  }

  X <- X[rownames(X) != "NA", colnames(X) != "NA", drop = FALSE]
  X <- X[!duplicated(rownames(X)), !duplicated(colnames(X)), drop = FALSE]
  return(X)
}


#' Exports individual columns of a selected dataset from depmap_datasets.h5 file. You can specify the row names either as ccle_names or depmap_ids.
#'
#' @param file : Please point out to the downloaded depmap_datasets.h5 file.
#' @param dataset : The dataset you want to export. You can check the available datasets with substr(setdiff(rhdf5::h5ls(file)$group, "/"),2,100)
#' @param feature_names : Particular column names to pull from the given dataset.
#' @param rownames_depmap_ids : Default TRUE, you can get the rownames as ccle_names by switching to FALSE.
#'
#' @return
#' @export
#'
#' @examples
read_features <- function(file = '/data/biomarker/current/depmap_datasets_public.h5',
                          dataset, feature_names = NULL,
                          rownames_depmap_ids = TRUE) {
  require(rhdf5)
  if(word(file, sep = fixed("://")) %in% c("s3", "http", "https")){
    s3 = TRUE
  } else{
    s3 = FALSE
  }

  row_meta <- h5read(file, name = paste0(dataset, "/row_meta"), s3 = s3)
  column_meta <- h5read(file, name = paste0(dataset, "/column_meta"), s3 = s3)

  if(!is.null(feature_names)){
    column_meta <- dplyr::filter(column_meta, column_name %in% feature_names)
  }

  if(nrow(column_meta) > 0){
    X <- h5read(file, name = paste0(dataset, "/mat"),
                index = list(NULL, column_meta$ix), s3 = s3)
    colnames(X) <- column_meta$column_name
    if(rownames_depmap_ids){
      rownames(X) <- row_meta$ModelID
    }else{
      rownames(X) <- row_meta$CCLEName
    }
    X <- X[rownames(X) != "NA", colnames(X) != "NA", drop = FALSE]
    X <- X[!duplicated(rownames(X)), !duplicated(colnames(X)), drop = FALSE]
    return(X)
  } else{
    return(NULL)
  }
}




#' Creates and returns two concatenated depmap feature matrices for given response and confounder matrices.
#' Returned matrices shares the same row names (depmap_id's). The rownames are limited to the intersection of the rownames of Y and W.
#' Additionally, any columns or rows that has more than 10% NA's are dropped.
#'
#' @param Y : Response matrix, rows are indexed by depmap_id's.
#' @param W : Confounder matrix, rows are indexed by depmap_id's.
#' @param file : File path for the compressed depmap datasets (.h5 file).
#'
#' @return Returns a list of two matrices X.DNA and X.RNA. X.DNA contains CopyNumber, Mutation, Fusion, Lineage, and OmicsSignatures features as columns.
#'         X.RNA additionally includes the RNA expression data as well.
#' @export
#'
#' @examples
RF_feature_sets <- function(Y, W = NULL, file = '/data/biomarker/current/depmap_datasets_public.h5') {
  require(magrittr)

  if(!is.null(W)){
    cl = intersect(rownames(Y), rownames(W))
    Y <- Y[cl, , drop = FALSE]
    W <- W[cl, , drop = FALSE]
    colnames(W) %<>% paste0("CONF_", .)
  }


  CN <- read_dataset(file = file, dataset = 'CopyNumber')
  colnames(CN) %<>% paste0("CN_", .)
  MUT <- read_dataset(file = file, dataset = 'Mutation')
  colnames(MUT) %<>% paste0("MUT_", .)
  FUS <- read_dataset(file = file, dataset = 'Fusion')
  colnames(FUS) %<>% paste0("FUS_", .)
  LIN <- read_dataset(file = file, dataset = 'Lineage')
  colnames(LIN) %<>% paste0("LIN_", .)
  SIG <- read_dataset(file = file, dataset = 'Signatures')
  colnames(SIG) %<>% paste0("SIG_", .)

  cl = intersect(rownames(MUT), rownames(CN)) %>%
    intersect(rownames(FUS)) %>%
    intersect(rownames(LIN)) %>%
    intersect(rownames(SIG)) %>%
    intersect(rownames(Y))

  X.DNA <- cbind(CN[cl, ], MUT[cl, ]) %>%
    cbind(FUS[cl, ]) %>%
    cbind(LIN[cl, ]) %>%
    cbind(SIG[cl, ])

  if(!is.null(W)){
    X.DNA <- cbind(X.DNA, W[cl, , drop = FALSE])
  }

  X.DNA <- X.DNA[rowMeans(is.finite(X.DNA)) > 0.9, ,drop = TRUE]

  rm(CN, MUT, FUS, LIN, SIG)
  EXP <- read_dataset(file = file, dataset = 'Expression')
  colnames(EXP) %<>% paste0("EXP_", .)


  cl = intersect(rownames(X.DNA), rownames(EXP))

  X.RNA <- EXP[cl, ]
  rm(EXP)

  X.RNA <- X.RNA[rowMeans(is.finite(X.RNA)) > 0.9, ,drop = TRUE]

  colnames(X.DNA) %<>% make.names()
  colnames(X.RNA) %<>% make.names()

  return(list(X.DNA = X.DNA,  X.RNA = X.RNA))
}


#' Fit two random forest models (DNA and DNA+RNA) for each column of Y and return the
#' summary results: Cross validated R2, PearsonScore, feature importances (mean, sd, stability)
#'
#' @param Y : Response matrix, rows are depmap_id's. 2 models are fit for each column. Columns should be named.
#' @param W : Confounder matrix. Should not contain and NA's, rows are indexed by depmap_ids, and aligned with Y. Default is NULL.
#' @param file : Path to the h5 file containing depmap datasets.
#' @param k : Number of cross validation folds. Default is 10.
#'
#' @return Returns a data frame containing the feature importances and model properties for each column of Y (y) and model (DNA or DNA+RNA)
#' @export
#'
#' @examples
multivariate_biomarker_table <- function(Y, W = NULL, file = '/data/biomarker/current/depmap_datasets_public.h5', k = 10) {
  X <- RF_feature_sets(Y, W = W, file = file)
  cl <- intersect(rownames(X$X.DNA), rownames(X$X.RNA))

  rf.DNA <- list(); rf.RNA <- list()
  for(ix in 1:dim(Y)[2]){

    y = Y[,ix]; y = y[is.finite(y)]
    cl_ = intersect(cl, names(y))

    rf_DNA <- random_forest(X$X.DNA[cl_, ], y[cl_], k = k)
    rf_RNA <- random_forest(cbind(X$X.DNA[cl_, ], X$X.RNA[cl_, ]), y[cl_], folds = rf_DNA$folds, k = k)

    rf.DNA[[ix]] <- rf_DNA$model_table %>%
      dplyr::mutate(y = colnames(Y)[ix],
                    model = "DNA") %>%
      dplyr::filter(RF.imp.stability > 0.5)

    rf.RNA[[ix]] <- rf_RNA$model_table %>%
      dplyr::mutate(y = colnames(Y)[ix],
                    model = "DNA+RNA") %>%
      dplyr::filter(RF.imp.stability > 0.5)

    print(paste0(ix, " out of ", dim(Y)[2]))
  }

  RF <- dplyr::bind_rows(dplyr::bind_rows(rf.DNA), dplyr::bind_rows(rf.RNA))

  return(RF)
}








#' Function that creates the multivariate biomarker results from the input file.
#'
#' @param in_path : Path for the input file.
#' @param out_path : Path to write the output file.
#' @param input_file_name : Name of the input file.
#' @param output_file_name : Name of the output file.
#' @param cell_line_column : Index for of cell lines. Input file should have this column. The strong default is depmap_id.
#' @param depmap_file : File path for the .h5 files for the depmap datasets.
#' @param treatment_columns : Set of columns of the input file to identify each profile.
#' @param response_column : Column of the input file that has the respose file.
#' @param aggregate_function : How to aggregate if the treatment columns and cell line column doesn't uniquely identifies the response value.
#' @param transform_function : Should the response variable be transformed before the analysis. Default is identity, a common choice is log2.
#'
#' @return
#' @export
#'
#' @examples
#' create_multivariate_biomarker_table( in_path = "data/MTS_SEQ003_KF/",
#'                                      treatment_columns = c("x_project_id", "pert_plate", "pert_name", "pert_dose"),
#'                                      depmap_file = "~/Downloads/mts_sequencing_analysis/data/depmap_datasets_internal_24q4.h5")
#'
#' create_multivariate_biomarker_table( in_path = "data/MTS_SEQ003_KF/", input_file_name = "drc.csv",
#'                                      output_file_name = "l2auc_multivariate_biomarkers",
#'                                      treatment_columns = c("x_project_id", "pert_plate", "pert_name"),
#'                                      response_column = "auc", transform_function = log2,
#'                                      depmap_file = "~/Downloads/mts_sequencing_analysis/data/depmap_datasets_internal_24q4.h5")
#'
create_multivariate_biomarker_table <- function(in_path, out_path = NULL,
                                                output_file_name = "l2fc_multivariate_biomarkers",
                                                depmap_file,
                                                treatment_columns = c("pert_id", "x_project_id", "pert_name", "pert_plate", "pert_dose"),
                                                response_column = "median_l2fc", aggregate_function = median, transform_function = function(x){x}) {
  require(data.table)
  require(tidyverse)
  require(rlang)

  # check if the input path exists
  if(!file.exists(in_path)){
    stop("Input path does not exist!")
  }

  # load the input file
  input_file <- data.table::fread(in_path)

  # check if the input table populated
  if(nrow(input_file) == 0){
    stop(paste0("Input file for multivariate biomarker is empty!"))
  }

  # check for missing critical columns
  necessary_columns <- unique(c(response_column, treatment_columns, "depmap_id"))

  if(any(!necessary_columns %in% colnames(input_file))){
    stop(paste0(paste0(setdiff(necessary_columns, colnames(input_file)), collapse = ", "),
                " are missing from multivariate biomarker input file!"))
  }

  # cast the responses into a matrix
  M <- input_file %>%
    dplyr::filter(is.finite(.data[[response_column]])) %>%
    tidyr::unite(cn, treatment_columns, sep = "::") %>%
    reshape2::acast(depmap_id ~ cn, value.var = response_column, fun.aggregate = aggregate_function) %>%
    transform_function()

  # generate the biomarker table
  multivariate_biomarker_table <- multivariate_biomarker_table(Y = M, file = depmap_file, k = 10)

  # Export the biomarker table -----
  print(paste0("Writing the multivariate output file to ", paste0(out_path, "/", output_file_name)))
  input_file %>%
    dplyr::select_at(vars(treatment_columns)) %>%
    dplyr::distinct() %>%
    tidyr::unite(y, treatment_columns, sep = "::", remove = FALSE) %>%
    dplyr::inner_join(multivariate_biomarker_table) %>%
    dplyr::select(-y) %>%
    write_csv(paste0(out_path, "/", output_file_name))
}



#' Function that creates the univariate biomarker results from the input file.
#'
#' @param in_path : Path for the input file.
#' @param out_path : Path to write the output file.
#' @param input_file_name : Name of the input file.
#' @param output_file_name : Name of the output file.
#' @param cell_line_column : Index for of cell lines. Input file should have this column. The strong default is depmap_id.
#' @param depmap_file : File path for the .h5 files for the depmap datasets.
#' @param treatment_columns : Set of columns of the input file to identify each profile.
#' @param response_column : Column of the input file that has the respose file.
#' @param aggregate_function : How to aggregate if the treatment columns and cell line column doesn't uniquely identifies the response value.
#' @param transform_function : Should the response variable be transformed before the analysis. Default is identity, a common choice is log2.
#' @param features : Set of depmap datasets to be used in the analysis. The default NULL makes sure all the datasets in the dpemap_file are used.
#' @param min_sample_size : Minimum sample size for a feature to be included in the analysis. The default is 100.
#' @param parallel : Should the analysis be parallelized to multiple cores.
#'
#' @return
#' @export
#'
#' @examples
#' create_univariate_biomarker_table( in_path = "data/MTS_SEQ003_KF/", input_file_name = "drc.csv",
#'                                    output_file_name = "l2auc_univariate_biomarkers",
#'                                    treatment_columns = c("x_project_id", "pert_plate", "pert_name"),
#'                                    response_column = "auc", transform_function = log2,
#'                                    depmap_file = "~/Downloads/mts_sequencing_analysis/data/depmap_datasets_internal_24q4.h5",
#'                                    parallel = TRUE)
#'
#' create_univariate_biomarker_table( in_path = "data/MTS_SEQ003_KF/",
#'                                    treatment_columns = c("x_project_id", "pert_plate", "pert_name", "pert_dose"),
#'                                    depmap_file = "~/Downloads/mts_sequencing_analysis/data/depmap_datasets_internal_24q4.h5",
#'                                    parallel = TRUE)
create_univariate_biomarker_table <- function(in_path, out_path,
                                              output_file_name = "l2fc_univariate_biomarkers", depmap_file,
                                                treatment_columns = c("pert_id", "x_project_id", "pert_name", "pert_plate", "pert_dose"),
                                                response_column = "median_l2fc", aggregate_function = median, transform_function = function(x){x},
                                              features = NULL, min_sample_size = 100, parallel = FALSE) {
  require(data.table)
  require(tidyverse)
  require(rlang)

  # check if the input path exists
  if(!file.exists(in_path)){
    stop("Input path does not exist!")
  }

  # load the input file
  input_file <- data.table::fread(in_path)

  # check if the input table populated
  if(nrow(input_file) == 0){
    stop(paste0("Input file for univariate biomarker is empty!"))
  }

  # check for missing critical columns
  necessary_columns <- unique(c(response_column, treatment_columns, "depmap_id"))
  print("treatment_columns:")
  print(treatment_columns)
  print("necessary_columns:")
  print(necessary_columns)
  print("input_columns:")
  print(colnames(input_file))
  print("response_column:")
  print(response_column)

  if(any(!necessary_columns %in% colnames(input_file))){
    stop(paste0(paste0(setdiff(necessary_columns, colnames(input_file)), collapse = ", "),
                " are missing in univariate biomarker input file!"))
  }

  # cast the responses into a matrix
  M <- input_file %>%
    dplyr::filter(is.finite(.data[[response_column]])) %>%
    tidyr::unite(cn, treatment_columns, sep = "::") %>%
    reshape2::acast(depmap_id ~ cn, value.var = response_column, fun.aggregate = aggregate_function) %>%
    transform_function()

  # generate the biomarker table
  univariate_biomarker_table <- univariate_biomarker_table(Y = M, features = features, file = depmap_file,  n.X.min = min_sample_size, parallel =  parallel)

  # Export the biomarker table -----
  print(paste0("Writing the univariate output file to ", paste0(out_path, "/", output_file_name)))
  input_file %>%
    dplyr::select_at(vars(treatment_columns)) %>%
    dplyr::distinct() %>%
    tidyr::unite(y, treatment_columns, sep = "::", remove = FALSE) %>%
    dplyr::inner_join(univariate_biomarker_table) %>%
    dplyr::select(-y) %>%
    write_csv(paste0(out_path, "/", output_file_name))
}