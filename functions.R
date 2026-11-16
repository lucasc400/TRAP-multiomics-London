# =============================================================================
# functions.R
#
# Helper functions and settings sourced by the analysis scripts. Packages are
# loaded by the calling scripts: lme4 (Multilmer) and nlme (run_mvn);
# parallel and FactoMineR are called with `::`.
#
#   Settings         significance thresholds, exposures, number of cores
#   Pre-processing   Multilmer()
#   MWAS             run_mvn(), mvn_by_exposure()
#   Clustering       GroupingSummary()
# =============================================================================


# Settings ----------------------------------------------------------------

# Per-test significance levels (see the Methods section of the paper)
alpha_metab <- 0.05 / 6040 # metabolomics, Bonferroni over the 6,040 features
alpha_metab_ent <- 0.05 / 284 # metabolomics, effective number of tests
alpha_mrna <- 0.05 / 202 # transcriptomics, effective number of tests

# TRAP exposures, as named in the data (CBLK = black carbon, BC in the paper)
ListExpo <- c("NO2", "PCNT", "PM25", "PM10", "CBLK")

# Number of cores used to run the MWAS in parallel
n_cores <- as.integer(Sys.getenv("N_CORES", unset = "5"))


# Pre-processing ----------------------------------------------------------

#' Linear mixed models fitted to each column of a data matrix
#'
#' Fits `y ~ 1 + fixed + (1 | random)` for each column of `ydata`. Random
#' effects with a variance below `tol` are dropped. Used to remove technical
#' variation ("denoising"): the residuals returned are the observations minus
#' the random intercepts.
#'
#' @param ydata matrix of outcomes (samples x features).
#' @param xdata data frame with the fixed and random effect variables.
#' @param fixed names of the fixed effect variables (optional).
#' @param random names of the random intercept variables.
#' @param var name of a fixed effect to test (continuous or binary only).
#' @param tol tolerance below which a random effect variance is set to zero.
#' @param progressbar show a progress bar.
#' @param residuals return the denoised data.
#' @return list with the number of observations, the coefficients for `var`
#'   (if provided), the random intercepts, the denoised data (if
#'   `residuals = TRUE`) and the intra-class correlation of each random effect.
Multilmer <- function(ydata, xdata, fixed = NULL, random, var = NULL, tol = 1e-5, progressbar = TRUE, residuals = FALSE) {
  if (!is.null(var)) {
    coefficients <- matrix(NA, nrow = ncol(ydata), ncol = 3)
    rownames(coefficients) <- colnames(ydata)
    colnames(coefficients) <- c("coef", "coef.se", "pval")
  }
  nobs <- rep(NA, ncol(ydata))
  names(nobs) <- colnames(ydata)

  if (residuals) {
    Residuals <- matrix(NA, nrow = nrow(ydata), ncol = ncol(ydata))
    rownames(Residuals) <- rownames(ydata)
    colnames(Residuals) <- colnames(ydata)
  }

  # Random intercepts: one matrix (levels x features) per random effect
  random_intercepts <- list()
  for (r in 1:length(random)) {
    mylevels <- unique(as.character(xdata[, random[r]]))
    ri <- matrix(NA, ncol = ncol(ydata), nrow = length(mylevels))
    rownames(ri) <- mylevels
    colnames(ri) <- colnames(ydata)
    random_intercepts[[r]] <- ri
  }

  if (progressbar) {
    progbar <- txtProgressBar(style = 3)
  }

  # Intra-class correlation
  icc <- matrix(NA, nrow = ncol(ydata), ncol = length(random))
  rownames(icc) <- colnames(ydata)
  colnames(icc) <- random

  for (k in 1:ncol(ydata)) {
    if (progressbar) {
      setTxtProgressBar(progbar, k / ncol(ydata))
    }
    x <- ydata[, k]
    nobs[k] <- nrow(na.exclude(cbind(x, xdata[c(fixed, random)])))
    if (!is.null(fixed)) {
      fe <- paste(fixed, collapse = "+")
    }
    re <- paste(paste0("(1|", random, ")"), collapse = "+")
    if (!is.null(fixed)) {
      f <- paste0("x~1+", fe, "+", re)
    } else {
      f <- paste0("x~1+", re)
    }
    model <- lmer(as.formula(f),
      data = xdata, REML = FALSE,
      control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = tol))
    )
    vcor <- data.frame(VarCorr(model))
    rownames(vcor) <- vcor[, 1]
    vcov.sum <- sum(vcor[, "vcov"])
    vcor <- vcor[random, "vcov"]
    icc[k, ] <- vcor / vcov.sum

    newrandom <- random[vcor > tol]
    if (length(newrandom) == 0) {
      # No random effect left: linear model
      if (!is.null(fixed)) {
        fe <- paste(fixed, collapse = "+")
        f <- paste0("x~1+", fe)
        model <- lm(as.formula(f), data = xdata)
        if (residuals) {
          resids <- residuals(model) + data.matrix(xdata[, fixed]) %*% coef(model)[fixed] + coef(model)[1]
          Residuals[, k] <- resids
        }
      } else {
        if (residuals) {
          Residuals[, k] <- x
        }
      }
      if (!is.null(var)) {
        newfixed <- fixed[!fixed %in% var]
        if (length(newfixed) == 0) {
          fe <- 1
        } else {
          fe <- paste(newfixed, collapse = "+")
        }
        f <- paste0("x~1+", fe)
        model0 <- lm(as.formula(f), data = xdata)
        myanova <- anova(model, model0)
        coefficients[k, 3] <- myanova$`Pr(>F)`[2]
      }
    } else {
      if (any(vcor < tol)) {
        # Some but not all random effect terms are smaller than tolerance
        re <- paste(paste0("(1|", newrandom, ")"), collapse = "+")
        if (!is.null(fixed)) {
          f <- paste0("x~1+", fe, "+", re)
        } else {
          f <- paste0("x~1+", re)
        }
        model <- lmer(as.formula(f),
          data = xdata, REML = FALSE,
          control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = tol))
        )
      }
      for (r in 1:length(newrandom)) {
        absolute_rid <- which(random == newrandom[r])
        rintercepts <- ranef(model)[[newrandom[r]]]
        random_intercepts[[absolute_rid]][rownames(rintercepts), k] <- unlist(rintercepts)
      }
      if (!is.null(var)) {
        newfixed <- fixed[!fixed %in% var]
        fe <- paste(newfixed, collapse = "+")
        re <- paste(paste0("(1|", newrandom, ")"), collapse = "+")
        f <- paste0("x~1+", fe, "+", re)
        model0 <- lmer(as.formula(f),
          data = xdata, REML = FALSE,
          control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = tol))
        )
        myanova <- anova(model, model0)
        if (sum(grepl(var, rownames(summary(model)$coefficients))) > 1) {
          stop("This code only works for continuous or binary 'var'. Please create dummy variables.")
        }
        coefficients[k, 3] <- myanova$`Pr(>Chisq)`[2]
      }
      if (residuals) {
        resids <- residuals(model) + data.matrix(xdata[, fixed]) %*% fixef(model)[fixed] + fixef(model)[1]
        Residuals[, k] <- resids
      }
    }
    if (!is.null(var)) {
      coefficients[k, 1:2] <- summary(model)$coefficients[grep(var, rownames(summary(model)$coefficients)), 1:2]
    }
  }
  names(random_intercepts) <- random

  out <- list(nobs = nobs)
  if (!is.null(var)) {
    out$coefficients <- as.data.frame(coefficients, stringsAsFactors = FALSE)
  }
  out$random_intercepts <- random_intercepts
  if (residuals) {
    out$residuals <- Residuals
  }
  out$icc <- icc
  return(out)
}


# MWAS --------------------------------------------------------------------

#' Multivariate normal (MVN) models for repeated measures
#'
#' For each feature (column of `ydata`), fits a generalised least squares
#' model `feature ~ rhs` with an unstructured within-subject correlation
#' (`corSymm` by `ID`) and a separate residual variance for each level of
#' `var_group` (by default `wave`, the site x time point combination). Each
#' model is fitted on the observations with complete data for the variables
#' it uses.
#'
#' @param ydata matrix of outcomes (samples x features), in the same row
#'   order as `xdata`.
#' @param xdata data frame with `ID`, `var_group` and the variables in `rhs`.
#' @param rhs right-hand side of the model formula, for example
#'   `"Expo.delta + Age + sex + BMI + Group"`.
#' @param var_group variable defining the groups of observations with their
#'   own residual variance.
#' @param on_failure what to do if a model cannot be fitted: `"stop"`, `"NA"`
#'   (missing estimates if the model gives an error) or `"NA_warning"`
#'   (missing estimates if the model gives an error or a warning).
#' @return list of three data frames, `pval`, `beta` and `tval`, with one row
#'   per feature (`OMIC_ID`) and one column per model term.
run_mvn <- function(ydata, xdata, rhs, var_group = "wave", on_failure = c("stop", "NA", "NA_warning")) {
  on_failure <- match.arg(on_failure)
  model_formula <- as.formula(paste("Omic ~", rhs))
  weights_formula <- as.formula(paste("~ 1 |", var_group))
  model_vars <- unique(c(all.vars(model_formula), "ID", var_group))

  fit_gls <- function(data) {
    gls(model_formula,
      data = data,
      correlation = corSymm(form = ~ 1 | ID),
      weights = varIdent(form = weights_formula),
      na.action = na.pass,
      control = list(singular.ok = TRUE)
    )
  }

  tables <- vector("list", ncol(ydata))
  for (k in seq_len(ncol(ydata))) {
    xdata$Omic <- ydata[, k]
    data_k <- xdata[complete.cases(xdata[, model_vars]), ]
    model <- switch(on_failure,
      "stop" = fit_gls(data_k),
      "NA" = tryCatch(fit_gls(data_k), error = function(e) NULL),
      "NA_warning" = tryCatch(fit_gls(data_k), error = function(e) NULL, warning = function(w) NULL)
    )
    if (!is.null(model)) {
      tables[[k]] <- summary(model)$tTable
    }
  }

  fitted <- which(!sapply(tables, is.null))
  if (length(fitted) == 0) {
    stop("None of the models could be fitted.")
  }
  if (length(fitted) < ncol(ydata)) {
    message(ncol(ydata) - length(fitted), " model(s) could not be fitted and are set to NA.")
  }
  model_terms <- rownames(tables[[fitted[1]]])

  extract <- function(column) {
    values <- vapply(tables, function(tt) {
      if (is.null(tt)) {
        return(rep(NA_real_, length(model_terms)))
      }
      tt[match(model_terms, rownames(tt)), column]
    }, numeric(length(model_terms)))
    values <- matrix(values, ncol = length(model_terms), byrow = TRUE, dimnames = list(colnames(ydata), model_terms))
    data.frame(OMIC_ID = colnames(ydata), values, check.names = FALSE)
  }
  return(list(pval = extract("p-value"), beta = extract("Value"), tval = extract("t-value")))
}

#' MVN models for each TRAP exposure, run in parallel
#'
#' For each exposure, sets `Expo.delta` to the change in that exposure
#' (column `<exposure>_delta` of `xdata`), runs run_mvn() and saves the
#' results to `sprintf(out_file, <exposure>)`.
#'
#' @param ydata,xdata,rhs passed to run_mvn(); `rhs` should include `Expo.delta`.
#' @param out_file file name template, for example
#'   `"Results/MVN/metab/%s_metab_mvn.rds"`.
#' @param ... further arguments passed to run_mvn().
#' @param exposures exposures to run.
mvn_by_exposure <- function(ydata, xdata, rhs, out_file, ..., exposures = ListExpo) {
  cl <- parallel::makeCluster(min(n_cores, length(exposures)))
  on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, library(nlme))
  parallel::clusterExport(cl, c("run_mvn", "mvn_one_exposure"))
  .mvn_y <- ydata
  .mvn_x <- xdata
  parallel::clusterExport(cl, c(".mvn_y", ".mvn_x"), envir = environment())
  invisible(parallel::parLapply(cl, exposures, mvn_one_exposure, rhs = rhs, out_file = out_file, ...))
}

# Worker for mvn_by_exposure(): uses the data exported to the worker
mvn_one_exposure <- function(Expo, rhs, out_file, ...) {
  .mvn_x$Expo.delta <- .mvn_x[[paste0(Expo, "_delta")]]
  saveRDS(run_mvn(.mvn_y, .mvn_x, rhs = rhs, ...), sprintf(out_file, Expo))
}


# Clustering --------------------------------------------------------------

#' Summarise groups of features by their first principal component
#'
#' @param x data matrix (samples x features).
#' @param group named vector of group memberships (names are features).
#' @param summary "centroid" (PC1 scores) or "medoid" (feature contributing
#'   most to PC1).
#' @param verbose show a progress bar.
#' @return list with `xsummarised` (samples x groups, scaled) and
#'   `description` (one row per feature: group, range of pairwise
#'   correlations within the group, variance explained by PC1, loadings and
#'   contributions).
GroupingSummary <- function(x, group, summary = "centroid", verbose = TRUE) {
  # Getting data transpose
  x <- t(x)

  # Number of groups
  k <- max(group)

  # Initialisation of the objects
  xsummarised <- matrix(NA, nrow = k, ncol = ncol(x))
  myrownames <- rep(NA, length = k)
  mytable <- data.frame(
    feature = names(group), group = unname(group),
    min_cor = rep(1, length(group)),
    max_cor = rep(1, length(group)),
    ev = rep(NA, length(group)),
    loadings = rep(NA, length(group)),
    contribution = rep(NA, length(group))
  )
  rownames(mytable) <- mytable[, 1]

  # Loop over the groups
  if (verbose) {
    pb <- utils::txtProgressBar(style = 3)
  }
  for (i in 1:k) {
    # Extracting data in the group
    xcluster <- x[names(group)[which(group == i)], , drop = FALSE]

    # PCA on the group
    mypca <- FactoMineR::PCA(t(xcluster), graph = FALSE, ncp = 1)

    # Proportion of explained variance by PC1
    mytable[mytable$group == i, "ev"] <- mypca$eig[1, 2] / 100

    # Loadings coefficients on PC1 (sign chosen so that most loadings are positive)
    if ((sum(mypca$var$coord > 0) / length(mypca$var$coord)) > 0.5) {
      mytable[rownames(xcluster), "loadings"] <- mypca$var$coord / sqrt(mypca$eig[1, 1])
    } else {
      mytable[rownames(xcluster), "loadings"] <- -mypca$var$coord / sqrt(mypca$eig[1, 1])
    }

    # Contribution of original variables to PC1
    mytable[rownames(xcluster), "contribution"] <- mypca$var$contrib / 100

    # Range in pairwise correlations
    if (nrow(xcluster) > 1) {
      mycor <- stats::cor(t(xcluster))
      mycor[upper.tri(mycor, diag = TRUE)] <- NA
      mytable[rownames(xcluster), "min_cor"] <- min(mycor, na.rm = TRUE)
      mytable[rownames(xcluster), "max_cor"] <- max(mycor, na.rm = TRUE)
    }

    # Using PC1-scores to summarise the group
    if (summary == "centroid") {
      if ((sum(mypca$var$coord > 0) / length(mypca$var$coord)) > 0.5) {
        xsummarised[i, ] <- mypca$ind$coord[, 1]
      } else {
        xsummarised[i, ] <- -mypca$ind$coord[, 1]
      }
    }

    # Using the variable most contributing to PC1-scores for summary. In a
    # group of two features both contribute equally, so the medoid is chosen
    # by numerical precision and can differ between platforms.
    if (summary == "medoid") {
      xsummarised[i, ] <- xcluster[which.max(abs(mypca$var$coord)), ]
      myrownames[i] <- rownames(xcluster)[which.max(abs(mypca$var$coord))]
    }

    if (verbose) {
      utils::setTxtProgressBar(pb, i / k)
    }
  }

  # Re-ordering the description table (multivariate groups first)
  mytable <- mytable[sort.list(mytable$group, decreasing = TRUE), ]

  # Returning the data with summarised features as columns
  colnames(xsummarised) <- colnames(x)
  if (summary == "centroid") {
    rownames(xsummarised) <- paste0("F", 1:k)
  } else {
    rownames(xsummarised) <- myrownames
  }
  xsummarised <- t(xsummarised)
  xsummarised <- scale(xsummarised)

  return(list(xsummarised = xsummarised, description = mytable))
}
