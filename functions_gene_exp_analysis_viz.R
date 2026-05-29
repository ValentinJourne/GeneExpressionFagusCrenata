#####################################
#Gene expression analysis in Fagus crenata
#all functions create and used
#####################################

#creator : Valentin Journ\'e; Kyushu University
#contact : journe.valentin@gmail.com
#file created on 2024-09-26
#updated on 2026 - 06

#####################################
#loading package
#####################################
#*Note, I used GPT to clean the header of the functions

# Detach all non-base R packages from the current session.
# This utility function removes all loaded packages except the core
# base R packages, helping to avoid namespace conflicts and ensuring
# a clean working environment before running analyses.
detachAllPackages <- function() {
  basic.packages <- c(
    "package:stats",
    "package:graphics",
    "package:grDevices",
    "package:utils",
    "package:datasets",
    "package:methods",
    "package:base"
  )

  package.list <- search()[ifelse(
    unlist(gregexpr("package:", search())) == 1,
    TRUE,
    FALSE
  )]

  package.list <- setdiff(package.list, basic.packages)

  if (length(package.list) > 0)
    for (package in package.list) detach(package, character.only = TRUE)
}


# Utility function to load R packages.
# Packages that are not already installed are automatically installed
# and then loaded into the current session.
using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require, character.only = TRUE))
  need <- libs[req == FALSE]
  if (length(need) > 0) {
    install.packages(need)
    lapply(need, require, character.only = TRUE)
  }
}

# Calculate the coefficient of variation (CV = SD / Mean)
# Used to quantify inter-annual variability in gene expression or flowering.
# Returns NA when the mean is zero or missing.
my_cv_fun <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(NA_real_)
  s / m
}

# Calculate inter-annual variability in gene expression.
# For each gene and individual, monthly expression values are aggregated
# within years (sum of expression across sampling months), and the
# coefficient of variation (CV) is calculated across years.
# Gene-level variability is then obtained by averaging CV values across
# individuals. The transformed coefficient of variation (kCV) is also
# computed following Lobry et al 2023 MEE, providing a bounded measure
# of variability ranging from 0 to 1.
functioncalc_gene_variability_cv = function(data) {
  gene_year_expr <- data %>%
    group_by(qseqid, TreeID, year) %>%
    summarise(
      year_expr = sum(level.exp, na.rm = TRUE), # or sum / window-based
      .groups = "drop"
    )
  cv_gene <- gene_year_expr %>%
    group_by(qseqid, TreeID) %>%
    summarise(
      CVi = my_cv_fun(year_expr),
      .groups = "drop"
    ) %>%
    group_by(qseqid) %>%
    summarise(
      CV = mean(CVi, na.rm = TRUE),
      CV_median = median(CVi, na.rm = TRUE),
      nobs = sum(!is.na(CVi)),
      .groups = "drop"
    ) %>%
    mutate(kCV = sqrt(CV^2 / (1 + CV^2)))
  return(cv_gene)
}

# Calculate synchrony in gene expression among individuals.
# For each gene, expression values are reshaped into a time-by-individual matrix,
# and pairwise correlations among individuals are calculated across sampling dates.
# Synchrony is defined as the mean pairwise correlation among individuals.
# The percentage of missing pairwise correlations is also reported (value would be the same for all genes)
functioncalc_gene_synchrony = function(
  data,
  method = c("pearson", "spearman")
) {
  method <- match.arg(method) #make sure only one matching

  synchrony_gene <- data %>%
    group_by(qseqid) %>%
    nest() %>%
    mutate(
      cor_mat = map(
        data,
        ~ {
          df <- .x %>%
            dplyr::select(TreeID, dateformat, level.exp) %>%
            pivot_wider(names_from = TreeID, values_from = level.exp)
          #cor matrix and then caluclate only for pairs with lower trim
          mat <- cor(
            df[, -1],
            use = "pairwise.complete.obs",
            method = method
          )
          mean_cor <- mean(mat[lower.tri(mat)], na.rm = TRUE)

          percna = sum(is.na(mat)) / prod(dim(mat))

          tibble(mean_synchrony = mean_cor, percentage.na = percna)
        }
      )
    ) %>%
    unnest(cor_mat) %>%
    dplyr::select(qseqid, mean_synchrony, percentage.na)
}


# Format daily climate data and extract temporal variables.
# Converts date strings into POSIXct format and derives day, month, year,
# calendar date, and day-of-year (DOY) information for subsequent analyses.
format_climate_temperature = function(data) {
  data %>%
    mutate(date.bis = as.POSIXct(Date, format = "%Y/%m/%d")) %>%
    mutate(
      day = day(date.bis),
      month = month(date.bis),
      year = year(date.bis)
    ) %>%
    mutate(
      datetime = as.Date(paste(year, month, day, sep = "-")),
      DOY = yday(datetime)
    )
}

# Standardize a numeric vector using z-score transformation.
# Values are centered on the mean and scaled by the standard deviation,
# resulting in a variable with mean = 0 and standard deviation = 1.
# Missing values are ignored when calculating summary statistics.
#I am using this for the heatmaps
cal_z_score <- function(x) {
  (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
}

# Replace a specified character string with a new value.
# Intended for cleaning imported datasets by substituting predefined
replace_string <- function(x) {
  if (is.character(x)) {
    x <- ifelse(x == string_to_replace, new_value, x)
  }
  return(x)
}

# Calculate the standard error of the mean (SEM).
std.error <- function(x, na.rm = TRUE) sd(x, na.rm = TRUE) / sqrt(length(x))

# Reshape Boruta importance history into a data frame for plotting or downstream filtering.
# The function extracts finite importance values from a Boruta object, optionally includes
# selected shadow attributes, assigns colors according to Boruta decisions
# (Confirmed, Tentative, Rejected, Shadow), orders features by median importance,
# and returns a wide data frame of importance values.
# If requested, the output can be filtered to retain only selected features.
#the function is based on the original one
reshape_the_Boruta_data <- function(
  x,
  whichShadow = c(TRUE, TRUE, TRUE),
  colCode = c('green', 'yellow', 'red', 'blue'),
  col = NULL,
  filter.output = F,
  select.features = c("a", "b")
) {
  if (is.null(x$ImpHistory))
    stop('Importance history was not stored during the Boruta run.')

  #Removal of -Infs and conversion to a list
  lz <- lapply(
    1:ncol(x$ImpHistory),
    function(i) x$ImpHistory[is.finite(x$ImpHistory[, i]), i]
  )
  colnames(x$ImpHistory) -> names(lz)

  #Selection of shadow meta-attributes
  numShadow <- sum(whichShadow)
  lz[c(rep(TRUE, length(x$finalDecision)), whichShadow)] -> lz

  generateCol <- function(x, colCode, col, numShadow) {
    #Checking arguments
    if (is.null(col) & length(colCode) != 4)
      stop('colCode should have 4 elements.')
    #Generating col
    if (is.null(col)) {
      rep(colCode[4], length(x$finalDecision) + numShadow) -> cc
      cc[c(x$finalDecision == 'Confirmed', rep(FALSE, numShadow))] <- colCode[1]
      cc[c(x$finalDecision == 'Tentative', rep(FALSE, numShadow))] <- colCode[2]
      cc[c(x$finalDecision == 'Rejected', rep(FALSE, numShadow))] <- colCode[3]
      col = cc
    }
    return(col)
  }

  #Generating color vector
  col <- generateCol(x, colCode, col, numShadow)

  #Ordering boxes due to attribute median importance
  ii <- order(sapply(lz, stats::median))
  lz[ii] -> lz
  col <- col[ii]
  lz_df <- do.call(rbind.data.frame, lz)
  df <- as.data.frame(t(lz_df))
  names(df) <- names(lz)
  rownames(df) <- NULL
  # NEW: filter output if requested
  if (filter.output && !is.null(select.features)) {
    select.features <- intersect(select.features, colnames(df))
    df <- df[, select.features, drop = FALSE]
  }
  return(df)
}

# Transform proportional data for beta regression.
# Applies the Smithson & Verkuilen (2006) adjustment to move values
# bounded at 0 and 1 into the open interval (0,1), which is required
# for standard beta regression models.
# The transformation depends on the sample size and preserves the
# relative ordering of observations.
y.transf.betareg <- function(y) {
  n.obs <- sum(!is.na(y))
  (y * (n.obs - 1) + 0.5) / n.obs
}

# Compute gene expression summaries across all possible monthly windows.
# For each selected gene, individual tree, and year, this function identifies
# all available combinations of sampled months and calculates both the mean
# and cumulative log2 expression within each window.
# The output contains one row per TreeID, year, gene, and month-window
# combination, and is used as input for logistic or beta regression analyses.
compute_window_summary <- function(gene_ids, data) {
  data %>%
    filter(qseqid %in% gene_ids) %>%
    group_by(TreeID, year, qseqid) %>%
    arrange(month) %>%
    summarise(data = list(pick(everything())), .groups = "drop") %>%
    mutate(
      window_data = map(
        data,
        ~ {
          df <- .x
          months_available <- unique(df$month)
          all_combos <- unlist(
            lapply(
              1:length(months_available),
              function(k) combn(months_available, k, simplify = FALSE)
            ),
            recursive = FALSE
          )

          windows <- map_dfr(all_combos, function(mset) {
            window_df <- df %>% filter(month %in% mset)
            tibble(
              window_months = paste(sort(mset), collapse = "-"),
              mean_expr = mean(window_df$level.exp.log2, na.rm = TRUE),
              cumsum_expr = sum(window_df$level.exp.log2, na.rm = TRUE)
            )
          })

          windows
        }
      )
    ) %>%
    dplyr::select(-data) %>% #dplyr here !
    unnest(window_data)
}

# Extract model summaries, confidence intervals, goodness-of-fit metrics, and AUC values.
# This function takes a data frame containing fitted models in a list-column,
# removes rows with NULL models, and extracts tidy coefficient tables,
# confidence intervals, global model summaries, and model performance metrics.
# For logistic regression models fitted with glm, the area under the ROC curve
# (AUC) is additionally calculated from predicted probabilities and observed
# responses. Model objects and data columns are removed from the final output
extract_model_info <- function(model_df, model_col = "model_cumsum") {
  model_sym <- rlang::sym(model_col)

  model_df %>%
    filter(!map_lgl(!!model_sym, is.null)) %>% # Filter out NULL models early
    mutate(
      tidied_mod = map(
        !!model_sym,
        ~ tryCatch(broom::tidy(.x), error = function(e) tibble())
      ),
      fixefmod_mod = map(
        !!model_sym,
        ~ tryCatch(as.data.frame(confint(.x)), error = function(e) tibble())
      ),
      glanced_mod = map(
        !!model_sym,
        ~ tryCatch(broom::glance(.x), error = function(e) tibble())
      ),
      GoF_mod = map(
        !!model_sym,
        ~ {
          model <- .x
          perf <- tryCatch(
            performance::model_performance(model),
            error = function(e) tibble()
          )
          perf <- perf %>% dplyr::select(-any_of(c("AIC", "BIC")))

          if (inherits(model, "glm")) {
            pred <- tryCatch(
              predict(model, type = "response"),
              error = function(e) rep(NA, nobs(model))
            )
            actual <- tryCatch(model$y, error = function(e) NA)
            auc <- tryCatch(
              as.numeric(pROC::auc(pROC::roc(actual, pred))),
              error = function(e) NA_real_
            )
            perf$AUC <- auc
          } else {
            perf$AUC <- NA_real_
          }
          perf
        }
      )
    ) %>%
    unnest(
      c(tidied_mod, fixefmod_mod, glanced_mod, GoF_mod),
      keep_empty = TRUE
    ) %>%
    dplyr::select(-matches("model"), -data)
}

# Fit beta regression models relating flowering intensity to gene expression.
# For each gene and seasonal window, a beta regression model is fitted using
# flowering intensity (transformed to the open interval (0,1)) as the response
# and gene expression metrics (e.g., cumulative expression) as predictors.
# Models that fail to converge or produce errors are safely returned as NULL.
fit.beta.regression = function(
  window_data,
  formula.fit.here = as.formula("flowering.percentage.trans ~ cumsum_expr")
) {
  window_data %>%
    group_by(qseqid, window_months) %>%
    nest() %>%
    mutate(
      model_cumsum = map(
        data,
        possibly(
          #need the possible here because sometime models are not workign
          ~ betareg(formula.fit.here, data = .x),
          otherwise = NULL
        )
      )
    )
}

# Fit logistic regression models relating flowering occurrence to gene expression.
# For each gene and seasonal window, a binomial generalized linear model (logit link)
# is fitted using flowering presence/absence as the response and gene expression
# metrics (e.g., cumulative expression) as predictors.
# One model is fitted per gene and seasonal window combination.
fit.logistic.regression = function(
  window_data,
  formula.fit.here = as.formula("fac.mastONOFF ~ cumsum_expr")
) {
  window_data %>%
    group_by(qseqid, window_months) %>%
    nest() %>%
    mutate(
      model_cumsum = map(
        data,
        ~ glm(formula.fit.here, data = .x, family = "binomial")
      )
    )
}

# Run regression analysis for one batch of gene-window expression summaries.
# The function reads a precomputed batch file, joins flowering intensity data
# by TreeID and year, fits either beta regression or logistic regression models
# for each gene and seasonal window, extracts model summaries and performance
# metrics, and saves the results as a .qs file. Existing output files are
# skipped to avoid rerunning completed batches.
#NOTE to MYSELF, at the time I created the code, I was using qs
#but now this packaage is unstable, so I might need to change and update with qs2 package
#or simply create CSV, but it takes too much space
run_regression_batch <- function(
  batch_index,
  input_dir,
  output_dir,
  response_formula = as.formula("flowering.percentage.trans ~ cumsum_expr"),
  average.flo.intensity.individual, # fruiting intensity dataset
  modeltype = "betareg"
) {
  library(qs)
  library(dplyr)
  library(purrr)
  library(betareg)
  library(performance)
  library(broom)
  library(rlang)

  # Define output file early to skip if exists
  output_file <- file.path(
    output_dir,
    paste0("statistic_summary_batch", batch_index, ".qs")
  )
  if (file.exists(output_file)) {
    message("Skipped batch ", batch_index, " (already exists)")
    return(NULL)
  }

  files <- list.files(
    path = input_dir,
    pattern = "\\.qs$",
    full.names = TRUE
  )

  if (batch_index > length(files)) {
    warning("Batch index ", batch_index, " out of bounds.")
    return(NULL)
  }

  window_data <- qs::qread(files[[batch_index]]) %>%
    left_join(average.flo.intensity.individual, by = c("TreeID", "year"))

  # Fit regression and extractinformatuoj
  if (modeltype == "betareg") {
    best_models_windows <- fit.beta.regression(
      window_data,
      formula.fit.here = response_formula
    )
  } else {
    best_models_windows <- fit.logistic.regression(
      window_data,
      formula.fit.here = response_formula
    )
  }

  summary_cumsum_gene <- extract_model_info(
    best_models_windows,
    model_col = "model_cumsum"
  )

  # Save output with low qs# Updata MAY 2026, maybe use qss2 now
  qs::qsave(summary_cumsum_gene, output_file)
  message("Saved batch ", batch_index, " → ", output_file)
}

# Generate a reverse day calendar for climate-window analyses.
# For each reference year, the function builds a calendar counting backwards
# from a specified reference day of year (refday), over a user-defined number
# of previous days. This is useful for aligning climate data to biological
# reference dates and extracting climate windows spanning one or more years.
# The output contains calendar dates, year, month, day-of-year, and reversed
# day indices for each reference year.
generate_reverse_day_calendar <- function(
  refday = 274,
  lastdays = 1095,
  yearback = 3,
  start_year = 1946
) {
  library(dplyr)
  library(lubridate)

  if (lastdays > 365 & yearback < 2) stop("For >365 days, set yearback >= 2")
  if (lastdays > 730 & yearback < 3) stop("For >730 days, set yearback >= 3")

  end_year <- start_year + yearback + 1

  DATE <- seq(
    as.Date(paste0(start_year, "-01-01")),
    as.Date(paste0(end_year, "-01-01")),
    by = "day"
  )
  dfata <- data.frame(
    DATE = DATE,
    YEAR = year(DATE),
    MONTHab = month.abb[month(DATE)],
    DOY = yday(DATE)
  )

  yearperiod <- start_year:(end_year - 1)
  vectotemp <- vector("list", length = length(yearperiod) - yearback)

  for (k in seq_along(vectotemp)) {
    yearsref <- yearperiod[k + yearback - 1]
    yearrefminusOne <- yearsref - yearback

    tt <- dfata %>%
      filter(YEAR >= yearrefminusOne & YEAR <= yearsref) %>%
      mutate(
        referenceFin = case_when(
          YEAR == yearsref & DOY == refday ~ 1,
          YEAR == yearsref & DOY > refday ~ NA_real_,
          TRUE ~ 0
        )
      ) %>%
      filter(!is.na(referenceFin)) %>%
      arrange(desc(DATE)) %>%
      mutate(days.reversed = 1:n())

    vectotemp[[k]] <- tt %>%
      filter(days.reversed < lastdays) %>%
      arrange(days.reversed) %>%
      mutate(YEAR = yearsref)
  }

  bind_rows(vectotemp)
}

# Extract yearly climate values from selected climate windows.
# For each best climate-window model, this function identifies the corresponding
# window of days using a reverse-day calendar, extracts the selected climate
# variable from the daily climate data for each year, and calculates the mean
# climate value within that window. Model information such as response,
# climate variable, statistic, function, Delta AICc, and window boundaries is
# retained in the output.
extract_climate_windows_by_year <- function(
  best_clim_mode,
  calendar_climate,
  temperature_df
) {
  library(dplyr)
  library(purrr)
  library(lubridate)

  result <- best_clim_mode %>%
    mutate(row_id = row_number()) %>%
    group_split(row_id) %>%
    map_dfr(function(row) {
      win_open <- row$WindowOpen
      win_close <- row$WindowClose
      clim_var <- as.character(row$climate)

      # Extract years present in temperature data
      years_in_data <- unique(temperature_df$year)

      # For each year: check if leap year and adjust win_open/win_close
      #i sued gpt to double check and improve that leap year
      window_by_year <- purrr::map_dfr(years_in_data, function(y) {
        is_leap <- lubridate::leap_year(y)
        win_open_adj <- if (is_leap) win_open + 1 else win_open
        win_close_adj <- if (is_leap) win_close + 1 else win_close

        # Get calendar DOYs for the adjusted window
        win_dates <- calendar_climate %>%
          dplyr::filter(
            days.reversed >= win_close,
            days.reversed <= win_open
          ) %>%
          dplyr::pull(DOY)

        # Adjust for leap year: add 1 to DOYs after Feb 28
        win_dates_adj <- if (is_leap) {
          ifelse(win_dates >= 60, win_dates + 1, win_dates)
        } else {
          win_dates
        }

        temp_merged <- temperature_df %>%
          dplyr::filter(year == y, DOY %in% win_dates_adj)

        if (nrow(temp_merged) == 0) return(NULL)

        temp_merged %>%
          dplyr::summarise(
            year = y,
            mean_value = mean(.data[[clim_var]], na.rm = TRUE),
            .groups = "drop"
          )
      })

      # Attach model info
      if (nrow(window_by_year) > 0) {
        window_by_year <- window_by_year %>%
          dplyr::mutate(
            response = row$response,
            climate = clim_var,
            type = row$type,
            stat = row$stat,
            func = row$func,
            DeltaAICc = row$DeltaAICc,
            WindowOpen = win_open,
            WindowClose = win_close
          )
      }

      window_by_year
    })

  return(result)
}


# Plot the predicted effect of one predictor from a fitted GAM.
# The function creates a prediction grid in which the focal predictor varies
# across its observed range while all other predictors are held at their mean.
# Predicted responses and 95% confidence intervals are then plotted as a smooth
# effect curve, optionally with the original observed data points overlaid.
plot_gam_effect <- function(
  model,
  data,
  predictor,
  response = "level.exp.log2",
  n_points = 100,
  pointdata = T
) {
  all_vars <- attr(terms(model), "term.labels")

  #seq grid
  grid <- data.frame(matrix(nrow = n_points, ncol = length(all_vars)))
  names(grid) <- all_vars

  for (v in all_vars) {
    if (v == predictor) {
      grid[[v]] <- seq(
        min(data[[v]], na.rm = TRUE),
        max(data[[v]], na.rm = TRUE),
        length.out = n_points
      )
    } else {
      grid[[v]] <- mean(data[[v]], na.rm = TRUE)
    }
  }

  # Predict
  preds <- predict(model, newdata = grid, se.fit = TRUE, type = "response")
  #intercept <- attr(preds, "constant"), if terms
  #grid$fit <- intercept + preds$fit[, paste0("s(", predictor, ")")]
  #grid$lower <- grid$fit - 1.96 * preds$se.fit[, paste0("s(", predictor, ")")]
  #grid$upper <- grid$fit + 1.96 * preds$se.fit[, paste0("s(", predictor, ")")]

  grid$fit <- preds$fit
  grid$lower <- preds$fit - 1.96 * preds$se.fit
  grid$upper <- preds$fit + 1.96 * preds$se.fit

  if (pointdata == T) {
    ggplot(grid, aes_string(x = predictor, y = "fit")) +
      geom_point(
        data = data,
        aes_string(x = predictor, y = response),
        alpha = 0.4
      ) +
      geom_line(color = "blue", linewidth = 0.8) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
      labs(title = paste("Effect of", predictor, "on", response)) +
      theme_minimal()
  } else {
    ggplot(grid, aes_string(x = predictor, y = "fit")) +
      geom_line(color = "blue", linewidth = 0.8) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
      labs(title = paste("Effect of", predictor, "on", response)) +
      theme_minimal()
  }
}

# Fit generalized additive models (GAMs) separately for each gene.
# For each gene, monthly expression and environmental/resource variables are
# first aggregated at the gene, year, TreeID level. A GAM is then fitted to
# model cumulative log2 gene expression as a smooth function of climate
# variables and resource traits, including temperature, rainfall, sunshine,
# nitrogen, NSC, and carbon. For each successfully fitted model, the function
# stores model summaries, smooth-effect plots, and derivatives of smooth terms
# to evaluate how gene expression changes along environmental gradients.
fit_gam_to_genes <- function(data, k = 4) {
  # List to store all model outputs
  all_results <- list()
  all_plots <- list()
  all_derivative <- list()
  unique_genes <- unique(data$qseqid)

  for (gene in unique_genes) {
    message("Fitting model for: ", gene)

    test <- data %>%
      filter(qseqid == gene) %>%
      group_by(qseqid, sseqid, year, TreeID) %>%
      summarise(
        level.exp.log2 = sum(level.exp.log2, na.rm = TRUE),
        rainfall = sum(rainfall, na.rm = TRUE),
        minTemp = mean(minTemp, na.rm = TRUE),
        maxTemp = mean(maxTemp, na.rm = TRUE),
        sunshine = mean(sunshine, na.rm = TRUE),
        NSC = mean(NSC, na.rm = TRUE),
        Nitrogen = mean(Nitrogen, na.rm = TRUE),
        Carbon = mean(Carbon, na.rm = TRUE),
        .groups = "drop"
      )

    # Skip if not enough data
    if (nrow(test) < 10) next

    # Fit GAM
    mod_gam <- tryCatch(
      {
        gam(
          level.exp.log2 ~
            s(minTemp, k = k) +
              s(rainfall, k = k) +
              s(sunshine, k = k) +
              s(Nitrogen, k = k) +
              s(NSC, k = k) +
              s(Carbon, k = k),
          #family = Gamma(link = "log"), #could not work efficiently
          data = test
        )
      },
      error = function(e) {
        warning("Model failed for: ", gene)
        return(NULL)
      }
    )

    if (!is.null(mod_gam)) {
      # Model summary
      summary_df <- cbind(
        broom::glance(mod_gam),
        broom::tidy(mod_gam)
      ) %>%
        mutate(qseqid = gene)

      all_results[[gene]] <- summary_df

      # Save plot
      p <- gratia::draw(mod_gam, residuals = TRUE)
      p[[1]][["labels"]][["title"]] <- paste0(gene, " s(minTemp)")

      all_plots[[gene]] <- p

      #get slope over gradient
      derivative <- gratia::derivatives(mod_gam, type = "central") %>%
        mutate(qseqid = gene)
      all_derivative[[gene]] = derivative
    }
  }

  # Combine results
  final_results <- bind_rows(all_results)
  final_derivative = bind_rows(all_derivative)
  return(list(
    results_table = final_results,
    results_derivative = final_derivative,
    plots = all_plots
  ))
}

# Summarize resource and nutrient variables across observations.
# The function aggregates carbohydrate, nutrient, and stoichiometric variables
# using either the mean or cumulative sum, depending on the selected method.
# Flowering status (fac.mastONOFF) is retained from the first observation.
# Intended for generating annual or seasonal summaries of resource availability.
summarise_elements_resources <- function(data, method = c("mean", "cumsum")) {
  method <- match.arg(method)

  summarise_fun <- switch(
    method,
    mean = ~ mean(.x, na.rm = TRUE),
    cumsum = ~ sum(.x, na.rm = TRUE)
  )

  data %>%
    summarise(
      across(
        c(
          starch,
          Fructose,
          Glucose,
          Maltose,
          Mannose,
          `Myo-Inositol`,
          Raffinose,
          Sucrose,
          SS,
          NSC,
          Carbon,
          Nitrogen,
          `C/N`
        ),
        summarise_fun
      ),
      fac.mastONOFF = first(fac.mastONOFF),
      .groups = "drop"
    )
}

# Format amino acid concentration data from a wide spreadsheet into long format.
# The function reconstructs column names from a two-row header, converts dates
# and concentration values to appropriate formats, and reshapes the data into
# a tidy long table. Compound names, sample identifiers, and measurement units
# are extracted from the original headers, producing one row per observation.
format.amino.acid.long = function(data) {
  header_rows <- data[1:2, ]
  data_rows <- data[-(1:2), ]
  header_matrix <- as.matrix(header_rows)
  # Fill missing compound names from the left
  for (j in 2:ncol(header_matrix)) {
    if (is.na(header_matrix[1, j]) || header_matrix[1, j] == "") {
      header_matrix[1, j] <- header_matrix[1, j - 1]
    }
  }
  headers <- paste0(header_matrix[1, ], "_", header_matrix[2, ])
  headers[1:2] <- c("Date", "Month")

  names(data_rows) <- headers
  data_rows <- data_rows %>%
    mutate(Date = as.Date(Date))
  data_rows <- data_rows %>%
    mutate(across(-c(Date, Month), as.numeric))

  long <- data_rows %>%
    pivot_longer(
      cols = -c(Date, Month),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    separate(
      Variable,
      into = c("Compound", "Sample"),
      sep = "_",
      remove = TRUE
    ) %>%
    drop_na(Value)

  long <- long %>%
    mutate(
      Unit = str_extract(Compound, "(?<=\\().+(?=\\))"), # extract text inside parentheses
      Compound = str_trim(str_remove(Compound, "\\s*\\(.*\\)")) # remove "(unit)" from name
    ) %>%
    dplyr::select(Date, Month, Compound, Unit, Sample, Value)
}


# Perform Gene Ontology (GO) enrichment analysis using GOstats.
# The function builds a gene set collection from a custom GO annotation table,
# defines a universe of annotated genes, and tests whether GO terms are
# over-represented in a user-defined gene list. Enrichment can be performed
# for Biological Process (BP) or Cellular Component (CC) ontology terms.
# The output includes raw and BH-adjusted p-values, odds ratios, expected
# counts, observed gene counts, universe counts, significance labels, and
# enrichment ratios.
#adapted from Kudoh et al 2025, Elife
GO_analysis <- function(
  list.genes,
  GO.gene.list,
  ontology = c("BP", "CC"),
  pvalueCutoff = 0.05
) {
  #load databse
  library(GOstats)
  library("GSEABase")
  library("qvalue")
  library(org.At.tair.db)
  library(AnnotationDbi)
  library(GO.db)

  goframeData <- data.frame(
    go_id = GO.gene.list$GO_ID,
    evidence = "IEA", # or omit if not available
    GeneID = GO.gene.list$qseqid
  )

  goFrame <- GOFrame(goframeData)
  goAllFrame <- GOAllFrame(goFrame)
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection()) #GSEABase

  all.Gene = unique(GO.gene.list$qseqid) #full list of GO genes

  #specify my list of genes of interest as vector
  gene_ls = list.genes$qseqid

  #now the function
  p <- GSEAGOHyperGParams(
    name = "Paramaters",
    geneSetCollection = gsc,
    geneIds = gene_ls,
    universeGeneIds = all.Gene,
    ontology = ontology, #or change either CC or BP
    pvalueCutoff = pvalueCutoff,
    conditional = TRUE,
    testDirection = "over"
  )

  result <- hyperGTest(p)

  # summarize
  raw_pval = pvalues(result)
  adj_pval = p.adjust(raw_pval, method = "BH") #to get the adjusted p-value
  odds = oddsRatios(result)
  expected_count = expectedCounts(result)
  gene_count = geneCounts(result)
  gene_univ_count = universeCounts(result)

  summary_table = data.frame(
    rank = 1:length(raw_pval),
    go_id = names(raw_pval),
    raw.pvalues = raw_pval,
    adjusted.pvalues = adj_pval,
    oddsRatios = odds,
    expectedCounts = expected_count,
    geneCounts = gene_count,
    universeCounts = gene_univ_count
  ) %>%
    mutate(
      significace = case_when(
        adjusted.pvalues < 0.001 ~ "***",
        adjusted.pvalues < 0.01 ~ "**",
        adjusted.pvalues < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      ratio = geneCounts / universeCounts
    )
  return(summary_table)
}


# Perform a paired Wilcoxon signed-rank test between years.
# The function compares matched observations from 2021 and 2022,
# pairing samples by their ordering after sorting by sample identity.
# If a column named 'concentration' is present, it is automatically
# renamed to 'Value' for compatibility. The Wilcoxon test statistic
# and associated p-value are returned.
#NOTE HERE YOU NEED 2021 and 2022
do.stat.element_gsh_sulfate_totalS = function(data) {
  summary_stats <- data %>%
    group_by(compound, Year) %>%
    summarise(
      n = n(),
      mean = mean(concentration, na.rm = TRUE),
      median = median(concentration, na.rm = TRUE),
      sd = sd(concentration, na.rm = TRUE),
      se = sd / sqrt(n),
      min = min(concentration, na.rm = TRUE),
      max = max(concentration, na.rm = TRUE),
      .groups = "drop"
    )

  print(summary_stats)

  #test on element

  compounds <- levels(data$compound)

  wilcox_results <- lapply(compounds, function(cmp) {
    g2021 <- data %>%
      filter(compound == cmp, Year == 2021) %>%
      arrange(Sample) %>% # CRITICAL part : same order for pairing
      pull(concentration)

    g2022 <- data %>%
      filter(compound == cmp, Year == 2022) %>%
      arrange(Sample) %>%
      pull(concentration)

    test <- wilcox.test(g2021, g2022, paired = TRUE, exact = FALSE)

    direction <- ifelse(
      mean(g2022) > mean(g2021),
      "HIGHER in 2022",
      "LOWER in 2022"
    )

    cat("\nCompound:", cmp, "\n")
    cat("  Mean 2021:", round(mean(g2021), 2), "\n")
    cat("  Mean 2022:", round(mean(g2022), 2), "\n")
    cat("  Direction:", direction, "\n")
    cat(
      "  V =",
      test$statistic,
      "| p-value =",
      round(test$p.value, 4),
      "|",
      ifelse(test$p.value < 0.05, "*** SIGNIFICANT", "not significant"),
      "\n"
    )

    data.frame(
      compound = cmp,
      p_value = test$p.value,
      direction = direction,
      mean_2021 = mean(g2021, na.rm = TRUE),
      mean_2022 = mean(g2022, na.rm = TRUE)
    )
  }) %>%
    bind_rows()

  sig_labels <- wilcox_results %>%
    mutate(
      sig_label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::select(compound, sig_label, p_value)

  sig_positions <- summary_stats %>%
    group_by(compound) %>%
    summarise(y_pos = max(mean + se) * 1.12, .groups = "drop")

  sig_labels <- left_join(sig_labels, sig_positions, by = "compound")

  return(list(
    summary_stats = summary_stats,
    wilcox_results = wilcox_results,
    sig_labels = sig_labels
  ))
}

# Perform a paired Wilcoxon signed-rank test between years.
# The function compares matched observations from 2021 and 2022,
# pairing samples by their ordering after sorting by sample identity.
# If a column named 'concentration' is present, it is automatically
# renamed to 'Value' for compatibility. The Wilcoxon test statistic
# and associated p-value are returned.
#NOTE HERE YOU NEED 2021 and 2022
#it is the same as before, except here I wanted less detaisl, but idea is the same
#just did not want all other stats
wilcox.short.test = function(data) {
  if (any(names(data) == 'concentration')) {
    data$Value = data$concentration
  }
  g2021 <- data %>%
    filter(Year == 2021) %>%
    arrange(Sample) %>% # should be same
    pull(Value)

  g2022 <- data %>%
    filter(Year == 2022) %>%
    arrange(Sample) %>%
    pull(Value)

  test <- wilcox.test(g2021, g2022, paired = TRUE, exact = FALSE)
  return(test)
}

# extract_climate_windows_by_year <- function(
#   best_clim_mode,
#   calendar_climate,
#   temperature_df
# ) {
#   library(dplyr)
#   library(purrr)
#
#   result <- best_clim_mode %>%
#     mutate(row_id = row_number()) %>%
#     group_split(row_id) %>%
#     map_dfr(function(row) {
#       # No extraction, row is already a 1-row tibble
#       win_open <- 53 #row$WindowOpen
#       win_close <- 37 #row$WindowClose
#       clim_var <- "rainfall" #as.character(row$climate)
#
#       # Filter calendar within window
#       win_dates <- calendar_climate %>%
#         dplyr::select(DOY, days.reversed) %>%
#         dplyr::filter(days.reversed >= win_close, days.reversed <= win_open)
#
#       # Join temperature with calendar to get year
#       temp_merged <- temperature_df %>%
#         inner_join(win_dates, by = "DOY") # Assumes DOY is in both
#
#       # For each year, compute the mean of the desired variable
#       yearly_clim <- temp_merged %>%
#         group_by(year) %>%
#         summarise(
#           mean_value = mean(.data[[clim_var]], na.rm = TRUE),
#           .groups = "drop"
#         ) %>%
#         mutate(
#           response = row$response,
#           climate = clim_var,
#           type = row$type,
#           stat = row$stat,
#           func = row$func,
#           DeltaAICc = row$DeltaAICc,
#           WindowOpen = win_open,
#           WindowClose = win_close
#         )
#
#       yearly_clim
#     })
#
#   return(result)
# }

#go analysis function from Kudoh
# my_GO_analysis = function(gsc, gene_ls, all.Gene, pvalueCutoff = .05) {
#   # parameter for GOstats
#   p <- GSEAGOHyperGParams(
#     name = "Paramaters",
#     geneSetCollection = gsc,
#     geneIds = gene_ls,
#     universeGeneIds = all.Gene,
#     ontology = "BP",
#     pvalueCutoff = pvalueCutoff,
#     conditional = TRUE,
#     testDirection = "over"
#   )
#
#   # Fisher's exact test
#   result <- hyperGTest(p)
#
#   # summarize
#   raw_pval = pvalues(result)
#   adj_pval = p.adjust(raw_pval, method = "BH")
#   odds = oddsRatios(result)
#   expected_count = expectedCounts(result)
#   gene_count = geneCounts(result)
#   gene_univ_count = universeCounts(result)
#
#   summary_table = data.frame(
#     rank = 1:length(raw_pval),
#     go_id = names(raw_pval),
#     raw.pvalues = raw_pval,
#     adjusted.pvalues = adj_pval,
#     oddsRatios = odds,
#     expectedCounts = expected_count,
#     geneCounts = gene_count,
#     universeCounts = gene_univ_count,
#     siginificance = ""
#   )
#
#   summary_table$siginificance[summary_table$adjusted.pvalues < 0.001] = "***"
#   summary_table$siginificance[summary_table$adjusted.pvalues < 0.01] = "**"
#   summary_table$siginificance[
#     (summary_table$adjusted.pvalues >= 0.01) &
#       (summary_table$adjusted.pvalues < 0.05)
#   ] = "*"
#   summary_table$Genes = ""
#
#   target_GOdatabase = GOdatabase[GOdatabase$Orthogroup_ID %in% gene_ls, ]
#
#   for (k in 1:dim(summary_table)[1]) {
#     GOid_k = summary_table$go_id[k]
#     OGs = unique(target_GOdatabase$Orthogroup_ID[
#       target_GOdatabase$go_id == GOid_k
#     ])
#     summary_table$Genes[k] = paste(OGs, collapse = ", ")
#   }
#
#   row.names(summary_table) = 1:dim(summary_table)[1]
#   summary_table = left_join(summary_table, go_termname, by = "go_id")
#
#   return(summary_table)
# }
