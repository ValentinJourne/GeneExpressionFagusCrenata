using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require, character.only = TRUE))
  need <- libs[req == FALSE]
  if (length(need) > 0) {
    install.packages(need)
    lapply(need, require, character.only = TRUE)
  }
}

#to calculate the CV
#i used the yearly CV (and not seasonal ones winthing year)
my_cv_fun <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(NA_real_)
  s / m
}

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

#now a simialr appraoch with sync accrsos indivdidual gene experesision
functioncalc_gene_synchrony = function(
  data,
  method = c("pearson", "spearman")
) {
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


#format climate file
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


cal_z_score <- function(x) {
  (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
}

replace_string <- function(x) {
  if (is.character(x)) {
    x <- ifelse(x == string_to_replace, new_value, x)
  }
  return(x)
}

std.error <- function(x, na.rm = TRUE) sd(x, na.rm = TRUE) / sqrt(length(x))

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

#transform 0-1 data
y.transf.betareg <- function(y) {
  n.obs <- sum(!is.na(y))
  (y * (n.obs - 1) + 0.5) / n.obs
}

#compute the window for log/beta reg analysis
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
    dplyr::select(-data) %>%
    unnest(window_data)
}

#extract informaiton from beta regression of fruit relation to gene expression
#remove column with model inside the name
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
#fit beta regression
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
          ~ betareg(formula.fit.here, data = .x),
          otherwise = NULL
        )
      )
    )
}

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

  # Save output with low qs
  qs::qsave(summary_cumsum_gene, output_file)
  message("Saved batch ", batch_index, " → ", output_file)
}

#generate a calendar for weather (afer climwin analysis)
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

#plot gam effect
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

#fit each gam model to gene
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

#we got more genes on average for low CV but counted more max GO terms with high CV
# make gene set collection for GOstats - Copy from Kudoh et al
#list.genes.cv is the list of genes
#GO.gene.list for GO list
#adapted from Kudoh 2025
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


#do wilcoxon test on GSH and sulfate total S
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
      arrange(Sample) %>% # CRITICAL: same order for pairing
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
