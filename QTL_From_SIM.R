# ============================================================================
# PART 1: evaluate_QTL Function with Performance Metrics
# ============================================================================

evaluate_QTL <- function(xlist, Ylist, Miss, correct_values, qtl_indices = NULL,
                         output_dir = "./Processed", hyper_template, 
                         run, lod_threshold = 3, tolerance = 5, verbose = FALSE,
                         simulate_pheno = FALSE, qtl_effect_size = 1, error_var = 1) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  X <- unlist(xlist)
  y <- unlist(Ylist)
  I.m <- unlist(Miss)

  # Handle cases where indices or dimensions may not match
  if (is.null(X) || is.null(y) || is.null(I.m)) {
    warning(paste("Missing data for run", run))
    return(NULL)
  }

  out6 <<- X
  out7 <<- I.m
  
  zz <- X
  zz[I.m,] <- NA
  zz[zz == 1] <- 2
  zz[zz == 0] <- 1
  zzy <- y

  # Option to simulate phenotype with QTL effects
  if (simulate_pheno && !is.null(qtl_indices)) {
    if (verbose) cat("Simulating phenotype with QTL effects...\n")

    # Simulate phenotype based on genotypes
    pheno_sim <- simulate_phenotype(
      xlist = X,
      qtl_indices = qtl_indices,
      qtl_effects = rep(qtl_effect_size, length(qtl_indices)),
      error_var = error_var
    )

    # Create phenotype vector (one value per individual)
    simulated_phenotype <- pheno_sim$phenotype

    if (verbose) {
      cat("  Heritability:", round(pheno_sim$heritability, 3), "\n")
      cat("  Phenotype range:", round(range(simulated_phenotype), 2), "\n")
    }
  } else {
    simulated_phenotype <- NULL
  }

  # Assuming a structure with 30 columns
  breaks <- list(1:6, 7:12, 13:18, 19:24, 25:30)
  result_list_all <- list()

  for (i in 1:5) {
    data_matrix <- matrix(zz[, breaks[[i]]], nrow = nrow(zz), ncol = length(breaks[[i]]))
    colnames <- paste0("Val", breaks[[i]])
    dimnames(data_matrix) <- list(1:nrow(zz), colnames)
    map_vector <- zzy[breaks[[i]]]
    names(map_vector) <- colnames
    result_list <- list(data = data_matrix, map = map_vector)
    class(result_list) <- "A"
    result_list_all[[i]] <- result_list
  }
  
  # Use provided EMSData template
  EMSData <- hyper_template
  for (i in 1:5) {
    EMSData$geno[[as.character(i)]] <- result_list_all[[i]]
  }
  EMSData$geno <- EMSData$geno[1:5]
  
  # REPLACE entire pheno with simulated data if requested, otherwise use zzy
  if (!is.null(simulated_phenotype)) {
    # Get original pheno structure info
    original_colnames <- colnames(EMSData$pheno)
    original_rownames <- rownames(EMSData$pheno)
    
    # Create new pheno data frame with simulated phenotype
    n_indiv <- length(simulated_phenotype)
    EMSData$pheno <- as.data.frame(simulated_phenotype)
    
    # Preserve column name if there was one
    if (!is.null(original_colnames) && length(original_colnames) > 0) {
      colnames(EMSData$pheno) <- original_colnames[1]
    } else {
      colnames(EMSData$pheno) <- "phenotype"
    }
    
    # Preserve row names if they existed and match length
    if (!is.null(original_rownames) && length(original_rownames) >= n_indiv) {
      rownames(EMSData$pheno) <- original_rownames[1:n_indiv]
    }
    
    if (verbose) {
      cat("  Completely replaced pheno with simulated data\n")
      cat("  New pheno dimensions:", nrow(EMSData$pheno), "x", ncol(EMSData$pheno), "\n")
    }
  } else {
    # Replace pheno with zzy data
    original_colnames <- colnames(EMSData$pheno)
    original_rownames <- rownames(EMSData$pheno)
    
    n_indiv <- length(zzy)
    EMSData$pheno <- as.data.frame(zzy)
    
    # Preserve column name if there was one
    if (!is.null(original_colnames) && length(original_colnames) > 0) {
      colnames(EMSData$pheno) <- original_colnames[1]
    } else {
      colnames(EMSData$pheno) <- "phenotype"
    }
    
    # Preserve row names if they existed and match length
    if (!is.null(original_rownames) && length(original_rownames) >= n_indiv) {
      rownames(EMSData$pheno) <- original_rownames[1:n_indiv]
    }
  }
  
  out1 <<- EMSData
  out2 <<- EMSData
  out3 <<- zzy
  out4 <<- EMSData$pheno

  something <- EMSData
  
  # Drop markers with all missing data
  something <- drop.nullmarkers(something)
  
  # Calculate genotype probabilities
  something <- calc.genoprob(something, step = 1, error.prob = 0.01)
  
  # Perform genome scans (without permutations for QTL detection)
  out.em <- scanone(something)
  out.hk <- scanone(something, method = "hk")
  
  # Simulate genotypes and run multiple imputation scan
  something <- sim.geno(something, step = 2, n.draws = 16, error.prob = 0.01)
  out.imp <- scanone(something, method = "imp")
  
  # Get max LOD scores to help set threshold
  max_lod_em <- max(out.em[,3], na.rm = TRUE)
  max_lod_hk <- max(out.hk[,3], na.rm = TRUE)
  max_lod_imp <- max(out.imp[,3], na.rm = TRUE)
  
  if (verbose) {
    cat("Run", run, "- Max LOD scores: EM=", round(max_lod_em, 2), 
        "HK=", round(max_lod_hk, 2), "IMP=", round(max_lod_imp, 2), "\n")
  }
  
  # Extract detected QTL positions for each method
  detected_em <- summary(out.em, threshold = lod_threshold, format = "allpeaks")
  detected_hk <- summary(out.hk, threshold = lod_threshold, format = "allpeaks")
  detected_imp <- summary(out.imp, threshold = lod_threshold, format = "allpeaks")
  
  # This code extracts detected Quantitative Trait Loci (QTL) from three different mapping approaches: 
  # expectation-maximization (EM), Haley-Knott (HK) regression, and multiple imputation methods. 
  # The summary() function applies a LOD score threshold to each method's scan results to identify statistically significant QTL peaks, 
  # with the format = "allpeaks" argument ensuring that all significant peaks across the genome are reported for each approach.
  # This process generates three separate datasets containing the genomic positions and statistical metrics of detected QTL,
  # facilitating direct comparison of mapping performance and detection consistency across the EM, HK, and multiple imputation methods.
  
  # Ensure all are data frames (handle different qtl package return types)
  convert_to_df <- function(x) {
    if (is.null(x) || length(x) == 0) {
      return(data.frame(chr = character(0), pos = numeric(0), lod = numeric(0)))
    }
    if (is.data.frame(x)) {
      return(x)
    }
    # Handle case where it's a named vector (single QTL)
    if (is.numeric(x) && !is.null(names(x))) {
      return(data.frame(
        chr = as.character(x["chr"]),
        pos = as.numeric(x["pos"]),
        lod = as.numeric(x["lod"]),
        stringsAsFactors = FALSE
      ))
    }
    # If nothing matched, return empty
    return(data.frame(chr = character(0), pos = numeric(0), lod = numeric(0)))
  }
  
  detected_em <- convert_to_df(detected_em)
  detected_hk <- convert_to_df(detected_hk)
  detected_imp <- convert_to_df(detected_imp)
  
  if (verbose && nrow(detected_em) > 0) {
    cat("Run", run, "- Detected", nrow(detected_em), "QTL(s) with EM at positions:", 
        paste(round(detected_em$pos, 1), collapse = ", "), "\n")
  }
  
  # Get correct QTL positions from correct_values
  correct_counts <- table(unlist(correct_values))
  correct_positions <- as.numeric(names(correct_counts))
  
  # Function to calculate performance metrics
  calculate_metrics <- function(detected_qtl, correct_pos, tol = tolerance) {
    # Extract detected positions - handle different return types
    if (is.null(detected_qtl)) {
      detected_pos <- numeric(0)
    } else if (is.data.frame(detected_qtl)) {
      # Data frame with multiple QTL detected
      if (nrow(detected_qtl) == 0) {
        detected_pos <- numeric(0)
      } else {
        detected_pos <- detected_qtl$pos
      }
    } else if (is.numeric(detected_qtl) && length(detected_qtl) > 0) {
      # Single QTL returned as named vector
      detected_pos <- as.numeric(detected_qtl["pos"])
    } else {
      detected_pos <- numeric(0)
    }
    
    # Calculate True Positives, False Positives, False Negatives
    TP <- 0
    matched_correct <- rep(FALSE, length(correct_pos))
    matched_detected <- rep(FALSE, length(detected_pos))
    
    # For each detected QTL, check if it's close to a correct position
    for (i in seq_along(detected_pos)) {
      for (j in seq_along(correct_pos)) {
        if (!matched_correct[j] && abs(detected_pos[i] - correct_pos[j]) <= tol) {
          TP <- TP + 1
          matched_correct[j] <- TRUE
          matched_detected[i] <- TRUE
          break
        }
      }
    }
    
    FP <- sum(!matched_detected)
    FN <- sum(!matched_correct)
    
    # Calculate Precision, Recall, F1
    Precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
    Recall <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
    F1_Score <- ifelse(Precision + Recall > 0, 
                       2 * (Precision * Recall) / (Precision + Recall), 0)
    
    return(list(
      True_Positives = TP,
      False_Positives = FP,
      False_Negatives = FN,
      Precision = Precision,
      Recall = Recall,
      F1_Score = F1_Score
    ))
  }
  
  # Calculate metrics for each method
  metrics_em <- calculate_metrics(detected_em, correct_positions)
  metrics_hk <- calculate_metrics(detected_hk, correct_positions)
  metrics_imp <- calculate_metrics(detected_imp, correct_positions)
  
  # Create summary data frame
  metrics_df <- data.frame(
    Method = c("QTL EM", "QTL HK", "QTL IMP"),
    True_Positives = c(metrics_em$True_Positives, 
                       metrics_hk$True_Positives, 
                       metrics_imp$True_Positives),
    False_Positives = c(metrics_em$False_Positives, 
                        metrics_hk$False_Positives, 
                        metrics_imp$False_Positives),
    False_Negatives = c(metrics_em$False_Negatives, 
                        metrics_hk$False_Negatives, 
                        metrics_imp$False_Negatives),
    Precision = c(metrics_em$Precision, 
                  metrics_hk$Precision, 
                  metrics_imp$Precision),
    Recall = c(metrics_em$Recall, 
               metrics_hk$Recall, 
               metrics_imp$Recall),
    F1_Score = c(metrics_em$F1_Score, 
                 metrics_hk$F1_Score, 
                 metrics_imp$F1_Score)
  )
  
  # Compile results
  results <- list(
    run = run,
    metrics_summary = metrics_df,
    detected_qtl = list(
      em = detected_em,
      hk = detected_hk,
      imp = detected_imp
    ),
    max_lod_scores = c(EM = max_lod_em, HK = max_lod_hk, IMP = max_lod_imp),
    correct_positions = correct_positions,
    lod_threshold_used = lod_threshold
  )
  
  
  return(results)
}


# ============================================================================
# PART 2: QTL_Summary Function
# ============================================================================

# QTL_Summary <- function(rdaFile, correct_values = NULL, use_model_sel = TRUE,
#                         lod_threshold = 3, tolerance = 5, verbose = FALSE,
#                         simulate_pheno = FALSE, qtl_effect_size = 1, error_var = 1,
#                         save_plots = TRUE, plot_dir = "Plots_QTL") {
#   library(qtl)
#   library(ggplot2)
#   library(gridExtra)
#   data(hyper)
# 
#   xlist <- list()
#   Ylist <- list()
#   I.m <- list()
#   FullMiss <- list()
#   correct_values_list <- list()
#   qtl_indices_list <- list()
# 
#   cat("Loading", length(rdaFile), "RDA files...\n")
# 
#   # Load all files
#   for (i in seq_along(rdaFile)) {
#     file_path <- rdaFile[i]
# 
#     tmp_env <- new.env()
#     loaded_objs <- tryCatch(load(file_path, envir = tmp_env), error = function(e) {
#       warning(paste("Error loading", file_path, ":", e$message))
#       return(NULL)
#     })
# 
#     if (is.null(loaded_objs) || length(loaded_objs) == 0) {
#       warning(paste("Failed to load file:", file_path))
#       next
#     }
# 
#     obj <- tmp_env[[loaded_objs[1]]]
# 
#     # Extract data from AllResults object
#     if (is.list(obj) && "xlist" %in% names(obj)) {
#     xlist[[i]] <- obj$xlist
#     Ylist[[i]] <- obj$Ylist
#     I.m[[i]] <- obj$I.m
#     FullMiss[[i]] <- obj$FullMiss
#       # Extract correct QTL positions for this dataset
#       if (is.null(correct_values)) {
#         if (use_model_sel && !is.null(obj$model.sel)) {
#           correct_values_list[[i]] <- obj$Ylist[obj$model.sel]
#           qtl_indices_list[[i]] <- obj$model.sel
#         } else if (!is.null(obj$EMS) && is.numeric(obj$EMS)) {
#           correct_values_list[[i]] <- obj$Ylist[obj$EMS]
#           qtl_indices_list[[i]] <- obj$EMS
#         }
#       } else {
#         correct_values_list[[i]] <- correct_values
#         qtl_indices_list[[i]] <- NULL
#       }
#     } else {
#       warning(paste("Could not extract data from file", i))
#     }
#   }
# 
#   matched_files[1]
#   # Remove NULL entries
#   valid_indices <- !sapply(xlist, is.null)
#   xlist <- xlist[valid_indices]
#   Ylist <- Ylist[valid_indices]
#   I.m <- I.m[valid_indices]
#   FullMiss <- FullMiss[valid_indices]
#   correct_values_list <- correct_values_list[valid_indices]
#   qtl_indices_list <- qtl_indices_list[valid_indices]
# 
#   n_datasets <- length(xlist)
#   cat("\nValid datasets to process:", n_datasets, "\n")
# 
#   if (n_datasets == 0) {
#     stop("No valid datasets found!")
#   }
# 
#   # Show example of true QTL positions
#   if (length(correct_values_list) > 0 && !is.null(correct_values_list[[1]])) {
#     cat("Example true QTL positions (dataset 1):",
#         paste(round(correct_values_list[[1]], 3), collapse = ", "), "\n\n")
#   }
# 
#   # Process each dataset
#   QTL_df <- vector("list", n_datasets)
# 
#   for (i in 1:n_datasets) {
#     # cat("Evaluating dataset", i, "of", n_datasets, "...\n")
# 
#     tryCatch({
#       # miss_data <- if (!is.null(FullMiss[[i]])) FullMiss[[i]] else I.m[[i]]
#       # Handle missing data properly
#       if (!is.null(FullMiss[[i]]) && length(FullMiss[[i]]) == length(xlist[[i]])) {
#         miss_data <- FullMiss[[i]]
#       } else if (!is.null(I.m[[i]]) && length(I.m[[i]]) == length(xlist[[i]])) {
#         miss_data <- I.m[[i]]
#       } else {
#         # No missing data - create a FALSE matrix matching xlist dimensions
#         cat("  Creating empty missing data indicator (no missingness)\n")
#         miss_data <- lapply(xlist[[i]], function(x) matrix(FALSE, nrow = nrow(x), ncol = ncol(x)))
#       }
# 
#       # Debug: Check data dimensions
#       cat("  xlist length:", length(xlist[[i]]), "\n")
#       cat("  Ylist length:", length(Ylist[[i]]), "\n")
#       cat("  Miss length:", length(miss_data), "\n")
# 
#       QTL_df[[i]] <- evaluate_QTL(
#         xlist = xlist[[i]],
#         Ylist = Ylist[[i]],
#         Miss = miss_data,
#         correct_values = correct_values_list[[i]],
#         qtl_indices = qtl_indices_list[[i]],
#         output_dir = "Processed",
#         hyper_template = hyper,
#         run = i,
#         lod_threshold = lod_threshold,
#         tolerance = tolerance,
#         verbose = verbose,
#         simulate_pheno = simulate_pheno,
#         qtl_effect_size = qtl_effect_size,
#         error_var = error_var
#       )
#       cat("  Done\n")
#     }, error = function(e) {
#       cat("  ERROR:", e$message, "\n")
#       cat("  Traceback:\n")
#       print(traceback())
#       QTL_df[[i]] <- NULL
#     })
#   }
# 
#   # Remove failed evaluations
#   QTL_df <- QTL_df[!sapply(QTL_df, is.null)]
# 
#   if (length(QTL_df) == 0) {
#     stop("No datasets were successfully evaluated.")
#   }
# 
#   cat("\nSuccessfully evaluated", length(QTL_df), "datasets\n")
# 
#   # Compile all metrics
#   all_metrics <- do.call(rbind, lapply(seq_along(QTL_df), function(i) {
#     if (!is.null(QTL_df[[i]]$metrics_summary)) {
#       df <- QTL_df[[i]]$metrics_summary
#       df$Run <- i
#       return(df)
#     }
#     return(NULL)
#   }))
# 
#   # Calculate average metrics by method
#   metrics_by_method <- aggregate(
#     cbind(True_Positives, False_Positives, False_Negatives,
#           Precision, Recall, F1_Score) ~ Method,
#     data = all_metrics,
#     FUN = mean,
#     na.rm = TRUE
#   )
# 
#   # Calculate standard deviations
#   metrics_sd <- aggregate(
#     cbind(True_Positives, False_Positives, False_Negatives,
#           Precision, Recall, F1_Score) ~ Method,
#     data = all_metrics,
#     FUN = sd,
#     na.rm = TRUE
#   )
#   names(metrics_sd)[-1] <- paste0(names(metrics_sd)[-1], "_SD")
# 
#   metrics_summary <- merge(metrics_by_method, metrics_sd, by = "Method")
# 
#   cat("\n=== Average Performance Across All Runs ===\n")
#   print(metrics_by_method, digits = 3)
# 
#   # Report max LOD scores
#   all_max_lods <- sapply(QTL_df, function(x) x$max_lod_scores)
#   cat("\n=== Average Max LOD Scores ===\n")
#   cat("EM:  ", round(mean(all_max_lods["EM",]), 3), "\n")
#   cat("HK:  ", round(mean(all_max_lods["HK",]), 3), "\n")
#   cat("IMP: ", round(mean(all_max_lods["IMP",]), 3), "\n")
#   cat("\nLOD threshold used:", lod_threshold, "\n")
#   cat("Position tolerance used:", tolerance, "\n")
# 
#   # ============================================================================
#   # CREATE BOXPLOTS
#   # ============================================================================
#   
#   cat("\n=== Generating Boxplots ===\n")
#   
#   # Create plot directory if saving
#   if (save_plots && !dir.exists(plot_dir)) {
#     dir.create(plot_dir, recursive = TRUE)
#     cat("Created plot directory:", plot_dir, "\n")
#   }
#   
#   # Define metrics to plot
#   metrics_to_plot <- c("True_Positives", "False_Positives", "False_Negatives",
#                        "Precision", "Recall", "F1_Score")
#   
#   # Create a list to store all plots
#   plot_list <- list()
#   
#   # Generate boxplot for each metric
#   for (metric in metrics_to_plot) {
#     p <- ggplot(all_metrics, aes(x = Method, y = .data[[metric]], fill = Method)) +
#       geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
#       geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
#       theme_bw(base_size = 12) +
#       theme(
#         legend.position = "none",
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         plot.title = element_text(hjust = 0.5, face = "bold")
#       ) +
#       labs(
#         title = gsub("_", " ", metric),
#         x = "Method",
#         y = metric
#       ) +
#       scale_fill_brewer(palette = "Set2")
#     
#     plot_list[[metric]] <- p
#     
#     # Save individual plot if requested
#     if (save_plots) {
#       filename <- file.path(plot_dir, paste0("boxplot_", metric, ".png"))
#       ggsave(filename, plot = p, width = 6, height = 5, dpi = 300)
#       cat("  Saved:", filename, "\n")
#     }
#   }
#   
#   # Create combined plot with all metrics
#   combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)
#   
#   # Save combined plot
#   if (save_plots) {
#     filename <- file.path(plot_dir, "boxplot_all_metrics.png")
#     ggsave(filename, plot = combined_plot, width = 15, height = 10, dpi = 300)
#     cat("  Saved combined plot:", filename, "\n")
#   }
#   
#   # Create LOD scores boxplot
#   lod_data <- data.frame(
#     Method = rep(rownames(all_max_lods), ncol(all_max_lods)),
#     LOD_Score = as.vector(all_max_lods)
#   )
#   
#   lod_plot <- ggplot(lod_data, aes(x = Method, y = LOD_Score, fill = Method)) +
#     geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
#     geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
#     geom_hline(yintercept = lod_threshold, linetype = "dashed", 
#                color = "red", linewidth = 1) +
#     annotate("text", x = 0.7, y = lod_threshold + 0.5, 
#              label = paste("Threshold =", lod_threshold), color = "red") +
#     theme_bw(base_size = 12) +
#     theme(
#       legend.position = "none",
#       plot.title = element_text(hjust = 0.5, face = "bold")
#     ) +
#     labs(
#       title = "Maximum LOD Scores by Method",
#       x = "Method",
#       y = "Max LOD Score"
#     ) +
#     scale_fill_brewer(palette = "Set2")
#   
#   plot_list[["LOD_Scores"]] <- lod_plot
#   
#   if (save_plots) {
#     filename <- file.path(plot_dir, "boxplot_LOD_scores.png")
#     ggsave(filename, plot = lod_plot, width = 6, height = 5, dpi = 300)
#     cat("  Saved:", filename, "\n")
#   }
#   
#   cat("\nBoxplot generation complete!\n")
#   
#   return(list(
#     # individual_results = QTL_df,
#     all_metrics = all_metrics,
#     average_metrics = metrics_by_method,
#     metrics_with_sd = metrics_summary,
#     plots = plot_list  # Return plots for later use
#     # n_datasets = length(QTL_df)
#   ))
# }

QTL_Summary <- function(rdaFile, correct_values = NULL, use_model_sel = TRUE,
                        lod_threshold = 3, tolerance = 5, verbose = FALSE,
                        simulate_pheno = FALSE, qtl_effect_size = 1, error_var = 1,
                        save_plots = TRUE, plot_dir = "Plots", save_name = rdaFile) {
  library(qtl)
  library(ggplot2)
  library(gridExtra)
  data(hyper)
  
  xlist <- list()
  Ylist <- list()
  I.m <- list()
  FullMiss <- list()
  correct_values_list <- list()
  qtl_indices_list <- list()
  
  cat("Loading", length(rdaFile), "RDA files...\n")
  
  # Load all files
  for (i in seq_along(rdaFile)) {
    file_path <- rdaFile[i]
    
    tmp_env <- new.env()
    loaded_objs <- tryCatch(load(file_path, envir = tmp_env), error = function(e) {
      warning(paste("Error loading", file_path, ":", e$message))
      return(NULL)
    })
    
    if (is.null(loaded_objs) || length(loaded_objs) == 0) {
      warning(paste("Failed to load file:", file_path))
      next
    }
    
    obj <- tmp_env[[loaded_objs[1]]]
    
    # Extract data from AllResults object
    if (is.list(obj) && "xlist" %in% names(obj)) {
      xlist[[i]] <- obj$xlist
      Ylist[[i]] <- obj$Ylist
      I.m[[i]] <- obj$I.m
      FullMiss[[i]] <- obj$FullMiss
      # Extract correct QTL positions for this dataset
      if (is.null(correct_values)) {
        if (use_model_sel && !is.null(obj$model.sel)) {
          correct_values_list[[i]] <- obj$Ylist[obj$model.sel]
          qtl_indices_list[[i]] <- obj$model.sel
        } else if (!is.null(obj$EMS) && is.numeric(obj$EMS)) {
          correct_values_list[[i]] <- obj$Ylist[obj$EMS]
          qtl_indices_list[[i]] <- obj$EMS
        }
      } else {
        correct_values_list[[i]] <- correct_values
        qtl_indices_list[[i]] <- NULL
      }
    } else {
      warning(paste("Could not extract data from file", i))
    }
  }
  
  matched_files[1]
  # Remove NULL entries
  valid_indices <- !sapply(xlist, is.null)
  xlist <- xlist[valid_indices]
  Ylist <- Ylist[valid_indices]
  I.m <- I.m[valid_indices]
  FullMiss <- FullMiss[valid_indices]
  correct_values_list <- correct_values_list[valid_indices]
  qtl_indices_list <- qtl_indices_list[valid_indices]
  
  n_datasets <- length(xlist)
  cat("\nValid datasets to process:", n_datasets, "\n")
  
  if (n_datasets == 0) {
    stop("No valid datasets found!")
  }
  
  # Show example of true QTL positions
  if (length(correct_values_list) > 0 && !is.null(correct_values_list[[1]])) {
    cat("Example true QTL positions (dataset 1):",
        paste(round(correct_values_list[[1]], 3), collapse = ", "), "\n\n")
  }
  
  # Process each dataset
  QTL_df <- vector("list", n_datasets)
  
  for (i in 1:n_datasets) {
    # cat("Evaluating dataset", i, "of", n_datasets, "...\n")
    
    tryCatch({
      # miss_data <- if (!is.null(FullMiss[[i]])) FullMiss[[i]] else I.m[[i]]
      # Handle missing data properly
      if (!is.null(FullMiss[[i]]) && length(FullMiss[[i]]) == length(xlist[[i]])) {
        miss_data <- FullMiss[[i]]
      } else if (!is.null(I.m[[i]]) && length(I.m[[i]]) == length(xlist[[i]])) {
        miss_data <- I.m[[i]]
      } else {
        # No missing data - create a FALSE matrix matching xlist dimensions
        cat("  Creating empty missing data indicator (no missingness)\n")
        miss_data <- lapply(xlist[[i]], function(x) matrix(FALSE, nrow = nrow(x), ncol = ncol(x)))
      }
      
      # Debug: Check data dimensions
      cat("  xlist length:", length(xlist[[i]]), "\n")
      cat("  Ylist length:", length(Ylist[[i]]), "\n")
      cat("  Miss length:", length(miss_data), "\n")
      
      QTL_df[[i]] <- evaluate_QTL(
        xlist = xlist[[i]],
        Ylist = Ylist[[i]],
        Miss = miss_data,
        correct_values = correct_values_list[[i]],
        qtl_indices = qtl_indices_list[[i]],
        output_dir = "Processed",
        hyper_template = hyper,
        run = i,
        lod_threshold = lod_threshold,
        tolerance = tolerance,
        verbose = verbose,
        simulate_pheno = simulate_pheno,
        qtl_effect_size = qtl_effect_size,
        error_var = error_var
      )
      cat("  Done\n")
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
      cat("  Traceback:\n")
      print(traceback())
      QTL_df[[i]] <- NULL
    })
  }
  
  # Remove failed evaluations
  QTL_df <- QTL_df[!sapply(QTL_df, is.null)]
  
  if (length(QTL_df) == 0) {
    stop("No datasets were successfully evaluated.")
  }
  
  cat("\nSuccessfully evaluated", length(QTL_df), "datasets\n")
  
  # Compile all metrics
  all_metrics <- do.call(rbind, lapply(seq_along(QTL_df), function(i) {
    if (!is.null(QTL_df[[i]]$metrics_summary)) {
      df <- QTL_df[[i]]$metrics_summary
      df$Run <- i
      df$Dataset <- tools::file_path_sans_ext(basename(rdaFile[i]))
      return(df)
    }
    return(NULL)
  }))
  
  # Calculate average metrics by method
  metrics_by_method <- aggregate(
    cbind(True_Positives, False_Positives, False_Negatives,
          Precision, Recall, F1_Score) ~ Method,
    data = all_metrics,
    FUN = mean,
    na.rm = TRUE
  )
  
  # Calculate standard deviations
  metrics_sd <- aggregate(
    cbind(True_Positives, False_Positives, False_Negatives,
          Precision, Recall, F1_Score) ~ Method,
    data = all_metrics,
    FUN = sd,
    na.rm = TRUE
  )
  names(metrics_sd)[-1] <- paste0(names(metrics_sd)[-1], "_SD")
  
  metrics_summary <- merge(metrics_by_method, metrics_sd, by = "Method")
  
  cat("\n=== Average Performance Across All Runs ===\n")
  print(metrics_by_method, digits = 3)
  
  # Report max LOD scores
  all_max_lods <- sapply(QTL_df, function(x) x$max_lod_scores)
  cat("\n=== Average Max LOD Scores ===\n")
  cat("EM:  ", round(mean(all_max_lods["EM",]), 3), "\n")
  cat("HK:  ", round(mean(all_max_lods["HK",]), 3), "\n")
  cat("IMP: ", round(mean(all_max_lods["IMP",]), 3), "\n")
  cat("\nLOD threshold used:", lod_threshold, "\n")
  cat("Position tolerance used:", tolerance, "\n")
  
  # ============================================================================
  # CREATE BOXPLOTS
  # ============================================================================
  
  cat("\n=== Generating Boxplots ===\n")
  
  # Create plot directory if saving
  if (save_plots && !dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    cat("Created plot directory:", plot_dir, "\n")
  }
  
  # Use save_name for file suffix if provided
  file_suffix <- if (!is.null(save_name)) paste0("_", save_name) else ""
  
  # Define metrics to plot
  metrics_to_plot <- c("True_Positives", "False_Positives", "False_Negatives",
                       "Precision", "Recall", "F1_Score")
  
  # Create a list to store all plots
  plot_list <- list()
  
  # Generate boxplot for each metric
  for (metric in metrics_to_plot) {
    p <- ggplot(all_metrics, aes(x = Method, y = .data[[metric]], fill = Method)) +
      geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
      geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      labs(
        title = gsub("_", " ", metric),
        x = "Method",
        y = metric
      ) +
      scale_fill_brewer(palette = "Set2")
    
    plot_list[[metric]] <- p
    
    # Save individual plot if requested
    if (save_plots) {
      filename <- file.path(plot_dir, paste0("boxplot_", metric, file_suffix, ".png"))
      ggsave(filename, plot = p, width = 6, height = 5, dpi = 300)
      cat("  Saved:", filename, "\n")
    }
  }
  
  # Create combined plot with all metrics
  combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)
  
  # Save combined plot
  if (save_plots) {
    filename <- file.path(plot_dir, paste0("boxplot_all_metrics", file_suffix, ".png"))
    ggsave(filename, plot = combined_plot, width = 15, height = 10, dpi = 300)
    cat("  Saved combined plot:", filename, "\n")
  }
  
  # Create LOD scores boxplot
  lod_data <- data.frame(
    Method = rep(rownames(all_max_lods), ncol(all_max_lods)),
    LOD_Score = as.vector(all_max_lods)
  )
  
  lod_plot <- ggplot(lod_data, aes(x = Method, y = LOD_Score, fill = Method)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
    geom_hline(yintercept = lod_threshold, linetype = "dashed", 
               color = "red", linewidth = 1) +
    annotate("text", x = 0.7, y = lod_threshold + 0.5, 
             label = paste("Threshold =", lod_threshold), color = "red") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = "Maximum LOD Scores by Method",
      x = "Method",
      y = "Max LOD Score"
    ) +
    scale_fill_brewer(palette = "Set2")
  
  plot_list[["LOD_Scores"]] <- lod_plot
  
  if (save_plots) {
    filename <- file.path(plot_dir, paste0("boxplot_LOD_scores", file_suffix, ".png"))
    ggsave(filename, plot = lod_plot, width = 6, height = 5, dpi = 300)
    cat("  Saved:", filename, "\n")
  }
  
  cat("\nBoxplot generation complete!\n")
  
  return(list(
    # individual_results = QTL_df,
    all_metrics = all_metrics,
    average_metrics = metrics_by_method,
    metrics_with_sd = metrics_summary,
    plots = plot_list  # Return plots for later use
    # n_datasets = length(QTL_df)
  ))
}


# ============================================================================
# USAGE EXAMPLE
# ============================================================================

# Quick diagnostic function
check_data_structure <- function(rdaFile) {
  cat("Checking first file:", rdaFile[1], "\n\n")
  
  tmp_env <- new.env()
  load(rdaFile[1], envir = tmp_env)
  obj <- tmp_env[[ls(tmp_env)[1]]]
  
  cat("Object class:", class(obj), "\n")
  cat("Object names:", paste(names(obj), collapse = ", "), "\n\n")
  
  if ("xlist" %in% names(obj)) {
    cat("xlist structure:\n")
    cat("  Type:", class(obj$xlist), "\n")
    cat("  Length:", length(obj$xlist), "\n")
    if (length(obj$xlist) > 0) {
      cat("  First element type:", class(obj$xlist[[1]]), "\n")
      if (is.matrix(obj$xlist[[1]]) || is.data.frame(obj$xlist[[1]])) {
        cat("  First element dims:", paste(dim(obj$xlist[[1]]), collapse = " x "), "\n")
      }
    }
    cat("\n")
  }
  
  if ("Ylist" %in% names(obj)) {
    cat("Ylist structure:\n")
    cat("  Type:", class(obj$Ylist), "\n")
    cat("  Length:", length(obj$Ylist), "\n")
    cat("  Values (first 10):", paste(head(obj$Ylist, 10), collapse = ", "), "\n")
    cat("  Range:", paste(range(obj$Ylist), collapse = " to "), "\n")
    cat("\n")
  }
  
  if ("FullMiss" %in% names(obj)) {
    cat("FullMiss structure:\n")
    cat("  Type:", class(obj$FullMiss), "\n")
    cat("  Length:", length(obj$FullMiss), "\n")
    cat("\n")
  }
  
  # Check if there are actual QTL positions stored
  if ("model.sel" %in% names(obj)) {
    cat("model.sel found - structure:\n")
    print(str(obj$model.sel, max.level = 2))
  }
  
  if ("EMS" %in% names(obj)) {
    cat("\nEMS found - structure:\n")
    print(str(obj$EMS, max.level = 2))
  }
  
  return(obj)
}

# Function to extract true QTL positions from your simulation
extract_true_qtls <- function(obj) {
  # Check different possible locations where true QTLs might be stored
  true_qtls <- NULL
  
  if (!is.null(obj$model.sel)) {
    cat("Checking model.sel for QTL positions...\n")
    print(obj$model.sel)
  }
  
  if (!is.null(obj$EMS)) {
    cat("\nChecking EMS for QTL positions...\n")
    if (is.list(obj$EMS)) {
      print(names(obj$EMS))
      if ("beta" %in% names(obj$EMS)) {
        cat("Beta coefficients:\n")
        print(obj$EMS$beta)
      }
      if ("selected" %in% names(obj$EMS)) {
        cat("Selected markers:\n")
        print(obj$EMS$selected)
      }
    } else {
      print(obj$EMS)
    }
  }
  
  # The Ylist contains genetic map positions
  cat("\nGenetic map positions (Ylist):\n")
  cat("Total markers:", length(obj$Ylist), "\n")
  
  # Convert marker indices to actual genetic map positions
  if (!is.null(obj$model.sel)) {
    true_positions <- obj$Ylist[obj$model.sel]
    cat("\nTrue QTL marker indices:", paste(obj$model.sel, collapse = ", "), "\n")
    cat("True QTL genetic positions:", paste(round(true_positions, 3), collapse = ", "), "\n")
  }
  
  if (!is.null(obj$EMS) && is.numeric(obj$EMS)) {
    ems_positions <- obj$Ylist[obj$EMS]
    cat("\nEMS marker indices:", paste(obj$EMS, collapse = ", "), "\n")
    cat("EMS genetic positions:", paste(round(ems_positions, 3), collapse = ", "), "\n")
  }
  
  return(list(
    model_sel_indices = obj$model.sel,
    model_sel_positions = if(!is.null(obj$model.sel)) obj$Ylist[obj$model.sel] else NULL,
    EMS_indices = obj$EMS,
    EMS_positions = if(!is.null(obj$EMS) && is.numeric(obj$EMS)) obj$Ylist[obj$EMS] else NULL,
    all_map_positions = obj$Ylist
  ))
}

# Function to simulate phenotype with QTL effects
simulate_phenotype <- function(xlist, qtl_indices, qtl_effects = NULL, 
                               error_var = 1, baseline = 0) {
  # xlist: genotype matrix (individuals x markers)
  # qtl_indices: which markers are causal QTLs
  # qtl_effects: effect sizes for each QTL (defaults to equal effects)
  # error_var: residual variance
  
  X <- as.matrix(xlist)
  n <- nrow(X)
  
  # Default effect sizes
  if (is.null(qtl_effects)) {
    qtl_effects <- rep(1, length(qtl_indices))
  }
  
  # Calculate genetic effects
  genetic_value <- baseline
  for (i in seq_along(qtl_indices)) {
    # Recode genotypes: 0,1 -> -1,1 for additive effects
    geno <- X[, qtl_indices[i]]
    geno_centered <- ifelse(geno == 1, 1, -1)
    genetic_value <- genetic_value + qtl_effects[i] * geno_centered
  }
  
  # Add noise
  phenotype <- genetic_value + rnorm(n, 0, sqrt(error_var))
  
  return(list(
    phenotype = phenotype,
    genetic_value = genetic_value,
    heritability = var(genetic_value) / var(phenotype)
  ))
}

# Helper function to get correct values from all files
get_all_correct_values <- function(rdaFiles, use_model_sel = TRUE) {
  cat("Extracting true QTL positions from", length(rdaFiles), "files...\n")
  
  all_positions <- list()
  
  for (i in seq_along(rdaFiles)) {
    tmp_env <- new.env()
    load(rdaFiles[i], envir = tmp_env)
    obj <- tmp_env[[ls(tmp_env)[1]]]
    
    if (use_model_sel && !is.null(obj$model.sel)) {
      # Use model.sel indices
      all_positions[[i]] <- obj$Ylist[obj$model.sel]
    } else if (!is.null(obj$EMS) && is.numeric(obj$EMS)) {
      # Use EMS indices
      all_positions[[i]] <- obj$Ylist[obj$EMS]
    }
  }
  
  cat("Extracted positions for", length(all_positions), "datasets\n")
  return(all_positions)
}

# Define your correct QTL positions
# correct_values <- c(5, 10, 15, 20, 25)  # Example positions

# STEP 1: First check your data structure
# obj <- check_data_structure(matched_files$sim_n100_miss0.05_noise111111)

# STEP 2: Then run the summary
# results <- QTL_Summary(
#   rdaFile = matched_files$sim_n100_miss0.05_noise111111,
#   correct_values = correct_values
# )

# Access results:
# results$average_metrics          # Mean performance by method
# results$metrics_with_sd          # With standard deviations
# results$all_metrics              # All individual runs
# results$individual_results[[1]]  # First run details

# # Run with simulated phenotype
# results <- QTL_Summary(
#   rdaFile = matched_files$sim_n100_miss0.05_noise111111,
#   use_model_sel = TRUE,        # Use model.sel as true QTLs
#   simulate_pheno = TRUE,        # Simulate phenotype with QTL effects
#   qtl_effect_size = 0.5,        # Effect size for each QTL
#   error_var = 1,                # Residual variance
#   lod_threshold = 2,            # LOD threshold
#   tolerance = 0.5,              # Position tolerance
#   verbose = FALSE
# )

# results <- QTL_Summary(
#   rdaFile = matched_files$sim_n100_miss-1_noise111111,
#   use_model_sel = TRUE,        # Use model.sel as true QTLs
#   simulate_pheno = TRUE,        # Simulate phenotype with QTL effects
#   qtl_effect_size = 0.5,        # Effect size for each QTL
#   error_var = 1,                # Residual variance
#   lod_threshold = 2,            # LOD threshold
#   tolerance = 0.5,              # Position tolerance
#   verbose = FALSE
# )

# # View results
# results$average_metrics
# 
# # Try different effect sizes to see performance
# results_strong <- QTL_Summary(
#   rdaFile = matched_files$sim_n100_miss0.05_noise111111,
#   simulate_pheno = TRUE,
#   qtl_effect_size = 1.0,        # Stronger effects
#   error_var = 1,
#   lod_threshold = 3,
#   verbose = FALSE
# )

# # Try different effect sizes to see performance
# results_weak <- QTL_Summary(
#   rdaFile = matched_files$sim_n100_miss0_noise111111,
#   simulate_pheno = FALSE,
#   qtl_effect_size = 1.0,        # Stronger effects
#   error_var = 1,
#   lod_threshold = 3,
#   verbose = FALSE
# )
# results_weak <- QTL_Summary(
#   rdaFile = matched_files$sim_n100_miss0.01_noise111111,
#   correct_values = correct_values
# )
# results_weak2 <- QTL_Summary(
#   rdaFile = matched_files$sim_n100_miss0_noise111111,
#   correct_values = correct_values
# )

# results_strong$average_metrics
# results_weak$average_metrics
# 
# plot(out1)
# plot(out2)
# hist(out1$pheno$bp)
# hist(out2$pheno$bp)

# # # Run QTL_Summary for all 18 matched files
# # # Initialize with named list - use setNames to ensure all names are preserved
# # all_results <- setNames(vector("list", length(matched_files)), names(matched_files))
# # 
# # for (i in 1:1) {
# #   # file_name <- names(matched_files)[i]
# #   file_path <- matched_files[[i]]
# # }
# # 
# # for (i in 1:length(matched_files)) {
# #   file_name <- names(matched_files)[i]
# #   file_path <- matched_files[[i]]
# # 
# #   cat("Processing file", i, "of", length(matched_files), ":", file_name, "\n")
# # 
# #   tryCatch({
# #     results <- QTL_Summary(
# #       rdaFile = file_path,  # Pass as single element (will be wrapped in vector internally)
# #       simulate_pheno = FALSE,
# #       qtl_effect_size = 1.0,
# #       error_var = 1,
# #       lod_threshold = 3,
# #       verbose = FALSE  # Set to TRUE if you want detailed output
# #     )
# # 
# #     all_results[[i]] <- results  # Use index instead of name
# #     cat("  Success! Processed", results$n_datasets, "datasets\n")
# # 
# #   }, error = function(e) {
# #     cat("  Error:", e$message, "\n")
# #     all_results[[i]] <<- NULL  # Use index instead of name
# #   })
# #   save <- "QTL_Summary_Results"
# #   save(all_results, file = paste0("Processed/", save, ".rda"))
# # }
# qwe

# # Summary of results
# cat("\n=== Summary ===\n")
# cat("Total files to process:", length(matched_files), "\n")
# 
# # Count successes and failures by checking each element
# n_success <- 0
# n_failed <- 0
# n_missing <- 0
# 
# for (i in 1:length(matched_files)) {
#   if (i <= length(all_results)) {
#     if (!is.null(all_results[[i]])) {
#       n_success <- n_success + 1
#     } else {
#       n_failed <- n_failed + 1
#     }
#   } else {
#     n_missing <- n_missing + 1
#   }
# }
# 
# cat("Successful:", n_success, "\n")
# cat("Failed (with errors):", n_failed, "\n")
# cat("Not processed:", n_missing, "\n")
# 
# # View which files succeeded/failed
# status_vec <- character(length(matched_files))
# datasets_processed <- integer(length(matched_files))
# 
# for (i in 1:length(matched_files)) {
#   if (i <= length(all_results) && !is.null(all_results[[i]])) {
#     status_vec[i] <- "Success"
#     datasets_processed[i] <- all_results[[i]]$n_datasets
#   } else if (i <= length(all_results) && is.null(all_results[[i]])) {
#     status_vec[i] <- "Failed"
#     datasets_processed[i] <- 0
#   } else {
#     status_vec[i] <- "Not Processed"
#     datasets_processed[i] <- 0
#   }
# }
# 
# success_status <- data.frame(
#   file = names(matched_files),
#   status = status_vec,
#   n_datasets = datasets_processed,
#   stringsAsFactors = FALSE
# )
# print(success_status)
# 
# cat("\n=== Failed Files ===\n")
# print(success_status[success_status$status != "Success", ])
# files<-success_status[success_status$status != "Success", ][1]
# 
# # Get failed files
# # Get failed files
# # Get failed files
# cat("\n=== Failed Files ===\n")
# failed_info <- success_status[success_status$status != "Success", ]
# print(failed_info)
# 
# # Extract the file names
# failed_file_names <- failed_info$file
# 
# cat("\n=== Re-running", length(failed_file_names), "failed files ===\n")
# 
# # Create a new results list for just these files
# failed_results <- setNames(vector("list", length(failed_file_names)), failed_file_names)
# 
# for (i in 1:length(failed_file_names)) {
#   file_name <- failed_file_names[i]
#   file_path <- matched_files[[file_name]]
# 
#   cat("\nProcessing", i, "of", length(failed_file_names), ":", file_name, "\n")
# 
#   tryCatch({
#     results <- QTL_Summary(
#       rdaFile = file_path,
#       simulate_pheno = TRUE,
#       qtl_effect_size = 1.0,
#       error_var = 1,
#       lod_threshold = 3,
#       verbose = TRUE  # Turn on verbose to see what's happening
#     )
# 
#     failed_results[[file_name]] <- results
#     cat("  SUCCESS! Processed", results$n_datasets, "datasets\n")
# 
#   }, error = function(e) {
#     cat("  ERROR:", e$message, "\n")
#     failed_results[[file_name]] <<- NULL
#   })
# }
# 
# # Summary of re-run
# cat("\n=== Re-run Summary ===\n")
# cat("Files attempted:", length(failed_file_names), "\n")
# cat("Now successful:", sum(!sapply(failed_results, is.null)), "\n")
# cat("Still failing:", sum(sapply(failed_results, is.null)), "\n")
