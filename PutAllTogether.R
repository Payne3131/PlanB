# #change file loaction
# 
# # EMS outputs
# 


rm(list = ls())
start.time <- Sys.time()
library('plyr')
source("EMS_Functions.R")

# === SETUP PARAMETERS === #
sample_sizes <- c(100, 50)
noise_sets <- list(
  c(1,1,1,1,1,1),
  c(0.1,0.1,0.1,0.1,0.1,0.1),
  c(1,0.1,0.5,0.4,0.3,0.9)
)
missing_vals <- c(0.05, 0.01, 0)

# Create or reset error log
error_log <- "simulation_errors.log"
if (file.exists(error_log)) file.remove(error_log)

# Make sure the save directory exists
save_dir <- "Sims"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# Function to create unique save names
create_save_name <- function(n, missing, noise, seed) {
  paste0("sim_n", n, "_miss", missing,
         "_noise", paste(noise, collapse = ""),
         "_seed", seed)
}

# === MAIN LOOP === #
required_successes <- 50

for (n in sample_sizes) {
  for (noise in noise_sets) {
    for (missing in missing_vals) {
      
      success_count <- 0
      attempt_seed <- 1
      
      # Create a file name pattern to match existing saved files for this combination
      file_pattern <- paste0("sim_n", n, "_miss", missing, "_noise", noise, "_seed\\d+\\.rda")
      existing_files <- list.files(save_dir, pattern = file_pattern)
      
      # Count how many matching files exist already
      success_count <- length(existing_files)
      print(success_count)
      
      # Skip simulation entirely if we already have enough
      if (success_count >= required_successes) {
        cat("Already have", success_count, "simulations for: n =", n, 
            "missing =", missing, "noise =", noise, "- Skipping.\n")
        next
      }
      
      while (success_count < required_successes && attempt_seed < 300) {
        
        save_name <- create_save_name(n, missing, noise, attempt_seed)
        save_path <- file.path(save_dir, paste0(save_name, ".rda"))  # Adjust extension as needed
        # file.exists("EMS_11_2025/Sims/sim_n100_miss0.05_noise111111_seed1.rda")
        if (file.exists(save_path)) {
          cat("File already exists, skipping simulation for:", save_name, "\n")
          success_count <- success_count + 1
          attempt_seed <- attempt_seed + 1
          print(success_count)
          next
        }
        
        cat("Attempting:", save_name, "\n")
        
        result <- tryCatch({
          generate_EMS_Simulation(
            n,
            chr = 5,
            marker.chr = 6,
            missing,
            noise,
            seed = attempt_seed,
            save = save_name
          )
          print(success_count)
          TRUE  # Success
        }, error = function(e) {
          error_message <- paste("Failed:", save_name, "- Error:", e$message)
          cat(error_message, "\n")
          write(paste(Sys.time(), error_message), file = error_log, append = TRUE)
          FALSE  # Failure
        })
        
        if (result) {
          success_count <- success_count + 1
        }
        
        attempt_seed <- attempt_seed + 1
      }
      
      if (success_count < required_successes) {
        warning(paste("Only", success_count, "successful runs for n =", n,
                      "missing =", missing, "noise =", paste(noise, collapse = ",")))
      }
    }
  }
}

end.time <- Sys.time()
print(paste("Total time:", round(difftime(end.time, start.time, units = "mins"), 2)))

#############Summaries#########################
#############Summaries#########################
#############Summaries#########################
# dir_path <- "EMS_11_2025/Sims"C
library(dplyr)
dir_path <- "Sims"
file_count <- length(list.files(path = dir_path, recursive = FALSE))
print(file_count)

save_names <- c(
  "sim_n100_miss0.05_noise111111",
  "sim_n100_miss0.05_noise0.10.10.10.10.10.1",
  "sim_n100_miss0.05_noise10.10.50.40.30.9",
  "sim_n100_miss0.01_noise111111",
  "sim_n100_miss0.01_noise0.10.10.10.10.10.1",
  "sim_n100_miss0.01_noise10.10.50.40.30.9",
  "sim_n100_miss0_noise111111",
  "sim_n100_miss0_noise0.10.10.10.10.10.1",
  "sim_n100_miss0_noise10.10.50.40.30.9",
  "sim_n50_miss0.05_noise111111",
  "sim_n50_miss0.05_noise0.10.10.10.10.10.1",
  "sim_n50_miss0.05_noise10.10.50.40.30.9",
  "sim_n50_miss0.01_noise111111",
  "sim_n50_miss0.01_noise0.10.10.10.10.10.1",
  "sim_n50_miss0.01_noise10.10.50.40.30.9",
  "sim_n50_miss0_noise111111",
  "sim_n50_miss0_noise0.10.10.10.10.10.1",
  "sim_n50_miss0_noise10.10.50.40.30.9"
)

rda_files <- list.files(dir_path, pattern = "^sim_n.*\\.rda$", full.names = TRUE)
# ???? <- rda_files[grepl(save_names[1], rda_files)]

# # Rename files that contain "miss-1"
# for (file in rda_files) {
#   if (grepl("miss-1_", file)) {
#     new_file <- gsub("miss-1_", "miss0_", file)
#     file.rename(from = file, to = new_file)
#     cat("Renamed:", basename(file), "->", basename(new_file), "\n")
#   }
# }
# 
# # Also update the save_names vector
# save_names <- gsub("miss-1_", "miss0_", save_names)

# Assuming save_names and rda_files are defined as shown
matched_files <- lapply(save_names, function(prefix) {
  rda_files[grepl(prefix, rda_files)]
})

# Optional: name the list elements with the corresponding save_name for clarity
names(matched_files) <- save_names

matched_files[1]


# Define correct values
correct_values <- c(1, 2, 3, 7, 8, 15)

source("EMS_Functions.R")  # Load your simulation function

# # Define the directory where .rda files are stored
dir_path <- "EMS_11_2025/Sims"

### worked 12/16 but overwrites saved plots
# EMS_BIC_Summaries <- function(rdaFile, save_plots = TRUE, plot_dir = "Plots_EMS_BIC") {
#   top_selection <- list()
#   bottom_selection <- list()
#   xlist <- list()
#   Ylist <- list()
#   I.m <- list()
#   FullMiss <- list()
#   
#   library(qtl)
#   library(ggplot2)
#   library(gridExtra)
#   library(dplyr)
#   data(hyper)
#   
#   for (file_path in rdaFile) {
#     tmp_env <- new.env()
#     loaded_objs <- tryCatch(load(file_path, envir = tmp_env), error = function(e) NULL)
#     
#     if (is.null(loaded_objs) || length(loaded_objs) == 0) {
#       warning(paste("Failed to load or no objects in file:", file_path))
#       next
#     }
#     
#     obj <- tmp_env[[loaded_objs[1]]]
#     file_label <- tools::file_path_sans_ext(basename(file_path))
#     
#     if (is.list(obj) && all(c("model.sel", "EMS") %in% names(obj))) {
#       top_selection[[file_label]] <- obj$model.sel
#       bottom_selection[[file_label]] <- obj$EMS
#       xlist[[file_label]] <- obj$xlist
#       Ylist[[file_label]] <- obj$Ylist
#       I.m[[file_label]] <- obj$I.m
#       FullMiss[[file_label]] <- obj$FullMiss
#     } else {
#       warning(paste("model.sel or EMS not found in object from", file_label))
#     }
#   }
#   
#   if (!exists("evaluate_methods")) {
#     stop("Function 'evaluate_methods' is not defined.")
#   }
#   
#   # Containers for evaluation results
#   top_results_df <- data.frame()
#   bottom_results_df <- data.frame()
#   
#   # Loop through all datasets
#   for (key in names(top_selection)) {
#     res <- evaluate_methods(
#       top_selection[[key]],
#       bottom_selection[[key]],
#       correct_values,
#       top_name = "Model Selection",
#       bottom_name = "EMS"
#     )
#     
#     top_results_df <- rbind(top_results_df, res %>% filter(Method == "Model Selection"))
#     bottom_results_df <- rbind(bottom_results_df, res %>% filter(Method == "EMS"))
#   }
#   
#   # Compute summary for top/bottom
#   top_summary <- colMeans(top_results_df[sapply(top_results_df, is.numeric)], na.rm = TRUE)
#   bottom_summary <- colMeans(bottom_results_df[sapply(bottom_results_df, is.numeric)], na.rm = TRUE)
#   
#   final_summary_df <- rbind(
#     `BIC` = top_summary,
#     `EMS` = bottom_summary
#   )
#   
#   
#   # ------------------------------
#   # 📌 FIXED DATA SUMMARY — TRUE PROPORTION MISSING
#   # ------------------------------
#   
#   Missing_Proportion <- sapply(FullMiss, function(fm) {
#     mean(is.na(fm), na.rm = FALSE)
#   })
#   
#   data_summary_df <- data.frame(
#     Mean_Markers              = mean(sapply(xlist, ncol), na.rm = TRUE),
#     Mean_Samples              = mean(sapply(xlist, nrow), na.rm = TRUE),
#     Mean_Missing_Proportion   = mean(Missing_Proportion, na.rm = TRUE),
#     SD_Missing_Proportion     = sd(Missing_Proportion, na.rm = TRUE),
#     Min_Missing_Proportion    = min(Missing_Proportion, na.rm = TRUE),
#     Max_Missing_Proportion    = max(Missing_Proportion, na.rm = TRUE),
#     Mean_Y_Variance           = mean(sapply(Ylist, var), na.rm = TRUE)
#   )
#   
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
#   # Combine results for plotting
#   all_results <- rbind(
#     top_results_df %>% mutate(Method = "BIC"),
#     bottom_results_df %>% mutate(Method = "EMS")
#   )
#   
#   # Add run numbers
#   all_results$Run <- rep(1:nrow(top_results_df), 2)
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
#     if (metric %in% names(all_results)) {
#       p <- ggplot(all_results, aes(x = Method, y = .data[[metric]], fill = Method)) +
#         geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
#         geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
#         theme_bw(base_size = 12) +
#         theme(
#           legend.position = "none",
#           plot.title = element_text(hjust = 0.5, face = "bold")
#         ) +
#         labs(
#           title = gsub("_", " ", metric),
#           x = "Method",
#           y = metric
#         ) +
#         scale_fill_manual(values = c("BIC" = "#66C2A5", "EMS" = "#FC8D62"))
#       
#       plot_list[[metric]] <- p
#       
#       # Save individual plot if requested
#       if (save_plots) {
#         filename <- file.path(plot_dir, paste0("boxplot_", metric, ".png"))
#         ggsave(filename, plot = p, width = 6, height = 5, dpi = 300)
#         cat("  Saved:", filename, "\n")
#       }
#     }
#   }
#   
#   # Create combined plot with all metrics
#   if (length(plot_list) > 0) {
#     combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)
#     
#     # Save combined plot
#     if (save_plots) {
#       filename <- file.path(plot_dir, "boxplot_all_metrics.png")
#       ggsave(filename, plot = combined_plot, width = 15, height = 10, dpi = 300)
#       cat("  Saved combined plot:", filename, "\n")
#     }
#   }
#   
#   # Create comparison plot (BIC vs EMS side by side)
#   comparison_data <- all_results %>%
#     select(Method, Run, all_of(metrics_to_plot[metrics_to_plot %in% names(all_results)])) %>%
#     tidyr::pivot_longer(cols = -c(Method, Run), names_to = "Metric", values_to = "Value")
#   
#   comparison_plot <- ggplot(comparison_data, aes(x = Method, y = Value, fill = Method)) +
#     geom_boxplot(alpha = 0.7) +
#     geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
#     facet_wrap(~Metric, scales = "free_y", ncol = 3) +
#     theme_bw(base_size = 11) +
#     theme(
#       legend.position = "bottom",
#       strip.background = element_rect(fill = "gray90"),
#       strip.text = element_text(face = "bold")
#     ) +
#     labs(
#       title = "BIC vs EMS: Performance Comparison Across All Metrics",
#       x = "Method",
#       y = "Value"
#     ) +
#     scale_fill_manual(values = c("BIC" = "#66C2A5", "EMS" = "#FC8D62"))
#   
#   plot_list[["Comparison"]] <- comparison_plot
#   
#   if (save_plots) {
#     filename <- file.path(plot_dir, "boxplot_comparison.png")
#     ggsave(filename, plot = comparison_plot, width = 14, height = 10, dpi = 300)
#     cat("  Saved comparison plot:", filename, "\n")
#   }
#   
#   # Create data characteristics boxplots
#   data_chars <- data.frame(
#     Dataset = names(Missing_Proportion),
#     Missing_Proportion = Missing_Proportion,
#     N_Markers = sapply(xlist, ncol),
#     N_Samples = sapply(xlist, nrow),
#     Y_Variance = sapply(Ylist, var)
#   )
#   
#   # Missing proportion plot
#   missing_plot <- ggplot(data_chars, aes(x = 1, y = Missing_Proportion)) +
#     geom_boxplot(fill = "#8DA0CB", alpha = 0.7, width = 0.5) +
#     geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
#     theme_bw(base_size = 12) +
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       plot.title = element_text(hjust = 0.5, face = "bold")
#     ) +
#     labs(
#       title = "Distribution of Missing Data Proportion",
#       x = "",
#       y = "Missing Proportion"
#     )
#   
#   plot_list[["Missing_Proportion"]] <- missing_plot
#   
#   # Markers and Samples plot
#   ns_data <- data_chars %>%
#     select(Dataset, N_Markers, N_Samples) %>%
#     tidyr::pivot_longer(cols = c(N_Markers, N_Samples), 
#                         names_to = "Type", values_to = "Count")
#   
#   ns_plot <- ggplot(ns_data, aes(x = Type, y = Count, fill = Type)) +
#     geom_boxplot(alpha = 0.7) +
#     geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
#     theme_bw(base_size = 12) +
#     theme(
#       legend.position = "none",
#       plot.title = element_text(hjust = 0.5, face = "bold")
#     ) +
#     labs(
#       title = "Distribution of Dataset Dimensions",
#       x = "",
#       y = "Count"
#     ) +
#     scale_fill_brewer(palette = "Pastel1") +
#     scale_x_discrete(labels = c("N_Markers" = "Markers", "N_Samples" = "Samples"))
#   
#   plot_list[["Dataset_Dimensions"]] <- ns_plot
#   
#   # Phenotype variance plot
#   var_plot <- ggplot(data_chars, aes(x = 1, y = Y_Variance)) +
#     geom_boxplot(fill = "#FFD92F", alpha = 0.7, width = 0.5) +
#     geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
#     theme_bw(base_size = 12) +
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       plot.title = element_text(hjust = 0.5, face = "bold")
#     ) +
#     labs(
#       title = "Distribution of Phenotype Variance",
#       x = "",
#       y = "Variance"
#     )
#   
#   plot_list[["Phenotype_Variance"]] <- var_plot
#   
#   # Combined data characteristics plot
#   data_char_combined <- grid.arrange(missing_plot, ns_plot, var_plot, ncol = 3)
#   
#   if (save_plots) {
#     filename <- file.path(plot_dir, "boxplot_missing_proportion.png")
#     ggsave(filename, plot = missing_plot, width = 6, height = 5, dpi = 300)
#     
#     filename <- file.path(plot_dir, "boxplot_dataset_dimensions.png")
#     ggsave(filename, plot = ns_plot, width = 6, height = 5, dpi = 300)
#     
#     filename <- file.path(plot_dir, "boxplot_phenotype_variance.png")
#     ggsave(filename, plot = var_plot, width = 6, height = 5, dpi = 300)
#     
#     filename <- file.path(plot_dir, "boxplot_data_characteristics.png")
#     ggsave(filename, plot = data_char_combined, width = 15, height = 5, dpi = 300)
#     cat("  Saved data characteristics plots\n")
#   }
#   
#   cat("\nBoxplot generation complete!\n")
#   
#   # Return both summaries and plots
#   return(list(
#     Method_Summary = final_summary_df,
#     Data_Summary   = data_summary_df,
#     all_results    = all_results,
#     data_characteristics = data_chars,
#     plots          = plot_list  # Return plots for later use
#   ))
# }

EMS_BIC_Summaries <- function(rdaFile, save_plots = TRUE, plot_dir = "Plots_EMS_BIC", 
                              save_name = rdaFile) {
  top_selection <- list()
  bottom_selection <- list()
  xlist <- list()
  Ylist <- list()
  I.m <- list()
  FullMiss <- list()
  
  library(qtl)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  data(hyper)
  
  for (file_path in rdaFile) {
    tmp_env <- new.env()
    loaded_objs <- tryCatch(load(file_path, envir = tmp_env), error = function(e) NULL)
    
    if (is.null(loaded_objs) || length(loaded_objs) == 0) {
      warning(paste("Failed to load or no objects in file:", file_path))
      next
    }
    
    obj <- tmp_env[[loaded_objs[1]]]
    file_label <- tools::file_path_sans_ext(basename(file_path))
    
    if (is.list(obj) && all(c("model.sel", "EMS") %in% names(obj))) {
      top_selection[[file_label]] <- obj$model.sel
      bottom_selection[[file_label]] <- obj$EMS
      xlist[[file_label]] <- obj$xlist
      Ylist[[file_label]] <- obj$Ylist
      I.m[[file_label]] <- obj$I.m
      FullMiss[[file_label]] <- obj$FullMiss
    } else {
      warning(paste("model.sel or EMS not found in object from", file_label))
    }
  }
  
  if (!exists("evaluate_methods")) {
    stop("Function 'evaluate_methods' is not defined.")
  }
  
  # Containers for evaluation results
  top_results_df <- data.frame()
  bottom_results_df <- data.frame()
  
  # Loop through all datasets
  for (key in names(top_selection)) {
    res <- evaluate_methods(
      top_selection[[key]],
      bottom_selection[[key]],
      correct_values,
      top_name = "Model Selection",
      bottom_name = "EMS"
    )
    
    top_results_df <- rbind(top_results_df, res %>% filter(Method == "Model Selection"))
    bottom_results_df <- rbind(bottom_results_df, res %>% filter(Method == "EMS"))
  }
  
  # Compute summary for top/bottom
  top_summary <- colMeans(top_results_df[sapply(top_results_df, is.numeric)], na.rm = TRUE)
  bottom_summary <- colMeans(bottom_results_df[sapply(bottom_results_df, is.numeric)], na.rm = TRUE)
  
  final_summary_df <- rbind(
    `BIC` = top_summary,
    `EMS` = bottom_summary
  )
  
  
  # ------------------------------
  # 📌 FIXED DATA SUMMARY — TRUE PROPORTION MISSING
  # ------------------------------
  
  Missing_Proportion <- sapply(FullMiss, function(fm) {
    mean(is.na(fm), na.rm = FALSE)
  })
  
  data_summary_df <- data.frame(
    Mean_Markers              = mean(sapply(xlist, ncol), na.rm = TRUE),
    Mean_Samples              = mean(sapply(xlist, nrow), na.rm = TRUE),
    Mean_Missing_Proportion   = mean(Missing_Proportion, na.rm = TRUE),
    SD_Missing_Proportion     = sd(Missing_Proportion, na.rm = TRUE),
    Min_Missing_Proportion    = min(Missing_Proportion, na.rm = TRUE),
    Max_Missing_Proportion    = max(Missing_Proportion, na.rm = TRUE),
    Mean_Y_Variance           = mean(sapply(Ylist, var), na.rm = TRUE)
  )
  
  
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
  
  # Combine results for plotting
  all_results <- rbind(
    top_results_df %>% mutate(Method = "BIC"),
    bottom_results_df %>% mutate(Method = "EMS")
  )
  
  # Add run numbers
  all_results$Run <- rep(1:nrow(top_results_df), 2)
  
  # Define metrics to plot
  metrics_to_plot <- c("True_Positives", "False_Positives", "False_Negatives",
                       "Precision", "Recall", "F1_Score")
  
  # Create a list to store all plots
  plot_list <- list()
  
  # Generate boxplot for each metric
  for (metric in metrics_to_plot) {
    if (metric %in% names(all_results)) {
      p <- ggplot(all_results, aes(x = Method, y = .data[[metric]], fill = Method)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
        geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
        theme_bw(base_size = 12) +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")
        ) +
        labs(
          title = gsub("_", " ", metric),
          x = "Method",
          y = metric
        ) +
        scale_fill_manual(values = c("BIC" = "#66C2A5", "EMS" = "#FC8D62"))
      
      plot_list[[metric]] <- p
      
      # Save individual plot if requested
      if (save_plots) {
        filename <- file.path(plot_dir, paste0("boxplot_", metric, file_suffix, ".png"))
        ggsave(filename, plot = p, width = 6, height = 5, dpi = 300)
        cat("  Saved:", filename, "\n")
      }
    }
  }
  
  # Create combined plot with all metrics
  if (length(plot_list) > 0) {
    combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)
    
    # Save combined plot
    if (save_plots) {
      filename <- file.path(plot_dir, paste0("boxplot_all_metrics", file_suffix, ".png"))
      ggsave(filename, plot = combined_plot, width = 15, height = 10, dpi = 300)
      cat("  Saved combined plot:", filename, "\n")
    }
  }
  
  # Create comparison plot (BIC vs EMS side by side)
  comparison_data <- all_results %>%
    select(Method, Run, all_of(metrics_to_plot[metrics_to_plot %in% names(all_results)])) %>%
    tidyr::pivot_longer(cols = -c(Method, Run), names_to = "Metric", values_to = "Value")
  
  comparison_plot <- ggplot(comparison_data, aes(x = Method, y = Value, fill = Method)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    facet_wrap(~Metric, scales = "free_y", ncol = 3) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = "BIC vs EMS: Performance Comparison Across All Metrics",
      x = "Method",
      y = "Value"
    ) +
    scale_fill_manual(values = c("BIC" = "#66C2A5", "EMS" = "#FC8D62"))
  
  plot_list[["Comparison"]] <- comparison_plot
  
  if (save_plots) {
    filename <- file.path(plot_dir, paste0("boxplot_comparison", file_suffix, ".png"))
    ggsave(filename, plot = comparison_plot, width = 14, height = 10, dpi = 300)
    cat("  Saved comparison plot:", filename, "\n")
  }
  
  # Create data characteristics boxplots
  data_chars <- data.frame(
    Dataset = names(Missing_Proportion),
    Missing_Proportion = Missing_Proportion,
    N_Markers = sapply(xlist, ncol),
    N_Samples = sapply(xlist, nrow),
    Y_Variance = sapply(Ylist, var)
  )
  
  # Missing proportion plot
  missing_plot <- ggplot(data_chars, aes(x = 1, y = Missing_Proportion)) +
    geom_boxplot(fill = "#8DA0CB", alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = "Distribution of Missing Data Proportion",
      x = "",
      y = "Missing Proportion"
    )
  
  plot_list[["Missing_Proportion"]] <- missing_plot
  
  # Markers and Samples plot
  ns_data <- data_chars %>%
    select(Dataset, N_Markers, N_Samples) %>%
    tidyr::pivot_longer(cols = c(N_Markers, N_Samples), 
                        names_to = "Type", values_to = "Count")
  
  ns_plot <- ggplot(ns_data, aes(x = Type, y = Count, fill = Type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = "Distribution of Dataset Dimensions",
      x = "",
      y = "Count"
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_x_discrete(labels = c("N_Markers" = "Markers", "N_Samples" = "Samples"))
  
  plot_list[["Dataset_Dimensions"]] <- ns_plot
  
  # Phenotype variance plot
  var_plot <- ggplot(data_chars, aes(x = 1, y = Y_Variance)) +
    geom_boxplot(fill = "#FFD92F", alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = "Distribution of Phenotype Variance",
      x = "",
      y = "Variance"
    )
  
  plot_list[["Phenotype_Variance"]] <- var_plot
  
  # Combined data characteristics plot
  data_char_combined <- grid.arrange(missing_plot, ns_plot, var_plot, ncol = 3)
  
  if (save_plots) {
    filename <- file.path(plot_dir, paste0("boxplot_missing_proportion", file_suffix, ".png"))
    ggsave(filename, plot = missing_plot, width = 6, height = 5, dpi = 300)
    
    filename <- file.path(plot_dir, paste0("boxplot_dataset_dimensions", file_suffix, ".png"))
    ggsave(filename, plot = ns_plot, width = 6, height = 5, dpi = 300)
    
    filename <- file.path(plot_dir, paste0("boxplot_phenotype_variance", file_suffix, ".png"))
    ggsave(filename, plot = var_plot, width = 6, height = 5, dpi = 300)
    
    filename <- file.path(plot_dir, paste0("boxplot_data_characteristics", file_suffix, ".png"))
    ggsave(filename, plot = data_char_combined, width = 15, height = 5, dpi = 300)
    cat("  Saved data characteristics plots\n")
  }
  
  cat("\nBoxplot generation complete!\n")
  
  # Return both summaries and plots
  return(list(
    Method_Summary = final_summary_df,
    Data_Summary   = data_summary_df,
    all_results    = all_results,
    data_characteristics = data_chars,
    plots          = plot_list  # Return plots for later use
  ))
}

source("QTL_From_SIM.R")  # Load your QTL functions
library(QTL)

# Run QTL_Summary for all 18 matched files
# Initialize with named list - use setNames to ensure all names are preserved

QTL_results <- setNames(vector("list", length(matched_files)), names(matched_files))

for (i in 1:length(matched_files)) {
  file_name <- names(matched_files)[i]
  file_path <- matched_files[[i]]
  
  cat("Processing file", i, "of", length(matched_files), ":", file_name, "\n")
  
  tryCatch({
    results <- QTL_Summary(
      rdaFile = file_path,  # Pass as single element (will be wrapped in vector internally)
      correct_values = correct_values
    #   simulate_pheno = FALSE,
    #   qtl_effect_size = 1.0,
    #   error_var = 1,
    #   lod_threshold = 3,
    #   verbose = FALSE  # Set to TRUE if you want detailed output
    )
    
    QTL_results[[i]] <- results  # Use index instead of name
    cat("  Success! Processed", results$n_datasets, "datasets\n")
    
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    QTL_results[[i]] <<- NULL  # Use index instead of name
  })
}




EMS_results <- lapply(matched_files, EMS_BIC_Summaries)
EMS_results

save(QTL_results, file = paste0("Processed/QTL_Summary_Results.rda"))
save(EMS_results, file = paste0("Processed/EMS_Summary_Results.rda"))

out <- rbind(QTL_results$sim_n100_miss0.05_noise111111$plots$True_Positives,
EMS_results$sim_n100_miss0.05_noise111111$plots$True_Positives)
out
#################
#################
#################

load("Processed/EMS_Summary_Results.rda")
load("Processed/QTL_Summary_Results.rda")

############ Testing #########
library(tibble)
library(dplyr)
library(tidyr)

# extract_ems <- function(sim) {
#   m <- EMS_results[[sim]]$Method_Summary
#   df <- as.data.frame(m)
#   df$Method <- rownames(df)
#   df
# }

# head(QTL_results[[names(QTL_results)[1]]]$average_metrics)
# 
# extract_qtl <- function(sim) {
#   QTL_results[[sim]]$average_metrics
# }
# 
# sim <- names(EMS_results)[1]
# 
# ems <- extract_ems(sim)
# qtl <- extract_qtl(sim)
# 
# combined <- bind_rows(ems, qtl)
# 
# combined_long <- combined |>
#   pivot_longer(-Method, names_to="Metric", values_to="Value")

library(ggplot2)

# plot_sim <- function(sim) {
#   
#   ems <- extract_ems(sim)
#   qtl <- extract_qtl(sim)
#   
#   df <- bind_rows(ems, qtl) |>
#     pivot_longer(-Method, names_to="Metric", values_to="Value")
#   
#   ggplot(df, aes(Method, Value, fill=Method)) +
#     geom_col(position="dodge") +
#     facet_wrap(~Metric, ncol=3, scales="free_y") +
#     theme_bw() +
#     labs(title=sim, x=NULL, y=NULL) +
#     theme(legend.position="none",
#           strip.text=element_text(size=11),
#           plot.title=element_text(face="bold"))
# }
# 
# for(sim in names(EMS_results)) {
#   ggsave(
#     paste0("SuppFig_", sim, ".pdf"),
#     plot_sim(sim),
#     width=12,
#     height=8
#   )
# }

# ems <- EMS_results[[sim]]$all_results
# qtl <- QTL_results[[sim]]$all_metrics
# 
# combined <- bind_rows(ems, qtl)

clean_ems <- function(sim) {
  EMS_results[[sim]]$all_results |>
    rename(
      TP = True_Positives,
      FP = False_Positives,
      FN = False_Negatives,
      F1 = F1_Score
    ) |>
    select(Method, TP, FP, FN, Precision, Recall, F1, Run)
}

clean_qtl <- function(sim) {
  QTL_results[[sim]]$all_metrics |>
    rename(
      TP = True_Positives,
      FP = False_Positives,
      FN = False_Negatives,
      F1 = F1_Score
    ) |>
    select(Method, TP, FP, FN, Precision, Recall, F1, Run)
}

# combined <- bind_rows(
#   clean_ems(sim),
#   clean_qtl(sim)
# )

# plot_sim <- function(sim) {
#   
#   df <- bind_rows(
#     clean_ems(sim),
#     clean_qtl(sim)
#   ) |>
#     pivot_longer(c(TP, FP, FN, Precision, Recall, F1),
#                  names_to="Metric", values_to="Value")
#   
#   ggplot(df, aes(Method, Value, fill=Method)) +
#     geom_boxplot(outlier.shape=NA, alpha=.7) +
#     geom_jitter(width=.15, alpha=.4, size=1) +
#     facet_wrap(~Metric, ncol=3, scales="free_y") +
#     theme_bw() +
#     labs(title=sim, x=NULL, y=NULL) +
#     theme(
#       legend.position="none",
#       strip.text = element_text(size=11),
#       plot.title = element_text(face="bold")
#     )
# }

plot_sim <- function(sim) {
  
  # ---- Parse simulation name ----
  parts <- strsplit(sim, "_")[[1]]
  
  n_val     <- gsub("n", "", parts[2])
  miss_val  <- as.numeric(gsub("miss", "", parts[3]))
  noise_raw <- gsub("noise", "", parts[4])
  
  # ---- Convert missing to percent ----
  miss_percent <- paste0(miss_val * 100, "%")
  
  # ---- Map noise patterns to readable labels ----
  noise_label <- dplyr::case_when(
    noise_raw == "111111" ~ "High Noise",
    noise_raw == "0.10.10.10.10.10.1" ~ "Low Noise",
    noise_raw == "10.10.50.40.30.9" ~ "Random Noise",
    TRUE ~ noise_raw
  )
  
  pretty_title <- paste0(
    "Sample Size: ", n_val,
    "  |  Missing: ", miss_percent,
    "  |  ", noise_label
  )
  
  # ---- Build data ----
  df <- bind_rows(
    clean_ems(sim),
    clean_qtl(sim)
  ) |>
    pivot_longer(c(TP, FP, FN, Precision, Recall, F1),
                 names_to="Metric", values_to="Value")
  
  # ---- Plot ----
  ggplot(df, aes(Method, Value, fill=Method)) +
    geom_boxplot(outlier.shape=NA, alpha=.7) +
    geom_jitter(width=.15, alpha=.4, size=1) +
    facet_wrap(~Metric, ncol=3, scales="free_y") +
    theme_bw() +
    labs(title=pretty_title, x=NULL, y=NULL) +
    theme(
      legend.position="none",
      strip.text = element_text(size=11),
      plot.title = element_text(face="bold"),
      axis.text.x = element_text(face = "bold")
    )
}

sim_names <- c()
for(sim in names(EMS_results)) {
  sim_names <- c(sim_names, sim)
  ggsave(
    paste0("FinalPaper/SuppFig_", sim, ".JPG"),
    plot_sim(sim),
    width=12,
    height=8
  )
}
sim_names

# 
# table(combined$Method)
# 
# plot_sim <- function(sim) {
#   
#   ems <- EMS_results[[sim]]$all_results
#   qtl <- QTL_results[[sim]]$all_metrics
#   
#   df <- bind_rows(ems, qtl) |>
#     pivot_longer(c(TP, FP, FN, Precision, Recall, F1),
#                  names_to="Metric", values_to="Value")
#   
#   ggplot(df, aes(Method, Value, fill=Method)) +
#     geom_boxplot(outlier.shape=NA, alpha=.7) +
#     geom_jitter(width=.15, alpha=.4, size=1) +
#     facet_wrap(~Metric, ncol=3, scales="free_y") +
#     theme_bw() +
#     labs(title=sim, x=NULL, y=NULL) +
#     theme(legend.position="none",
#           strip.text=element_text(size=11),
#           plot.title=element_text(face="bold"))
# }
# 
# for(sim in names(EMS_results)) {
#   ggsave(
#     paste0("SuppFig_", sim, ".pdf"),
#     plot_sim(sim),
#     width=12,
#     height=8
#   )
# }

###HEATMAPS###
library(tidyverse)

all_sims <- map_dfr(names(EMS_results), function(sim) {
  
  bind_rows(
    clean_ems(sim),
    clean_qtl(sim)
  ) |>
    mutate(Simulation = sim)
})

heatmap_data <- all_sims |>
  group_by(Simulation, Method) |>
  summarise(mean_F1 = mean(F1, na.rm=TRUE), .groups="drop") |>
  mutate(Simulation = fct_reorder(Simulation, parse_number(Simulation))) |>
  mutate(Simulation =  str_remove(Simulation, "^sim_"))

library(ggplot2)

ggplot(heatmap_data, aes(Method, Simulation, fill=mean_F1)) +
  geom_tile(color="white") +
  scale_fill_viridis_c(name="Mean F1") +
  theme_minimal() +
  labs(
    title="Mean F1 score across 18 simulation conditions",
    x="Method",
    y="Simulation condition"
  ) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

heatmap <- ggplot(heatmap_data, aes(Method, Simulation, fill = mean_F1)) +
  geom_tile(color = "white", linewidth = 0.8) +
  
  # Better color scale: option="magma" or "plasma" has more contrast than default viridis
  scale_fill_viridis_c(
    name = "Mean F1",
    option = "plasma",
    direction = 1,          # light = high F1 (intuitive: bright = good)
    limits = c(0, 1),        # anchor to full F1 range for honest comparison
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = c("0", "0.25", "0.50", "0.75", "1.0")
  ) +
  
  # Add the actual values in each cell
  geom_text(aes(label = round(mean_F1, 2),
                color = mean_F1 > 0.5),   # dark text on light tiles, light on dark
            size = 2.8, fontface = "bold") +
  scale_color_manual(values = c("white", "grey20"), guide = "none") +
  
  # Reorder rows so n100 and n50 are grouped, and missingness increases top-to-bottom
  # (add this upstream with: mutate(Simulation = fct_reorder(Simulation, ...)) if needed)
  
  theme_minimal(base_size = 11) +
  labs(
    title = "Mean F1 Score Across 18 Simulation Conditions",
    subtitle = "Higher = better; anchored to full [0, 1] range",
    x = "Method",
    y = NULL                  # y-axis label adds clutter; title is self-explanatory
  ) +
  theme(
    plot.title          = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.title.position = "plot",   # <-- anchors title to full plot width, not just panel
    plot.subtitle    = element_text(color = "grey40", size = 9, margin = margin(b = 8), hjust = 0.5),
    axis.text.x      = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y      = element_text(size = 8.5),
    panel.grid       = element_blank(),   # grid lines are noise on a heatmap
    legend.position  = "right",
    legend.key.height = unit(1.8, "cm"),  # taller legend bar
    plot.margin      = margin(10, 15, 10, 10)
  )

ggsave("SuppFig_heatmap.pdf", plot = heatmap, 
       width = 14, height = 9, dpi = 600)
#Win loss#
wins <- all_sims |>
  filter(Method %in% c("EMS", "BIC", "QTL EM", "QTL HK", "QTL IMP")) |>
  select(Simulation, Run, Method, F1) |>
  pivot_wider(names_from=Method, values_from=F1)

win_table <- wins |>
  summarise(
    EMS_beats_BIC  = mean(EMS > `BIC`,  na.rm=TRUE),
    EMS_beats_QTL_EM  = mean(EMS > `QTL EM`,  na.rm=TRUE),
    EMS_beats_QTL_HK  = mean(EMS > `QTL HK`,  na.rm=TRUE),
    EMS_beats_QTL_IMP = mean(EMS > `QTL IMP`, na.rm=TRUE)
  )

round(100 * win_table, 1)

# Make sure you have tidyverse loaded
library(dplyr)
library(tidyr)

# Filter to only the methods of interest and pivot wider already done
# wins <- all_sims |> ... (your previous code)

# Create a function to do paired Wilcoxon tests comparing EMS to other methods
wilcox_test_results <- wins %>%
  summarise(
    EMS_vs_BIC   = wilcox.test(EMS, `BIC`, paired = TRUE, alternative = "greater")$p.value,
    EMS_vs_QTL_EM = wilcox.test(EMS, `QTL EM`, paired = TRUE, alternative = "greater")$p.value,
    EMS_vs_QTL_HK = wilcox.test(EMS, `QTL HK`, paired = TRUE, alternative = "greater")$p.value,
    EMS_vs_QTL_IMP = wilcox.test(EMS, `QTL IMP`, paired = TRUE, alternative = "greater")$p.value
  )

# Round the p-values for easier reading
round(wilcox_test_results, 4)

results_table <- bind_cols(
  win_table %>% rename_with(~ paste0("Percent_Win_", .)),
  wilcox_test_results %>% rename_with(~ paste0("pval_", .))
)

round(results_table, 3)

results_table_display <- results_table %>%
  mutate(across(starts_with("pval"), ~ ifelse(. < 1e-4, "<0.0001", round(., 4))),
         across(starts_with("Percent_Win"), ~ round(100 * ., 1)))

results_table_display

effect_sizes <- wins %>%
  summarise(
    median_diff_BIC   = median(EMS - BIC, na.rm=TRUE),
    median_diff_QTL_EM = median(EMS - `QTL EM`, na.rm=TRUE),
    median_diff_QTL_HK = median(EMS - `QTL HK`, na.rm=TRUE),
    median_diff_QTL_IMP = median(EMS - `QTL IMP`, na.rm=TRUE)
  )

variability <- all_sims %>%
  filter(Method %in% c("EMS","BIC","QTL EM","QTL HK","QTL IMP")) %>%
  group_by(Simulation, Method) %>%
  summarise(
    mean_F1 = mean(F1, na.rm=TRUE),
    sd_F1   = sd(F1, na.rm=TRUE),
    IQR_F1  = IQR(F1, na.rm=TRUE),
    .groups = "drop"
  )

ranks <- all_sims %>%
  filter(Method %in% c("EMS","BIC","QTL EM","QTL HK","QTL IMP")) %>%
  group_by(Simulation, Run) %>%
  mutate(Rank = rank(-F1, ties.method = "average")) %>%  # higher F1 = better
  summarise(mean_rank = mean(Rank, na.rm=TRUE), .groups="drop")

# “Across all 18 simulation conditions, EMS achieved higher F1 scores than QTL-EM in 83% of Monte-Carlo runs, QTL-HK in 91% of runs, and QTL-IMP in 90% of runs (Table S1). Mean performance across conditions is visualized in Figure 3.”

#########maybe############

library(dplyr)
library(tidyr)
library(purrr)

# Select the methods of interest
methods_of_interest <- c("EMS", "BIC", "QTL EM", "QTL HK", "QTL IMP")

# Pivot data wider to have one column per method per run
wins <- all_sims %>%
  filter(Method %in% methods_of_interest) %>%
  select(Simulation, Run, Method, F1) %>%
  pivot_wider(names_from = Method, values_from = F1)

# Function to calculate %wins, median delta F1, Wilcoxon p-value per comparison vs EMS
summarize_vs_EMS <- function(df, comparator) {
  tibble(
    Comparison = paste("EMS vs", comparator),
    Percent_Win = mean(df$EMS > df[[comparator]], na.rm=TRUE) * 100,
    Median_Delta_F1 = median(df$EMS - df[[comparator]], na.rm=TRUE),
    Wilcox_p = wilcox.test(df$EMS, df[[comparator]], paired=TRUE, alternative="greater")$p.value
  )
}

# Apply across all comparators
comparisons <- methods_of_interest[methods_of_interest != "EMS"]
summary_table <- map_dfr(comparisons, ~summarize_vs_EMS(wins, .x))

# Optional: format p-values and round numbers for publication
summary_table <- summary_table %>%
  mutate(
    Wilcox_p = ifelse(Wilcox_p < 1e-4, "<0.0001", signif(Wilcox_p, 4)),
    Percent_Win = round(Percent_Win, 1),
    Median_Delta_F1 = round(Median_Delta_F1, 3)
  )

# Optional: add SD and IQR of F1 for EMS and comparator
stats_table <- all_sims %>%
  filter(Method %in% methods_of_interest) %>%
  group_by(Method) %>%
  summarise(
    Mean_F1 = mean(F1, na.rm=TRUE),
    SD_F1 = sd(F1, na.rm=TRUE),
    IQR_F1 = IQR(F1, na.rm=TRUE),
    .groups = "drop"
  )

# Merge summary table with stats if you want
# This is optional depending on how much detail you want
summary_table
stats_table

library(dplyr)
library(tidyr)
library(purrr)

# Methods of interest
methods_of_interest <- c("EMS", "BIC", "QTL EM", "QTL HK", "QTL IMP")

# Pivot F1 data wider per run
wins <- all_sims %>%
  filter(Method %in% methods_of_interest) %>%
  select(Simulation, Run, Method, F1) %>%
  pivot_wider(names_from = Method, values_from = F1)

# Function to calculate summary per comparator vs EMS
summarize_vs_EMS <- function(df, comparator) {
  
  # % Wins and median difference
  percent_win <- mean(df$EMS > df[[comparator]], na.rm=TRUE) * 100
  median_delta <- median(df$EMS - df[[comparator]], na.rm=TRUE)
  
  # Paired Wilcoxon p-value
  wilcox_p <- wilcox.test(df$EMS, df[[comparator]], paired=TRUE, alternative="greater")$p.value
  
  # F1 stats for EMS and comparator
  EMS_stats <- df$EMS
  comparator_stats <- df[[comparator]]
  
  tibble(
    Comparison = paste("EMS vs", comparator),
    Percent_Win = round(percent_win, 1),
    Median_Delta_F1 = round(median_delta, 3),
    Wilcox_p = ifelse(wilcox_p < 1e-4, "<0.0001", signif(wilcox_p, 4)),
    EMS_Mean = round(mean(EMS_stats, na.rm=TRUE), 3),
    EMS_SD = round(sd(EMS_stats, na.rm=TRUE), 3),
    EMS_IQR = round(IQR(EMS_stats, na.rm=TRUE), 3),
    Comparator_Mean = round(mean(comparator_stats, na.rm=TRUE), 3),
    Comparator_SD = round(sd(comparator_stats, na.rm=TRUE), 3),
    Comparator_IQR = round(IQR(comparator_stats, na.rm=TRUE), 3)
  )
}

# Apply across all comparators
comparisons <- methods_of_interest[methods_of_interest != "EMS"]
supp_table <- map_dfr(comparisons, ~summarize_vs_EMS(wins, .x))

# View or save table
supp_table

# Optional: save to CSV for supplement
# write.csv(supp_table, "Supplementary_Table_EMS_vs_Others.csv", row.names = FALSE)



############ Testing #########


save_names <- c(
  "sim_n100_miss0.05_noise111111",
  "sim_n100_miss0.05_noise0.10.10.10.10.10.1",
  "sim_n100_miss0.05_noise10.10.50.40.30.9",
  "sim_n100_miss0.01_noise111111",
  "sim_n100_miss0.01_noise0.10.10.10.10.10.1",
  "sim_n100_miss0.01_noise10.10.50.40.30.9",
  "sim_n100_miss0_noise111111",
  "sim_n100_miss0_noise0.10.10.10.10.10.1",
  "sim_n100_miss0_noise10.10.50.40.30.9",
  "sim_n50_miss0.05_noise111111",
  "sim_n50_miss0.05_noise0.10.10.10.10.10.1",
  "sim_n50_miss0.05_noise10.10.50.40.30.9",
  "sim_n50_miss0.01_noise111111",
  "sim_n50_miss0.01_noise0.10.10.10.10.10.1",
  "sim_n50_miss0.01_noise10.10.50.40.30.9",
  "sim_n50_miss0_noise111111",
  "sim_n50_miss0_noise0.10.10.10.10.10.1",
  "sim_n50_miss0_noise10.10.50.40.30.9"
)

EMS_results

QTL_results

str(QTL_results,max.level = 3)
str(EMS_results,max.level = 3)
QTL_results$sim_n100_miss0.05_noise111111

EMS_results$sim_n100_miss0.05_noise111111$Method_Summary

QTL_results$sim_n100_miss0.05_noise111111$average_metrics

combined <- rbind(
  transform(EMS_results$sim_n100_miss0.05_noise111111$Method_Summary,
            Method = rownames(EMS_results$sim_n100_miss0.05_noise111111$Method_Summary)),
  QTL_results$sim_n100_miss0.05_noise111111$average_metrics
)
rownames(combined) <- NULL

# initialize list to store all combined results
combined_all <- list()

for (sim in save_names) {
  
  # extract EMS summary and add rownames as a column
  ems <- EMS_results[[sim]]$Method_Summary
  ems <- transform(ems, Method = rownames(ems))
  
  # extract QTL summary
  qtl <- QTL_results[[sim]]$average_metrics
  
  # combine
  combined <- rbind(ems, qtl)
  rownames(combined) <- NULL
  
  # store in list using sim name as section name
  combined_all[[sim]] <- combined
}

combined_all$sim_n100_miss0.05_noise111111

# combined_all now contains:
#   $sim_n100_miss0.05_noise111111
#   $sim_n100_miss0.05_noise0.10.10.10.10.10.1
#   $sim_n100_miss0.05_noise10.10.50.40.30.9
#   ...
#   All simulations automatically processed.

library(dplyr)
library(tidyr)

library(tidyverse)
library(stringr)  # Explicitly load stringr
library(kableExtra)



# Combine all sims into one long table
combined_long <- bind_rows(
  lapply(names(combined_all), function(sim) {
    df <- combined_all[[sim]]
    df$Simulation <- sim
    return(df)
  })
)

tidy_table <- combined_long %>%
  pivot_longer(
    cols = -c(Method, Simulation),   # keep Method and Simulation as identifiers
    names_to = "Metric",             # new column for metric names
    values_to = "Value"              # new column for metric values
  )

wide_table <- tidy_table %>%
  pivot_wider(
    names_from = Metric,
    values_from = Value
  )

# Separate the Simulation column into components
wide_table_separated <- wide_table %>%
  separate(
    Simulation, 
    into = c("sim", "sample_size", "missingness", "noise"), 
    sep = "_", 
    remove = FALSE) %>%
  mutate(
    # Clean up the extracted columns
    sample_size = sub("n", "", sample_size),
    sample_size = as.numeric(sample_size),
    missingness = sub("miss", "", missingness),
    missingness = as.numeric(missingness),
    missingness = if_else(missingness == 0, 0.1, missingness),
    noise = sub("noise", "", noise),
    # Categorize noise patterns
    noise = case_when(
      noise == "111111" ~ "Uniform Strong",
      grepl("^0\\.10\\.1", noise) ~ "Uniform Weak",
      TRUE ~ "Random"
    )
   ) %>%
      select(-Simulation, -sim)  # Remove the "sim" column as it's redundant
wide_table_separated

library(ggplot2)

ggplot(wide_table_separated, aes(x = missingness, y = F1_Score, color = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_grid(noise ~ sample_size, labeller = label_both) +
  labs(
    x = "Missingness",
    y = "F1 Score",
    title = "F1 Score vs Missingness",
    color = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "lightgray")
  ) +
  scale_color_brewer(palette = "Set1")

library(ggplot2)
library(tidyr)

# Reshape data to long format
plot_data <- wide_table_separated %>%
  pivot_longer(
    cols = c(True_Positives, False_Positives, False_Negatives, Precision, Recall, F1_Score),
    names_to = "Metric",
    values_to = "Value"
  )

# Create plot for all metrics, faceted by Metric and noise
ggplot(plot_data, aes(x = missingness, y = Value, color = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_grid(Metric ~ noise, scales = "free_y", labeller = label_both) +
  labs(x = "Missingness", y = "Value", title = "All Metrics vs Missingness by Noise Level") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "lightgray")
  ) +
  scale_color_brewer(palette = "Set1")+

facet_grid(Metric ~ noise + sample_size, scales = "free_y", labeller = label_both)

# Create a fancy table
wide_table_separated %>%
  select(Method, sample_size, missingness, noise, 
         True_Positives, False_Positives, False_Negatives,
         Precision, Recall, F1_Score) %>%
  kbl(col.names = c("Method", "Sample Size", "Missingness", "Noise",
                    "True Positives", "False Positives", "False Negatives",
                    "Precision", "Recall", "F1 Score"),
      caption = "Performance Metrics Across Simulation Conditions",
      digits = 3,
      align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE,
                position = "center",
                font_size = 12) %>%
  column_spec(1, bold = TRUE) %>%
  row_spec(0, bold = TRUE, background = "#f0f0f0") %>%
  collapse_rows(columns = 2:4, valign = "top")

# Get unique sample sizes
sample_sizes <- sort(unique(wide_table_separated$sample_size))

# Create list to store tables
table_5meth_F1Score <- list()

# Create separate tables for each sample size

n <- 50
table <- wide_table_separated %>%
  filter(sample_size == n) %>%
  select(Method, missingness, noise, 
         True_Positives, False_Positives, False_Negatives,
         Precision, Recall, F1_Score) %>%
  kbl(col.names = c("Method", "Missingness", "Noise",
                    "True Positives", "False Positives", "False Negatives",
                    "Precision", "Recall", "F1 Score"),
      # caption = paste0("Performance Metrics for Sample Size n = ", n),
      digits = 3,
      align = c("l", "c", "c", "c", "c", "c", "c", "c", "c")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE,
                position = "center",
                font_size = 12) %>%
  column_spec(1, bold = TRUE) %>%
  row_spec(0, bold = TRUE, background = "#f0f0f0") %>%
  collapse_rows(columns = 1:3, valign = "top")

# Save to list with descriptive name
table_name <- paste0("table_n", n)
table_5meth_F1Score[[table_name]] <- table

# Print table
print(table)
cat("\n\n")

# Save as HTML file
save_kable(table, file = paste0("FinalPaper/table_", n, ".html"))

n <- 100
table <- wide_table_separated %>%
  filter(sample_size == n) %>%
  select(Method, missingness, noise, 
         True_Positives, False_Positives, False_Negatives,
         Precision, Recall, F1_Score) %>%
  kbl(col.names = c("Method", "Missingness", "Noise",
                    "True Positives", "False Positives", "False Negatives",
                    "Precision", "Recall", "F1 Score"),
      # caption = paste0("Performance Metrics for Sample Size n = ", n),
      digits = 3,
      align = c("l", "c", "c", "c", "c", "c", "c", "c", "c")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE,
                position = "center",
                font_size = 12) %>%
  column_spec(1, bold = TRUE) %>%
  row_spec(0, bold = TRUE, background = "#f0f0f0") %>%
  collapse_rows(columns = 1:3, valign = "top")

# Save to list with descriptive name
table_name <- paste0("table_n", n)
table_5meth_F1Score[[table_name]] <- table

# Print table
print(table)
cat("\n\n")

# Save as HTML file
save_kable(table, file = paste0("FinalPaper/table_", n, ".html"))



# Option 1: Overall Performance Comparison
# "Across all simulation conditions, EMS consistently demonstrated superior performance with the highest F1 scores (ranging from 0.90-0.93), achieving near-perfect recall while maintaining high precision (0.87-0.88). Model Selection showed moderate performance with high recall (0.64-1.0) but lower precision (0.32-0.66) due to elevated false positive rates. The QTL methods (EM, HK, and IMP) exhibited perfect precision but suffered from very low recall (0.15-0.24), indicating a conservative approach that missed many true positives."
# Option 2: Method-Specific Narrative
# "We compared five methods across varying sample sizes and noise conditions. EMS emerged as the most balanced approach, correctly identifying 94-100% of true positives while keeping false positives under 1.0 on average. In contrast, Model Selection demonstrated higher sensitivity but lower specificity, with false positive rates ranging from 3.7 to 8.6. The three QTL approaches showed nearly identical performance, maintaining perfect precision at the cost of detecting only 15-24% of true signals."
# Option 3: Condition-Specific Findings
# "Performance varied substantially across noise conditions. Under uniform strong noise (111111), all methods performed better, with EMS achieving F1 scores of 0.93. However, under uniform weak noise (0.10.10.10.10.10.1), method performance diverged more dramatically. Sample size also impacted results, with [describe patterns you observe across n=100, 200, etc.]."
# Option 4: Precision-Recall Tradeoff Focus
# "Our results reveal a clear precision-recall tradeoff among methods. EMS achieved optimal balance (precision: 0.87-0.88, recall: 0.95-1.0). Model Selection favored recall over precision, while QTL methods prioritized precision, accepting lower recall rates. This suggests researchers should select methods based on whether minimizing false positives or false negatives is more critical for their application."

table_5meth_F1Score$table_n50
table_5meth_F1Score$table_n100

# Alternative: Create a summary table (averaged by Method and Noise Type)
summary_table <- wide_table_separated %>%
  group_by(Method, noise) %>%
  summarise(
    n = n(),
    TP = mean(True_Positives),
    FP = mean(False_Positives),
    FN = mean(False_Negatives),
    Precision = mean(Precision),
    Recall = mean(Recall),
    F1 = mean(F1_Score),
    .groups = "drop"
  )

summary_table %>%
  kbl(col.names = c("Method", "Noise", "N", "True Pos.", "False Pos.", 
                    "False Neg.", "Precision", "Recall", "F1 Score"),
      caption = "Average Performance Metrics by Method and Noise Type",
      digits = 3,
      align = c("l", "l", "c", "c", "c", "c", "c", "c", "c")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE,
                position = "center",
                font_size = 12) %>%
  column_spec(1, bold = TRUE, width = "3cm") %>%
  column_spec(2, width = "3cm") %>%
  row_spec(0, bold = TRUE, background = "#e8e8e8") %>%
  collapse_rows(columns = 1, valign = "middle") %>%
  add_header_above(c(" " = 3, "Counts" = 3, "Metrics" = 3))
