
source("EMS/EMS_Functions.R")
library(MASS)
library('plyr')
library(abind)

generate_matrix <- function(from, to, rows, columns, seed) {
  # Generate values from 1 to 100
  set.seed(seed)
  values <- from:to
  
  # Repeat the values as needed to fill the 10x30 matrix
  values_repeated <- rep(values, length.out = rows * columns)
  random_values <- sample(from:to, rows * columns, replace = TRUE)
  
  # Create the matrix with 10 rows and 30 columns
  matrix_data <- matrix(random_values, nrow = rows, ncol = columns, byrow = TRUE)
  
  return(matrix_data)
}

generate_missing <- function(from, to, size, seed) {
  # Generate values from 1 to 100
  set.seed(seed)
  values <- from:to
  
  # # Repeat the values as needed to fill the 10x30 matrix
  # values_repeated <- rep(values, length.out = rows * columns)
  # random_values <- sample(from:to, rows * columns, replace = TRUE)
  
  # Create the matrix with 10 rows and 30 columns
  # matrix_data <- matrix(random_values, nrow = rows, ncol = columns, byrow = TRUE)
  
  indices <- sample(1:length(values), size = size)
  
  # Split data into 10% and 90% portions
  data_10 <- values[indices]  # 10% portion
  data_90 <- values[-indices]  # 90% portion
  if (size == 0) {
    data_10<- "null"
    data_90<- values
  }
  
  return(list(data_10, data_90))
}

# Function to generate a 10x30 matrix with values from 1 to 100
seed = 1
set.seed(seed)

# original data for EMS example Sunil Rao
n = 100
chr = 5
marker.chr = 6
p = chr*marker.chr # total dim of regression model
m.ind = 1:p # 
subj.ind = 1:n # subject indicator
coef = 10#.1 to 3 works for the seed 3
size = .05



# Call the function to generate and print the matrix
MM <- generate_matrix(1, n, 10, p, seed)

out <- generate_missing(1, n, size = size, seed)
I.m <- out[[1]]
I.o <- out[[2]]

X.index.1 <- Generate_X_Index(n,seed)
save(X.index, file="EMS/X.index_notwork.rda")

# load("EMS/M.n.100.rda")
# load("EMS/I.m.n.100.rda")
# load("EMS/I.o.n.100.rda")
# # load("EMS/S.ind.n.100.rda") #dont need to run?
# load("EMS/X.index.n.100.rda")
# # load("EMS/EMS.test.05.23.24.n.100.rda") # resukts

# I.mo <- c(I.m, I.o)


# len.m = length(0)
# len.o = length(I.mo) #checked to see if I needed this to be two parts or one
c.m =(marker.chr - 1)*chr

# I.index = numeric(n)


theta = .2


beta = numeric(p)
beta[1] = coef
beta[2] = coef
beta[3] = coef
beta[7] = coef
beta[8] = coef
beta[15] = coef

index = 1:p

start.time <- Sys.time()
result = vector("list", 1)
r=1
set.seed(seed)
X = t(sapply(rep(chr,n),inheritance,theta = theta, marker.chr = marker.chr))

len.m = length(I.m) # number of missing
len.o = length(I.o) # number of observed
I.index = numeric(n)
I.index[I.m] = 1

epsilon = 0.1*rnorm(n, mean = 0, sd =1)
y = X%*%beta + epsilon
index.c = index
count.prep = matrix(NA, nrow = c.m, ncol =2)
iter = 0
for(k in 1:chr){
  for(j in 1:(marker.chr -1)){
    iter = iter +1
    index.j = (k - 1)*marker.chr + j
    index.j.plus = index.j +1
    count.prep[iter,] = get.prep.theta.hat(index.j, index.j.plus)
  }
}
N.count = apply(count.prep,2,sum)
theta.pre = N.count[2]/N.count[1]
gamma.est = theta.pre/(1-theta.pre)	

eps = 1
iteration = 1	
while(eps > 10^(-4)& iteration <= 100){
  theta.i = rep(NA,n)
  for(i in 1:n){
    theta.i[i] = get.theta.comp(X.index[[i]],X[i,],gamma.est)
  }
  theta.hat = sum(theta.i)/(n*chr*(marker.chr-1))
  eps = abs(theta.hat - theta.pre)
  theta.pre = theta.hat
  gamma.est = theta.pre/(1-theta.pre)	
  iteration = iteration + 1
}



y.o = y[I.o] #Rows with no missing
X.sub = X[I.o,] #Rows with no missing
# if (X.sub)
beta.c = solve(t(X.sub)%*%X.sub)%*%t(X.sub)%*%y.o
y.o.hat = X.sub%*%beta.c
var.c = as.numeric(t(y.o - y.o.hat)%*%(y.o - y.o.hat))/(len.o - p)
check = 1
while(check != 0){
  model.dic = get.model.dic(index.c, gamma.est,var.c,beta.c)
  model.dic = unique(lapply(model.dic,sort))
  len.model = length(model.dic)
  
  BIC.v = sapply(model.dic,get.BIC,index.c=index.c,X.index=X.index,x.design=X,gamma.est=gamma.est,y=y,var.est=var.c,beta.est=beta.c)
  model.sel = model.dic[[which(BIC.v == min(BIC.v))]]
  var.c = get.var(index.c,model.sel,X.index,X,gamma.est,y,var.c,beta.c)
  beta.c = get.beta.hat(index.c,model.sel,X.index,X,gamma.est,y,var.c,beta.c)

  if(length(model.sel) == length(index.c)){
    if(sum(model.sel %in% index.c) == length(index.c)){
      check = 0
    }else{
      index.c = model.sel
      check = check + 1
    }
  }else{
    index.c = model.sel
    check = check + 1
  }
}

result.1 = model.sel
result.2 = FB.selection(p,X,y,29) #is the EMS method
result[[r]] = list(result.1,result.2)

result.1
result.2
save(result, file="EMS/Sim_EMS_Seed_5_CoefChange")
five <- as.character(seed)
save(result, file = paste0("EMS/Sim_EMS_", five, "_CoefChange2"))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

rm(list = ls())
start.time1 <- Sys.time()
source("EMS/EMS_Functions.R")
generate_EMS_Simulation(100, 5, 6, .05, c(1,1,1,1,1,1), seed =1, save = "testing_6_10")
generate_EMS_Simulation(50, 5, 6, .05, c(1,1,1,1,1,1), seed = 1, save = "testing50_noise.05")

test1 <- generate_EMS_Simulation(100, 5, 6, .05, c(1,1,1,1,1,1), seed = 1, save = "testing")
test2 <- generate_EMS_Simulation(100, 5, 6, -1, c(1,1,1,1,1,1), seed = 1, save = "testing2")

test1.1 <- generate_EMS_Simulation(100, 5, 6, -1, c(.1,.1,.1,.1,.1,.1), seed = 1, save = "testing.1")
test2.1 <- generate_EMS_Simulation(100, 5, 6, -1, c(3,3,3,3,3,3), seed = 1, save = "testing2.1")

test1.2 <- generate_EMS_Simulation(200, 5, 6, -1, c(1,1,1,1,1,1), seed = 1, save = "testing.2")
test2.2 <- generate_EMS_Simulation(50, 5, 6, -1, c(1,1,1,1,1,1), seed = 1, save = "testing2.2")

test1.3 <- generate_EMS_Simulation(200, 5, 6, -1, c(.5,.1,.9,1,.3,1), seed = 1, save = "testing.3")
test2.3 <- generate_EMS_Simulation(50, 5, 6, .01, c(.5,.5,.5,.5,.5,.5), seed = 1, save = "testing2.3")
end.time1 <- Sys.time()
time.taken <- end.time1 - start.time1
time.taken

start.time11 <- Sys.time()
fake1 <- generate_EMS_Simulation(100, 5, 6, .05, c(.5,.5,.5,.5,.5,.5), seed = 1, save = "fake1")
end.time11 <- Sys.time()

start.time12 <- Sys.time()
fake2 <- tryCatch({
  generate_EMS_Simulation(100, 5, 6, -1, c(1,1,1,1,1,1), seed = 1, save = "fake2")
}, error = function(e) {
  message("Error in fake2: ", e$message)
  NULL
})
end.time12 <- Sys.time()

start.time13 <- Sys.time()#when using .05 for missing, the coef will only work for 0:1, but using the orginal missing methood it will work for 0:1?
fake3 <- tryCatch({
  generate_EMS_Simulation(100, 5, 6, .05, c(1,1,1,1,1,1), seed = 1, save = "fake3")
}, error = function(e) {
  message("Error in fake3: ", e$message)
  NULL
})
end.time13 <- Sys.time()

end.time11 - start.time11
end.time12 - start.time12
end.time13 - start.time13

# Define sets
top_selection <- test1$model.sel
bottom_selection <- test1$EMS
correct_values <- c(1, 2, 3, 7, 8, 15)

# Precision, Recall, F1 for Top
tp_top <- length(intersect(top_selection, correct_values))
fp_top <- length(setdiff(top_selection, correct_values))
fn_top <- length(setdiff(correct_values, top_selection))
precision_top <- tp_top / (tp_top + fp_top)
recall_top <- tp_top / (tp_top + fn_top)
f1_top <- 2 * (precision_top * recall_top) / (precision_top + recall_top)

# Precision, Recall, F1 for Bottom
tp_bottom <- length(intersect(bottom_selection, correct_values))
fp_bottom <- length(setdiff(bottom_selection, correct_values))
fn_bottom <- length(setdiff(correct_values, bottom_selection))
precision_bottom <- tp_bottom / (tp_bottom + fp_bottom)
recall_bottom <- tp_bottom / (tp_bottom + fn_bottom)
f1_bottom <- 2 * (precision_bottom * recall_bottom) / (precision_bottom + recall_bottom)

# Results
data.frame(
  Method = c("Model Selection", "EMS"),
  Precision = c(precision_top, precision_bottom),
  Recall = c(recall_top, recall_bottom),
  F1_Score = c(f1_top, f1_bottom)
)

# It's calculating precision, recall, and F1 score for each selection method by:
# 
# Finding true positives (tp): values that are both in the selection and in the correct values
# Finding false positives (fp): values that are in the selection but not in the correct values
# Finding false negatives (fn): values that are in the correct values but not in the selection
# 
# 
# It's creating a data frame to compare the performance metrics of both methods side by side
# 
# These metrics are commonly used in classification tasks to evaluate how well a selection method performs:
#   
#   Precision measures accuracy of positive predictions (what percentage of selected items are actually correct)
# Recall measures completeness (what percentage of correct items were actually selected)
# F1 score is the harmonic mean of precision and recall, providing a balanced measure
# 
# The code is essentially evaluating which selection method ("Top" or "Bottom") better identifies the correct values defined in correct_values.

# Function to calculate precision, recall, and F1 metrics
calculate_metrics <- function(selections, correct_values) {
  # Ensure inputs are flattened vectors (in case of nested lists)
  selections <- unlist(selections)
  correct_values <- unlist(correct_values)
  
  # Calculate true positives, false positives, false negatives
  tp <- length(intersect(selections, correct_values))
  fp <- length(setdiff(selections, correct_values))
  fn <- length(setdiff(correct_values, selections))
  
  # Calculate metrics
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  # Return as named list
  return(list(
    true_positives = tp,
    false_positives = fp,
    false_negatives = fn,
    precision = precision,
    recall = recall,
    f1_score = f1
  ))
}

# Define correct values
correct_values <- c(1, 2, 3, 7, 8, 15)

# Define selections (can combine multiple data sources)
top_selection <- c(test1$model.sel, test1.1$model.sel, test1.2$model.se1l, test1.3$model.sel, test2.1$model.sel, test2.2$model.se1l, test2.3$model.sel)
bottom_selection <- c(test1$EMS, test1.1$EMS, test1.2$EMS, test1.3$EMS, test2.1$EMS, test2.2$EMS, test2.3$EMS)

# Calculate metrics for both methods
top_metrics <- calculate_metrics(top_selection, correct_values)
bottom_metrics <- calculate_metrics(bottom_selection, correct_values)

# Create comparison data frame
results <- data.frame(
  Method = c("Model Selection", "EMS"),
  True_Positives = c(top_metrics$true_positives, bottom_metrics$true_positives),
  False_Positives = c(top_metrics$false_positives, bottom_metrics$false_positives),
  False_Negatives = c(top_metrics$false_negatives, bottom_metrics$false_negatives),
  Precision = c(top_metrics$precision, bottom_metrics$precision),
  Recall = c(top_metrics$recall, bottom_metrics$recall),
  F1_Score = c(top_metrics$f1_score, bottom_metrics$f1_score)
)

# Display results
print(results)


# test1.1 <- generate_EMS_Simulation(100, 5, 6, 1, .1, seed = 3, save = "testing.1")
# rm(list = ls())
# source("EMS/EMS_Functions.R")
# test2 <- generate_EMS_Simulation(100, 5, 6, .3, c(1,1,1,1,1,1), seed = 1, save = "testing2")
# test3 <- generate_EMS_Simulation(100, 5, 6, 1, 1, seed = 7, save = "testing3")
# test4 <- generate_EMS_Simulation(200, 5, 6, 1, 1, seed = 7, save = "testing4")

# load("EMS/Sim_EMS_1_CoefChange2_testing.rda")
# test1 <- result
# load("EMS/Sim_EMS_1_CoefChange2_testing2.rda")
# test2 <- result

# load("EMS/Sim_EMS_Seed_1")
# seed1 <- result
# load("EMS/Sim_EMS_Seed_3")
# seed3 <- result
# load("EMS/Sim_EMS_Seed_5")
# seed5 <- result
# load("EMS/Sim_EMS_Seed_5")
# seed5_CoefChange <- result

# load("EMS/Sim_EMS_3_CoefChange2_testing.rda")
# test1 <- result
# load("EMS/Sim_EMS_3_CoefChange2_testing.1.rda")
# test1.1 <- result
# load("EMS/Sim_EMS_5_CoefChange2_testing2.rda")
# test2 <- result
# load("EMS/Sim_EMS_7_CoefChange2_testing3.rda")
# test3 <- result
# 
# 
# seed1[1]
# seed3[1]
# seed5[1]
# seed5_CoefChange[1]


# setup runs for simulations
# def()

# #################################################################################
# ### data to output to the qlt pipeline, plug thisn into the pipeline and run
# # zz <- X[I.m,]
# # I.m
# 
# # ?replace()
# zz = X
# zz[I.m,] = NA
# zz[zz==1] = 2
# zz[zz==0] = 1
# 
# zzy = y
# 
# zz
# zzy
# 
# # str(x)
# # str(y)
# # 
# # make something to build and separate the data by the marker values to match hyper
# 
# marker1 <- zz[,1:6]
# marker1y <- zzy[1:6]

# rm(list = ls())
# start.time <- Sys.time()
# source("C:/Users/Payne/OneDrive/Documents/EMS_Functions.R")  # Load your simulation function
# 
# # === SETUP PARAMETERS === #
# n_seeds <- 50
# seed_list <- seq(1, 100, by = 2) #sample(1:10000, n_seeds)  # Or just use 1:50 if you want predictable seeds
# 
# sample_sizes <- c(100, 50)  # Use different sample sizes
# noise_sets <- list(
#   c(1,1,1,1,1,1),
#   c(.1,.1,.1,.1,.1,.1),
#   c(1,.1,.5,.4,.3,.9)
# )
# 
# missing_vals <- c(0.05, 0.01, -1)  # Different missing data strategies
# 
# # === MAIN LOOP === #
# results <- list()
# counter <- 1
# 
# for (seed in seed_list) {
#   for (n in sample_sizes) {
#     for (noise in noise_sets) {
#       for (missing in missing_vals) {
#         
#         save_name <- paste0("sim_n", n, "_miss", missing,
#                             "_noise", paste(noise, collapse = ""),
#                             "_seed", seed)
#         
#         cat("Running simulation:", save_name, "\n")
#         
#         results[[counter]] <- tryCatch({
#           generate_EMS_Simulation(
#             Sample_size = n,
#             chr = 5,
#             marker.chr = 6,
#             missing = missing,
#             noise_coef = noise,
#             seed = seed,
#             save = save_name
#           )
#         }, error = function(e) {
#           warning(paste("Failed:", save_name, "-", e$message))
#           NULL
#         })
#         
#         counter <- counter + 1
#       }
#     }
#   }
# }
# 
# end.time <- Sys.time()
# print(paste("Total time:", round(difftime(end.time, start.time, units = "mins"), 2)))


# Possible Considerations:
# Error Handling in tryCatch:
#
# If generate_EMS_Simulation consistently fails with certain parameters, you might want to log the configurations that cause the failure to investigate later.
#
# Simulation Output:
#
# The save_name variable suggests that each simulation result is saved with a unique filename. Ensure that the function generate_EMS_Simulation supports saving results to the filesystem and that you have appropriate write permissions.
#
# Performance:
#
# The loops involve a high number of iterations: for n_seeds (50), sample_sizes (2), noise_sets (3), and missing_vals (3), this will lead to
# 50
# ×
# 2
# ×
# 3
# ×
# 3
# =
# 900
# 50×2×3×3=900 simulations. Depending on the complexity of the simulation, this might take a significant amount of time, and you may want to optimize or parallelize the process if it's taking too long.
#
# Would you like suggestions for optimizing this code or any additional information on any specific part?
#

# rm(list = ls())
# start.time <- Sys.time()
# source("EMS/EMS_Functions.R")  # Load your simulation function
# 
# # === SETUP PARAMETERS === #
# n_seeds <- 50
# seed_list <- 1#seq(1, 100, by = 2)  # Or just use 1:50 if you want predictable seeds
# 
# sample_sizes <- c(100, 50)  # Use different sample sizes
# noise_sets <- list(
#   c(1,1,1,1,1,1),
#   c(.1,.1,.1,.1,.1,.1),
#   c(1,.1,.5,.4,.3,.9)
# )
# 
# # missing_vals <- c(0.05, 0.01, -1, 0)  # Different missing data strategies
# missing_vals <- c( 0)  # Different missing data strategies
# 
# # Create a log file to track errors
# error_log <- "simulation_errors.log"
# if (file.exists(error_log)) file.remove(error_log)  # Remove old log if it exists
# 
# # === MAIN LOOP === #
# results <- list()
# counter <- 1
# 
# for (seed in seed_list) {
#   for (n in sample_sizes) {
#     for (noise in noise_sets) {
#       for (missing in missing_vals) {
#         
#         save_name <- paste0("sim_n", n, "_miss", missing,
#                             "_noise", paste(noise, collapse = ""),
#                             "_seed", seed)
#         
#         cat("Running simulation:", save_name, "\n")
#         
#         # results[[counter]]
#         results<- tryCatch({
#           # Run the simulation function
#           generate_EMS_Simulation(
#             n,
#             chr = 5,
#             marker.chr = 6,
#             missing,
#             noise,
#             seed = seed,
#             save = save_name
#           )
#         }, error = function(e) {
#           # Capture and log detailed error information
#           error_message <- paste("Failed:", save_name, "- Error:", e$message)
#           cat(error_message, "\n")  # Print the error to the console
#           
#           # Log the error to a file
#           write(paste(Sys.time(), error_message), file = error_log, append = TRUE)
#           
#           # Return NULL for this simulation so the loop continues
#           return(NULL)
#         })
#         
#         counter <- counter + 1
#       }
#     }
#   }
# }
# 
# end.time <- Sys.time()
# print(paste("Total time:", round(difftime(end.time, start.time, units = "mins"), 2)))

# Assuming test1, test1.1, ..., test2.3 are defined in your environment
# source("EMS/EMS_Functions.R")  # Load your simulation function
# 
# correct_values <- c(1, 2, 3, 7, 8, 15)
# 
# file_list <- list(
#   "sim_n100_miss0.05_noise111111_seed1",
#   "sim_n100_miss0.05_noise0.10.10.10.10.1_seed1",
#   "sim_n100_miss0.05_noise1.10.50.40.30.9_seed1",
#   "sim_n100_miss0.01_noise111111_seed1",
#   "sim_n100_miss0.01_noise0.10.10.10.10.1_seed1",
#   "sim_n100_miss0.01_noise1.10.50.40.30.9_seed1",
#   "sim_n100_miss-1_noise111111_seed1",
#   "sim_n100_miss-1_noise0.10.10.10.10.1_seed1",
#   "sim_n100_miss-1_noise1.10.50.40.30.9_seed1"#,
#   # "sim_n50_miss0.05_noise111111_seed1",
#   # "sim_n50_miss0.05_noise0.10.10.10.10.1_seed1",
#   # "sim_n50_miss0.05_noise1.10.50.40.30.9_seed1",
#   # "sim_n50_miss0.01_noise111111_seed1",
#   # "sim_n50_miss0.01_noise0.10.10.10.10.1_seed1",
#   # "sim_n50_miss0.01_noise1.10.50.40.30.9_seed1",
#   # "sim_n50_miss-1_noise111111_seed1",
#   # "sim_n50_miss-1_noise0.10.10.10.10.1_seed1",
#   # "sim_n50_miss-1_noise1.10.50.40.30.9_seed1"
# )

# # Directory where files are saved
# dir_path <- "EMS/"
# 
# # Initialize lists
# top_selection <- list()
# bottom_selection <- list()
# 
# # Loop over each file name
# for (file_name in file_list) {
#   file_path <- paste0(dir_path, file_name, ".rda")
#   
#   # Load the .rda file and capture the loaded object name
#   loaded_objs <- load(file_path)
#   
#   # Extract the first (or only) object from the file
#   obj <- get(loaded_objs[1])
#   
#   # Extract model selection and EMS components
#   top_selection[[file_name]] <- obj$model.sel
#   bottom_selection[[file_name]] <- obj$EMS
# }

# Evaluate using your function
# evaluate_methods(top_selection, bottom_selection, correct_values,
#                  top_name = "Model Selection", bottom_name = "EMS")

# 
# out2 <- load("EMS/sim_n100_miss0.05_noise111111_seed1.rda")
# 
# 
# source("EMS/EMS_Functions.R")  # Load your simulation function
# 
# # Define correct values
# correct_values <- c(1, 2, 3, 7, 8, 15)
# 
# # Define the directory where .rda files are stored
# dir_path <- "EMS/Sims"
# 
# # List all .rda files in the directory
# rda_files <- list.files(dir_path, pattern = "^sim_n.*\\.rda$", full.names = TRUE)
# rda_files_miss0 <- rda_files[grepl("miss-1_", rda_files)]
# 
# # Initialize output lists
# top_selection <- list()
# bottom_selection <- list()
# 
# testing <- function(rdaFile) {
#   # Loop over all .rda files
#   for (file_path in rdaFile) {
#     # Print file being processed (optional for debugging)
#     cat("Processing:", file_path, "\n")
#   
#     # Load the .rda file and get the loaded object name(s)
#     loaded_objs <- tryCatch(load(file_path), error = function(e) NULL)
#   
#     # If load failed or nothing loaded, skip
#     if (is.null(loaded_objs) || length(loaded_objs) == 0) {
#       warning(paste("Failed to load or no objects in file:", file_path))
#       next
#     }
#   
#     # Extract object (first loaded object)
#     obj <- get(loaded_objs[1], envir = .GlobalEnv)
#   
#     # File label for storing in list (cleaned from file name)
#     file_label <- tools::file_path_sans_ext(basename(file_path))
#   
#     # Try to extract components, if they exist
#     if (!is.null(obj$model.sel) && !is.null(obj$EMS)) {
#       top_selection[[file_label]] <- obj$model.sel
#       bottom_selection[[file_label]] <- obj$EMS
#     } else {
#       warning(paste("model.sel or EMS not found in object from", file_label))
#     }
#   }
#   
#   # Evaluate (assuming this function exists in your environment)
#   evaluate_methods(top_selection, bottom_selection, correct_values,
#                    top_name = "Model Selection", bottom_name = "EMS")
# }
# 
# testing <- function(rdaFile) {
#   top_selection <- list()
#   bottom_selection <- list()
#   
#   for (file_path in rdaFile) {
#     cat("Processing:", file_path, "\n")
#     
#     tmp_env <- new.env()
#     loaded_objs <- tryCatch(load(file_path, envir = tmp_env), error = function(e) NULL)
#     
#     if (is.null(loaded_objs) || length(loaded_objs) == 0) {
#       warning(paste("Failed to load or no objects in file:", file_path))
#       next
#     }
#     
#     obj <- tmp_env[[loaded_objs[1]]]  # safer than .GlobalEnv
#     
#     file_label <- tools::file_path_sans_ext(basename(file_path))
#     
#     if (is.list(obj) && all(c("model.sel", "EMS") %in% names(obj))) {
#       top_selection[[file_label]] <- obj$model.sel
#       bottom_selection[[file_label]] <- obj$EMS
#     } else {
#       warning(paste("model.sel or EMS not found in object from", file_label))
#     }
#   }
#   
#   if (!exists("evaluate_methods")) {
#     stop("Function 'evaluate_methods' is not defined.")
#   }
#   
#   evaluate_methods(top_selection, bottom_selection, correct_values,
#                    top_name = "Model Selection", bottom_name = "EMS")
# }
# 
# testing(rda_files_miss0)

