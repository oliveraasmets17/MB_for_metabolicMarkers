


# Command line arguments ----
#------------------------------------------------------------#
#                                                            #
#                SET UP COMMAND LINE ARGUMENTS               # 
#                                                            #
#------------------------------------------------------------#
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")

# Expecting:
# target timeframe model_input_list vector_of_OTUs output_file

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 

# Declare target variable
target <- args[1]

# Declare modelling timeframe
timeframe <- args[2]

# Declare model list used 
input_list <- args[3]

# Declare features of interest
variables_of_interest <- args[4:(length(args)-1)]

# Declare output file name 
output_file <- args[length(args)]








# Program setup phase ----
#------------------------------------------------------------#
#                                                            #
#                     PROGRAM SETUP PHASE                    # 
#                                                            #
#------------------------------------------------------------#

# Load packages
library("dplyr")
library("caret")
library("stringr")
library("readr")
library("doParallel")
library("DALEX")

# Load saved RF models
models <- readRDS(file = input_list)

# Define input dataset
if (timeframe == "18m"){
  master <- readRDS(file = "full_MLdataset_50_for_18m.rds")
} else if (timeframe == "48m"){
  master <- readRDS(file = "full_MLdataset_50_for_48m.rds")
}

# Define all targets 
baselines <- c("P_gl0", "P_gl120", "matsuda", "disposition", "secretion", "HbA1cP", "P_ins0", "P_ins120")
targets <- paste("Predict_", baselines, sep = "")


# Define additional covariates
covar_of_interest <- c("Age", "BMI", "F_famdmfamily", "F_elevgl")


# IDs for genera used
genera_of_interest <- readRDS(file = "genera50.rds")


# Clean dataset for machine learning purposes
master_clean <- master %>%
  dplyr::filter(is.na(Reimb_DM_date_122017) == TRUE & is.na(Death_Date_Feb2017) == TRUE) %>%
  dplyr::select(all_of(baselines), all_of(targets), all_of(covar_of_interest), all_of(genera_of_interest)) %>%
  dplyr::filter(complete.cases(.)) %>%
  as.data.frame()







# Explain results ----
#------------------------------------------------------------#
#                                                            #
#                 EXPLAIN MODELLING RESULTS                  # 
#                                                            #
#------------------------------------------------------------#

# Return model names
all_models <- names(models)

# Filter out models considered
considered_models <- all_models[str_detect(all_models, fixed(target))]

# Start variable explanations
ALE_output <- data.frame()

for (model_index in 1:length(considered_models)){
  
  # Select model 
  model_used = models[[considered_models[model_index]]]
  
  # Derive parameters
  model_parameters = unlist(str_split(considered_models[model_index], pattern = "_"))
  seed = model_parameters[length(model_parameters)] %>% parse_number()
  
  # Recreate data split
  proportion = 0.75
  
  set.seed(seed)
  trainIndex = sample(1:nrow(master_clean), size = round(proportion*nrow(master_clean)), replace = FALSE)
  
  # Split data
  test = master_clean[-trainIndex ,]
  test_y = test %>% dplyr::pull(target)
  test_x = test %>% dplyr::select(-one_of(target))
  
  # Build explainer
  explainer = DALEX::explain(model_used, label = "rf", 
                             data = test_x, 
                             y = test_y,
                             colorize = FALSE,
                             verbose = FALSE)
  
  # Acumulated local efects data
  ALE = ingredients::accumulated_dependency(explainer, variables = variables_of_interest)
  
  ALE_data = ALE %>% 
    as.data.frame() %>% 
    dplyr::mutate(seed = seed,
                  target = target)
  
  # Append data
  ALE_output = rbind(ALE_output, ALE_data)
}


print(head(ALE_output))
saveRDS(ALE_output, file = output_file)
