

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
# timeframe RF_output model_List_output

# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 

# Declare modelling timeframe
timeframe <- args[1]

# Declare outputs
RF_output_file <- args[2]

# Declare target variable
modelList_output_file <- args[3]





# Program setup phase ----
#------------------------------------------------------------#
#                                                            #
#                     PROGRAM SETUP PHASE                    # 
#                                                            #
#------------------------------------------------------------#

# Load packages
library("dplyr")
library("caret")


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








# Run RF models ----
#------------------------------------------------------------#
#                                                            #
#                   RUN RANDOMFOREST MODELS                  # 
#                                                            #
#------------------------------------------------------------#

run_RF <- TRUE
if (run_RF == TRUE){
  
  # Build the model
  RF_output <- data.frame()
  model_List <- list()
  
  for (s in 1:200){
    
    # Set random seed
    seed <- s
    
    # Split data
    proportion <- 0.75
    
    set.seed(seed)
    trainIndex <- sample(1:nrow(master_clean), size = round(proportion*nrow(master_clean)), replace = FALSE)
    
    # Split data
    train <- master_clean[trainIndex, ]
    test <- master_clean[-trainIndex, ]
    
    set.seed(seed)
    control <- trainControl(method = "cv", 
                            number = 10, 
                            savePredictions = "final",
                            search = "random", 
                            allowParallel = FALSE)
    
    for (i in targets){
      target_var <- train %>% dplyr::pull(i)
      predictors <- train %>% dplyr::select(-contains("Predict"))
      
      rf <- train(y = target_var,
                  x = predictors,
                  trControl = control,
                  method = "ranger",
                  tuneLength = 100, 
                  preProcess = c("center", "scale"),
                  importance = 'permutation')
      
      modelname = paste(i, "_", s, sep = "")
      model_List[[modelname]] = rf
      
      PredictorSet <- varImp(rf)$importance %>%
        tibble::rownames_to_column(var = "predictor_code") %>%
        dplyr::arrange(desc(Overall))
      
      run_result = PredictorSet %>%
        dplyr::mutate(target = i, 
                      seed = s,
                      RMSE = caret::RMSE(predict.train(rf, newdata = test), test[ ,i]), 
                      MAE = caret::MAE(predict.train(rf, newdata = test), test[ ,i]))
      
      RF_output = rbind(RF_output, run_result)
    }
  }
  
  # Add significant metadata
  RF_output = RF_output %>%
    dplyr::mutate(method = "RF", 
                  preprocess = "CLR", 
                  timeWindow = timeframe)
  
  # Save the outputs
  saveRDS(RF_output, file = RF_output_file)
  saveRDS(model_List, file = modelList_output_file)
  
}

