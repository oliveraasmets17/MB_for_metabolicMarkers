# MB_for_metabolicMarkers

This folder contains codes for the paper "Machine learning reveals time-varying microbial predictors with complex effects on glucose regulation". 

1) Build_RF_baseline_models.R - skript to run random forest models considering only clinical markers as predictors. Outputs a list consisting of a data frame of performance estimates and variable importances, and a list of random forest models for all target variables and random seeds. 
2) Build_RF_models.R - skript to run random forest models considering clinical markers and CLR-transformed microbial markers as predictors. Outputs a list consisting of a data frame of performance estimates and variable importances, and a list of random forest models for all target variables and random seeds. 
4) Analyze_modelling_results.R - skipt that compares the results of the baseline models with models including microbial predictors, selects top microbial markers and plots the "screeplot" for feature importance. Uses output of Build_RF_baseline_models.R and Build_RF_models.R.
3) Explain_results.R - skript that describes the effect of the predictors based on the random forest models generated by Build_RF_models.R and selected predictors. Outputs a dataframe of ALE results.  
5) Plot_ALE_results.R - skript that plots the ALE results over 200 runs based on the results produced by Explain_results.R. 
