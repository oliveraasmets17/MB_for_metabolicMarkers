

# Program setup phase ----
#------------------------------------------------------------#
#                                                            #
#                     PROGRAM SETUP PHASE                    # 
#                                                            #
#------------------------------------------------------------#

# Load packages
library("tidyverse")

# Load taxonomy data
taxonomy_genus_help <- readRDS(file = "C:/Users/oliver17/Desktop/Doktorantuur/METSIM/RData/taxonomy_genus_help.rds") %>% 
  dplyr::mutate(n = nchar(taxa),
                taxa = case_when(taxa == "Clostridiales vadinBB60 group / uncultured bacterium" ~ "Clostridiales vadinBB60 group\n uncultured bacterium",
                                 taxa == "Clostridiales vadinBB60 group / uncultured organism" ~ "Clostridiales vadinBB60 group\n uncultured organism",
                                 taxa == "Clostridiales vadinBB60 group / gut metagenome" ~ "Clostridiales vadinBB60 group\n gut metagenome",
                                 taxa == "Coriobacteriales Incertae Sedis / uncultured" ~ "Coriobacteriales Incertae Sedis\n uncultured",
                                 taxa == "uncultured Thermoanaerobacterales bacterium" ~ "Thermoanaerobacterales\n uncultured bacterium",
                                 taxa == "Mollicutes RF39 / uncultured bacterium" ~ "Mollicutes RF39\n uncultured bacterium",
                                 taxa == "Izimaplasmatales / uncultured organism" ~ "Izimaplasmatales\n uncultured organism",
                                 taxa == "Muribaculaceae / uncultured bacterium" ~ "Muribaculaceae\n uncultured bacterium", 
                                 TRUE ~ taxa))


# Define the most relevant taxa used
RF_18m_interesting <- readRDS(file = "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_ML_T2D/18m_features.rds")
RF_48m_interesting <- readRDS(file = "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_ML_T2D/48m_features.rds")

taxa_importance <- rbind(RF_18m_interesting %>% dplyr::mutate(timeframe = "18m"), 
                         RF_48m_interesting %>% dplyr::mutate(timeframe = "48m")) %>%
  dplyr::select(target, timeframe, predictor_code, mean_imp) %>% 
  dplyr::left_join(taxonomy_genus_help, by = c("predictor_code" = "OTU")) %>%
  dplyr::mutate(target = case_when(target == "Predict_HbA1cP" ~  "HbA1c",
                                   target == "Predict_secretion" ~  "Secretion index",
                                   target == "Predict_P_ins120" ~  "2h insulin",
                                   target == "Predict_P_ins0" ~  "Fasting insulin",
                                   TRUE ~ "something else"))  

find_taxa_order = function(p_target, p_timeframe){
  data_out = taxa_importance %>%
    dplyr::filter(target == p_target & timeframe == p_timeframe) %>%
    dplyr::arrange(desc(mean_imp)) %>%
    dplyr::pull(taxa)
  
  return(data_out)
}



# Load modelling results - accumulated local effects
dir <- "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_ML_T2D/Results"

all_files <- list.files(path = dir)
output <- data.frame()
for (file in all_files){
  
  metadata = unlist(str_split(file, pattern = "_"))
  
  data = readRDS(file = file.path(dir, file)) %>% 
    dplyr::mutate(timeframe = metadata[3])
  
  output = rbind(output, data)
}

output_final <- output %>% 
  dplyr::left_join(taxonomy_genus_help, by = c("_vname_" = "OTU"))



# Calculate percentails of the CLR abundances
master18 <- readRDS(file = "C:/Users/oliver17/Desktop/Doktorantuur/METSIM/RData/full_MLdataset_50_for_18m.rds")
master48 <- readRDS(file = "C:/Users/oliver17/Desktop/Doktorantuur/METSIM/RData/full_MLdataset_50_for_48m.rds")

CLR_CI <- master18[ ,c(1, 33:204)] %>%
  tidyr::gather(key, value, -sampleCode) %>%
  dplyr::group_by(key) %>%
  dplyr::summarise(q2.5 = quantile(value, 0.025, na.rm = T),
                   q50 = quantile(value, 0.5, na.rm = T), 
                   q97.5 = quantile(value, 0.975, na.rm = T))

# Merge the percentiles
output_final <- output_final %>%
  dplyr::left_join(CLR_CI, by = c("_vname_" = "key"))


# Aggregate data - produce mean trajectory
output_mean <- output %>%
  dplyr::left_join(CLR_CI, by = c("_vname_" = "key")) %>%
  dplyr::rename("OTU" = "_vname_") %>% 
  dplyr::filter(`_x_` >= q2.5 & `_x_` <= q97.5) %>%
  dplyr::mutate(x_agg = ((`_x_`*10) %/% 5)/2) %>% 
  dplyr::group_by(x_agg, OTU, timeframe, target) %>% 
  dplyr::summarise(y_agg = mean(`_yhat_`)) %>%
  dplyr::ungroup() %>% 
  dplyr::select(target, timeframe, OTU, x_agg, y_agg) %>% 
  dplyr::arrange(target, timeframe, OTU, x_agg, y_agg) %>% 
  as.data.frame() %>%
  dplyr::left_join(taxonomy_genus_help[ ,c("OTU", "taxa")], by = "OTU")




# Plotting 18 month predictors for secretion index
plot_help18 <- plot_importance1 %>%
  dplyr::select(target, predictor_code, label, timeframe) %>%
  dplyr::mutate(target = case_when(target == "HbA1c" ~ "Predict_HbA1cP",
                                   target == "Secretion index" ~ "Predict_secretion",
                                   target == "2h insulin" ~ "Predict_P_ins120"))  

plot_help48 <- plot_importance2 %>%
  dplyr::select(target, predictor_code, label, timeframe) %>%
  dplyr::mutate(target = case_when(target == "Fasting insulin" ~ "Predict_P_ins0",
                                   target == "Secretion index" ~ "Predict_secretion",
                                   target == "2h insulin" ~ "Predict_P_ins120"))

output_mean_help <- output_mean %>%
  dplyr::left_join(plot_help18, by = c("OTU" = "predictor_code", "target", "timeframe")) %>%
  dplyr::left_join(plot_help48, by = c("OTU" = "predictor_code", "target", "timeframe")) %>%
  dplyr::filter(label.x == 1 | label.y == 1) %>%
  dplyr::mutate(target = case_when(target == "Predict_HbA1cP" ~  "HbA1c",
                                   target == "Predict_secretion" ~  "Secretion index",
                                   target == "Predict_P_ins120" ~  "2h insulin",
                                   target == "Predict_P_ins0" ~  "Fasting insulin",
                                   TRUE ~ "something else"),
                target_new = paste(target, timeframe, sep = " "))

output_final_help <- output_final %>%
  dplyr::left_join(plot_help18, by = c("_vname_" = "predictor_code", "target", "timeframe")) %>%
  dplyr::left_join(plot_help48, by = c("_vname_" = "predictor_code", "target", "timeframe")) %>%
  dplyr::filter(label.x == 1 | label.y == 1) %>%
  dplyr::mutate(target = case_when(target == "Predict_HbA1cP" ~  "HbA1c",
                                   target == "Predict_secretion" ~  "Secretion index",
                                   target == "Predict_P_ins120" ~  "2h insulin",
                                   target == "Predict_P_ins0" ~  "Fasting insulin",
                                   TRUE ~ "something else"),
                target_new = paste(target, timeframe, sep = " "))



# Define function to draw the ALE plots
make_plot = function(p_target, label_x, label_y, p_ncol, p_tf){

  taxa_order <- find_taxa_order(p_target = p_target, p_timeframe = p_tf)
  
  if (is.na(label_y)){
    d1 <- output_final_help %>% dplyr::filter(target == p_target & label.x == label_x)
    d2 <- output_mean_help %>% dplyr::filter(target == p_target & label.x == label_x )
  } else{
    d1 <- output_final_help %>% dplyr::filter(target == p_target & label.y == label_y)
    d2 <- output_mean_help %>% dplyr::filter(target == p_target & label.y == label_y)
  }
  d1 <- d1 %>% dplyr::mutate(taxa = factor(taxa, levels = taxa_order))
  d2 <- d2 %>% dplyr::mutate(taxa = factor(taxa, levels = taxa_order))
  
  d1 <- d1 %>% dplyr::mutate(taxa = case_when(taxa == "Clostridiales vadinBB60 group / uncultured bacterium" ~
                                                "Clostridiales vadinBB60 \n group \n uncultured bacterium",
                                              taxa == "Clostridiales vadinBB60 group / gut metagenome" ~
                                                "Clostridiales vadinBB60 \n group \n gut metagenome",
                                              taxa == "Muribaculaceae / metagenome" ~
                                                "Muribaculaceae \n metagenome",
                                              taxa == "Rhodospirillales / uncultured" ~
                                                "Rhodospirillales \n uncultured",
                                              taxa == "Prevotellaceae / uncultured" ~ 
                                                "Prevotellaceae \n uncultured",
                                              TRUE ~ as.character(taxa)))
  
  d2 <- d2 %>% dplyr::mutate(taxa = case_when(taxa == "Clostridiales vadinBB60 group / uncultured bacterium" ~
                                                "Clostridiales vadinBB60 \n group \n uncultured bacterium",
                                              taxa == "Clostridiales vadinBB60 group / gut metagenome" ~
                                                "Clostridiales vadinBB60 \n group \n gut metagenome",
                                              taxa == "Muribaculaceae / metagenome" ~
                                                "Muribaculaceae \n metagenome",
                                              taxa == "Rhodospirillales / uncultured" ~
                                                "Rhodospirillales \n uncultured",
                                              taxa == "Prevotellaceae / uncultured" ~ 
                                                "Prevotellaceae \n uncultured",
                                              TRUE ~ as.character(taxa)))
  
  
  p <- ggplot() + 
    geom_path(data = d1, aes(x = `_x_`, y = `_yhat_`, group = seed), color = "#0892d0", size = 0.3, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 0, linetype = 3, size = 1) +
    geom_path(data = d2, aes(x = x_agg, y = y_agg, group = taxa), color = "darkorange", size = 1.5) + 
    facet_wrap(vars(taxa), scales = "free", nrow = p_ncol) + 
    xlab("") + 
    ylab("") +
    theme_bw() + 
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          strip.text = element_text(size = 7),
          axis.text = element_text(size = 7),
          strip.background = element_rect(size = 1))
  
  return(p)
}

p1 <- make_plot(p_target = "2h insulin", label_x = 1, label_y = NA, p_ncol = 1, p_tf = "18m")
p2 <- make_plot(p_target = "2h insulin", label_x = NA, label_y = 1, p_ncol = 1, p_tf = "48m")

p3 <- make_plot(p_target = "HbA1c", label_x = 1, label_y = NA, p_ncol = 1, p_tf = "18m")
p4 <- make_plot(p_target = "Fasting insulin", label_x = NA, label_y = 1, p_ncol = 1, p_tf = "48m")

p5 <- make_plot(p_target = "Secretion index", label_x = 1, label_y = NA, p_ncol = 1, p_tf = "18m")
p6 <- make_plot(p_target = "Secretion index", label_x = NA, label_y = 1, p_ncol = 1, p_tf = "48m")
