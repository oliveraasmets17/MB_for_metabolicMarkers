

# Program setup phase ----
#------------------------------------------------------------#
#                                                            #
#                     PROGRAM SETUP PHASE                    # 
#                                                            #
#------------------------------------------------------------#

# Load packages
library("tidyverse")
library("xtable")
library("xlsx")
library("ggrepel")



# Read metadata
RData_folder <- "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_METSIM/RData/"
physeqGenus <- readRDS(file = file.path(RData_folder, "physeqGenus.rds"))

# Read help data
taxonomy_genus_help <- readRDS(file = file.path(RData_folder, "taxonomy_genus_help.rds"))

taxonomy_data_help <- readRDS(file = file.path(RData_folder, "taxonomy_data_help.rds")) %>% 
  dplyr::mutate(phylum = substring(phylum, 6),
                class = substring(class, 6),
                order = substring(order, 6),
                family = substring(family, 6),
                genus = substring(genus, 6)) %>% 
  dplyr::select(-one_of("species"))

# Read precalculated modeling results
Results_folder <- "C:/Users/oliver17/Desktop/Doktorantuur/hpc/Results"
RF_18m_b1 <- readRDS(file = file.path(Results_folder, "RF_CLR_18m_sensitivity_output_b1.rds"))
RF_18m_b2 <- readRDS(file = file.path(Results_folder, "RF_CLR_18m_sensitivity_output_b2.rds"))

RF_48m_b1 <- readRDS(file = file.path(Results_folder, "RF_CLR_48m_sensitivity_output_b1.rds"))
RF_48m_b2 <- readRDS(file = file.path(Results_folder, "RF_CLR_48m_sensitivity_output_b2.rds"))

# Read in baseline results
RF_18m_bl <- readRDS(file = file.path(Results_folder, "RF_CLR_18m_sensitivity_baseline_output.rds"))
RF_48m_bl <- readRDS(file = file.path(Results_folder, "RF_CLR_48m_sensitivity_baseline_output.rds"))





# Analyze modelling results ----
#------------------------------------------------------------#
#                                                            #
#                 ANALYZE MODELLING RESULTS                  # 
#                                                            #
#------------------------------------------------------------#

# Gather performance data
perf_RF_18m <- rbind(RF_18m_b1, RF_18m_b2) %>% dplyr::distinct(target, seed, method, RMSE) %>% dplyr::rename("RMSE_18" = "RMSE")
perf_RF_48m <- rbind(RF_48m_b1, RF_48m_b2) %>% dplyr::distinct(target, seed, method, RMSE) %>% dplyr::rename("RMSE_48" = "RMSE")

perf_RF_18m_bl <- RF_18m_bl %>% dplyr::distinct(target, seed, method, RMSE) %>% dplyr::rename("RMSE_18_bl" = "RMSE")
perf_RF_48m_bl <- RF_48m_bl %>% dplyr::distinct(target, seed, method, RMSE) %>% dplyr::rename("RMSE_48_bl" = "RMSE")

perf_RF <- perf_RF_18m %>% 
  dplyr::left_join(perf_RF_48m, by = c("target", "seed", "method")) %>% 
  dplyr::left_join(perf_RF_18m_bl, by = c("target", "seed", "method")) %>%
  dplyr::left_join(perf_RF_48m_bl, by = c("target", "seed", "method"))  


# Summarize changes in RMSE and counts for winners
summary_table <- perf_RF %>% 
  dplyr::mutate(chg_18 = RMSE_18 - RMSE_18_bl,
                chg_48 = RMSE_48 - RMSE_48_bl) %>%
  dplyr::group_by(target) %>%
  dplyr::summarise(mean18 = mean(chg_18, na.rm = T),
                   sd18 = 1.96*sd(chg_18, na.rm = T),
                   n18 = sum(RMSE_18 < RMSE_18_bl),
                   mean48 = mean(chg_48, na.rm = T),
                   sd48 = 1.96*sd(chg_48, na.rm = T),
                   n48 = sum(RMSE_48 < RMSE_48_bl)) %>%
  dplyr::mutate(p18 = round(n18/200*100, 2),
                p48 = round(n48/200*100, 2),
                n18_out = paste(n18, " (", p18, "%)", sep = ""),
                n48_out = paste(n48, " (", p48, "%)", sep = ""),
                mean18_out  = paste(round(mean18, 3), " (",round(sd18, 4), ")", sep = ""),
                mean48_out = paste(round(mean48, 3), " (",round(sd48, 4), ")", sep = "")) 

summary_table_out <- summary_table%>% 
  dplyr::select(target, mean18_out, n18_out, mean48_out, n48_out)

print(xtable(summary_table_out), include.rownames = FALSE)
write.xlsx(summary_table_out, file = "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_ML_T2D/modelling_summary.xlsx")


# Bonferroni tests
bonferroni_output <- data.frame()
for (i in c(summary_table$n18, summary_table$n48)){
  
  test <- binom.test(i, 200, alternative = "greater")
  run <- data.frame(result = i, p_value = test$p.value) 
  
  bonferroni_output <- rbind(bonferroni_output, run)
}

bonferroni_output_BH <- bonferroni_output %>%
  dplyr::mutate(Bonferroni = p.adjust(p_value, method = "bonferroni"),
                Significance = ifelse(Bonferroni <= 0.05, "Significant", "NS"))


# Summarizing the results for feature selection - select by total average Overall-score
RF_18m_summarized <- rbind(RF_18m_b1, RF_18m_b2) %>% 
  dplyr::group_by(target, predictor_code) %>% 
  dplyr::summarise(n = n(), mean_imp = mean(Overall)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(str_detect(predictor_code, fixed("OTU."))) %>%
  dplyr::arrange(target, desc(mean_imp)) %>% 
  dplyr::left_join(taxonomy_data_help, by = c("predictor_code" = "OTU"))

RF_48m_summarized <- rbind(RF_48m_b1, RF_48m_b2) %>% 
  dplyr::group_by(target, predictor_code) %>% 
  dplyr::summarise(n = n(), mean_imp = mean(Overall), median_imp = median(Overall), sd_imp = sd(Overall), min_imp = min(Overall), max_imp = max(Overall)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(str_detect(predictor_code, fixed("OTU."))) %>%
  dplyr::arrange(target, desc(mean_imp)) %>% 
  dplyr::left_join(taxonomy_data_help, by = c("predictor_code" = "OTU"))




# Filtering relevant taxa
RF_18m_interesting <- RF_18m_summarized %>%
  dplyr::filter(target %in% c("Predict_HbA1cP", "Predict_secretion",  "Predict_P_ins120")) %>%
  dplyr::group_by(target) %>% 
  dplyr::top_n(50, wt = mean_imp) %>% 
  dplyr::select(target, predictor_code, phylum, order, class, family, genus, mean_imp) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()


RF_48m_interesting <- RF_48m_summarized %>%
  dplyr::filter(target %in% c("Predict_P_ins120", "Predict_P_ins0", "Predict_secretion")) %>%
  dplyr::group_by(target) %>% 
  dplyr::top_n(50, wt = mean_imp) %>% 
  dplyr::select(target, predictor_code, phylum, order, class, family, genus, mean_imp) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()


# Save the files
saveRDS(RF_18m_interesting, file = "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_ML_T2D/18m_features.rds")
saveRDS(RF_48m_interesting, file = "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_ML_T2D/48m_features.rds")







# Plot the results ----
#------------------------------------------------------------#
#                                                            #
#                      PLOTTING RESULTS                      # 
#                                                            #
#------------------------------------------------------------#

# Plot the peaks for feature importance scores
plot_importance1 <- RF_18m_interesting %>% 
  dplyr::mutate(help = 1) %>% 
  dplyr::left_join(taxonomy_genus_help, by = c("predictor_code" = "OTU")) %>%
  dplyr::group_by(target) %>%
  dplyr::mutate(x_new = cumsum(help)) %>%
  dplyr::ungroup() %>%
  # Choose arbitrarily genera to highlight
  dplyr::mutate(label = case_when(target == "Predict_HbA1cP" & x_new <= 5 ~ 1,
                                  target == "Predict_P_ins120" & x_new <= 5 ~ 1,
                                  target == "Predict_secretion" & x_new <= 3 ~ 1,
                                  TRUE ~ 0),
                target = case_when(target == "Predict_HbA1cP" ~  "HbA1c",
                                   target == "Predict_secretion" ~  "Secretion index",
                                   target == "Predict_P_ins120" ~  "2h insulin"),
                timeframe = "18m", 
                taxa = case_when(taxa == "Clostridiales vadinBB60 group / gut metagenome" ~ "Clostridiales vadinBB60 group \n gut metagenome",
                                 taxa == "Clostridiales vadinBB60 group / uncultured bacterium" ~ "Clostridiales vadinBB60 group \n uncultured bacterium",
                                 taxa == "Clostridiales vadinBB60 group / uncultured organism" ~ "Clostridiales vadinBB60 group \n uncultured organism",
                                 TRUE ~ taxa))

p1 <- ggplot(plot_importance1, aes(x = x_new, y = mean_imp, color = factor(label))) + 
  geom_point() +
  geom_label_repel(data = plot_importance1 %>% 
                     dplyr::filter(label == 1 & target == "HbA1c"), aes(x = x_new, y = mean_imp, label = taxa), 
                   size = 3.5, hjust = -0.3, nudge_x = 15, segment.size = 0.2,
                   direction = "y", force = 100, seed = 3, xlim = c(5, NA)) +
  geom_label_repel(data = plot_importance1 %>% 
                     dplyr::filter(label == 1 & target != "HbA1c"), aes(x = x_new, y = mean_imp, label = taxa), 
                   size = 3.5, hjust = -0.3, nudge_x = 15, segment.size = 0.2,
                   direction = "both", force = 2, xlim = c(5, NA)) +
  facet_wrap(vars(target), scales = "free") +
  ylab("Average permutation importance score") +
  xlab("") + 
  scale_color_manual(values = c("gray", "black"), guide = FALSE) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.background = element_rect(size = 1), 
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
p1

ggsave(plot = p1, filename = "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_ML_T2D/Feature_importance_screePlot_18m.png", width = 30, height = 10, units = "cm")
ggsave(plot = p1, filename = "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_ML_T2D/Feature_importance_screePlot_18m.svg", width = 30, height = 10, units = "cm")


# Plot the peaks for feature importance scores
plot_importance2 <- RF_48m_interesting %>% 
  dplyr::mutate(help = 1) %>% 
  dplyr::left_join(taxonomy_genus_help, by = c("predictor_code" = "OTU")) %>%
  dplyr::group_by(target) %>%
  dplyr::mutate(x_new = cumsum(help)) %>%
  dplyr::ungroup() %>%
  # Choose arbitrarily genera to highlight
  dplyr::mutate(label = case_when(target == "Predict_P_ins0" & x_new <= 3 ~ 1,
                                  target == "Predict_P_ins120" & x_new <= 5 ~ 1,
                                  target == "Predict_secretion" & x_new <= 3 ~ 1,
                                  TRUE ~ 0),
                target = case_when(target == "Predict_P_ins0" ~  "Fasting insulin",
                                   target == "Predict_secretion" ~  "Secretion index",
                                   target == "Predict_P_ins120" ~  "2h insulin"),
                timeframe = "48m")

p2 <- ggplot(plot_importance2, aes(x = x_new, y = mean_imp, color = factor(label))) + 
  geom_point() +
  geom_label_repel(data = plot_importance2 %>% 
                     dplyr::filter(label == 1), aes(x = x_new, y = mean_imp, label = taxa), 
                   size = 3.5, hjust = -0.3, nudge_y = 0.05, nudge_x = 15, segment.size = 0.2,
                   direction = "both", force = 2) +
  facet_wrap(vars(target), scales = "free") +
  ylab("Average permutation importance score") +
  xlab("") + 
  scale_color_manual(values = c("gray", "black"), guide = FALSE) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.background = element_rect(size = 1),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
p2

ggsave(plot = p2, filename = "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_ML_T2D/Feature_importance_screePlot_48m.png", width = 30, height = 10, units = "cm")
ggsave(plot = p2, filename = "C:/Users/oliver17/Desktop/Doktorantuur/Projekt_ML_T2D/Feature_importance_screePlot_48m.svg", width = 30, height = 10, units = "cm")
