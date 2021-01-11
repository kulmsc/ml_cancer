library(ggplot2)
library(cowplot)
library(openxlsx)
library(reshape2)
library(forestplot)
library(viridis)
theme_set(theme_cowplot())

cancer_to_comp <- c("bladder", "breast", "colorectum", "endometrium",
                    "kidney", "lung", "pancreas",
                    "stomach", "Lymphocytic_leukemia", "melanoma",
                    "Non-Hodgkins_Lymphoma", "prostate")
ntree <- 200
tort <- "ass"

#interesting interactions: PGS with X

get_new_labels <- function(old_labels, the_decoder = feat_decoder){
  old_labels <- as.character(old_labels)
  new_labels <- rep("", length(old_labels))
  if(!all(old_labels %in% the_decoder[,1])){
    print(old_labels[!(old_labels %in% the_decoder[,1])])
    print("missing feat decoder")
    exit()
  }
  for(i in 1:length(old_labels)){
    new_labels[i] <- the_decoder[the_decoder[,1] == old_labels[i], 2]
  }
  return(new_labels)
}

feat_decoder <- read.table("label_decoder.csv", stringsAsFactors = F, sep = ",")
feat_decoder[,1] <- trimws(feat_decoder[,1])
feat_decoder[,2] <- trimws(feat_decoder[,2])

cancer_decoder <- data.frame(old_names = cancer_to_comp,
                             new_names = c("Bladder", "Breast", "Colorectum", "Endometrium",
                                           "Kidney", "Lung", "Pancreas", "Stomach", "Lymphocytic Leukemia",
                                           "Melanoma", "Non-Hodgkins Lymphoma", "Prostate"),
                             stringsAsFactors = F)


###########################################################
get_base_stats <- FALSE
if(get_base_stats){
  all_data <- read.table("all_data.txt.gz", stringsAsFactors = F, header = T)
  all_columns <- as.data.frame(fread("all_data_columns.txt", sep = "\t", header = F))
  colnames(all_data) <- all_columns[,1]
  
  
  wanted_cols <- c("Year of Birth", "Sex", "weight", "Standing Height", "Time at current address",
   "Age completed full time education", "Number of days/week of moderate physical activity", 
   "Alcohol Frequency, Daily", "Alcohol Frequency, One to three times a month", "Alcohol Frequency, Never",
   "Smoking status, Never", "Smoking status, Previous", "Smoking status, Current", 
   "Had menopause", "Age when periods started (menarche)", "I10", "R10", "R07", "K29", "E78",
   "R69", "Age..All.usual.residents..Rural.Urban..Total..measures..Value")
  
  
  summary_data <- all_data[,colnames(all_data) %in% wanted_cols]
  summary_data$`Year of Birth` <- 2020 - summary_data$`Year of Birth` 
  summary_data <- summary_data[,order(colnames(summary_data))[rank(wanted_cols)]]
  
  colnames(summary_data) <- c("Age", "Sex", "Weight", "Height", "Time at Current Address", "Age Completed Education",
                              "No. Days per Week Exercise", "Alcohol Frequency - Daily", "Alcohol Frequency - 1-3 per Month",
                              "Alcohol Frequency - Never", "Smoking Status - Never", "Smoking Status - Previous",
                              "Smoking Status - Current", "Past Menopause", "Age of Menarche", "ICD - Hypertension",
                              "ICD - Abdominal/Pelvic Pain", "ICD - Pain in Throat", "ICD - Gastritis and duodenitis", 
                              "ICD - Disorders of lipoprotein metabolism", "ICD: Unknown causes of morbidity",
                              "Census - Age")
  summary_data$Age <- 2020 - summary_data$Age
  
  real_summary <- data.frame(matrix(0, nrow = ncol(summary_data), ncol = 5))
  for(i in 1:ncol(summary_data)){
    if(length(unique(summary_data[,i])) <= 2){
      real_summary[i,5] <- sum(summary_data[,i])
    } else {
      real_summary[i,1] <- signif(mean(summary_data[,i], na.rm=T),3)
      real_summary[i,2] <- signif(sd(summary_data[,i], na.rm=T),3)
      real_summary[i,3] <- signif(min(summary_data[,i], na.rm=T),3)
      real_summary[i,4] <- signif(max(summary_data[,i], na.rm=T),3)
    }
  }
  
  real_summary <- data.frame(colnames(summary_data), real_summary)
  write.table(real_summary, "summary_of_stats.txt", row.names = F, col.names = F, sep = "\t", quote = F)
}

############################################################3
## Feature Analysis

#HEAT MAP
list_feat <- list()
for(i in 1:length(cancer_to_comp)){
  list_feat[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/process_feat_imp.RDS"))
  list_feat[[i]]$pro_gain <- list_feat[[i]]$Gain/sum(list_feat[[i]]$Gain)
  list_feat[[i]]$cancer <- cancer_to_comp[i]
}
all_xgb_coef <- do.call("rbind", list_feat)
tab_feat_name <- table(unlist(lapply(list_feat, function(x) x[,1])))
good_names <- names(sort(tab_feat_name, decreasing = T)[1:10])

gain_holder <- matrix(NA, nrow = length(good_names), ncol = length(cancer_to_comp))
for(i in 1:length(cancer_to_comp)){
  for(j in 1:length(good_names)){
    if(good_names[j] %in% list_feat[[i]][,1]){
      gain_holder[j,i] <- list_feat[[i]]$pro_gain[list_feat[[i]][,1] == good_names[j]]
    }
  }
}

quick_fun <- function(x){x[is.na(x)] <- mean(x, na.rm = T); x}
new_gain_holder <- apply(gain_holder, 2, quick_fun)
row_order <- hclust(dist(new_gain_holder))$order
col_order <- hclust(dist(t(new_gain_holder)))$order

gain_holder <- as.data.frame(gain_holder)
colnames(gain_holder) <- cancer_to_comp
gain_holder$inter <- good_names

plot_df <- reshape2::melt(gain_holder, "inter")
plot_df$inter <- factor(plot_df$inter, levels = gain_holder$inter[row_order])
plot_df$variable <- factor(plot_df$variable, levels = colnames(gain_holder)[col_order])

the_plot <- ggplot(plot_df, aes(variable, inter, fill = value)) + geom_raster() +
 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
 scale_fill_viridis() +
 labs(x = "Cancer", y = "Feature", fill = "Gain\nFraction") +
 scale_y_discrete(labels = get_new_labels(levels(plot_df$inter))) +
 scale_x_discrete(labels = get_new_labels(levels(plot_df$variable), cancer_decoder))
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_feat_imp.png"), the_plot, "png", height = 7, width = 9)


#MEAN XGB COEFS
quick_fun <- function(x){x[is.na(x)] <- 0; y <- mean(x); y}
plot_df <- data.frame(feat = gain_holder$inter, 
                      mean_gain = apply(gain_holder[,1:length(cancer_to_comp)], 1, quick_fun ))
plot_df$feat <- factor(plot_df$feat, levels = plot_df$feat[order(plot_df$mean_gain)])
the_plot <- ggplot(plot_df, aes(mean_gain, feat)) + geom_point() +
 labs(x = "Mean Gain Proportion", y = "Feature") +
  scale_y_discrete(labels = get_new_labels(levels(plot_df$feat)))
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_mean_feat_imp.png"), the_plot, "png", height = 6, width = 7)


#LINEAR COEFS
all_linear_coefs <- matrix(0, nrow = length(good_names), ncol = length(cancer_to_comp))
for(i in 1:length(cancer_to_comp)){
  linear_coefs <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/linear_coef.RDS"))
  linear_coefs <- linear_coefs[!is.na(linear_coefs$Estimate),]
  linear_coefs$pro_est <- abs(linear_coefs$Estimate)/sum(abs(linear_coefs$Estimate))
  for(j in 1:length(good_names)){
    if(good_names[j] %in% linear_coefs$feat_imp_name){
      all_linear_coefs[j,i] <- linear_coefs$pro_est[linear_coefs$feat_imp_name == good_names[j]]
    }
  }
}

mean_linear_coef <- apply(all_linear_coefs, 1, mean)
plot_df$linear_mean_gain <- mean_linear_coef
plot_df <- reshape2::melt(plot_df, "feat")
the_plot <- ggplot(plot_df, aes(value, feat, color = variable)) + geom_point() +
 scale_color_discrete(labels = c("XGB", "Linear")) +
 labs(x = "Proportion of Quasi-Gain", y = "Feature", color = "Model") +
  scale_y_discrete(labels = get_new_labels(levels(plot_df$feat))) 
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_mean_quasi_imp.png"), the_plot, "png", height = 6, width = 8)


###########################################
###########################################
#Ten model features
list_ten_feat <- list()
for_table_ten_feat <- list()
for(i in 1:length(cancer_to_comp)){
  list_ten_feat[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/ten_feat_mod_coefs.RDS"))
  
  for_table_ten_feat[[i]] <- data.frame(get_new_labels(rownames(list_ten_feat[[i]][-1,])), list_ten_feat[[i]][-1,])
  for_table_ten_feat[[i]]$categ <- get_new_labels(rownames(list_ten_feat[[i]][-1,]), feat_decoder[,c(1,3)])
  for_table_ten_feat[[i]]$weight <- abs(for_table_ten_feat[[i]]$Estimate)/sum(abs(for_table_ten_feat[[i]]$Estimate))
  rownames(for_table_ten_feat[[i]]) <- NULL
  colnames(for_table_ten_feat[[i]]) <- c("Feature", "Coefficient", "Std. Error", "Z Value", "P Value", "Category", "Weight")
  for_table_ten_feat[[i]][,2] <- signif(for_table_ten_feat[[i]][,2], 3)
  for_table_ten_feat[[i]][,3] <- signif(for_table_ten_feat[[i]][,3], 3)
  for_table_ten_feat[[i]][,4] <- signif(for_table_ten_feat[[i]][,4], 3)
  for_table_ten_feat[[i]][,5] <- signif(for_table_ten_feat[[i]][,5], 3)
  for_table_ten_feat[[i]][,7] <- signif(for_table_ten_feat[[i]][,7], 3)
  for_table_ten_feat[[i]]$Cancer <- get_new_labels(cancer_to_comp[i], cancer_decoder)
}

tab_feat_name <- table(unlist(lapply(list_ten_feat, function(x) rownames(x[-1,]))))
good_names <- names(sort(tab_feat_name, decreasing = T)[1:10])

suppz_table_1 <- do.call("rbind", for_table_ten_feat)
suppz_table_1 <- suppz_table_1[,c(8,1,6,2,3,5,7)]
suppz_table_2 <- suppz_table_1[,c(1,2,4,5,6)]
suppz_table_3 <- suppz_table_1[,c(1,2,3,7)]
write.table(suppz_table_1, "supp_table_1.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(suppz_table_2, "supp_table_2.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(suppz_table_3, "supp_table_3.txt", row.names = F, col.names = T, quote = F, sep = "\t")

beta_holder <- matrix(NA, nrow = length(good_names), ncol = length(cancer_to_comp))
for(i in 1:length(cancer_to_comp)){
  for(j in 1:length(good_names)){
    if(good_names[j] %in% rownames(list_ten_feat[[i]])){
      beta_holder[j,i] <- list_ten_feat[[i]][rownames(list_ten_feat[[i]]) == good_names[j], 1]
    }
  }
}

for(i in 1:length(cancer_to_comp)){
  plot_df <- as.data.frame(list_ten_feat[[i]][-1,])
  colnames(plot_df) <- c("beta", "se", "Z", "p")
  plot_df$inter <- factor(rownames(plot_df), levels = rownames(plot_df)[order(plot_df$beta)])
  the_plot <- ggplot(plot_df, aes(beta, inter)) + geom_point() +
    scale_y_discrete(labels = get_new_labels(levels(plot_df$inter))) +
    geom_errorbarh(aes(xmin = beta - se, xmax = beta + se), height = 0.2) +
    geom_vline(aes(xintercept=0)) +
    labs(x = paste0("Coefficient: ", cancer_to_comp[i]), y = "Feature")
  plot(the_plot)
}

quick_fun <- function(x){x[is.na(x)] <- mean(x, na.rm = T); x}
new_beta_holder <- apply(beta_holder, 2, quick_fun)
row_order <- hclust(dist(new_beta_holder))$order
col_order <- hclust(dist(t(new_beta_holder)))$order

beta_holder <- as.data.frame(beta_holder)
colnames(beta_holder) <- cancer_to_comp
beta_holder$inter <- good_names

plot_df <- reshape2::melt(beta_holder, "inter")
plot_df$inter <- factor(plot_df$inter, levels = beta_holder$inter[row_order])
plot_df$variable <- factor(plot_df$variable, levels = colnames(beta_holder)[col_order])

the_plot <- ggplot(plot_df, aes(variable, inter, fill = value)) + geom_raster() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_viridis() +
  labs(x = "Cancer", y = "Feature", fill = "Coefficient") +
  scale_y_discrete(labels = get_new_labels(levels(plot_df$inter))) +
  scale_x_discrete(labels = get_new_labels(levels(plot_df$variable), cancer_decoder))
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_linear_feat_heat.png"), the_plot, "png", height = 6, width = 8)

######

list_ten_feat <- list()
for(i in 1:length(cancer_to_comp)){
  list_ten_feat[[i]] <- as.data.frame(readRDS(paste0(tort, "_results_", ntree, "/",
                                                     cancer_to_comp[i], "/ten_feat_mod_coefs.RDS")))[-1,]
  list_ten_feat[[i]]$cancer <- cancer_to_comp[i]
  list_ten_feat[[i]]$feat <- rownames(list_ten_feat[[i]])
  rownames(list_ten_feat[[i]]) <- NULL
}
plot_df <- do.call("rbind", list_ten_feat)
all_tenfeat_coef <- plot_df
plot_df$feat <- get_new_labels(plot_df$feat)
plot_df$type <- ""
colnames(plot_df)[1:4] <- c("beta", "se", "z", "p")

for(i in 1:nrow(plot_df)){
  plot_df$type[i] <- feat_decoder[feat_decoder[,2] == plot_df$feat[i],3]
}

plot_df$cancer <- factor(plot_df$cancer, levels = cancer_to_comp)
plot_df$type <- factor(plot_df$type, levels = unique(feat_decoder$V3))

the_plot <- ggplot(plot_df, aes(cancer, beta, color = type)) +
   geom_point(size = 2.1, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
  geom_line(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1:length(cancer_to_comp) + 0.5, color = "grey") +
  labs(x = "Cancer", y = "Coefficient", color = "Type") +
  scale_x_discrete(labels = get_new_labels(levels(plot_df$cancer), cancer_decoder)) +
  coord_flip()
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/linear_feat_fig1.png"), the_plot, "png", height = 8, width = 6)
  

plot_df$weight <- 0
for(x in cancer_to_comp){
  plot_df$weight[plot_df$cancer == x] <- abs(plot_df$beta[plot_df$cancer == x])/sum(abs(plot_df$beta[plot_df$cancer == x]))
}

feat_measure_prop_statz <- sum(plot_df$weight[plot_df$type == "Measure"])/sum(plot_df$weight)
feat_icd_prop_statz <- sum(plot_df$weight[plot_df$type == "ICD"])/sum(plot_df$weight)
feat_census_prop_statz <- sum(plot_df$weight[plot_df$type == "Census"])/sum(plot_df$weight)
feat_answer_prop_statz <- sum(plot_df$weight[plot_df$type == "Answer"])/sum(plot_df$weight)
feat_pgs_prop_statz <- sum(plot_df$weight[plot_df$type == "PGS"])/sum(plot_df$weight)
feat_biomarker_prop_statz <- sum(plot_df$weight[plot_df$type == "Biomarker"])/sum(plot_df$weight)

happy_own_health_total_statz <- sum(plot_df[plot_df$feat == "Happiness w/ Own Health",]$weight)/sum(plot_df$weight)
happy_own_health_pancreas_statz <- plot_df[plot_df$weight == max(plot_df[plot_df$feat == "Happiness w/ Own Health",]$weight),]

###
all_xgb_coef <- all_xgb_coef[,c("Interaction", "cancer", "pro_gain")]
all_tenfeat_coef <- all_tenfeat_coef[,c("feat", "cancer", "Estimate")]
for(cancer in cancer_to_comp){
  all_tenfeat_coef$Estimate[all_tenfeat_coef$cancer == cancer] <- 
    abs(all_tenfeat_coef$Estimate[all_tenfeat_coef$cancer == cancer])/
    sum(abs(all_tenfeat_coef$Estimate[all_tenfeat_coef$cancer == cancer]))
}
colnames(all_tenfeat_coef) <- c("Interaction", "cancer", "pro_gain")

all_tenfeat_coef$type <- "tenfeat"
all_xgb_coef$type <- "xgb"
all_coef <- rbind(all_tenfeat_coef, all_xgb_coef)
tab_coef <- table(all_coef$Interaction)
keep_names <- names(head(sort(tab_coef, decreasing = T), 10))
all_coef <- all_coef[all_coef$Interaction %in% keep_names,]
all_coef$Interaction <- as.factor(all_coef$Interaction)

the_plot <- ggplot(all_coef, aes(pro_gain, Interaction)) + #geom_boxplot() +
  geom_point(aes(color = type), position = position_jitterdodge(jitter.width = 0)) +
  labs(x = "Proportion of Effect on Model", y = "Feature", color = "Model") +
  scale_y_discrete(labels = get_new_labels(levels(all_coef$Interaction))) +
  scale_color_discrete(labels = c("Linear", "XGB"))
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/ogcomp_feat_fig1.png"), the_plot, "png", height = 5.5, width = 6.5)


##################################################################
## Interactions

list_pval <- list()
for(i in 1:length(cancer_to_comp)){
  list_pval[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/process_inter_coef.RDS"))
  list_pval[[i]]$cancer <- cancer_to_comp[i]
  list_pval[[i]] <- list_pval[[i]][list_pval[[i]]$p_val < 0.0001 & list_pval[[i]]$p_val != 0 & !is.na(list_pval[[i]]$p_val),]
}
inter_pval <- do.call("rbind", list_pval)
inter_pval <- inter_pval[,c(1,17,18,19,20:25)]
colnames(inter_pval)[6:9] <- c("beta", "se", "z", "p")
inter_pval <- inter_pval[order(inter_pval$p_val),]

inter_pval$cancer <- factor(inter_pval$cancer, levels = cancer_to_comp)
while(any(duplicated(inter_pval$Interaction))){
inter_pval$Interaction[duplicated(inter_pval$Interaction)] <- 
  paste0(inter_pval$Interaction[duplicated(inter_pval$Interaction)], "2")
}
inter_pval$Interaction <- factor(inter_pval$Interaction, levels = inter_pval$Interaction[order(inter_pval$p_val)])
the_plot <- ggplot(inter_pval, aes(Interaction, -log10(p_val), color = cancer)) + geom_point() +
  theme(axis.text.x=element_blank()) +
  labs(x = "Interaction", y = "-log10(P-Value)", color = "Cancer")
plot(the_plot)

inter_pval$read_feat_1 <- get_new_labels(inter_pval$feat1)
inter_pval$read_feat_2 <- get_new_labels(inter_pval$feat2)
suppz_table_5 <- inter_pval[,c(10,11,12,6,7,9)]
colnames(suppz_table_5) <- c("Cancer", "Feature - 1", "Feature - 2", "Coefficient", "Std. Error", "P Value")
suppz_table_5$Cancer <- get_new_labels(suppz_table_5$Cancer, cancer_decoder)
suppz_table_5[,4] <- signif(suppz_table_5[,4], 3)
suppz_table_5[,5] <- signif(suppz_table_5[,5], 3)
suppz_table_5[,6] <- signif(suppz_table_5[,6], 3)

write.table(suppz_table_5, "supp_table_5.txt", row.names = F, col.names = T, quote = F, sep = "\t")

exit()
###
#now with forest plot

inter_pval <- inter_pval[inter_pval$feat1 != inter_pval$feat2,]
plot_df <- inter_pval[abs(inter_pval$beta) > sort(abs(inter_pval$beta), decreasing = T)[28] &
                        inter_pval$p_val < sort(inter_pval$p_val)[28],]

feat1 <- c("Feature", as.list(paste0(get_new_labels(plot_df$feat1), ":", get_new_labels(plot_df$feat2))))
cancer <- c("Cancer", as.list(get_new_labels(as.character(plot_df$cancer), cancer_decoder)))
beta <- c("OR", as.list(signif(exp(plot_df$beta), 2)))
pval <- c("P-value", as.list(signif(plot_df$p, 2)))
plot_list <- list(feat1, cancer, beta, pval)
plot_mat <- cbind(c(NA, plot_df$beta - plot_df$se), c(NA, plot_df$beta), c(NA, plot_df$beta + plot_df$se))
bool_summ <- rep(FALSE, nrow(plot_df)+1)
bool_summ[1] <- TRUE
png(paste0(tort, "_results_", ntree, "/complinxgb_feat_fig1.png"), width = 1000, height = 400)
forestplot::forestplot(plot_list, plot_mat,
           zero = 0,
           boxsize = 0.3,
           is.summary = bool_summ,
           txt_gp = fpTxtGp(xlab = gpar(cex = 2.5), ticks = gpar(cex = 1.5), label = gpar(cex = 1.2)),
           lwd.xaxis = 2, 
           lwd.zero = 2,
           lwd.ci = 3,
           col = fpColors(zero = "black", lines = "black"),
           graphwidth = unit(8, "cm"),
           colgap = unit(0.01, "npc"))
dev.off()

####################
#Now looking at the simulations


list_pval <- list()
for(i in 1:length(cancer_to_comp)){
  list_pval[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/all_pval.RDS"))
  list_pval[[i]]$cancer <- cancer_to_comp[i]
}

plot_df <- do.call("rbind", list_pval)
plot_df$cancer <- as.factor(plot_df$cancer)
plot_df$value[plot_df$value < 1e-15] <- 1e-15
the_plot <- ggplot(plot_df, aes(cancer,-log10(value),  fill = variable)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(labels = get_new_labels(levels(plot_df$cancer), cancer_decoder)) +
  scale_fill_discrete(labels = c("Random", "From Top 100", "From Top 10", "From XGB")) +
  labs(x = "Cancer", y = "-log10(P-Value)", fill = "Features\nConsidered")
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_inter_pval.png"), the_plot, "png", height = 5, width = 7)

the_plot <- ggplot(plot_df, aes(variable, -log10(value))) + geom_boxplot(notch = T) +
  scale_x_discrete(labels = c("Random", "From Top 100", "From Top 10", "From XGB")) +
  labs(x = "Features Considered", y = "-log10(P-Value)")
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/compile_inter_pval.png"), the_plot, "png", height = 4, width = 6)

interaction_sim_statz <- wilcox.test(plot_df$value[plot_df$variable == "top_p"], plot_df$value[plot_df$variable == "my_p"])



##################################################################
##################################################################
##################################################################
## Accuracy Metrics

#AUC ###################################33
#Start with overal model comparisons
#These are all the ML models with all features
to_supp_table <- matrix(0, nrow = 12, ncol = 5)
list_auc <- list()
for(i in 1:length(cancer_to_comp)){
  list_auc[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/auc_all_ml.RDS"))
  list_auc[[i]]$cancer <- cancer_to_comp[i]
  to_supp_table[i,] <- paste(signif(list_auc[[i]]$auc,3), paste0("(", signif((list_auc[[i]]$auc - list_auc[[i]]$lo),3), ")"), sep = " ")
}
rownames(to_supp_table) <- get_new_labels(cancer_to_comp, cancer_decoder)
colnames(to_supp_table) <- list_auc[[1]]$names
write.table(to_supp_table, "supp_table_7.txt", col.names = T, row.names = T, quote = F, sep = "\t")


plot_df <- do.call("rbind", list_auc)
rank_cancer <- plot_df$cancer[plot_df$names == "XGB"][order(plot_df$auc[plot_df$names == "XGB"])]
plot_df$cancer <- factor(plot_df$cancer, levels = rank_cancer)

the_plot <- ggplot(plot_df, aes(auc, cancer, color = names)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
  geom_errorbarh(position = position_jitterdodge(jitter.width = 0, jitter.height = 0), aes(xmin = lo, xmax = hi), height = 0.2) +
  labs(x = "AUC", y = "Cancer") +
  scale_y_discrete(labels = get_new_labels(levels(plot_df$cancer), cancer_decoder)) +
  geom_hline(yintercept = 1:11+0.5)
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_ml_auc.png"), the_plot, "png", height = 4, width = 6)

the_plot <- ggplot(plot_df, aes(names, auc)) + geom_boxplot(notch = F) + geom_point() +
  labs(x = "Model", y = "AUC")
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_ml_auc_box.png"), the_plot, "png", height = 4, width = 6)

mean_knn_statz <- mean(plot_df$auc[plot_df$names == "KNN"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "KNN"])/sqrt(sum(plot_df$names == "KNN"))
mean_knn_statz <- signif(c(mean_knn_statz - temp_ci, mean_knn_statz, mean_knn_statz+temp_ci), 3)

mean_svm_statz <- mean(plot_df$auc[plot_df$names == "SVM"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "SVM"])/sqrt(sum(plot_df$names == "SVM"))
mean_svm_statz <- signif(c(mean_svm_statz - temp_ci, mean_svm_statz, mean_svm_statz+temp_ci), 3)


###################################################################################
# JUST 10 but still ML

list_auc <- list()
to_supp_table <- matrix(0, nrow = 12, ncol = 5)
for(i in 1:length(cancer_to_comp)){
  list_auc[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/auc_all_ten_ml.RDS"))
  list_auc[[i]]$cancer <- cancer_to_comp[i]
  to_supp_table[i,] <- paste(signif(list_auc[[i]]$auc,3), paste0("(", signif((list_auc[[i]]$auc - list_auc[[i]]$lo),3), ")"), sep = " ")
}
rownames(to_supp_table) <- get_new_labels(cancer_to_comp, cancer_decoder)
colnames(to_supp_table) <- list_auc[[1]]$names
write.table(to_supp_table, "supp_table_8.txt", col.names = T, row.names = T, quote = F, sep = "\t")

plot_df <- do.call("rbind", list_auc)
rank_cancer <- plot_df$cancer[plot_df$names == "XGB"][order(plot_df$auc[plot_df$names == "XGB"])]
plot_df$cancer <- factor(plot_df$cancer, levels = rank_cancer)

the_plot <- ggplot(plot_df, aes(auc, cancer, color = names)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
  geom_errorbarh(position = position_jitterdodge(jitter.width = 0, jitter.height = 0), aes(xmin = lo, xmax = hi), height = 0.2) +
  labs(x = "AUC", y = "Cancer") +
  scale_y_discrete(labels = get_new_labels(levels(plot_df$cancer), cancer_decoder)) +
  geom_hline(yintercept = 1:11+0.5)
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_ten_ml_auc.png"), the_plot, "png", height = 4, width = 6)

plot_df$names <- factor(as.character(plot_df$names), levels = c("KNN", "XGB", "SVM", "RF", "Linear: Ten Feat"))
the_plot <- ggplot(plot_df, aes(names, auc)) + geom_boxplot(notch = F) + geom_point() +
  labs(x = "Model", y = "AUC")
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_ten_ml_auc_box.png"), the_plot, "png", height = 4, width = 6)

mean_xgb_statz <- mean(plot_df$auc[plot_df$names == "XGB"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "XGB"])/sqrt(sum(plot_df$names == "XGB"))
mean_xgb_statz <- signif(c(mean_xgb_statz - temp_ci, mean_xgb_statz, mean_xgb_statz+temp_ci), 3)

mean_rf_statz <- mean(plot_df$auc[plot_df$names == "RF"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "RF"])/sqrt(sum(plot_df$names == "RF"))
mean_rf_statz <- signif(c(mean_rf_statz - temp_ci, mean_rf_statz, mean_rf_statz+temp_ci), 3)

mean_knn_statz <- mean(plot_df$auc[plot_df$names == "KNN"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "KNN"])/sqrt(sum(plot_df$names == "KNN"))
mean_knn_statz <- signif(c(mean_knn_statz - temp_ci, mean_knn_statz, mean_knn_statz+temp_ci), 3)

mean_svm_statz <- mean(plot_df$auc[plot_df$names == "SVM"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "SVM"])/sqrt(sum(plot_df$names == "SVM"))
mean_svm_statz <- signif(c(mean_svm_statz - temp_ci, mean_svm_statz, mean_svm_statz+temp_ci), 3)


#########################################################################
#######################$$$$$$$$$$$$$$$$$$$$$
#Includes Linear, Random Forest, XGB (kind of unnecessary now)

list_auc <- list()
for(i in 1:length(cancer_to_comp)){
  list_auc[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/auc_all.RDS"))
  list_auc[[i]]$cancer <- cancer_to_comp[i]
}

plot_df <- do.call("rbind", list_auc)
sub_plot_df <- plot_df[plot_df$names == "XGB",]
rf_aucs <- plot_df[plot_df$names == "Random\nForest",]
plot_df$cancer <- factor(plot_df$cancer, levels = sub_plot_df$cancer[order(sub_plot_df$auc)])

mean_rf_statz <- mean(plot_df$auc[plot_df$names == "Random\nForest"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "Random\nForest"])/sqrt(sum(plot_df$names == "Random\nForest"))
mean_rf_statz <- signif(c(mean_rf_statz - temp_ci, mean_rf_statz, mean_rf_statz+temp_ci), 3)

mean_xgb_statz <- mean(plot_df$auc[plot_df$names == "XGB"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "XGB"])/sqrt(sum(plot_df$names == "XGB"))
mean_xgb_statz <- signif(c(mean_xgb_statz - temp_ci, mean_xgb_statz, mean_xgb_statz+temp_ci), 3)

mean_linear_statz <- mean(plot_df$auc[plot_df$names == "Linear"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "Linear"])/sqrt(sum(plot_df$names == "Linear"))
mean_linear_statz <- signif(c(mean_linear_statz - temp_ci, mean_linear_statz, mean_linear_statz+temp_ci), 3)
                       

the_plot <- ggplot(plot_df, aes(auc, cancer, color = names)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
 geom_errorbarh(position = position_jitterdodge(jitter.width = 0, jitter.height = 0), aes(xmin = lo, xmax = hi), height = 0.2) +
 labs(x = "AUC", y = "Cancer") +
  scale_y_discrete(labels = get_new_labels(levels(plot_df$cancer), cancer_decoder)) +
  geom_hline(yintercept = 1:11+0.5)
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_auc.png"), the_plot, "png", height = 4, width = 6)

the_plot <- ggplot(plot_df, aes(names, auc)) + geom_boxplot(notch = F) + geom_point() +
  labs(x = "Model", y = "AUC")
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_auc_box.png"), the_plot, "png", height = 4, width = 6)


################################### Now doing finer AUC comparisons
#Between XGB with and without PGS

list_auc <- list()
for(i in 1:length(cancer_to_comp)){
  list_auc[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/auc_pgs.RDS"))
  list_auc[[i]]$cancer <- cancer_to_comp[i]
}

plot_df <- do.call("rbind", list_auc)
sub_plot_df <- plot_df[plot_df$names == "XGB: w/ PGS",]
plot_df$cancer <- factor(plot_df$cancer, levels = sub_plot_df$cancer[order(sub_plot_df$auc)])

the_plot <- ggplot(plot_df, aes(auc, cancer, color = names)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
  geom_errorbarh(position = position_jitterdodge(jitter.width = 0, jitter.height = 0), aes(xmin = lo, xmax = hi), height = 0.2) +
  labs(x = "AUC", y = "Cancer", color = "Model") +
  scale_y_discrete(labels = get_new_labels(levels(plot_df$cancer), cancer_decoder)) +
  geom_hline(yintercept = 1:11+0.5)
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/pgs_auc.png"), the_plot, "png", height = 4, width = 6)

############################################
#Linear with and without interactions

list_auc <- list()
for(i in 1:length(cancer_to_comp)){
  list_auc[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/auc_inter.RDS"))
  list_auc[[i]]$cancer <- cancer_to_comp[i]
}

plot_df <- do.call("rbind", list_auc)
plot_df <- plot_df[plot_df$names != "XGB",]
sub_plot_df <- plot_df[plot_df$names == "Linear: w/ interaction",]
plot_df$cancer <- factor(plot_df$cancer, levels = sub_plot_df$cancer[order(sub_plot_df$auc)])

the_plot <- ggplot(plot_df, aes(auc, cancer, color = names)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
  geom_errorbarh(position = position_jitterdodge(jitter.width = 0, jitter.height = 0), aes(xmin = lo, xmax = hi), height = 0.2) +
  labs(x = "AUC", y = "Cancer", color = "Model") +
  scale_y_discrete(labels = get_new_labels(levels(plot_df$cancer), cancer_decoder)) +
  geom_hline(yintercept = 1:11+0.5)
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/inter_auc.png"), the_plot, "png", height = 4, width = 6)


################################################3
# Model with all linear comparisons

list_auc <- list()
to_supp_table <- matrix(0, nrow = 12, ncol = 6)
for(i in 1:length(cancer_to_comp)){
  list_auc[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/auc_justlinear.RDS"))
  list_auc[[i]]$cancer <- cancer_to_comp[i]
  to_supp_table[i,] <- paste(signif(list_auc[[i]]$auc,3), paste0("(", signif((list_auc[[i]]$auc - list_auc[[i]]$lo),3), ")"), sep = " ")
}
rownames(to_supp_table) <- get_new_labels(cancer_to_comp, cancer_decoder)
colnames(to_supp_table) <- list_auc[[1]]$names
write.table(to_supp_table, "supp_table_9.txt", col.names = T, row.names = T, quote = F, sep = "\t")



plot_df <- do.call("rbind", list_auc)
plot_df <- plot_df[-which(plot_df$names %in% c("Inter: Five by 2")),]
sub_plot_df <- plot_df[plot_df$names == "Linear: All Feat",]
plot_df$cancer <- factor(plot_df$cancer, levels = sub_plot_df$cancer[order(sub_plot_df$auc)])



the_plot <- ggplot(plot_df, aes(auc, cancer, color = names)) + 
  geom_errorbarh(position = position_jitterdodge(jitter.width = 0, jitter.height = 0), aes(xmin = lo, xmax = hi), height = 0.2) +
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
  geom_hline(yintercept = 1:11+0.5) +
  labs(x = "AUC", y = "Cancer", color = "Model Type") +
  scale_y_discrete(labels =  get_new_labels(levels(plot_df$cancer), cancer_decoder))

plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/inter_justlinear_box.png"), the_plot, "png", height = 4, width = 6)

mean_tenfeat_statz <- mean(plot_df$auc[plot_df$names == "Linear: Ten Feat"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "Linear: Ten Feat"])/sqrt(sum(plot_df$names == "Linear: Ten Feat"))
mean_tenfeat_statz <- signif(c(mean_tenfeat_statz - temp_ci, mean_tenfeat_statz, mean_tenfeat_statz+temp_ci), 3)

mean_fivefeat_statz <- mean(plot_df$auc[plot_df$names == "Linear: Five Feat"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "Linear: Five Feat"])/sqrt(sum(plot_df$names == "Linear: Five Feat"))
mean_fivefeat_statz <- signif(c(mean_fivefeat_statz - temp_ci, mean_fivefeat_statz, mean_fivefeat_statz+temp_ci), 3)

mean_threefeat_statz <- mean(plot_df$auc[plot_df$names == "Linear: Three Feat"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "Linear: Three Feat"])/sqrt(sum(plot_df$names == "Linear: Three Feat"))
mean_threefeat_statz <- signif(c(mean_threefeat_statz - temp_ci, mean_threefeat_statz, mean_threefeat_statz+temp_ci), 3)

three_ten_p_statz <- wilcox.test(plot_df$auc[plot_df$names == "Linear: Three Feat"],
                           plot_df$auc[plot_df$names == "Linear: Ten Feat"], paired = T)$p.value

sub_linear_ten <- plot_df[plot_df$names == "Linear: Ten Feat",]
lin_ten_min_auc_statz <- sub_linear_ten[which.min(sub_linear_ten$auc),]
lin_ten_max_auc_statz <- sub_linear_ten[which.max(sub_linear_ten$auc),]
lin_ten_best_3_statz <- sub_linear_ten[sub_linear_ten$auc >= sort(sub_linear_ten$auc, decreasing = T)[3],]


###############################################################
#Now including the external predictions (along with Linear, XGB, Linear Ten Feat)

list_auc <- list()
list_auc2 <- list()
for(i in 1:length(cancer_to_comp)){
  list_auc[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/auc_w_external.RDS"))
  list_auc[[i]]$cancer <- cancer_to_comp[i]
  list_auc2[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/auc_justlinear.RDS"))
  list_auc2[[i]]$cancer <- cancer_to_comp[i]
}

plot_df <- rbind(do.call("rbind", list_auc), do.call("rbind", list_auc2))
plot_df$names <- as.character(plot_df$names)
plot_df$names[plot_df$names == "Q-Pred"] <- "QCancer"
plot_df <- plot_df[plot_df$names %in% c("Linear", "XGB", "QCancer", "Linear: Ten Feat"),]
sub_plot_df <- plot_df[plot_df$names == "Linear",]
plot_df$cancer <- factor(plot_df$cancer, levels = sub_plot_df$cancer[order(sub_plot_df$auc)])
plot_df$names <- factor(plot_df$names, levels = c("QCancer", "Linear", "Linear: Ten Feat", "XGB"))

the_plot <- ggplot(plot_df, aes(auc, cancer, color = names)) +
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2,
                 position = position_jitterdodge(jitter.width = 0, jitter.height = 0)) +
  labs(x = "AUC", y = "Cancer", color = "Model:") +
  scale_y_discrete(labels = get_new_labels(levels(plot_df$cancer), cancer_decoder)) +
  geom_hline(yintercept = 1:11+0.5, color = "grey") +
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.position = "bottom") +
  guides(color=guide_legend(ncol=2))
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/external_auc.png"), the_plot, "png", height = 8, width = 6)


the_plot <- ggplot(plot_df, aes(names, auc)) + geom_boxplot(notch = F) + geom_point() +
  labs(x = "Model", y = "AUC")
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/external_box_auc.png"), the_plot, "png", height = 4, width = 6)


mean_qpred_statz <- mean(plot_df$auc[plot_df$names == "QCancer"])
temp_ci <- 1.96*sd(plot_df$auc[plot_df$names == "QCancer"])/sqrt(sum(plot_df$names == "QCancer"))
mean_qpred_statz <- signif(c(mean_qpred_statz - temp_ci, mean_qpred_statz, mean_qpred_statz+temp_ci), 3)

ll_tenfeat_statz <- plot_df[plot_df$cancer == "Lymphocytic_leukemia" & plot_df$names == "Linear: Ten Feat",]
ll_qpred_statz <- plot_df[plot_df$cancer == "Lymphocytic_leukemia" & plot_df$names == "QCancer",]

temp <- plot_df[plot_df$cancer %in% plot_df$cancer[plot_df$names == "QCancer"],]
abstrat_p1_statz <- wilcox.test(temp$auc[temp$names == "XGB"], temp$auc[temp$names == "Linear: Ten Feat"], paired = T)$p.value
abstrat_p2_statz <- wilcox.test(temp$auc[temp$names == "QCancer"], temp$auc[temp$names == "Linear: Ten Feat"], paired = T)$p.value
abstrat_p3_statz <- wilcox.test(plot_df$auc[plot_df$names == "XGB"], plot_df$auc[plot_df$names == "Linear: Ten Feat"], paired = T)$p.value

for_table_df <- rbind(plot_df, rf_aucs)
final_table_df <- data.frame(matrix(0, nrow = length(unique(for_table_df$cancer)), ncol = 6))
colnames(final_table_df) <- c("cancer", as.character(unique(for_table_df$names)))
final_table_df[,1] <- as.character(unique(for_table_df$cancer))
for(i in 1:nrow(final_table_df)){
  for(j in 2:ncol(final_table_df)){
    auc_val <- for_table_df$auc[for_table_df$names == colnames(final_table_df)[j] & for_table_df$cancer == final_table_df[i,1]]
    ci_val <- for_table_df$hi[for_table_df$names == colnames(final_table_df)[j] & for_table_df$cancer == final_table_df[i,1]] - auc_val
    final_table_df[i,j] <- paste0(signif(auc_val, 3), " (", signif(ci_val, 3), ")")
  }
}

final_table_df$cancer <- get_new_labels(final_table_df$cancer, cancer_decoder)
write.table(final_table_df, "supp_table_6.txt", row.names = F, col.names = T, sep = "\t", quote = F)

###########################################################################
##################################### COMPARE VARIOUS ROC CURVES
i <- 1
roc_objs <- list()
roc_lin_to_xgb <- matrix(0, nrow = length(cancer_to_comp), ncol = 2)
roc_rf_to_xgb <- matrix(0, nrow = length(cancer_to_comp), ncol = 2)
for(name_cancer in cancer_to_comp){
  roc_objs[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", name_cancer, "/roc_objs.1.RDS"))
  temp <- roc.test(roc_objs[[i]]$linear, roc_objs[[i]]$xgb) #positive z-stat when the first auc > second auc
  roc_lin_to_xgb[i,1] <- temp$p.value
  roc_lin_to_xgb[i,2] <- unname(temp$statistic)
  temp <- roc.test(roc_objs[[i]]$rf, roc_objs[[i]]$xgb)
  roc_rf_to_xgb[i,1] <- temp$p.value
  roc_rf_to_xgb[i,2] <- unname(temp$statistic)
  i <- i + 1
}

roc_p_rf_xgb_statz <- sum(roc_rf_to_xgb[,1] < 0.05)
roc_p_lin_xgb_statz <- sum(roc_lin_to_xgb[,1] < 0.05)


i <- 1
roc_objs <- list()
other_objs <- list()
roc_lin_to_ten <- matrix(0, nrow = length(cancer_to_comp), ncol = 2)
roc_ten_to_five <- matrix(0, nrow = length(cancer_to_comp), ncol = 2)
roc_ten_to_three <- matrix(0, nrow = length(cancer_to_comp), ncol = 2)
roc_ten_to_qpred <- matrix(0, nrow = length(cancer_to_comp), ncol = 2)
for(name_cancer in cancer_to_comp){
  roc_objs[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", name_cancer, "/roc_objs.2.RDS"))
  if(file.exists(paste0(tort, "_results_", ntree, "/", name_cancer, "/rocobj_w_external.RDS"))){
    other_objs[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", name_cancer, "/rocobj_w_external.RDS"))
    external_good <- TRUE
  } else {
    other_objs[[i]] <- NULL
    external_good <- FALSE
  }
  temp <- roc.test(roc_objs[[i]]$linear, roc_objs[[i]]$ten)
  roc_lin_to_ten[i,1] <- temp$p.value
  roc_lin_to_ten[i,2] <- unname(temp$statistic)
  temp <- roc.test(roc_objs[[i]]$ten, roc_objs[[i]]$five)
  roc_ten_to_five[i,1] <- temp$p.value
  roc_ten_to_five[i,2] <- unname(temp$statistic)
  temp <- roc.test(roc_objs[[i]]$ten, roc_objs[[i]]$three)
  roc_ten_to_three[i,1] <- temp$p.value
  roc_ten_to_three[i,2] <- unname(temp$statistic)
  if(external_good){
  temp <- roc.test(roc_objs[[i]]$ten, other_objs[[i]])
  roc_ten_to_qpred[i,1] <- temp$p.value
  roc_ten_to_qpred[i,2] <- unname(temp$statistic)
  }
  i <- i + 1
}

roc_p_lin_ten_1_statz <- sum(roc_lin_to_ten[,1] > 0.05)
roc_p_lin_ten_2_statz <- sum(roc_lin_to_ten[,1] < 0.05 & roc_lin_to_ten[,2] > 0)
roc_p_lin_ten_3_statz <- sum(roc_lin_to_ten[,1] < 0.05 & roc_lin_to_ten[,2] < 0)
roc_p_ten_5_statz <- sum(roc_ten_to_five[,1] < 0.05)
roc_p_ten_3_statz <- sum(roc_ten_to_five[,1] < 0.05)


roc_comp_df <- data.frame(matrix(0, nrow = length(cancer_to_comp), ncol = 7))
roc_comp_df[,1] <- get_new_labels(cancer_to_comp, cancer_decoder)
roc_comp_df[,2] <- signif(roc_lin_to_xgb[,1] * sign(roc_lin_to_xgb[,2]), 3)
roc_comp_df[,3] <- signif(roc_rf_to_xgb[,1] * sign(roc_rf_to_xgb[,2]), 3)
roc_comp_df[,4] <- signif(roc_lin_to_ten[,1] * sign(roc_lin_to_ten[,2]), 3)
roc_comp_df[,5] <- signif(roc_ten_to_five[,1] * sign(roc_ten_to_five[,2]), 3)
roc_comp_df[,6] <- signif(roc_ten_to_three[,1] * sign(roc_ten_to_three[,2]), 3)
roc_comp_df[,7] <- signif(roc_ten_to_qpred[,1] * sign(roc_ten_to_qpred[,2]), 3)
colnames(roc_comp_df) <- c("Cancer", "Linear:XGB", "RF:XGB", "Linear:Ten", "Ten:Five", "Ten:Three", "Ten:Q-Pred")
write.table(roc_comp_df, "suppz_roc_comp_6.txt", row.names = F, col.names = T, quote = F, sep = "\t")

##################################################
i <- 1
all_acc <- matrix(0, nrow = length(cancer_to_comp), ncol = 2)
all_tpr <-matrix(0, nrow = length(cancer_to_comp), ncol = 2)
all_fpr <- matrix(0, nrow = length(cancer_to_comp), ncol = 2)

for(name_cancer in cancer_to_comp){
  if(file.exists(paste0(tort, "_results_", ntree, "/", name_cancer, "/qpred_stats.RDS"))){
    qstat <- readRDS(paste0(tort, "_results_", ntree, "/", name_cancer, "/qpred_stats.RDS"))
  } else {
    qstat <- rep(0, 4)
  }
  tenstat <- readRDS(paste0(tort, "_results_", ntree, "/", name_cancer, "/ten_feat_stats.RDS"))
  all_acc[i,] <- c(qstat[3], tenstat[3])
  all_tpr[i,] <- c(qstat[2], tenstat[2])
  all_fpr[i,] <- c(qstat[1], tenstat[1])
  i <- i +1
}

sub_acc <- all_acc[roc_ten_to_qpred[,1] < 0.05 & roc_ten_to_qpred[,1] != 0,]
acc_statz <- mean((sub_acc[,2] - sub_acc[1])/sub_acc[,1])

qstat_stat_df <- data.frame(get_new_labels(cancer_to_comp, cancer_decoder), 
                            signif(all_acc[,1], 3), signif(all_tpr[,1], 3), signif(all_fpr[,1], 3))
colnames(qstat_stat_df) <- c("Cancer", "Acc.", "TPR", "FPR")
write.table(qstat_stat_df, "supp_qstat_stat.txt", row.names = F, col.names = T, quote = F, sep = "\t")


# TPR, FPR, ACC #############################################################
list_stats <- list()
list_auc <- list()
for(i in 1:length(cancer_to_comp)){
  list_stats[[i]] <- as.data.frame(t(readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/ten_feat_stats.RDS"))))
  list_stats[[i]]$cancer <- cancer_to_comp[i]
  auc_val <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/auc_justlinear.RDS"))[4,]
  list_stats[[i]] <- cbind(list_stats[[i]], auc_val)
  
  valid_lab <- read.table(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/jupyter/",
                                 cancer_to_comp[i], "_valid_lab.txt"), stringsAsFactors = F, header = F) 
  train_lab <- read.table(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/jupyter/",
                                 cancer_to_comp[i], "_train_lab.txt"), stringsAsFactors = F, header = F)
  test_lab <- read.table(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/jupyter/",
                                cancer_to_comp[i], "_test_lab.txt"), stringsAsFactors = F, header = F)
  cc <- table(valid_lab) + table(train_lab) + table(test_lab)
  cc <- t(as.matrix(cc))
  colnames(cc) <- c("zero", "one")
  list_stats[[i]] <- cbind(list_stats[[i]], cc)
}

all_stats <- do.call("rbind", list_stats)
all_stats <- data.frame(cancer = all_stats$cancer, controls = all_stats$zero, cases = all_stats$one,
                        auc = signif(all_stats$auc, 3), auc_ci = signif(all_stats$auc - all_stats$lo, 3), acc = signif(all_stats$acc, 3),
                        fpr = signif(all_stats$fpr, 3), tpr = signif(all_stats$tpr, 3))
write.table(all_stats, "all_stats.txt", row.names = F, col.names = T, quote = F, sep = '\t')


# Odds Ratio ##########################################

list_or <- list()
for(i in 1:length(cancer_to_comp)){
 list_or[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/or_all.RDS"))
 list_or[[i]]$cancer <- cancer_to_comp[i]
}

all_or <- do.call("rbind", list_or)
all_or <- all_or[all_or$model %in% c("Linear", "XGB"),]
sub_or <- all_or[all_or$model == "XGB",]
all_or$cancer <- factor(all_or$cancer, levels = sub_or$cancer[order(sub_or$or)])

the_plot <- ggplot(all_or, aes(or, cancer, color = model)) +
geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
 geom_errorbarh(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5), aes(xmin = lo, xmax = hi), height = 0.2) +
 labs(x = "Odds Ratio", y = "Cancer", color = "Model\nType") +
  scale_y_discrete(labels = get_new_labels(levels(all_or$cancer), cancer_decoder)) +
  geom_hline(yintercept = 1:11+0.5)
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_ors.png"), the_plot, "png", height = 4, width = 6)


###

list_or <- list()
list_or2 <- list()
j <- 1
for(i in 1:length(cancer_to_comp)){
  if(file.exists(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/or_q_external.RDS"))){
    list_or[[j]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/or_q_external.RDS"))
    list_or[[j]]$cancer <- cancer_to_comp[i]
    j <- j + 1
  }
  list_or2[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/or_linear_ten.RDS"))
  list_or2[[i]]$cancer <- cancer_to_comp[i]
}

all_q_or <- do.call("rbind", list_or)
all_q_or$type <-"QCancer"
all_lin_or <- do.call("rbind", list_or2)
all_lin_or$type <- "Linear"
all_or <- rbind(all_q_or, all_lin_or)

or_write_df <- data.frame(matrix(0, nrow = length(unique(all_lin_or$cancer)), ncol = length(unique(all_lin_or$cut_off))+1))
for(i in 1:nrow(or_write_df)){
  for(j in 2:ncol(or_write_df)){
    or_write_df[i,j] <- signif(all_lin_or$or[all_lin_or$cancer == cancer_to_comp[i] & all_lin_or$cut_off == unique(all_lin_or$cut_off)[j-1]], 3)
  }
}
colnames(or_write_df) <- c("Cancer", paste("C", unique(all_lin_or$cut_off)))
or_write_df[,1] <- get_new_labels(cancer_to_comp, cancer_decoder)
write.table(or_write_df, "supp_ors_8.txt", row.names = F, col.names = T, quote = F, sep = "\t")


plot_spec <- function(cancer_name, title_name){
  spec_or <- all_or[all_or$cancer == cancer_name,]
  the_plot <- ggplot(spec_or, aes(cut_off, or, color = type)) + geom_point() + 
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
    labs(x = "Cut Off", y = "Odds Ratio", color = "Model", title = title_name)
  return(the_plot)
}

the_plot <- plot_spec("breast", "Breast")
ggsave(paste0(tort, "_results_", ntree, "/breast_ors.png"), the_plot, "png", height = 3.3, width = 5)
the_plot <- plot_spec("Non-Hodgkins_Lymphoma", "Non-Hodgkins Lymphoma")
ggsave(paste0(tort, "_results_", ntree, "/nhl_ors.png"), the_plot, "png", height = 3.3, width = 5)
the_plot <- plot_spec("prostate", "Prostate")
ggsave(paste0(tort, "_results_", ntree, "/prostate_ors.png"), the_plot, "png", height = 3.3, width = 5)


or_best_3_statz <- all_lin_or[all_lin_or$cut_off == 0.8 & all_lin_or$cancer %in% c("lung", "bladder", "Lymphocytic_leukemia"),]

sub_q_or <- all_q_or[all_lin_or$cut_off == 0.8,]
sub_q_or <- sub_q_or[!is.na(sub_q_or[,1]),]
sub_q_or <- sub_q_or[order(sub_q_or$cancer),]
sub_lin_or <- all_lin_or[all_lin_or$cut_off == 0.8,]
sub_lin_or <- sub_lin_or[sub_lin_or$cancer %in% sub_q_or$cancer,]
sub_lin_or <- sub_lin_or[order(sub_lin_or$cancer),]
temp <- sub_lin_or[,1] - sub_q_or[,1]
sub_lin_or[temp >= sort(temp, decreasing = T)[3],]

sub_q_or <- all_q_or[all_lin_or$cut_off == 0.99,]
sub_q_or <- sub_q_or[!is.na(sub_q_or[,1]),]
sub_q_or <- sub_q_or[order(sub_q_or$cancer),]
sub_q_or <- sub_q_or[sub_q_or$cancer %in% c("Non-Hodgkins_Lymphoma", "breast", "prostate"),]
sub_lin_or <- all_lin_or[all_lin_or$cut_off == 0.99,]
sub_lin_or <- sub_lin_or[sub_lin_or$cancer %in% sub_q_or$cancer,]
sub_lin_or <- sub_lin_or[order(sub_lin_or$cancer),]
sub_lin_or <- sub_lin_or[sub_lin_or$cancer %in% c("Non-Hodgkins_Lymphoma", "breast", "prostate"),]
fig_or_perc_diff_statz <- min((sub_lin_or$or - sub_q_or$or)/sub_q_or$or)
fig_or_min_statz <- min(sub_lin_or$or)

sub_lin_or <- all_lin_or[all_lin_or$cut_off == 0.99,]
sub_lin_auc <- plot_df[plot_df$names == "Linear: Ten Feat",]
sub_lin_or <- sub_lin_or[order(sub_lin_or$cancer),]
sub_lin_auc <- sub_lin_auc[order(sub_lin_auc$cancer),]
or_auc_corr_statz <- cor(sub_lin_or$or, sub_lin_auc$auc, method = "spearman")

# Cancer Time ############################################

list_time <- list()
case_list_time <- list()
for(i in 1:length(cancer_to_comp)){
  list_time[[i]] <- readRDS(paste0(tort, "_results_", ntree, "/", cancer_to_comp[i], "/cancer_time.RDS"))
  case_list_time[[i]] <- list_time[[i]][list_time[[i]]$pheno == 1,]
  
  list_time[[i]]$cancer <- cancer_to_comp[i]
  list_time[[i]]$time_group <- 0
  limit_1 <- quantile(list_time[[i]]$time, 1/3)
  limit_2 <- quantile(list_time[[i]]$time, 2/3)
  list_time[[i]]$time_group <- 3
  list_time[[i]]$time_group[list_time[[i]]$time < limit_1] <- 1
  list_time[[i]]$time_group[list_time[[i]]$time >= limit_1 & list_time[[i]]$time <= limit_2] <- 2
  
  ###################################
  
  case_list_time[[i]]$cancer <- cancer_to_comp[i]
  case_list_time[[i]]$time_group <- 0
  limit_1 <- quantile(case_list_time[[i]]$time, 1/3)
  limit_2 <- quantile(case_list_time[[i]]$time, 2/3)
  case_list_time[[i]]$time_group <- 3
  case_list_time[[i]]$time_group[case_list_time[[i]]$time < limit_1] <- 1
  case_list_time[[i]]$time_group[case_list_time[[i]]$time >= limit_1 & case_list_time[[i]]$time <= limit_2] <- 2
  
  case_list_time[[i]]$icd_time_group <- 0
  limit_1 <- quantile(case_list_time[[i]]$time_from_icd, 1/3)
  limit_2 <- quantile(case_list_time[[i]]$time_from_icd, 2/3)
  case_list_time[[i]]$icd_time_group <- 3
  case_list_time[[i]]$icd_time_group[case_list_time[[i]]$time_from_icd < limit_1] <- 1
  case_list_time[[i]]$icd_time_group[case_list_time[[i]]$time_from_icd >= limit_1 &
                                       case_list_time[[i]]$time_from_icd <= limit_2] <- 2
}

all_time <- do.call("rbind", list_time)
all_time$time_group <- as.factor(all_time$time_group)

case_all_time <- do.call("rbind", case_list_time)
case_all_time$time_group <- as.factor(case_all_time$time_group)
case_all_time$icd_time_group <- as.factor(case_all_time$icd_time_group)

#ggplot(case_all_time, aes(time/365, prob, color = cancer)) + geom_smooth(se=F) +
#  labs(x = "Study Time (Years)", y = "XGB Probability", color = "Cancer")
the_plot <- ggplot(all_time, aes(prob, cancer, fill = time_group)) + geom_boxplot() +
  labs(x = "XGB Probability", y = "Cancer", fill = "All Time\nGroup") +
  scale_y_discrete(labels = get_new_labels(unique(all_time$cancer), cancer_decoder))
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_prob_time.png"), the_plot, "png", height = 4, width = 6)


the_plot <- ggplot(case_all_time, aes(prob, cancer, fill = time_group)) + geom_boxplot() +
  labs(x = "XGB Probability", y = "Cancer", fill = "Cases Time\nIn Study\nGroup") +
  scale_y_discrete(labels = get_new_labels(unique(case_all_time$cancer), cancer_decoder))
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_prob_time_in_study.png"), the_plot, "png", height = 4, width = 6)


the_plot <- ggplot(case_all_time, aes(prob, cancer, fill = icd_time_group)) + geom_boxplot() +
  labs(x = "XGB Probability", y = "Cancer", fill = "Cases Time\nFrom ICD\nGroup") +
  scale_y_discrete(labels = get_new_labels(unique(case_all_time$cancer), cancer_decoder))
plot(the_plot)
ggsave(paste0(tort, "_results_", ntree, "/all_prob_time_from_ICD.png"), the_plot, "png", height = 4, width = 6)



wilcox.test(case_all_time$prob[case_all_time$cancer == "bladder" & case_all_time$time_group == 1],
       case_all_time$prob[case_all_time$cancer == "bladder" & case_all_time$time_group == 2])



##################
pval_holder <- matrix(0, nrow = length(cancer_to_comp), ncol = 2)
beta_holder <- matrix(0, nrow = length(cancer_to_comp), ncol = 2)

for(i in 1:length(cancer_to_comp)){
  sub_case <- case_all_time[case_all_time$cancer == cancer_to_comp[i],]
  mod <- lm(prob ~ time , data = sub_case)
  pval_holder[i,1] <- summary(mod)$coef[2,4]
  beta_holder[i,1] <- summary(mod)$coef[2,1]

  mod <- lm(prob ~ time_from_icd , data = sub_case)
  pval_holder[i,2] <- summary(mod)$coef[2,4]
  beta_holder[i,2] <- summary(mod)$coef[2,1]
}
