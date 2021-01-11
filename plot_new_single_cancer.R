library(ggplot2)
library(cowplot)
library(openxlsx)
library(stringr)
library(reshape2)
library(pROC)
library(epitools)
library(glmnet)
library(data.table)
theme_set(theme_cowplot())

simple_auc <- function(TPR, FPR){
 # inputs already sorted, best scores first 
 dFPR <- c(diff(FPR), 0)
 dTPR <- c(diff(TPR), 0)
 sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}



#################################################################


tort <- "ass" #either final or all
ml_metric <- "Expected.Gain"
ntree <- 200

##################################################################

single_roc_plot <- function(df, the_label, plot_label){
  roc_obj <- roc(the_label[,1] ~ df[,1])
  if(roc_obj$auc < 0.5){
   roc_obj <- roc(abs(the_label[,1]-1) ~ df[,1]) 
  }
  roc_df <- data.frame("fpr" = rev(1 - roc_obj$specificities),
                       "tpr" = rev(roc_obj$sensitivities))

  the_plot <- ggplot(roc_df, aes(fpr, tpr)) + geom_line(size=2) +
   geom_abline(intercept = 0, slope = 1) +
   labs(x = "False Positive Rate", y = "True Positive Rate",
        caption = paste0(plot_label, " AUC: ", signif(roc_obj$auc, 3)))
  

  
  cut_1 <- roc_df[which.min(abs(1-rowSums(roc_df))),]
  thresh_1 <- roc_obj$thresholds[which.min(abs(1-rowSums(roc_df)))]
  
  temp <- roc_df[,2] + 1 - roc_df[,1]
  cut_2 <- roc_df[which.max(temp[temp != 1]),]
  thresh_2 <- roc_obj$thresholds[which.max(temp[temp != 1])]
  
  my_binary <- rep(1, nrow(df))
  my_binary[df[,1] > thresh_1] <- 0
  acc1 <- sum(my_binary == the_label[,1])/nrow(df)

  my_binary <- rep(1, nrow(df))
  my_binary[df[,1] > thresh_2] <- 0
  acc2 <- sum(my_binary == the_label[,1])/nrow(df)
    
  if(acc1 > acc2){
    cut <- cut_1
    thresh <- thresh_1
    acc <- acc1
  } else {
    cut <- cut_2
    thresh <- thresh_2
    acc <- acc2
  }
  
  other_stats <- c(as.numeric(cut), acc, thresh)
  names(other_stats) <- c("fpr", "tpr", "acc", "thresh")
  
  return(list(the_plot, roc_df, roc_obj, other_stats))
}


multi_roc_plot <- function(df_list, label_list, plot_label_list){
  roc_df_list <- list()
  auc_phrase_list <- list()
  roc_obj_list <- list()
  for(i in 1:length(df_list)){
    roc_obj <- roc(label_list[[i]][,1] ~ df_list[[i]][,1])
    if(roc_obj$auc < 0.5){
     roc_obj <- roc(abs(label_list[[i]][,1]-1) ~ df_list[[i]][,1])
    }
    roc_obj_list[[i]] <- roc_obj
    roc_df_list[[i]] <- data.frame("fpr" = rev(1 - roc_obj$specificities),
                        "tpr" = rev(roc_obj$sensitivities),
                        "label" = plot_label_list[[i]])
    auc_phrase_list[[i]] <- paste0(plot_label_list[[i]], ": ", signif(roc_obj$auc, 3))
  }
  
  df <- do.call("rbind", roc_df_list)
  auc_phrase <- paste(auc_phrase_list, collapse = ", ")
  
  the_plot <- ggplot(df, aes(fpr, tpr, color = label)) + geom_line(size=2) +
   geom_abline(intercept = 0, slope = 1) +
   labs(x = "False Positive Rate", y = "True Positive Rate", color = "Model",
        caption = paste0("AUC: ", auc_phrase))
  plot(the_plot)
  
  names(roc_obj_list) <- plot_label_list
  return(roc_obj_list)
}

multi_auc_plot <- function(pred_list, label_list, plot_name_list){
 auc_plot_df <- data.frame(matrix(0, nrow = length(pred_list), ncol = 3))
 for(i in 1:length(pred_list)){
  roc_obj <- roc(label_list[[i]][,1] ~ pred_list[[i]][,1])
  auc_plot_df[i,] <- as.numeric(ci.auc(roc_obj))
 }
 auc_plot_df$names <- unlist(plot_name_list)
 colnames(auc_plot_df) <- c("lo", "auc", "hi", "names")
 auc_plot_df$names <- factor(auc_plot_df$names, levels = auc_plot_df$names[order(auc_plot_df$auc)])
 
 the_plot <- ggplot(auc_plot_df, aes(auc, names)) + geom_point() +
  labs(x = "AUC", y = "Prediction") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2)
 return(list(the_plot, auc_plot_df))
}

plot_single_prob <- function(prob, label, plot_label){
 prob_df <- data.frame(prob = prob[,1], label = as.factor(label[,1]))
 
 the_plot <- ggplot(prob_df, aes(prob, fill = label, color = label)) + geom_density(alpha = 0.5) +
  labs(y = "Density", x = paste("Probability by", plot_label, "model"), fill = "Phenotype") +
  guides(color = F) +
  scale_fill_discrete(labels = c("non-cancer", "cancer"))
 plot(the_plot)
}


plot_single_or <- function(pred, label, model_type){
  all_ors <- data.frame(matrix(0, nrow = 6, ncol = 3))
  colnames(all_ors) <- c("or", "lo", "hi")
  i <- 1
  for(cut_off in c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995)){
    hi_cut_off <- quantile(pred[,1], cut_off)
    lo_cut_off <- quantile(pred[,1], 0.5)
    the_table <- matrix(c(sum(pred[,1] > hi_cut_off & label[,1] == 1),
                          sum(pred[,1] > hi_cut_off & label[,1] == 0),
                          sum(pred[,1] < lo_cut_off & label[,1] == 1),
                          sum(pred[,1] < lo_cut_off & label[,1] == 0)),
                        nrow = 2)
    all_ors[i,] <- as.numeric(oddsratio.wald(the_table)$measure[2,])
    i <- i + 1
  }
  all_ors$cut_off <- as.factor(c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995))
  
  the_plot <- ggplot(all_ors, aes(cut_off, or)) + geom_point() + 
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
    labs(x = "Quantile Cut Off", y = paste("Odds Ratio of", model_type, "Model"))

  return(list(the_plot, all_ors))
}


plot_multiple_or <- function(pred_list, label_list, plot_name_list, cut_val = 0.9){
  all_ors <- data.frame(matrix(0, nrow = length(pred_list), ncol = 3))
  colnames(all_ors) <- c("or", "lo", "hi")
  for(i in 1:length(pred_list)){
    hi_cut_off <- quantile(pred_list[[i]][,1], cut_val)
    lo_cut_off <- quantile(pred_list[[i]][,1], 0.5)
    the_table <- matrix(c(sum(pred_list[[i]][,1] > hi_cut_off & label_list[[i]][,1] == 1),
                          sum(pred_list[[i]][,1] > hi_cut_off & label_list[[i]][,1] == 0),
                          sum(pred_list[[i]][,1] < lo_cut_off & label_list[[i]][,1] == 1),
                          sum(pred_list[[i]][,1] < lo_cut_off & label_list[[i]][,1] == 0)),
                        nrow = 2)
    all_ors[i,] <- as.numeric(oddsratio.wald(the_table)$measure[2,])
  }
  
  all_ors$model <- unlist(plot_name_list)
  all_ors <- all_ors[!(is.na(all_ors$or)),]
  print(all_ors)
  for(i in 1:nrow(all_ors)){
    if(all_ors[i,1] < 1){
      print(all_ors[i,])
      all_ors[i,1:3] <- exp(abs(log(all_ors[i,1:3])))
    }
  }
  all_ors$model <- factor(all_ors$model, levels = all_ors$model[order(all_ors$or)])
  the_plot <- ggplot(all_ors, aes(or, model)) + geom_point() + 
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    labs(x = "Odds Ratio", y = "Model")
  return(list(the_plot, all_ors))
}

get_new_labels <- function(old_labels){
  old_labels <- as.character(old_labels)
  new_labels <- rep("", length(old_labels))
  if(!all(old_labels %in% feat_decoder[,1] | old_labels == "")){
    print(old_labels[!(old_labels %in% feat_decoder[,1])])
    print("missing feat decoder")
    exit()
  }
  for(i in 1:length(old_labels)){
    if(old_labels[i] != ""){
      new_labels[i] <- feat_decoder[feat_decoder[,1] == old_labels[i], 2]
    }
  }
  return(new_labels)
}

feat_decoder <- read.table("label_decoder.csv", stringsAsFactors = F, sep = ",")
feat_decoder[,1] <- trimws(feat_decoder[,1])
feat_decoder[,2] <- trimws(feat_decoder[,2])

################################################



temp_func <- function(name_cancer){

  print("set up data")
  ########################################################################
  ######################### SET UP DATA ##################################
 ############################################################################
  
  #Doing my own feature positive/negative ranking

  valid_data <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                    name_cancer, "_valid.txt"), stringsAsFactors = F, header = F) 
  train_data <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                     name_cancer, "_train.txt"), stringsAsFactors = F, header = F)
  test_data <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                     name_cancer, "_test.txt"), stringsAsFactors = F, header = F)
  
  valid_cols <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                            name_cancer, "_valid_cols.txt"), stringsAsFactors = F, header = F, sep = "\t", quote = "$") 
  train_cols <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                            name_cancer, "_train_cols.txt"), stringsAsFactors = F, header = F, sep = "\t", quote = "$")
  test_cols <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                            name_cancer, "_test_cols.txt"), stringsAsFactors = F, header = F, sep = "\t", quote = "$")

  valid_lab <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                  name_cancer, "_valid_lab.txt"), stringsAsFactors = F, header = F) 
  train_lab <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                  name_cancer, "_train_lab.txt"), stringsAsFactors = F, header = F)
  test_lab <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                 name_cancer, "_test_lab.txt"), stringsAsFactors = F, header = F)

  valid_cols$V2 <- gsub("+", "", gsub("'", "", gsub("-", "", gsub(",", "", gsub("/", "", gsub("\\(", "",
                               gsub("\\)", "", gsub(" ", "_", valid_cols$V1))))))), fixed = T)
  train_cols$V2 <- gsub("+", "", gsub("'", "", gsub("-", "", gsub(",", "", gsub("/", "", gsub("\\(", "", gsub("\\)", "",
                               gsub(" ", "_", train_cols$V1))))))), fixed = T)
  test_cols$V2 <- gsub("+", "", gsub("'", "", gsub("-", "", gsub(",", "", gsub("/", "", gsub("\\(", "",
                               gsub("\\)", "", gsub(" ", "_", test_cols$V1))))))), fixed = T)
  
  colnames(valid_data) <- valid_cols$V2
  colnames(train_data) <- train_cols$V2
  colnames(test_data) <- test_cols$V2
  
  valid_data$cancer.phenotype <- valid_lab[,1]
  train_data$cancer.phenotype <- train_lab[,1]
  test_data$cancer.phenotype <- test_lab[,1]
  
  
  feat_imp <- read.xlsx(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                               name_cancer, "_feature_interactions.xlsx"), sheet = 1)
  feat_imp$Interaction <- gsub("+", "", gsub("'", "", gsub("-", "", gsub(",", "", gsub("/", "", gsub("\\(", "", gsub("\\)", "",
                                          gsub(" ", "_", feat_imp$Interaction))))))), fixed = T)


  inter_imp <- read.xlsx(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                name_cancer, "_feature_interactions.xlsx"), sheet = 2)
  inter_imp$feat1 <- str_split(inter_imp[,1], fixed("|"), simplify = T)[,1]
  inter_imp$feat2 <- str_split(inter_imp[,1], fixed("|"), simplify = T)[,2]
  inter_imp$feat1 <- gsub("+", "", gsub("'", "", gsub("-", "", gsub(",", "", gsub("/", "", gsub("\\(", "", gsub("\\)", "",
                                   gsub(" ", "_", inter_imp$feat1))))))), fixed = T)
  inter_imp$feat2 <- gsub("+", "", gsub("'", "", gsub("-", "", gsub(",", "", gsub("/", "", gsub("\\(", "", gsub("\\)", "",
                                   gsub(" ", "_", inter_imp$feat2))))))), fixed = T)

  # if(ml_metric != "Gain"){
  #   colnames(feat_imp)[colnames(feat_imp) == "Gain"] <- "old_gain"
  #   colnames(feat_imp)[colnames(feat_imp) == ml_metric] <- "Gain"
  #   feat_imp <- feat_imp[order(feat_imp$Gain, decreasing = T), ]
  # 
  #   colnames(inter_imp)[colnames(inter_imp) == "Gain"] <- "old_gain"
  #   colnames(inter_imp)[colnames(inter_imp) == ml_metric] <- "Gain"
  #   inter_imp <- inter_imp[order(inter_imp$Gain, decreasing = T), ]
  # }
  # 
  # 
  # 
  # add_to <- feat_imp$Interaction[1:15][!(feat_imp$Interaction[1:15] %in% trimws(feat_decoder[,1]))]
  # if(length(add_to) > 0){
  #   print("fix labels")
  #   add_to <- cbind(add_to, rep("", length(add_to)))
  #   colnames(add_to) <- colnames(feat_decoder)
  #   feat_decoder <- rbind(feat_decoder, add_to)
  #   write.table(feat_decoder, "check_names.csv", row.names = F, col.names = F, quote=T, sep = ",")
  #   exit()
  # }
  # 
  # 
  # feat_direc <- rep(0, nrow(feat_imp))
  # for(i in 1:length(feat_direc)){
  #  mod <- glm(as.formula(paste("cancer.phenotype ~ ", feat_imp$Interaction[i])),
  #      family = "binomial", data = train_data)
  #  feat_direc[i] <- sign(as.numeric(mod$coefficients[2]))
  # }
  # 
  # feat_imp$direc <- feat_direc
  # feat_imp$direc[feat_imp$direc == 0] <- 1
  # saveRDS(feat_imp, paste0(tort, "_results_", ntree, "/", name_cancer, "/process_feat_imp.RDS"))
  # 
  # 
  # 
  #   
  #   inter_direc <- rep(0, nrow(inter_imp))
  #   inter_uni_p <- rep(NA, nrow(inter_imp))
  #   inter_uni_coef <- matrix(NA, nrow = nrow(inter_imp), ncol = 4)
  #   for(i in 1:length(inter_direc)){
  #     if(inter_imp$feat1[i] == inter_imp$feat2[i]){
  #       mod <- glm(as.formula(paste("cancer.phenotype ~ ",
  #               inter_imp$feat1[i], "+poly(", inter_imp$feat1[i], ",2)")),
  #               family = "binomial", data = train_data)
  #       to_line <- 3
  #     } else {
  #       mod <- glm(as.formula(paste("cancer.phenotype ~ ",
  #                                   inter_imp$feat1[i], "*", inter_imp$feat2[i])),
  #                  family = "binomial", data = train_data)
  #       to_line <- 4
  #     }
  #    
  #     if(nrow(summary(mod)$coef) == to_line){
  #       print("good")
  #       inter_direc[i] <- mod$coefficients[to_line]
  #       inter_uni_p[i] <- summary(mod)$coef[to_line,4]
  #       inter_uni_coef[i,] <- summary(mod)$coef[to_line,]
  #     }
  #   }
  #   
  #   inter_imp$direc <- sign(inter_direc)
  #   inter_imp$direc[inter_imp$direc == 0] <- 1 #when NAs appear in the model
  #   inter_imp$p_val <- inter_uni_p
  #   saveRDS(inter_imp, paste0(tort, "_results_", ntree, "/", name_cancer, "/process_inter_imp.RDS"))
  #   saveRDS(cbind(inter_imp, inter_uni_coef), paste0(tort, "_results_", ntree, "/",
  #                                                    name_cancer, "/process_inter_coef.RDS"))
  #   
  #   
  #   print("get feature/info plots")
  #   ############################################################################
  #   ################## GET FEATURE INFO/PLOTS ##################################
  #   ############################################################################
  #   
  #   #Now make the feature plots
  #   
  #   
  #   feat_imp$pro_gain <- feat_imp$Gain/sum(feat_imp$Gain)
  #   extent_feat_imp <- min(10, nrow(feat_imp))
  #   sub_feat_imp <- feat_imp[1:extent_feat_imp,]
  #   sub_feat_imp$Interaction <- factor(sub_feat_imp$Interaction,
  #                           levels = sub_feat_imp$Interaction[order(sub_feat_imp$pro_gain)])
  #   
  #   
  #   the_plot <- ggplot(sub_feat_imp, aes(pro_gain, Interaction, fill = as.factor(direc))) + geom_col() +
  #    labs(x = "Fraction of Gain", y = "Feature", fill = "Direction") +
  #    scale_y_discrete(labels = get_new_labels(levels(sub_feat_imp$Interaction)))
  #   plot(the_plot)
  #   ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/feat_imp.png"),
  #          the_plot, "png", height=6, width=7)
  #   
  #   
  #     inter_imp$pro_gain <- inter_imp$Gain/sum(inter_imp$Gain)
  #     sub_inter_imp <- inter_imp[1:min(10, nrow(inter_imp)),]
  #     sub_inter_imp$Interaction <- paste0(get_new_labels(sub_inter_imp$feat1), "|", get_new_labels(sub_inter_imp$feat2))
  #     sub_inter_imp$Interaction <- factor(sub_inter_imp$Interaction,
  #                           levels = sub_inter_imp$Interaction[order(sub_inter_imp$pro_gain)])
  #   
  #     the_plot <- ggplot(sub_inter_imp, aes(pro_gain, Interaction, fill = as.factor(direc))) + geom_col() +
  #      labs(x = "Fraction of Gain", y = "Feature", fill = "Direction") +
  #       theme(axis.text.y = element_text(size=10))
  #     plot(the_plot)
  #     ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/inter_imp.png"),
  #            the_plot, "png", height=6, width=10)
  #   
  #   
  #   print("linear plots")
  #   ##########################################################
  #   ##########   Now make my own LINEAR FEATURE plots  ##################
  #   ##########################################################
  #   
  #   mod <- glm(cancer.phenotype ~ ., data = train_data[,colnames(train_data) %in%
  #                 c("cancer.phenotype", feat_imp$Interaction)], family = "binomial")
  #   
  #   linear_coef <- as.data.frame(summary(mod)$coef)[-1,]
  #   linear_coef$feat_imp_name <- rownames(linear_coef)
  #   linear_coef <- linear_coef[linear_coef$feat_imp_name %in% feat_imp$Interaction,]
  #   linear_coef <- linear_coef[order(linear_coef$feat_imp_name)[rank(feat_imp$Interaction)],]
  #   #linear_coef$feat_imp_name[nchar(linear_coef$feat_imp_name) > 40] <-
  #   # paste0(substr(linear_coef$feat_imp_name[nchar(linear_coef$feat_imp_name) > 40], 1, 40), "...")
  #   
  #   
  #   #MAKE PLOT HERE WITH LINEAR COEF
  #   colnames(linear_coef)[2] <- "stderr"
  #   plot_df <- linear_coef[1:nrow(sub_feat_imp),]
  #   plot_df$feat_imp_name <- factor(plot_df$feat_imp_name, levels = rev(plot_df$feat_imp_name))
  #   the_plot <- ggplot(plot_df, aes(exp(Estimate), feat_imp_name)) +
  #    geom_point() + geom_vline(aes(xintercept = 1)) +
  #    geom_errorbarh(aes(xmin = exp(Estimate - stderr), xmax = exp(Estimate + stderr)), height = 0.2) +
  #    labs(x = "Odds Ratio", y = "Feature") +
  #     scale_y_discrete(labels = get_new_labels(levels(plot_df$feat_imp_name)))
  #   
  #   plot(the_plot)
  #   ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/all_linear_feat_terms.png"),
  #                    the_plot, "png", height=6, width=8)
  #   saveRDS(linear_coef, paste0(tort, "_results_", ntree, "/", name_cancer, "/linear_coef.RDS"))
  #   
  #   
  #   #try to main quasi-gain from effect
  #   linear_coef$Gain <- abs(linear_coef$Estimate)/sum(abs(linear_coef$Estimate))
  #   plot_df <- data.frame(xgb_gain = sub_feat_imp$pro_gain,
  #                         linear_gain = linear_coef$Gain[1:nrow(sub_feat_imp)],
  #                        feat = as.character(sub_feat_imp$Interaction))
  #   plot_df <- reshape2::melt(plot_df, "feat")
  #   plot_df$feat <- factor(plot_df$feat,
  #                              levels = sub_feat_imp$Interaction[order(sub_feat_imp$pro_gain)])
  #   the_plot <- ggplot(plot_df, aes(value, feat, fill = variable)) +
  #     geom_bar(position="dodge",stat="identity") +
  #    labs(x = "Quasi-Gain", y = "Feature", fill = "Model\nType") +
  #     scale_y_discrete(labels = get_new_labels(levels(plot_df$feat))) +
  #     scale_fill_discrete(labels = c("XGB Gain", "Linear Gain"))
  #   plot(the_plot)
  #   ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/linear_and_xgb_gain.png"),
  #          the_plot, "png", height=6, width=8)
  #   
  #   print("interaction plots")
  #   ##########################################################
  #   ##########   Now make my own INTERACTION FEATURE plots  ##################
  #   ##########################################################
  # 
  #   #prepare the formula - 100 linear terms, 50 interaction terms
  #   #the_formula <- paste0("cancer.phenotype ~ ", paste(colnames(train_data)[
  #   # -which(colnames(train_data) == "cancer.phenotype")], collapse = "+"))
    the_formula <- paste0("cancer.phenotype ~ ", paste(feat_imp$Interaction, collapse = "+"))
    inter_terms <- rep("", min(c(50, nrow(inter_imp))))
    other_possible_terms <- rep("", min(c(50, nrow(inter_imp))))
    for(i in 1:min(c(50, nrow(inter_imp)))){
     if(inter_imp$feat1[i] == inter_imp$feat2[i]){
       inter_terms[i] <- paste0("I(", inter_imp$feat1[i], "^2)")
     } else {
       inter_terms[i] <- paste0(inter_imp$feat1[i], ":", inter_imp$feat2[i])
       other_possible_terms[i] <- paste0(inter_imp$feat2[i], ":", inter_imp$feat1[i])
     }
    }

  #   # all_between_inter_p <- rep(0, length(inter_terms))
  #   # for(i in 1:length(inter_terms)){
  #   #   the_formula <- paste0(the_formula, "+", inter_terms[i])
  #   #   mod <- glm(as.formula(the_formula), data = train_data, family = "binomial")
  #   #   all_between_inter_p[i] <- summary(mod)$coef[nrow(summary(mod)$coef),4]
  #   # }
  # 
    the_formula <- paste0(the_formula, "+", paste(inter_terms, collapse = "+"))
    mod <- glm(as.formula(the_formula), data = train_data, family = "binomial")

    linear_inter_pred_valid <- predict(mod, valid_data)
    complex_valid_roc_obj <- roc(valid_data$cancer.phenotype ~ predict(mod, valid_data))
    complex_valid_roc <- data.frame("X2" = rev(1 - complex_valid_roc_obj$specificities),
                              "X1" = rev(complex_valid_roc_obj$sensitivities))

    linear_inter_pred_test <- predict(mod, test_data)
    complex_test_roc_obj <- roc(valid_data$cancer.phenotype ~ predict(mod, valid_data))
    complex_test_roc <- data.frame("X2" = rev(1 - complex_test_roc_obj$specificities),
                                    "X1" = rev(complex_test_roc_obj$sensitivities))
  # 
  #   #ODDS RATIOS - of coefficients
  #   inter_coef <- as.data.frame(summary(mod)$coef)
  #   inter_coef <- inter_coef[grepl(":", rownames(inter_coef)) | grepl("^2", rownames(inter_coef), fixed = T),]
  #   inter_coef$inter_imp_name <- rownames(inter_coef)
  # 
  #   new_inter_terms <- inter_terms[(inter_terms %in% inter_coef$inter_imp_name) |
  #                                (other_possible_terms %in% inter_coef$inter_imp_name)]
  #   new_other_possible_terms <- other_possible_terms[(inter_terms %in% inter_coef$inter_imp_name) |
  #                                (other_possible_terms %in% inter_coef$inter_imp_name)]
  #   new_inter_terms[!(new_inter_terms %in% inter_coef$inter_imp_name)] <-
  #    new_other_possible_terms[!(new_inter_terms %in% inter_coef$inter_imp_name)]
  #   inter_coef <- inter_coef[order(inter_coef$inter_imp_name)[rank(new_inter_terms)],]
  # 
  #   colnames(inter_coef)[2] <- "stderr"
  #   plot_df <- inter_coef[1:nrow(sub_inter_imp),]
  #   plot_df$inter_imp_name <- factor(plot_df$inter_imp_name, levels = rev(plot_df$inter_imp_name))
  # 
  #   temp_feat1 <- str_split(levels(plot_df$inter_imp_name), ":", simplify = T)[,1]
  #   temp_feat2 <- str_split(levels(plot_df$inter_imp_name), ":", simplify = T)[,2]
  #   square_bool <- grepl("^2", temp_feat1, fixed = T)
  #   if(sum(square_bool) > 0){
  #     temp_feat1[grepl("^2", temp_feat1, fixed = T)] <- unlist(strsplit(unlist(lapply(
  #       strsplit(grep("^2", temp_feat1, value = T, fixed = T),
  #                "(", fixed = T), function(x) x[2])), "^2)", fixed = T))
  #   }
  #   temp_feat1 <- get_new_labels(temp_feat1)
  #   temp_feat2 <- get_new_labels(temp_feat2)
  #   new_feat <- paste0(temp_feat1, ":", temp_feat2)
  #   if(sum(square_bool) > 0){
  #     new_feat[square_bool] <- paste0(substr(new_feat, 1, nchar(new_feat)-1)[square_bool], "^2")
  #   }
  # 
  #   the_plot <- ggplot(plot_df, aes(exp(Estimate), inter_imp_name)) +
  #    geom_point() + geom_vline(aes(xintercept = 1)) +
  #    geom_errorbarh(aes(xmin = exp(Estimate - stderr), xmax = exp(Estimate + stderr)), height = 0.2) +
  #    labs(x = "Odds Ratio", y = "Interaction") +
  #     scale_y_discrete(labels = new_feat) +
  #     theme(axis.text.y = element_text(size=10))
  # 
  #   plot(the_plot)
  #   ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/all_linear_inter_terms.png"),
  #          the_plot, "png", height=6, width=9)
  # 
  #   #P-VALUES - of coefficients
  #   colnames(plot_df)[4] <- "P"
  #   the_plot <- ggplot(plot_df, aes(-log10(P), inter_imp_name)) +
  #    geom_point() + geom_vline(aes(xintercept = -log10(0.05))) +
  #    labs(x = "-log10(P-Value)", y = "Interaction") +
  #     scale_y_discrete(labels = new_feat) +
  #     theme(axis.text.y = element_text(size=10))
  # 
  #   plot(the_plot)
  #   ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/all_linear_inter_pval.png"),
  #          the_plot, "png", height=6, width=9)
  #   saveRDS(inter_coef, paste0(tort, "_results_", ntree, "/", name_cancer, "/interaction_things.RDS"))
  # 
  #   #different way of calculating p-values - from the more univariate models
  #   new_plot_df <- inter_imp[1:10,]
  #   ggplot(new_plot_df, aes(-log10(p_val), Interaction)) + geom_point() +
  #     scale_y_discrete(labels = new_feat) +
  #     labs(x = "-log10(P-Value)", y = "Interaction")
  # 
  # 
  #   # #############################################################################
  #   # #CHECK FOR INTERACTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #   # random_pval <- rep(NA, length(feat_imp$Interaction))
  #   # random_names_1 <- sample(colnames(train_data)[-which(colnames(train_data) == "cancer.phenotype")],
  #   #                          length(feat_imp$Interaction))
  #   # random_names_2 <- sample(colnames(train_data)[-which(colnames(train_data) == "cancer.phenotype")],
  #   #                          length(feat_imp$Interaction))
  #   # for(k in 1:length(feat_imp$Interaction)){
  #   #   mod <- glm(as.formula(paste0("cancer.phenotype ~ ", random_names_1[k],"*",random_names_2[k])),
  #   #              data = train_data, family = "binomial")
  #   #   if(nrow(summary(mod)$coef) == 4){
  #   #     random_pval[k] <- summary(mod)$coef[4,4]
  #   #   }
  #   # }
  #   # 
  #   # top_feat_pval <- rep(NA, length(feat_imp$Interaction))
  #   # random_names_1 <- sample(feat_imp$Interaction)
  #   # random_names_2 <- sample(feat_imp$Interaction)
  #   # for(k in 1:length(feat_imp$Interaction)){
  #   #   mod <- glm(as.formula(paste0("cancer.phenotype ~ ", random_names_1[k],"*",random_names_2[k])),
  #   #              data = train_data, family = "binomial")
  #   #   if(nrow(summary(mod)$coef) == 4){
  #   #     top_feat_pval[k] <- summary(mod)$coef[4,4]
  #   #   }
  #   # }
  #   # 
  #   # very_top_feat_pval <- rep(NA, length(feat_imp$Interaction))
  #   # random_names_1 <- sample(feat_imp$Interaction[order(feat_imp$pro_gain, decreasing = T)][1:10], 1600, replace = T)
  #   # random_names_2 <- sample(feat_imp$Interaction[order(feat_imp$pro_gain, decreasing = T)][1:10], 1600, replace = T)
  #   # both_random <- paste0(random_names_1, random_names_2)
  #   # random_names_1 <- random_names_1[!duplicated(both_random)]
  #   # random_names_2 <- random_names_2[!duplicated(both_random)]
  #   # for(k in 1:length(feat_imp$Interaction)){
  #   #   mod <- glm(as.formula(paste0("cancer.phenotype ~ ", random_names_1[k],"*",random_names_2[k])),
  #   #              data = train_data, family = "binomial")
  #   #   if(nrow(summary(mod)$coef) == 4){
  #   #     very_top_feat_pval[k] <- summary(mod)$coef[4,4]
  #   #   }
  #   # }
  #   # 
  #   # pval_plot_df <- data.frame(random_p = random_pval, top_p = top_feat_pval,
  #   #                            best_p = very_top_feat_pval, my_p = inter_imp$p_val[1:length(random_pval)])
  #   # pval_plot_df <- reshape2::melt(pval_plot_df)
  #   # the_plot <- ggplot(pval_plot_df, aes(variable, -log10(value))) + geom_boxplot(notch = T) +
  #   #   labs(x = "Features", y = "-log10(P-Value)") +
  #   #   scale_x_discrete(labels = c("Random", "Top 100", "Top 10", "XGB Inter."))
  #   # plot(the_plot)
  #   # ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/all_pval_comp.png"),
  #   #        the_plot, "png", height=5, width=6)
  #   # saveRDS(pval_plot_df, paste0(tort, "_results_", ntree, "/", name_cancer, "/all_pval.RDS"))


    print("read in data")
    ##########################################################
    ##########   READ IN DATA  ##################
    ##########################################################
    
  valid_lab <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                   name_cancer, "_valid_lab.txt"), stringsAsFactors = F, header = F) 
  train_lab <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                   name_cancer, "_train_lab.txt"), stringsAsFactors = F, header = F)
  test_lab <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                  name_cancer, "_test_lab.txt"), stringsAsFactors = F, header = F)
  
  pgs_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                              name_cancer, "_pgs_test_prob.txt"))
  pgs_valid_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                              name_cancer, "_pgs_valid_prob.txt"))
    
    
  linear_nopgs_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                              name_cancer, ".linear_nopgs_test_prob.txt"))
  linear_nopgs_valid_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                              name_cancer, ".linear_nopgs_valid_prob.txt"))
  linear_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                              name_cancer, ".linear_test_prob.txt"))
  linear_valid_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                               name_cancer, ".linear_valid_prob.txt"))
  
  rf_nopgs_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                              name_cancer, ".rf_nopgs_test_prob.txt"))
  rf_nopgs_valid_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                               name_cancer, ".rf_nopgs_valid_prob.txt"))
  rf_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                        name_cancer, ".rf_test_prob.txt"))
  rf_valid_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                         name_cancer, ".rf_valid_prob.txt"))
  
  svm_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                     name_cancer, ".svm_test_prob.txt"))
  knn_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                     name_cancer, ".knn_test_prob.txt"))
  
  xgb_nopgs_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                              name_cancer, ".xgb_nopgs_test_prob.txt"))
  xgb_nopgs_valid_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                               name_cancer, ".xgb_nopgs_valid_prob.txt"))
  xgb_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                        name_cancer, ".xgb_test_prob.txt"))
  xgb_valid_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                         name_cancer, ".xgb_valid_prob.txt"))
  
  xgb_ten_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                     name_cancer, ".xgb_test_prob.ten.txt"))
  rf_ten_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                         name_cancer, ".rf_test_prob.ten.txt"))
  svm_ten_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                         name_cancer, ".svm_test_prob.ten.txt"))
  knn_ten_test_prob <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
                                         name_cancer, ".knn_test_prob.ten.txt"))
  
  
  ##########################################
  # Making all of the models #
  ##########################################
  #make simple model
  mod <- glm(paste0("cancer.phenotype ~",
             paste(feat_imp$Interaction[1:5], collapse = "+")),
      data = train_data, family = "binomial")
  five_feat_test_pred <- data.frame(predict(mod, test_data))
  five_feat_valid_pred <- data.frame(predict(mod, valid_data))
  
  mod <- glm(paste0("cancer.phenotype ~",
                    paste(feat_imp$Interaction[1:3], collapse = "+")),
             data = train_data, family = "binomial")
  three_feat_test_pred <- data.frame(predict(mod, test_data))
  three_feat_valid_pred <- data.frame(predict(mod, valid_data))
  
  
  #make ten normal feature model
  mod <- glm(paste0("cancer.phenotype ~",
                    paste(feat_imp$Interaction[1:10], collapse = "+")),
             data = train_data, family = "binomial")
  realten_feat_test_pred <- data.frame(predict(mod, test_data))
  realten_feat_valid_pred <- data.frame(predict(mod, valid_data))
  saveRDS(summary(mod)$coefficients, paste0(tort, "_results_", ntree, "/", name_cancer, "/ten_feat_mod_coefs.RDS"))
  
  
  #make simple interaction model
  mod <- glm(paste0("cancer.phenotype ~",
                    paste(feat_imp$Interaction[1:5], collapse = "+"),
                    "+", paste(inter_terms[1:5], collapse = "+")),
             data = train_data, family = "binomial")
  ten_feat_test_pred <- data.frame(predict(mod, test_data))
  ten_feat_valid_pred <- data.frame(predict(mod, valid_data))
  
  #glmnet time
  zero_glmnet <- function(x){
    net_mod <- glmnet(y = as.matrix(train_data$cancer.phenotype),
                    x = as.matrix(train_data[,-which(colnames(train_data) == "cancer.phenotype")]),
                    lambda = x)
    y <- sum(as.numeric(net_mod$beta) != 0) - 10
    return(y)
  }
  
  the_root <- uniroot(zero_glmnet, lower = 5e-8, upper = 1)
  net_mod <- glmnet(y = as.matrix(train_data$cancer.phenotype),
                    x = as.matrix(train_data[,-which(colnames(train_data) == "cancer.phenotype")]),
                    lambda = the_root$root)
  beta_mat <- as.matrix(net_mod$beta)
  mod <- glm(paste0("cancer.phenotype ~",
                    paste(rownames(beta_mat)[beta_mat[,1] != 0], collapse = "+")),
             data = train_data, family = "binomial")
  enet_feat_test_pred <- data.frame(predict(mod, test_data))
  enet_feat_valid_pred <- data.frame(predict(mod, valid_data))
  saveRDS(summary(mod)$coefficients, paste0(tort, "_results_", ntree, "/", name_cancer, "/ten_feat_enet_coefs.RDS"))
  

  print("roc plot")
  ##########################################################
  ##########   MAKE ROC PLOT  ##################
  ##########################################################
  
  #comparisons I want to make: pgs, no pgs; linear, xgb, rf; valid, test; interaction, none

  single_roc <- single_roc_plot(xgb_test_prob, test_lab, "my linear")
  single_roc <- single_roc_plot(realten_feat_test_pred, test_lab, "my linear")
  saveRDS(single_roc[[4]], paste0(tort, "_results_", ntree, "/", name_cancer, "/ten_feat_stats.RDS"))


  #multi_roc_plot(list(xgb_test_prob, xgb_valid_prob),
  #               list(test_lab, valid_lab),
  #               list("pgs_test", "pgs_valid"))
  roc_objs <- multi_roc_plot(list(linear_test_prob, rf_test_prob, xgb_test_prob, svm_test_prob, knn_test_prob),
                 list(test_lab, test_lab, test_lab, test_lab, test_lab),
                 list("linear", "rf", "xgb", "svm", "knn"))
  saveRDS(roc_objs, paste0(tort, "_results_", ntree, "/", name_cancer, "/roc_objs.1.RDS"))
  
  #multi_roc_plot(list(xgb_nopgs_test_prob, xgb_test_prob),
  #               list(test_lab, test_lab),
  #               list("xgb_no_pgs", "xgb"))
  #multi_roc_plot(list(linear_test_prob, as.matrix(linear_inter_pred_test), xgb_test_prob),
  #               list(test_lab, as.matrix(test_data$cancer.phenotype), test_lab),
  #               list("linear", "linear_w_interaction", "xgb"))
  
  roc_objs <- multi_roc_plot(list(linear_test_prob, five_feat_test_pred, ten_feat_test_pred, three_feat_test_pred),
                 list(test_lab, test_lab, test_lab, test_lab),
                 list("linear", "five", "ten", "three"))
  saveRDS(roc_objs, paste0(tort, "_results_", ntree, "/", name_cancer, "/roc_objs.2.RDS"))
  
  roc_objs <- multi_roc_plot(list(ten_feat_test_pred, rf_ten_test_prob, xgb_ten_test_prob, svm_ten_test_prob, knn_ten_test_prob),
                             list(test_lab, test_lab, test_lab, test_lab, test_lab),
                             list("linear", "rf", "xgb", "svm", "knn"))
  saveRDS(roc_objs, paste0(tort, "_results_", ntree, "/", name_cancer, "/roc_objs.3.RDS"))
  
  
  print("auc plot")
  ##########################################################
  ##########   MAKE AUC PLOTS  ##################
  ##########################################################
  
  
  auc_items <- multi_auc_plot(list(xgb_test_prob, xgb_valid_prob),
                 list(test_lab, valid_lab),
                 list("PGS - test set", "PGS - valid set"))
  the_plot <- auc_items[[1]]
  saveRDS(auc_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_test_valid.RDS"))
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_test_valid.png"),
         the_plot, "png", height=4, width=5)
  
  auc_items <- multi_auc_plot(list(linear_test_prob, rf_test_prob, xgb_test_prob),
                 list(test_lab, test_lab, test_lab),
                 list("Linear", "Random\nForest", "XGB"))
  the_plot <- auc_items[[1]]
  saveRDS(auc_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_all.RDS"))
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_all_models.png"),
         the_plot, "png", height=4, width=5)
  
  auc_items <- multi_auc_plot(list(xgb_nopgs_test_prob, xgb_test_prob),
                 list(test_lab, test_lab),
                 list("XGB: No PGS", "XGB: w/ PGS"))
  the_plot <- auc_items[[1]]
  saveRDS(auc_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_pgs.RDS"))
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_pgs_comp.png"),
         the_plot, "png", height=4, width=5)
  
  auc_items <- multi_auc_plot(list(linear_test_prob, as.matrix(linear_inter_pred_test), xgb_test_prob),
                 list(test_lab, as.matrix(test_data$cancer.phenotype), test_lab),
                 list("Linear: w/out interaction", "Linear: w/ interaction", "XGB"))
  the_plot <- auc_items[[1]]
  saveRDS(auc_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_inter.RDS"))
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_inter_comp.png"),
         the_plot, "png", height=4, width=5)
  
  auc_items <- multi_auc_plot(list(linear_test_prob, five_feat_test_pred, ten_feat_test_pred),
                 list(test_lab, test_lab, test_lab),
                 list("Linear: All Feat", "Linear: Five Feat", "Linear: Ten Feat"))
  the_plot <- auc_items[[1]]
  saveRDS(auc_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_few_justlinear.RDS"))
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_few_justlinear.png"),
         the_plot, "png", height=4, width=5)
  
  
  auc_items <- multi_auc_plot(list(linear_test_prob, five_feat_test_pred, ten_feat_test_pred,
                                   realten_feat_test_pred, enet_feat_test_pred, three_feat_test_pred),
                              list(test_lab, test_lab, test_lab, test_lab, test_lab, test_lab),
                              list("Linear: All Feat", "Linear: Five Feat", "Inter: Five by 2",
                                   "Linear: Ten Feat", "Enet: Five Feat", "Linear: Three Feat"))
  the_plot <- auc_items[[1]]
  saveRDS(auc_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_justlinear.RDS"))
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_justlinear.png"),
         the_plot, "png", height=4, width=5)
  
  auc_items <- multi_auc_plot(list(xgb_test_prob, svm_test_prob, rf_test_prob, knn_test_prob, linear_test_prob),
                              list(test_lab, test_lab, test_lab, test_lab, test_lab),
                              list("XGB", "SVM", "RF", "KNN", "Linear"))
  the_plot <- auc_items[[1]]
  saveRDS(auc_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_all_ml.RDS"))
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_all_ml.png"),
         the_plot, "png", height=4, width=5)
  
  auc_items <- multi_auc_plot(list(xgb_ten_test_prob, rf_ten_test_prob, svm_ten_test_prob, knn_ten_test_prob, ten_feat_test_pred),
                              list(test_lab, test_lab, test_lab, test_lab, test_lab),
                              list("XGB", "RF", "SVM", "KNN", "Linear: Ten Feat"))
  the_plot <- auc_items[[1]]
  saveRDS(auc_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_all_ten_ml.RDS"))
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_all_ten_ml.png"),
         the_plot, "png", height=4, width=5)
  return(1)
  
  print("make probability plot")
  #########################################################################
  ##################33  Make probability plot #############################
  #########################################################################
  

  
  the_plot <- plot_single_prob(linear_test_prob, test_lab, "linear")
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/prob_density_linear.png"),
         the_plot, "png", height=5, width=6)

  the_plot <- plot_single_prob(xgb_test_prob, test_lab, "xgb")
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/prob_density_xgb.png"),
         the_plot, "png", height=5, width=6)

  
  print("get oods ratio")
  ########################################################################
  #####################  Get odds ratios #################################
  #########################################################################
  
  or_items <- plot_single_or(linear_test_prob, test_lab, "Linear")
  plot(or_items[[1]])
  or_items <- plot_single_or(xgb_test_prob, test_lab, "XGB")
  plot(or_items[[1]])
  or_items <- plot_single_or(realten_feat_test_pred, test_lab, "Linear: Ten")
  plot(or_items[[1]])
  saveRDS(or_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/or_linear_ten.RDS"))
  
  or_items <- plot_multiple_or(list(linear_test_prob, as.matrix(linear_inter_pred_test),
                                    xgb_nopgs_test_prob, xgb_valid_prob, xgb_test_prob),
                             list(test_lab, as.matrix(test_data$cancer.phenotype), test_lab,
                                  valid_lab, test_lab),
                             list("Linear: w/out Interaction", "Linear: w/ Interaction",
                                  "XGB: w/out PGS", "XGB: valid", "XGB"))
  plot(or_items[[1]])
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/prob_or_all.png"),
         or_items[[1]], "png", height=5, width=6)
  
  
  or_items <- plot_multiple_or(list(linear_test_prob,
                                    xgb_nopgs_test_prob, xgb_valid_prob, xgb_test_prob),
                               list(test_lab, test_lab,
                                    valid_lab, test_lab),
                               list("Linear", "XGB: w/out PGS", "XGB: valid", "XGB"))
  plot(or_items[[1]])
  saveRDS(or_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/or_all.RDS"))
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/prob_or_all2.png"),
         or_items[[1]], "png", height=5, width=6)
  
  
  or_items <- plot_multiple_or(list(linear_test_prob, realten_feat_test_pred, xgb_test_prob),
                               list(test_lab, test_lab, test_lab),
                               list("Linear", "Linear: Ten", "XGB"))
  plot(or_items[[1]])
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/prob_or_all2.png"),
         or_items[[1]], "png", height=5, width=6)
  
  or_items <- plot_multiple_or(list(xgb_test_prob, svm_test_prob, rf_test_prob, knn_test_prob, ten_feat_test_pred),
                              list(test_lab, test_lab, test_lab, test_lab, test_lab),
                              list("XGB", "SVM", "RF", "KNN", "Linear: Ten Feat"))
  saveRDS(or_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/prob_or_all_ml.RDS"))
  ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/prob_or_all_ml.png"),
         or_items[[1]], "png", height=4, width=5)
  
  
  ########################################################################
  #####################  Get time to #################################
  #########################################################################

  # eid <- read.table(paste0(tort, "_results_", ntree, "/", name_cancer, "/jupyter/",
  #                                    name_cancer, "_valid_eid.txt"))
  # 
  # #get date of birth
  # birth_info <- read.table("back_results_200/icd_dates/birth_info", sep=",", header=T)
  # birth_info <- birth_info[birth_info[,1] %in% eid[,1],]
  # birth_info <- birth_info[order(birth_info[,1])[rank(eid[,1])],]
  # birth_info$birth <- as.Date(paste0("15-", birth_info[,3], "-", birth_info[,2]), "%d-%m-%Y")
  # 
  # #get date of cancer
  # decoder <- read.table("back_results_200/icd_dates/decoder", stringsAsFactors = F)
  # date_cancer <- as.data.frame(fread(paste0("back_results_200/icd_dates/dates.",
  #                                  which(decoder[,1] == name_cancer) - 1, ".txt")))
  # date_cancer[date_cancer[,2] == 0, 2] <- ""
  # date_cancer <- date_cancer[date_cancer[,1] %in% eid[,1],]
  # date_cancer <- date_cancer[order(date_cancer[,1])[rank(eid[,1])],]
  # colnames(date_cancer) <- c("eid", "end_date")
  # 
  # #Get date of attending the testing center
  # attend_date <- read.table("date_attend_center.txt", stringsAsFactors=F)
  # attend_date <- attend_date[attend_date[,1] %in% eid[,1],]
  # attend_date <- attend_date[order(attend_date[,1])[rank(eid[,1])],]
  # attend_date[,2] <- as.Date(as.character(attend_date[,2]), "%Y%m%d")
  # 
  # date_cancer$pheno <- 0
  # date_cancer$pheno[date_cancer[,2] != ""] <- 1
  # date_cancer[date_cancer[,2] == "", 2] <- "09/30/2020"
  # date_cancer[,2] <- as.Date(date_cancer[,2], "%m/%d/%Y")
  # 
  # death <- read.table("back_results_200/icd_dates/death.txt", stringsAsFactors = F, header = T)
  # death <- death[!duplicated(death$eid),]
  # death$date_of_death <- as.Date(death$date_of_death, "%d/%m/%Y")
  # death <- death[death$eid %in% date_cancer$eid[date_cancer$pheno == 0],]
  # death <- death[order(death$eid)[rank(date_cancer$eid[date_cancer$eid %in% death$eid])],]
  # date_cancer$end_date[date_cancer$eid %in% death$eid] <- death$date_of_death
  # date_cancer$death <- 0
  # date_cancer$death[date_cancer$eid %in% death$eid] <- 1
  # 
  # icd_dates <- read.table("back_results_200/icd_dates/cancer_last.txt.gz", stringsAsFactors = F)
  # icd_eid <- read.table("back_results_200/icd_dates/cancer_eids.txt.gz", stringsAsFactors = F)
  # icd_dates <- icd_dates[seq(1, ncol(icd_dates), 2)]
  # icd_dates <- data.frame(eid = icd_eid[,1], last_date = icd_dates[, which(decoder[,1] == name_cancer)], stringsAsFactors = F)
  # icd_dates <- icd_dates[icd_dates$eid %in% date_cancer$eid,]
  # add_on <- data.frame(eid = date_cancer$eid[!(date_cancer$eid %in% icd_dates$eid)],
  #                      last_date = rep("2021-01-01", nrow(date_cancer)-nrow(icd_dates)),
  #                      stringsAsFactors = F)
  # icd_dates <- rbind(icd_dates, add_on)
  # icd_dates <- icd_dates[order(icd_dates$eid)[rank(date_cancer$eid)],]
  # date_cancer$last_icd <- as.Date(icd_dates$last_date)
  # 
  # 
  # #set up the final, usable data frame
  # date_cancer$prob <- xgb_test_prob[,1]
  # date_cancer$start_date <- attend_date[,2]
  # date_cancer$birth_date <- birth_info$birth
  # date_cancer$age <- as.numeric(Sys.Date() - date_cancer$birth_date)
  # date_cancer$time <- as.numeric(date_cancer$end_date - date_cancer$start)
  # date_cancer$last_icd[date_cancer$last_icd == "2021-01-01"] <- date_cancer$start_date[date_cancer$last_icd == "2021-01-01"]
  # date_cancer$time_from_icd <- NA
  # date_cancer$time_from_icd[date_cancer$pheno == 1] <- date_cancer$end_date[date_cancer$pheno == 1] -
  #   date_cancer$last_icd[date_cancer$pheno == 1]
  # date_cancer <- date_cancer[order(date_cancer$prob),]
  # date_cancer$pheno <- as.factor(date_cancer$pheno)
  # saveRDS(date_cancer, paste0(tort, "_results_", ntree, "/", name_cancer, "/cancer_time.RDS"))
  # 
  # 
  # easy_date_cancer <- date_cancer[c(sample(which(date_cancer$pheno == 0), 1000),
  #                                   which(date_cancer$pheno == 1)),]
  # the_plot <- ggplot(easy_date_cancer, aes(time/365, prob, color = pheno)) +
  #   geom_point() + labs(x = "Years of Study", y = "XGB Probability", color = "Status") +
  #   scale_color_discrete(breaks = c(1, 0), labels = c("Cancer", "Non-Cancer"))
  # plot(the_plot)
  # 
  # 
  # 
  # the_plot <- ggplot(date_cancer[date_cancer$pheno == 1,], aes(time/365, prob)) +
  #   geom_point(aes(color = time_from_icd/365)) + scale_color_viridis() +
  #   labs(x = "Years of Study", y = "XGB Probability", caption = "Cases Only", color = "Time From\nLast ICD")
  # ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/study_time_prob.png"),
  #        the_plot, "png", height=5, width=6)
  # 
  # ################## AGE PLOTS #########################################
  # 
  # the_plot <- ggplot(easy_date_cancer, aes(age/365, prob, color = pheno)) +
  #   geom_point() + labs(x = "Age", y = "XGB Probability", color = "Status") +
  #   scale_color_discrete(labels = c("Cancer", "Non-Cancer"))
  # plot(the_plot)
  # ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/age_time_prob.png"),
  #        the_plot, "png", height=5, width=6)
  # 
  # the_plot <- ggplot(date_cancer[date_cancer$pheno == 1,], aes(age/365, prob)) +
  #   geom_point() + labs(x = "Age", y = "XGB Probability", caption = "Cases Only")
  # plot(the_plot)
  # 
  # 
  # ##################### STRATIFY BY TIME TO LAST ICD ########################
  # early_group <- date_cancer$eid[date_cancer$time_from_icd < quantile(date_cancer$time_from_icd, 0.33, na.rm = T) &
  #                                  !is.na(date_cancer$time_from_icd)]
  # middle_group <- date_cancer$eid[date_cancer$time_from_icd > quantile(date_cancer$time_from_icd, 0.33, na.rm = T) &
  #                                   date_cancer$time_from_icd < quantile(date_cancer$time_from_icd, 0.66, na.rm = T) &
  #                                   !is.na(date_cancer$time_from_icd)]
  # late_group <- date_cancer$eid[date_cancer$time_from_icd > quantile(date_cancer$time_from_icd, 0.66, na.rm = T) &
  #                                 !is.na(date_cancer$time_from_icd)]
  # 
  # early_group <- c(early_group, date_cancer$eid[date_cancer$pheno == 0])
  # middle_group <- c(middle_group, date_cancer$eid[date_cancer$pheno == 0])
  # late_group <- c(late_group, date_cancer$eid[date_cancer$pheno == 0])
  # 
  # auc_items <- multi_auc_plot(list(xgb_test_prob[date_cancer$eid %in% early_group,1, drop=F],
  #                                  xgb_test_prob[date_cancer$eid %in% middle_group,1, drop=F],
  #                                  xgb_test_prob[date_cancer$eid %in% late_group,1, drop=F]),
  #                             list(test_lab[date_cancer$eid %in% early_group,1, drop=F], 
  #                                  test_lab[date_cancer$eid %in% middle_group,1, drop=F],
  #                                  test_lab[date_cancer$eid %in% late_group,1, drop=F]),
  #                             list("XGB: early ICD", "XGB: middle ICD", "XGB: late ICD"))
  # ggsave(paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_group_by_icd.png"),
  #        auc_items[[1]], "png", height=5, width=6)
  # saveRDS(auc_items[[2]], paste0(tort, "_results_", ntree, "/", name_cancer, "/auc_by_icd_group.RDS"))
  
}


#STILL NEED TO DO STOMACH
#"bladder", "breast", "colorectum", "endometrium",
#"kidney", "lung", "pancreas",
#"stomach",

for(x in cancer_to_comp <- c( "bladder", "breast", "colorectum", "endometrium",
                              "kidney", "lung", "pancreas",
                              "stomach","Lymphocytic_leukemia", "melanoma",
                             "Non-Hodgkins_Lymphoma", "prostate")){
 temp_func(x)
}




