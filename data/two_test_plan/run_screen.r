library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(devtools)
load_all("/dcs04/scharpf/data/annapragada/useful.stuff.aa")
library(caret)
library(recipes)
library(pROC)
library(data.table)
delfi<-fread("delfi_features.csv") %>% select(-starts_with("cov"))
meta<-fread("Clinical_Metadata_spreadsheet_8_11.csv")
meta<-meta %>% mutate(type=if_else(HCCStatus=="Yes","cancer","healthy"))
meta$grp_id<-sapply(strsplit(meta$id,"P"),"[",1)

d_c<-meta %>% filter(forveryclean==1 & Disease == "Non-cancer")
d_c$cat<-"Danish Cocos"
#cirr<-meta %>% filter(forveryclean==1 & Disease == "Cirrhosis")
#cirr$cat<-"Cirrhosis"
#hbv<-meta %>% filter(forveryclean==1 & Disease == "HBV")
#hbv$cat<-"HBV"
#hcv<-meta %>% filter(forveryclean==1 & Disease == "HCV")# & `HIV status`==0)
#hcv$cat<-"HCV"
hcc_naive<-meta %>% filter(forveryclean==1 & Disease == "HCC" & `Treatment before blood draw`==0)
hcc_naive$cat<-"HCC Tx Naive"
hcc_tx<-meta %>% filter(forveryclean==1 & Disease == "HCC" & `Treatment before blood draw`==1)
hcc_tx$cat<-"HCC Prior Tx"

d1<-rbind(hcc_naive,hcc_tx,d_c)

d1_delfi<-inner_join(delfi,d1 %>% select(id,type),by="id")
library(data.table)
tfbs<-fread("cohort_rel_cov2.txt")
d1_tf<-inner_join(tfbs,d1 %>% select(id,type),by="id")
#PCs from the TFs have individual AUCs around .55
#175 Chipseq individually with AUC greater than .7


ref<-fread("TF_names_ref.csv") %>% select(-V1)
ref<-ref %>% filter(liv != "other")
#liver_TFs<-read_csv("../data/Liver_COSMIC.csv")
#liver_TFs$perc<-liver_TFs$`Mutated samples`/liver_TFs$`Samples tested`
#TF_sel<-liver_TFs %>% filter(perc > .1)
#TF_sel$name<-sapply(str_split(TF_sel$`Gene name`,"_"),"[",1)
#TF_sel<-ref %>% filter(TF %in% TF_sel$name)
d1_tf_sel<-d1_tf %>% select(id,type,ref$stat)


d1_delfi<-inner_join(d1_delfi,d1_tf_sel %>% select(-type))


#### MODELING
recipe_seq <- recipe(type ~ ., data=d1_delfi) %>%
    update_role(id, new_role = "ID") %>%
    step_pca(starts_with("ratio"), prefix = "ratio_pc_",threshold=.9)     %>%
    step_corr(all_predictors(), threshold=0.95) %>%
    step_nzv(all_predictors())


glmnetGrid <- expand.grid(
    alpha = 1,
    lambda = 10^seq(-5, -1, length.out = 100))
#### Train models
set.seed(1234)
ctrl_all <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 10,
                     verboseIter = TRUE,
                     savePredictions="final",
                     classProbs=TRUE,
                     index=createMultiFolds(d1_delfi$type, 5, 10),
                     summaryFunction = twoClassSummary)

set.seed(1234)
model_delfi <- caret::train(recipe_seq,
                          data = d1_delfi,
                          method = "glmnet",
                          tuneGrid = glmnetGrid,
                          trControl = ctrl_all)
set.seed(1234)
model_delfi_gbm<-caret::train(recipe_seq,
                          data = d1_delfi,
                          method = "gbm",
                          trControl = ctrl_all)

#### RESULTS
pred.all <- model_delfi$pred
pred.all <- pred.all %>% dplyr::group_by(rowIndex) %>% dplyr::summarize(score.delfi.lr = mean(cancer))
all_data <- d1_delfi %>% dplyr::mutate(rowIndex = 1:n())
labels <- all_data %>% dplyr::select(rowIndex,id,type)
preds_all <- dplyr::inner_join(labels, pred.all, by="rowIndex")
labels <- d1 %>% dplyr::select(id, BCLC,Disease,grp_id,`Draw_number`,cat,`Treatment before blood draw`)
preds_all <- dplyr::inner_join(preds_all,labels, by="id")

pred.all <- model_delfi_gbm$pred
pred.all <- pred.all %>% dplyr::group_by(rowIndex) %>% dplyr::summarize(score.delfi.gbm = mean(cancer))
all_data <- d1_delfi %>% dplyr::mutate(rowIndex = 1:n())
labels <- all_data %>% dplyr::select(rowIndex,id,type)
preds_all2 <- dplyr::inner_join(labels, pred.all, by="rowIndex")
labels <- d1 %>% dplyr::select(id)
preds_all2 <- dplyr::inner_join(preds_all2,labels, by="id")

preds_all<-inner_join(preds_all,preds_all2 %>% select(-type,-rowIndex),by="id")

preds_final <- preds_all %>% select(-rowIndex)

write.csv(preds_final,"delfi-tf-results_screen.csv")
saveRDS(model_delfi,"delfi-tf-LR_screen.rds")
saveRDS(model_delfi_gbm,"delfi-tf-GBM_screen.rds")







