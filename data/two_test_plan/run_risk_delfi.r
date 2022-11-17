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

#d_c<-meta %>% filter(forveryclean==1 & Disease == "Non-cancer")
#d_c$cat<-"Danish Cocos"
cirr<-meta %>% filter(forveryclean==1 & Disease == "Cirrhosis")
cirr$cat<-"Cirrhosis"
hbv<-meta %>% filter(forveryclean==1 & Disease == "HBV")
hbv$cat<-"HBV"
hcv<-meta %>% filter(forveryclean==1 & Disease == "HCV")# & `HIV status`==0)
hcv$cat<-"HCV"
hcc_naive<-meta %>% filter(forveryclean==1 & Disease == "HCC" & `Treatment before blood draw`==0)
hcc_naive$cat<-"HCC Tx Naive"
hcc_tx<-meta %>% filter(forveryclean==1 & Disease == "HCC" & `Treatment before blood draw`==1)
hcc_tx$cat<-"HCC Prior Tx"

#d1<-rbind(hcc_naive,hcc_tx,cirr,hbv,d_c,hcv)
d1<-rbind(hcc_naive,hcc_tx,cirr,hbv,hcv)

#d1<-rbind(hcc_naive,hcc_tx,cirr,hbv,d_c)

d1_delfi<-inner_join(delfi,d1 %>% select(id,type),by="id")
library(data.table)


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

write.csv(preds_final,"delfi-results_risk.csv")
saveRDS(model_delfi,"delfi-LR_risk.rds")
saveRDS(model_delfi_gbm,"delfi-GBM_risk.rds")







