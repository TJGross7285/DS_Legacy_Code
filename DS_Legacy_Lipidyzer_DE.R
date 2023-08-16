##################### (10/17/2019) DS Legacy Lipidyzer DE Analysis 

####Set up DS Lipidyzer Data for Analysis
FA_conc<-read.csv("Fatty_Acid_Concentration.csv",check.names=FALSE,header=FALSE,stringsAsFactors=FALSE,na.strings=c(".","0"),sep=",",strip.white=TRUE)
colnames(FA_conc)<-c("UCI ID","Georgetown ID",as.character(FA_conc[4,3:dim(FA_conc)[2]]))
FA_conc<-FA_conc[-c(1:4,151:152),]

LC_conc<-read.csv("Lipid_Class_Concentration.csv",check.names=FALSE,header=FALSE,stringsAsFactors=FALSE,na.strings=c(".","0"),sep=",",strip.white=TRUE)
colnames(LC_conc)<-c("UCI ID","Georgetown ID",as.character(LC_conc[4,3:dim(LC_conc)[2]]))
LC_conc<-LC_conc[-c(1:4,151:152),]

LS_conc<-read.csv("Lipid_Species_Concentrations.csv",check.names=FALSE,header=FALSE,stringsAsFactors=FALSE,na.strings=c(".","0"),sep=",",strip.white=TRUE)
colnames(LS_conc)<-c("UCI ID","Georgetown ID",as.character(LS_conc[4,3:787]))
LS_conc<-LS_conc[-c(1:4,151:152),3:787]

total_fatty_acid<-read.csv("Total_Fatty_Acid.csv",check.names=FALSE,header=FALSE,stringsAsFactors=FALSE,na.strings=c(".","0"),sep=",",strip.white=TRUE)
colnames(total_fatty_acid)<-c("UCI ID","Georgetown ID",as.character(total_fatty_acid[4,3:dim(total_fatty_acid)[2]]))
total_fatty_acid<-total_fatty_acid[-c(1:4,151:152),]

####Check that SampleIDs Match 
all.equal(as.character(FA_conc[,1]),as.character(LC_conc[,1]),as.character(LS_conc[,1]),as.character(total_fatty_acid[,1]))

####Write futher intermediate files 
total<-cbind(FA_conc[,1:2],FA_conc[,-c(1:2)],LC_conc[,-c(1:2)],LS_conc[,-c(1:2)],total_fatty_acid[,-c(1:2)])
total$`UCI ID`<-as.character(total$`UCI ID`)

pheno<-read.csv("pheno_table.csv",check.names=FALSE)
pheno$`UCI ID`<-as.character(pheno$`UCI ID`)

joined<-dplyr::inner_join(pheno,total,by="UCI ID")
na<-apply(joined,2,is.na)
na_index<-apply(na,2,sum)/dim(joined)[1]
joined<-joined[,na_index<.333]

write.csv(joined,file="intermediate.csv")

####Read intermediate file
library(lubridate)
Y<-read.csv("intermediate.csv",check.names=FALSE)
sample_age<-(today()-dmy(Y$`Frozen Date`))/dyears()
Y<-cbind(as.data.frame(sample_age),Y)
pheno<-Y$`Dementia Status`
levels(pheno)<-list("DS-NAD"=c("Non-AD"),"DS-AD"=c("AD"))

####Pull concentration data, log2 transform, and impute NA/Inf values 
library(impute)
abunds<-log2(Y[,-c(1:28)])
T_abunds<-t(abunds)
edata<-impute.knn(T_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
rownames(edata)<-colnames(abunds)

####Calculate SVs 
library(sva)
mod<-model.matrix(~as.factor(pheno))
mod0<-model.matrix(~1,data=pheno)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)
design_table<-cbind(as.data.frame(pheno),as.data.frame(svobj$sv))
colnames(design_table)<-c("pheno","SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8","SV9","SV10","SV11","SV12")

####Save data object for modeling
save(svobj,edata,mod,design_table,file="Legacy_for_fsva_CLUSTER.RData")

####Set up accessory objects for DE incorporating pre-post timepoints SV1+SV2+SV3+SV4+SV5+SV6+SV7+SV8+SV9+SV10+SV11+SV12
design1<-as.data.frame(model.matrix(~0+pheno,data=design_table))
colnames(design1)[1:2]<-c("Control","MCIAD")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cm <- makeContrasts(
	`DSAD-DSNAD` = MCIAD-Control,
	levels=design1)
fit1_F <- contrasts.fit(fit1, cm)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
write.csv(T, file="DS_DE_LegacNoSV.csv")

####Write volcano plot to file 
library(EnhancedVolcano)
colnames(T)[2]<-"Log2FC"
pdf(file="Legacy_Lipidyzer_Volcano_10_17_19.pdf")
EnhancedVolcano(T,lab=T$Feature,x="Log2FC",y="P.Value",xlim=c(-.5,.5),
				ylim=c(0,10),FCcutoff=.2,pLabellingCutoff=.025)
dev.off()



















































####Vocano plot of lipidyzer (ALL METABOLITES)
pdf(file="Volcano_Plot_Lipidyzer_DS.pdf")
plot(table$Log2FC,-log10(table$P.Value), main="All Lipidyzer Metabolites: DS-AD Relative to DS-ND",xlab="Log2FC", ylab="-log10(p-value)")
dev.off()
#######################
#######################

####Filter DE deatures to include FDR<.05
modeling_features<-table%>%filter(adj.P.Val<.05)%>%arrange()
modeling_index<-rownames(expr) %in% T$FeatureNames 
modeling_abunds<-t(expr[modeling_index==TRUE,])

#### Set up Caret control object and outcome vector for classification
library(caret)
library(pROC)
DiseaseState<-pheno
levels(DiseaseState)<-c("Control","Case")
cctrl1 <- trainControl(method = "repeatedcv", repeats=10, number = 10, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)

####Train random forest; 10-fold CV repeated 10 times
set.seed(122)
RandomForest<-train(abunds,DiseaseState, method="rf", trControl = cctrl1, metric = "ROC")

####Generate ROC object for BoostedLogistic performance 
library(pROC)
indexRF<-RandomForest$pred$mtry==2
RF_roc_resampled<-roc(RandomForest$pred$obs[indexRF], RandomForest$pred$Case[indexRF], ci=TRUE, auc=TRUE)

####Train boosted logistic regression; 10-fold CV repeated 10 times
boostGrid<-expand.grid(nIter=28)
set.seed(122)
BoostedLogistic<-train(modeling_abunds,DiseaseState, method="LogitBoost", trControl = cctrl1, metric = "ROC",tuneGrid=boostGrid)

####Generate ROC object for BoostedLogistic performance 
library(pROC)
indexBoosted<-BoostedLogistic$pred$nIter==28
Boosted_roc_resampled<-roc(BoostedLogistic$pred$obs[indexBoosted], BoostedLogistic$pred$Case[indexBoosted], ci=TRUE, auc=TRUE)

####Generate pdf of results 
pdf(file="DS_Models_Lipidyzer.pdf")
###Random Forest performance over resampling
rf_resampled<-plot.roc(RF_roc_resampled, ylab="True Positive Rate", xlab="False Positive Rate", main="Random Forest Across Resampling",  ci=TRUE, print.auc=TRUE, legacy.axes=TRUE, print.auc.x=ifelse(RF_roc_resampled$percent, 50, .45),
print.auc.y=ifelse(RF_roc_resampled$percent, 50, .45))
rf_res_ciobj <- ci.se(rf_resampled, specificities= seq(0, 1, .01))
plot(rf_res_ciobj, type="shape", col="cornflowerblue")
###Plot variable importance
plot(varImp(RandomForest),top=25)

###LogitBoost performance over resampling
logitboost_resampled<-plot.roc(Boosted_roc_resampled, ylab="True Positive Rate", xlab="False Positive Rate", main="Boosted Logistic Regression Across Resampling",  ci=TRUE, print.auc=TRUE, legacy.axes=TRUE, print.auc.x=ifelse(Boosted_roc_resampled$percent, 50, .45),
print.auc.y=ifelse(Boosted_roc_resampled$percent, 50, .45))
logitb_res_ciobj <- ci.se(logitboost_resampled, specificities= seq(0, 1, .01))
plot(logitb_res_ciobj, type="shape", col="cornflowerblue")
###Plot variable importance
plot(varImp(BoostedLogistic),top=25)
dev.off()

###Write boosted logistic regression predictors to file
write.csv(predictors(BoostedLogistic),file="Predictors_LogitBoost.csv")



#####################################
#####################################

####Conduct correlations against clinical variables using unadjusted data
colnames(Y)[6]<-"Sex"
subject_var<-Y[,c(1,6,7,10:27)]
subject_var$Sex<-as.numeric(subject_var$Sex)
index<-caret::nearZeroVar(subject_var)
subject_var<-subject_var[,-index]
subject_var<-subject_var[,c(1,3:7)]
subject_var$`Level of ID`<-as.numeric(subject_var$`Level of ID`)
subject_var$`Level of ID`[108]<-NA
correlate_clin<-stats::cor(t(expr),subject_var,method="spearman",use="complete.obs")

pdf(file="Correlation_Plots.pdf")
gplots::heatmap.2(correlate_clin,trace="none")
dev.off()

write.csv(correlate_clin,file="DS_CORR.csv")







library(caret)
library(pROC)

pdf(file="Candidate_Plots.pdf")
A<-roc(response=pheno,predictor=modeling_abunds$`PC(18:1/22:5)`,ci=TRUE,auc=TRUE)
B<-roc(response=pheno,predictor=modeling_abunds$`PC(18:1/20:5)`,ci=TRUE,auc=TRUE)

cctrl1 <- trainControl(method = "repeatedcv", repeats=10, number = 10, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
DiseaseState<-pheno
levels(DiseaseState)<-c("Control","Case")
set.seed(122)
CandidateRF<-train(modeling_abunds[colnames(modeling_abunds)=="PC(18:1/20:5)"|colnames(modeling_abunds)=="PC(18:1/22:5)"],DiseaseState, method="rf", trControl = cctrl1, metric = "ROC")
indexRF<-CandidateRF$pred$mtry==2
RF_roc_resampled<-roc(CandidateRF$pred$obs[indexRF], CandidateRF$pred$Case[indexRF], ci=TRUE, auc=TRUE)

plot.roc(A, ylab="True Positive Rate", xlab="False Positive Rate", main="PC(18:1/22:5)",  ci=TRUE, print.auc=TRUE, legacy.axes=TRUE, print.auc.x=ifelse(A$percent, 50, .45),
print.auc.y=ifelse(A$percent, 50, .45))
logitb_res_ciobj <- ci.se(A, specificities= seq(0, 1, .01))
plot(logitb_res_ciobj, type="shape", col="cornflowerblue")

plot.roc(B, ylab="True Positive Rate", xlab="False Positive Rate", main="PC(18:1/20:5)",  ci=TRUE, print.auc=TRUE, legacy.axes=TRUE, print.auc.x=ifelse(B$percent, 50, .45),
print.auc.y=ifelse(B$percent, 50, .45))
logitb_res_ciobj <- ci.se(B, specificities= seq(0, 1, .01))
plot(logitb_res_ciobj, type="shape", col="cornflowerblue")

plot.roc(RF_roc_resampled, ylab="True Positive Rate", xlab="False Positive Rate", main="PC(18:1/22:5) and PC(18:1/20:5)",  ci=TRUE, print.auc=TRUE, legacy.axes=TRUE, print.auc.x=ifelse(RF_roc_resampled$percent, 50, .45),
print.auc.y=ifelse(RF_roc_resampled$percent, 50, .45))
logitb_res_ciobj <- ci.se(RF_roc_resampled, specificities= seq(0, 1, .01))
plot(logitb_res_ciobj, type="shape", col="cornflowerblue")

dev.off()




R<-read.csv("Worked.csv",check.names=FALSE)
up<-R%>%filter(log2FoldChange>.75 & padj<.05)%>%select(gene)
down<-R%>%filter(log2FoldChange< -.75 & padj<.05)%>%select(gene)
write.table(up, file = "up_RAS.txt", sep = "\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(down, file = "down_RAS.txt", sep = "\t",row.names=FALSE,col.names=FALSE,quote=FALSE)




library(dplyr)
library(readr)
library(tidytext)
#######library(cba)

R_down<-read_lines("Down.txt")
text_df_down<- tibble(line = 1:length(R_down), text = R_down)
F_down<-substr(text_df_down$text,1,25)
tidy_books_down<- text_df_down %>% unnest_tokens(word, text)
down<-tidy_books_down %>% count(word, sort = TRUE)
write.csv(down%>%filter(n>1),file="Down_Terms.csv")

R_up<-read_lines("up_fixed.txt")
text_df_up <- tibble(line = 1:length(R_up), text = R_up)
F_up<-substr(text_df_up$text,1,25)
tidy_books_up<- text_df_up %>% unnest_tokens(line,word, text)
up<-tidy_books_up %>% count(word, sort = TRUE)
write.csv(up%>%filter(n>1),file="Up_Terms.csv")



# define a nice color palette
pal <- brewer.pal(8,"Dark2")

# plot the 50 most common words
tokens_clean %>% 
  with(wordcloud(word, n, random.order = FALSE, max.words = 50, colors=pal))















up_dtm<- down%>%
  cast_dtm(line, word, n)

index<-caret::nearZeroVar(dtm)
clean_dtm<-dtm[,-index]
clean_dtm<-as.matrix(dtm)

set.seed(122)
clusterM<-rockCluster(clean_dtm[1:150,],n=10)
fitted<-fitted(clusterM)




