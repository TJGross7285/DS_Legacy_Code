library(impute)
library(limma)
library(sva)
library(dplyr)
library(GridOnClusters)
library(bnlearn)

####Read in abundance data and subset to final table of features 
read_in_neg_lip<-read.csv("joined_lipo_neg.csv",check.names=FALSE, na.strings=c("NA",".","","0"))[,-c(1:26)]
read_in_neg_met<-read.csv("joined_metab_neg.csv",check.names=FALSE, na.strings=c("NA",".","","0"))[,-c(1:26)]
read_in_pos_lip<-read.csv("joined_lipo_pos.csv",check.names=FALSE, na.strings=c("NA",".","","0"))[,-c(1:26)]
read_in_pos_met<-read.csv("joined_metab_pos.csv",check.names=FALSE, na.strings=c("NA",".","","0"))[,-c(1:26)]
all.equal(read_in_neg_lip$SampleID,read_in_neg_met$SampleID,
	      read_in_pos_lip$SampleID,read_in_pos_met$SampleID)

legacy_pos<-cbind(read_in_pos_lip,read_in_pos_met)
legacy_pos_mz<-gsub("_.+$","",colnames(legacy_pos))
legacy_pos_rt<-gsub("^.+_","",colnames(legacy_pos))
legacy_neg<-cbind(read_in_neg_lip,read_in_neg_met)
legacy_neg_mz<-gsub("_.+$","",colnames(legacy_neg))
legacy_neg_rt<-gsub("^.+_","",colnames(legacy_neg))

Chemistry<-as.factor(c(rep("Lipidomics",dim(read_in_pos_lip)[2]),
	                 rep("Metabolomics",dim(read_in_pos_met)[2]),
	                 rep("Lipidomics",dim(read_in_neg_lip)[2]),
	                 rep("Metabolomics",dim(read_in_neg_met)[2])))         

meta_pheno<-read.csv("joined_lipo_neg.csv",check.names=FALSE)[,c(1:26)]


meta_abunds<-cbind(as.data.frame(Chemistry),
	               rbind(cbind(legacy_pos_mz,legacy_pos_rt),cbind(legacy_neg_mz,legacy_neg_rt)))
abunds<-cbind(legacy_pos,legacy_neg)
abunds[abunds==0]<-NA
Mode<-as.factor(c(rep("positive",dim(legacy_pos)[2]),rep("negative",dim(legacy_neg)[2])))

####Impute missing data with KNN imputation
library(impute)
T_abunds<-t(abunds)
edata<-impute.knn(T_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
edata<-log2(apply(edata,2,as.numeric))
rownames(edata)<-rownames(T_abunds)

feature_meta<-cbind(colnames(abunds),as.data.frame(Mode),as.data.frame(meta_abunds))
feature_meta_use<-cbind(rownames(T_abunds),feature_meta)
colnames(feature_meta_use)[1]<-"Feature"


####Conduct SVA over samples 
pheno_table<-meta_pheno[,c(3,8)]
colnames(pheno_table)<-c("StudyID","DiseaseState")
mod<-model.matrix(~as.factor(DiseaseState),data=pheno_table)
mod0<-model.matrix(~1,data=pheno_table)
n.sv_BE<-num.sv(edata,mod,seed=122)
n.sv_LEEK<-num.sv(edata,mod,seed=122,method="leek")
svobj_LEEK<-sva(edata,mod,mod0,n.sv=n.sv_LEEK)$sv
svobj_BE<-sva(edata,mod,mod0,n.sv=n.sv_BE)$sv
total_LEEK<-cbind(as.data.frame(pheno_table),as.data.frame(svobj_LEEK))
total_BE<-cbind(as.data.frame(pheno_table),as.data.frame(svobj_BE))
colnames(total_LEEK)[1:2]<-c("SampleID","Main")
colnames(total_BE)[1:2]<-c("SampleID","Main")
design_table_BE<-cbind(model.matrix(~0+as.factor(Main),data=total_BE),as.data.frame(svobj_BE))
colnames(design_table_BE)[1:2]<-c("AD","Non_AD")

design_table_LEEK<-cbind(model.matrix(~0+as.factor(Main),data=total_LEEK),as.data.frame(svobj_LEEK))
colnames(design_table_LEEK)[1:2]<-c("AD","Non_AD")


#####`MCIAD-Converter`= MCIAD-Converter,
					#`Super-Control`= Super-Control,
					#`Super-MCIAD`= Super-MCIAD,
					#`Super-Converter`= Super-Converter,
					#`MCIAD-Control`= MCIAD-Control,
####Fit linear model for pairwise contrasts 
arrayw<-arrayWeights(edata, design=design_table_BE)  ########design_table_LEEK
fit1<-lmFit(edata,design_table_BE,weights=arrayw)
cm1 <- makeContrasts(`AD-Non_AD`= AD-`Non_AD`,
					 levels=design_table_BE)

fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-c("Feature")
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(5,6)]<-c("MZ","RT")

write.table(joinT%>%filter(Mode=="positive")%>%select("MZ","RT","P.Value","AD.Non_AD"),
			sep="\t",file="DS_Legacy_POS_BE.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="negative")%>%select("MZ","RT","P.Value","AD.Non_AD"),
			sep="\t",file="DS_Legacy_NEG_BE.txt",row.names=FALSE)

cd /Volumes/NO\ NAME/Active_Projects_4_9_2021/UCI\ Legacy  
mummichog -f DS_Legacy_POS_BE.txt -o DS_Legacy_POS_BE -m positive 
mummichog -f DS_Legacy_NEG_BE.txt -o DS_Legacy_NEG_BE -m negative 

###########PIUMet BE
write.table(joinT%>%filter(Chemistry=="Lipidomics")%>%mutate(`Prize`= -log10(P.Value))%>%select("MZ","Mode","Prize"),
			sep="\t",file="DS_Legacy_BE_PIUmet_LIP.txt",row.names=FALSE)
write.table(joinT%>%filter(Chemistry=="Metabolomics")%>%mutate(`Prize`= -log10(P.Value))%>%select("MZ","Mode","Prize"),
			sep="\t",file="DS_Legacy_BE_PIUmet_METAB.txt",row.names=FALSE)




arrayw<-arrayWeights(edata, design=design_table_LEEK)  ########design_table_LEEK
fit1<-lmFit(edata,design_table_LEEK,weights=arrayw)
cm1 <- makeContrasts(`AD-Non_AD`= AD-`Non_AD`,
					 levels=design_table_LEEK)

fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(5,6)]<-c("MZ","RT")

write.table(joinT%>%filter(Mode=="positive")%>%select("MZ","RT","P.Value","AD.Non_AD"),
			sep="\t",file="DS_Legacy_POS_LEEK.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="negative")%>%select("MZ","RT","P.Value","AD.Non_AD"),
			sep="\t",file="DS_Legacy_NEG_LEEK.txt",row.names=FALSE)

cd /Volumes/NO\ NAME/Active_Projects_4_9_2021/UCI\ Legacy  
mummichog -f DS_Legacy_POS_LEEK.txt -o DS_Legacy_POS_LEEK -m positive 
mummichog -f DS_Legacy_NEG_LEEK.txt -o DS_Legacy_NEG_LEEK -m negative 

###########PIUMet LEEK
write.table(joinT%>%filter(Chemistry=="Lipidomics")%>%mutate(`Prize`= -log10(P.Value))%>%select("MZ","Mode","Prize"),
			sep="\t",file="DS_Legacy_LEEK_PIUmet_LIP.txt",row.names=FALSE)
write.table(joinT%>%filter(Chemistry=="Metabolomics")%>%mutate(`Prize`= -log10(P.Value))%>%select("MZ","Mode","Prize"),
			sep="\t",file="DS_Legacy_LEEK_PIUmet_METAB.txt",row.names=FALSE)





#######BayesNet
meta_pheno<-read.csv("joined_lipo_neg.csv",check.names=FALSE)[,c(1:26)]
colnames(meta_pheno)[c(5)]<- c("Sex")
apoePheno<-read.csv("Mapstone Lipidomics Data2 10.21.18 with ApoE 7.2.21.csv",check.names=FALSE)
colnames(apoePheno)[c(1)]<- c("UCI ID")
joined<-dplyr::inner_join(apoePheno,meta_pheno,by="UCI ID")[,-c(6)]%>%filter(`Data Source`!= "Missing" &`Data Source`!= "NEED DATA CHECK")

LegLeek<-read.csv("Legacy_LEEK.csv",check.names=FALSE)
LegBE<-read.csv("Legacy_BE.csv",check.names=FALSE)
joined_LEEK<-dplyr::inner_join(LegLeek,joined,by="SampleID")
joined_BE<-dplyr::inner_join(LegBE,joined,by="SampleID")


pheno_factors<-apply(joined_BE[,c(21,25,26,33:46)],2,as.factor)
pheno_cont_BE<-apply(joined_BE[,c(4:18,24,29:32)],2,as.numeric)
pheno_cont_LEEK<-apply(joined_LEEK[,c(4,10,15:18)],2,as.numeric)
discrete_BE<-cbind(as.data.frame(discretize.jointly(pheno_cont_BE,k=c(2:50))$D),as.data.frame(pheno_factors))
discrete_BE<-discrete_BE[,-caret::nearZeroVar(discrete_BE)]
discrete_LEEK<-cbind(as.data.frame(discretize.jointly(pheno_cont_LEEK,k=c(2:50))$D),as.data.frame(pheno_factors))
discrete_LEEK<-discrete_LEEK[,-caret::nearZeroVar(discrete_LEEK)]
discrete_BE[,1:dim(discrete_BE)[2]]<-lapply(discrete_BE[,1:dim(discrete_BE)[2]],as.factor)
discrete_LEEK[,1:dim(discrete_LEEK)[2]]<-lapply(discrete_LEEK[,1:dim(discrete_LEEK)[2]],as.factor)


G<-h2pc(discrete_BE)
H<-arc.strength(G,discrete_BE)
strength.plot(G,H)


G<-h2pc(discrete_LEEK)
H<-arc.strength(G,discrete_LEEK)
strength.plot(G,H)



























####Carry out SVA
edata<-t(abunds)
mod<-model.matrix(~as.factor(DiseaseState), data=meta)
mod0<-model.matrix(~1,data=meta)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)
modSv<-cbind(mod,svobj$sv)
mod0Sv<-cbind(mod0,svobj$sv)
pValuesSv<-f.pvalue(edata,modSv,mod0Sv)
qValuesSv<-p.adjust(pValuesSv,method="fdr")

median_case<-apply(as.data.frame(edata[,meta$DiseaseState=="AD"]),1,median)
median_control<-apply(edata[,meta$DiseaseState=="Non-AD"],1,median)
log2FC<-log2(median_case/median_control)
mz<-gsub("_.+$","",rownames(edata))
rt<-gsub("^.+_","",rownames(edata))
table<-cbind(as.data.frame(mz),as.data.frame(rt),as.data.frame(pValuesSv),as.data.frame(qValuesSv),as.data.frame(log2FC))
colnames(table)<-c("m/z","retention time","PValue","QValue","Log2FC")
write.csv(table, file="SVA_DE.csv")