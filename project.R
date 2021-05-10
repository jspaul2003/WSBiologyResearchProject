#GLOBAL OPTIONS:

#This is done to get rid of unimportant warnings
options(warn=-1)
#To turn it back on run the following:
#options(warn=0)


#PACKAGES:

#install.packages("plyr")
library(plyr)
#install.packages("fitdistrplus")
library(fitdistrplus)
#install.packages("utils")
library(utils)
#install.packages("TSstudio")
library(TSstudio)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("chron")
library(chron)
#install.packages("expss")
library(expss)
#install("ggplot2")
library(ggplot2)
#install.packages("reshape2")
library(reshape2)
#install.packages("stats4")
library(stats4)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("latticeExtra")
library(latticeExtra)
#install.packages("Hmisc")     
library(Hmisc)
#install.packages("sem")
library(sem)
#install.packages("rgl")
library(rgl)
#install.packages("multcomp")
library(multcomp)
#install.packages("leaps")
library(leaps)
#install.packages("aplpack")
library(aplpack)
#install.packages("Rcmdr")
#library(Rcmdr)
#We will use the Regclass package on Rstudio cloud since
#Rcmdr wont work for some reason- probably to do with the
#GUI interface support
#install.packages("regclass")
library(regclass)
#install.packages("MASS")
library(MASS)
#install.packages("car")
library(car)
#install.packages("quantmod")
library(quantmod)
#install.packages("nnet")
library(nnet)
#install.packages("neuralnet")
library(neuralnet)
#install.packages("glmnet")
library(glmnet)
#install.packages("miscTools")
library(miscTools)
#install.packages("sandwich")
library(sandwich)
#install.packages("actuar")
library(actuar)
#install.packages("DescTools")
library(DescTools)

#SETTING SEED FOR CONSISTENT RESULTS
set.seed(2.71828)

#Converting Categorical Variables to Numerical

METABRIC_RNA_Mutation$Numeric_Breast_Surgery_Type <- as.numeric(factor(
  METABRIC_RNA_Mutation$type_of_breast_surgery,
  levels = c("MASTECTOMY","BREAST CONSERVING"))) - 1

METABRIC_RNA_Mutation$Numeric_Cancer_Type <- as.numeric(factor(
  METABRIC_RNA_Mutation$cancer_type,
  levels = c("Breast Cancer","Breast Sarcoma"))) - 1

METABRIC_RNA_Mutation$Numeric_Cancer_Type_Detailed <- as.numeric(factor(
  METABRIC_RNA_Mutation$cancer_type_detailed,
  levels = c("Breast Invasive Lobular Carcinoma","Breast Invasive Ductal Carcinoma",
             "Breast Mixed Ductal and Lobular Carcinoma","Breast","Metaplastic Breast Cancer"))) - 1

METABRIC_RNA_Mutation$Numeric_Cellularity <- as.numeric(factor(
  METABRIC_RNA_Mutation$cellularity,
  levels = c("High","Moderate","Low"))) - 1

METABRIC_RNA_Mutation$Numeric_Pam50_Claudin_Low_Subtype <- as.numeric(factor(
  METABRIC_RNA_Mutation$`pam50_+_claudin-low_subtype`,
  levels = c("claudin-low","LumA","LumB","Her2","Basal","Normal"))) - 1

METABRIC_RNA_Mutation$Numeric_Er_Status_Measured_By_ihc <- as.numeric(factor(
  METABRIC_RNA_Mutation$er_status_measured_by_ihc,
  levels = c("Positive","Negative"))) - 1

METABRIC_RNA_Mutation$Numeric_Er_Status <- as.numeric(factor(
  METABRIC_RNA_Mutation$er_status,
  levels = c("Positive","Negative"))) - 1

METABRIC_RNA_Mutation$Numeric_Her2_Status_Measured_By_snp6 <- as.numeric(factor(
  METABRIC_RNA_Mutation$her2_status_measured_by_snp6,
  levels = c("NEUTRAL","LOSS","GAIN","UNDEF"))) - 1

METABRIC_RNA_Mutation$Numeric_Her2_Status <- as.numeric(factor(
  METABRIC_RNA_Mutation$her2_status,
  levels = c("Positive","Negative"))) - 1

METABRIC_RNA_Mutation$Numeric_Tumor_Other_Histologic_Subtype <- as.numeric(factor(
  METABRIC_RNA_Mutation$tumor_other_histologic_subtype,
  levels = c("Ductal/NST","Mixed","Lobular","Tubular/ cribriform","Mucinous",
             "Medullary","Other","Metaplastic"))) - 1

METABRIC_RNA_Mutation$Numeric_Inferred_Menopausal_State <- as.numeric(factor(
  METABRIC_RNA_Mutation$inferred_menopausal_state,
  levels = c("Post","Pre"))) - 1

METABRIC_RNA_Mutation$Numeric_Primary_Tumor_Laterality <- as.numeric(factor(
  METABRIC_RNA_Mutation$primary_tumor_laterality,
  levels = c("Left","Right"))) - 1

METABRIC_RNA_Mutation$Numeric_Oncotree_Code <- as.numeric(factor(
  METABRIC_RNA_Mutation$oncotree_code,
  levels = c("IDC","BREAST","IMMC","MDLC","ILC"))) - 1

METABRIC_RNA_Mutation$Numeric_Pr_Status <- as.numeric(factor(
  METABRIC_RNA_Mutation$pr_status,
  levels = c("Positive","Negative"))) - 1

METABRIC_RNA_Mutation$Numeric_3_Gene_Classifier_Subtype <- as.numeric(factor(
  METABRIC_RNA_Mutation$`3-gene_classifier_subtype`,
  levels = c("ER+/HER2- High Prolif","ER-/HER2-","HER+","ER+/HER2- Low Prolif","HER2+"))) - 1

METABRIC_RNA_Mutation$Numeric_Death_from_Cancer <- as.numeric(factor(
  METABRIC_RNA_Mutation$death_from_cancer,
  levels = c("Living","Died of Disease","Died of Other Causes"))) - 1

#vif
#Finds Variance Inflation Factors. For use in Rstudio Cloud
#since it does not support RCmdr. Function is calling regclass
#VIF function. Being used to avoid rewriting a lot of code below.
vif=function(fit){
  VIF(fit)
}

#dropvif1:
#For use in model building section. Takes linear model
#Finds largest vif. If greater than 10, it drops the relevant
#variable. Returns data with dropped variable

dropvif1=function(fit,train.data){
  vifs=vif(fit)
  maxvif=which.max(vifs)
  print(maxvif)
  if(vifs[[maxvif[[1]]]]>10){
    print(paste("Dropped",colnames(train.data)[maxvif[[1]]]))
    train.data <- subset(train.data, select=-c(maxvif[[1]]))
  }
  return(train.data)
}

#normalize:
#Normalizes data for neural net models
normalize=function(x){
  return((x-min(x))/(max(x) -min(x) ))
}

#testlassoridge1:
#runs diagnostics for normality of residuals and plots residual
#fitted values plots. 
testlassoridge1=function(fits){
  residuals=train.data$overall_survival_months-fits
  plot(fits[1:length(fits)]~residuals[1:length(residuals)])
  print("Shapiro Wilk for Normality of residuals")
  print(shapiro.test(residuals[1:length(residuals)]))
}


METABRIC_RNA_Mutation <- read_csv("METABRIC_RNA_Mutation.csv")

nums <- unlist(lapply(METABRIC_RNA_Mutation, is.numeric))  
temp=METABRIC_RNA_Mutation[,nums]
temp$mutation_nottingham_prognostic_index=I(log(temp$nottingham_prognostic_index))
temp$tumor_size=I(log(temp$tumor_size))
temp$lymph_nodes_examined_positive=I(sqrt(temp$lymph_nodes_examined_positive))
temp$age_at_diagnosis=I(log(temp$age_at_diagnosis))
temp$mutation_count=I(log(temp$mutation_count))
temp=subset(temp,select =-c(patient_id))

temp= temp %>% select(-overall_survival_months, everything())

#boxplot(subset(temp,select=-c(overall_survival_months)))

sample <-sample.int(n =nrow(temp),size =floor(.80*nrow(temp)), replace = F)
train.data <- temp[sample, ]
test.data <-temp[-sample,]
train.data[sapply(train.data, is.infinite)] <- NA
train.data=drop_na(train.data)
test.data=drop_na(test.data)

#model 1
fit=lm(overall_survival_months~.,data=train.data)
fit1=fit
summary(fit1)
plot(fit1)
ncvTest(fit1)

#We have linearly related variables- bad

attributes(alias(fit)$Complete)$dimnames[[1]]

train.data=subset(train.data,select=-c(hras_mut,siah1_mut,smarcb1_mut, rasgef1b_mut))

#run model again
fit=lm(overall_survival_months~.,data=train.data)
fit1=fit
#vif(fit)

#deal with multicollinearity
length=1
while(length!=ncol(train.data)){
  length=ncol(train.data)
  train.data=dropvif1(fit,train.data)
  fit=lm(overall_survival_months~. ,data=train.data)
}

#Model 2
fit2=lm(overall_survival_months~. ,data=train.data); summary(fit2)
plot(fit2)
ncvTest(fit2)

#Model 3
fit3=stepwise(fit2,direction='forward/backward',criterion='AIC',trace='false'); summary(fit3)
plot(fit3)
ncvTest(fit3)

#Model 4
fit4=stepwise(fit2,direction='forward/backward',criterion='BIC',trace='false'); summary(fit4)
plot(fit4)
ncvTest(fit4)

#model 5
fit5=lm(overall_survival_months~.^2,data=train.data); summary(fit5)
plot(fit5)
ncvTest(fit5)

#model 6
fit6=stepwise(fit4,direction='forward/backward',criterion='AIC',trace='false'); summary(fit6)
plot(fit6)
ncvTest(fit6)

#model 7
fit7=stepwise(fit4,direction='forward/backward',criterion='BIC',trace='false'); summary(fit7)
plot(fit7)
ncvTest(fit7)

#model 8, LASSO Model
X = model.matrix(overall_survival_months~.,train.data )[,-1]
Y=train.data$overall_survival_months
cv=cv.glmnet(X,Y,alpha =1)
model=glmnet(X,Y,alpha =1, lambda=cv$lambda.min)

X1=cbind(1,X)
testlassoridge1(X1%*%coef(model));

#model 9, Ridge Model
cv=cv.glmnet(X,Y,alpha =0)
model2=glmnet(X,Y,alpha =0, lambda=cv$lambda.min)

testlassoridge1(X1%*%coef(model2))

#model 10, neural net
#we need to normalize the data
ntrain.data=as.data.frame(lapply(train.data,normalize))
nn=neuralnet(overall_survival_months~.-hla.g , data=ntrain.data, hidden=c(16,25), linear.output =T, threshold=0.5)

#Checking our models
#we also need to remove the columns in testing data (for testing lasso and ridge)
#that we removed in training
test.data=test.data[which(colnames(test.data)%in%colnames(train.data))]
sse1=sum(test.data$overall_survival_months-predict(fit,new=test.data))^2
sse2=sum(test.data$overall_survival_months-predict(fit2,new=test.data))^2
sse3=sum(test.data$overall_survival_months-predict(fit3,new=test.data))^2
sse4=sum(test.data$overall_survival_months-predict(fit4,new=test.data))^2
sse5=sum(test.data$overall_survival_months-predict(fit5,new=test.data))^2
sse6=sum(test.data$overall_survival_months-predict(fit6,new=test.data))^2

X.test = model.matrix(overall_survival_months~.,test.data )
fits = X.test%*%coef(model)
sse8=sum(test.data$overall_survival_months-fits)^2

fits2 =X.test%*%coef(model2)
sse9=sum(test.data$overall_survival_months-fits2)^2

nn.results <- compute(nn, test.data)
sse10=sum(test.data$overall_survival_months-nn.results$net.result)^2

