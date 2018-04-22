# init----
resp="KWHSPC"
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load(paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/final.Rdata'))

# include----
library(corrplot)
library(gam)
library(ModelMetrics)
library(randomForest)
library(e1071)
library(rpart)
library(rpart.plot)
# functions----
completeness<-function(dat)
{
  dat<-as.data.frame(sapply(dat,as.character))
  dat[dat=="."]<-NA
  a=nrow(dat)
  b=length(dat)
  c=sum(is.na(dat))/(a*b)
  return(round(100*(1-c),digits=2))
}
crossValidate<-function(cvtype,folds,dataset,model,resp)
{
  df<-dataset
  l <- vector("list", 2)
  if (cvtype=="kfold")
  {
    df$knum<-sample(1:folds,nrow(df),replace = TRUE)
    rmse_kfold<-0
    for (i in 1:folds)
    {
      df.test<-df[df$knum==i,]
      df.train<-df[!df$knum==i,]
      pred<-predict(model,df.test)
      pred[is.na(pred)]<-mean(pred,na.rm = T)
      rmse_kfold<-cbind(rmse_kfold,rmse(df.test[,resp],pred))
    }
    l[[1]]<-rmse_kfold[,-1]
    l[[2]]<-mean(rmse_kfold[,-1])
    return (l)
  }
  else if (cvtype=="LOOCV"||cvtype=="loocv")
  {
    rmse_loocv<-0
    for (i in 1:nrow(df))
    {
      df.test<-df[i,]
      df.train<-df[-i,]
      pred<-predict(model,df.test)
      pred[is.na(pred)]<-mean(df.train[,resp])
      rmse_loocv<-cbind(rmse_loocv,rmse(df.test[,resp],pred))
    }
    l[[1]]<-rmse_loocv[,-1]
    l[[2]]<-mean(rmse_loocv[,-1])
    return(l)
  }
}
loadBart<-function()
{
  .libPaths(c(.libPaths(), "/home/sonin/Rlibs"))
  options(java.parameters="-Xmx100g")
  library('rJava')
  library('bartMachine')
  set_bart_machine_num_cores(20)
}

# dfcreate----
master<-read.csv("inp/recs2009_public.csv")
master.az<-master[which(master$REPORTABLE_DOMAIN==24),]
ucol<-read.csv('inp/usecol.csv')
ucol.x<-ucol$Variable.Name
ucol.x<-as.character(ucol.x)
ucol.y<-ucol[which(ucol$Keep == 'R'),'Variable.Name']
ucol.y<-as.character(ucol.y)
ucol.of<-ucol[which(ucol$FT=='of'),"Variable.Name"]
ucol.of<-as.character(ucol.of)
ucol.f<-ucol[which(ucol$FT=='f'),"Variable.Name"]
ucol.f<-as.character(ucol.f)
ucol.nf<-ucol[which(ucol$FT=='nf'),"Variable.Name"]
ucol.nf<-as.character(ucol.nf)
master.az<-master.az[,ucol.x]
master.az$KWHSPC<-rowSums(master.az[,ucol.y])
rownames(master.az)<-NULL
master.az<-master.az[,!names(master.az) %in% ucol.y]
master.az[,ucol.of]<-lapply(master.az[ucol.of],factor,ordered=TRUE)
master.az[,ucol.f]<-lapply(master.az[ucol.f],factor)
master<-master.az # Subset with only data for AZ and useful columns
set.seed(9)
rows<-sample(1:nrow(master),0.80*nrow(master),replace = F)
master.train<-master[rows,]
master.test<-master[-rows,]
cat("\014")
rm(master.az,ucol,ucol.x,ucol.y,rows)
rmse<-data.frame(matrix(ncol=3,nrow=0))
colnames(rmse)<-c('Model','RMSE.IS','RMSE.OS')
# EDA----
lT2<-rapply(master.train,function(x)length(unique(x)))>2 # All factors with less than 2 unique values were dropped for the purposes of EDA
lT2<-names(lT2[which(lT2==FALSE)])
dnt<-append(lT2,resp)
form.lin<-as.formula(paste0(resp,"~", paste0(colnames(master.train[,!names(master.train) %in% dnt]),collapse="+")))
form.nonlin<-as.formula(paste0(resp,"~", paste0("s(", colnames(master.train[,!names(master.train) %in% c(ucol.f,dnt)]),", d=4)" ,collapse="+"),
                               '+', paste0(colnames(master.train[,!names(master.train) %in% c(ucol.nf,ucol.of,dnt)]),collapse="+")))

eda.gam.lin<-gam(formula=form.lin,data=master.train)
eda.gam.nonlin<-gam(formula=form.nonlin,data=master.train)
rm(form.lin,form.nonlin,dnt,lT2,ucol.f,ucol.nf,ucol.of,ucol.x,ucol.y)
keep<-read.csv('inp/keep.csv')
keep<-keep$out
keep<-as.character(keep)
keep<-append(keep,resp)
master.train<-master.train[,keep] #subsetting only the significant columns
master.test<-master.test[,keep]
rm(keep)

png(filename="plots/scatterplot.png",width=100,height=100,units="in",res=300)
pairs(KWHSPC~.,data=master.train)
dev.off()
cat("\014")
df.train<-master.train
df.test<-master.test
# BART1----
bart1<-bartMachine(X=df.train[,!names(df.train) %in% resp],y=df.train[,resp],
                   serialize = TRUE,run_in_sample = TRUE)
rmse <- rbind(rmse,data.frame('Model'='BART1','RMSE.IS'=rmse(bart1$y,bart1$y_hat_train),
                         'RMSE.OS'=rmse(predict(bart1,df.test[,!names(df.test) %in% resp]),df.test[,resp])))
png(filename = "plots/1.bart1_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=bart1$y,y=bart1$y_hat_train,xlab="Actual",ylab="Predicted",
     main="Y vs. Y-hat")
dev.off()
png(filename = "plots/2.bart1_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(bart1$residuals,main="Normality Plot for Residuals")
qqline(bart1$residuals)
dev.off()
png(filename = "plots/3.bart1_yyhat_credible.png",width=10,height=10,units = 'in',
    res=300)
plot_y_vs_yhat(bart1, credible_intervals = TRUE)
dev.off()
png(filename = "plots/4.bart1_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=bart1$residuals,x=bart1$y,xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()
png(filename = "plots/5.bart1_var_imp.png",width=10,height=10,units = 'in',
    res=300)
investigate_var_importance(bart1, num_replicates_for_avg = 20)
dev.off()
# BART1.CV----
bart1.cv<-bartMachineCV(X=df.train[,!names(df.train) %in% resp],y=df.train[,resp],
                        serialize = TRUE,k_folds = 10)
rmse <- rbind(rmse,data.frame('Model'='BART1.CV','RMSE.IS'=rmse(bart1.cv$y,bart1.cv$y_hat_train),
                              'RMSE.OS'=rmse(predict(bart1.cv,df.test[,!names(df.test) %in% resp]),df.test[,resp])))
colnames(rmse)<-c('Model','RMSE.IS','RMSE.OS')
png(filename = "plots/6.bart1.cv_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=bart1.cv$y,y=bart1.cv$y_hat_train,xlab="Actual",ylab="Predicted",
     main="Y vs. Y-hat")
dev.off()
png(filename = "plots/7.bart1.cv_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(bart1.cv$residuals,main="Normality Plot for Residuals")
qqline(bart1.cv$residuals)
dev.off()
png(filename = "plots/8.bart1.cv_yyhat_credible.png",width=10,height=10,units = 'in',
    res=300)
plot_y_vs_yhat(bart1.cv, credible_intervals = TRUE)
dev.off()
png(filename = "plots/9.bart1.cv_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=bart1.cv$residuals,x=bart1.cv$y,xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()
png(filename = "plots/10.bart1.cv_var_imp.png",width=10,height=10,units = 'in',
    res=300)
investigate_var_importance(bart1.cv, num_replicates_for_avg = 20)
dev.off()
# RF1----
rf1<-randomForest(KWHSPC~. ,data=df.train)
rmse <- rbind(rmse,data.frame('Model'="RF1",'RMSE.IS'=rmse(rf1$y,rf1$predicted),
                              'RMSE.OS'=rmse(predict(rf1,df.test),df.test$KWHSPC),stringsAsFactors = F))
png(filename = "plots/11.rf1_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=rf1$y,y=rf1$predicted,xlab="Actual",ylab="Predicted",
     main="Y vs. Y-hat")
dev.off()
png(filename = "plots/12.rf11_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(rf1$predicted-rf1$y,main="Normality Plot for Residuals")
qqline(rf1$predicted-rf1$y)
dev.off()
png(filename = "plots/13.rf11_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(rf1$predicted-rf1$y,main="Normality Plot for Residuals")
qqline(rf1$predicted-rf1$y)
dev.off()
png(filename = "plots/14.rf11_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=rf1$predicted-rf1$y,x=rf1$predicted,xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()
png(filename = "plots/15.rf11_resid_v_predictors_matrix.png",width=100,height=100,units = 'in',
    res=300)
pairs(rf1$predicted-rf1$y~.,data=df.train)
dev.off()
# SVM----
svm1<-svm(KWHSPC~.,data=df.train)
rmse <- rbind(rmse,data.frame('Model'="SVM1",'RMSE.IS'=rmse(df.train$KWHSPC,svm1$fitted),
                              'RMSE.OS'=rmse(predict(svm1,df.test),df.test$KWHSPC),stringsAsFactors = F))
png(filename = "plots/16.svm1_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=df.train$KWHSPC,y=svm1$fitted,xlab="Actual",ylab="Predicted",
     main="Y vs. Y-hat")
dev.off()
png(filename = "plots/17.svm1_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(svm1$residuals,main="Normality Plot for Residuals")
qqline(svm1$residuals)
dev.off()
png(filename = "plots/18.svm1_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=svm1$residuals,x=svm1$fitted,xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()
# RPart----
rpart1<-rpart(KWHSPC~.,data=df.train)
rmse <- rbind(rmse,data.frame('Model'="RPART1",'RMSE.IS'=rmse(rpart1$y,predict(rpart1,df.train)),
                              'RMSE.OS'=rmse(predict(rpart1,df.test),df.test$KWHSPC),stringsAsFactors = F))
png(filename = "plots/19.rpart1_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=df.train$KWHSPC,y=predict(rpart1,df.train),xlab="Actual",ylab="Predicted",
     main="Y vs. Y-hat")
dev.off()
png(filename = "plots/20.rpart1_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(predict(rpart1,df.train)-rpart1$y,main="Normality Plot for Residuals")
qqline(predict(rpart1,df.train)-rpart1$y)
dev.off()
png(filename = "plots/21.rpart1_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=predict(rpart1,df.train)-rpart1$y,x=predict(rpart1,df.train),xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()

save(list=ls(all=T),file='final.RData')
