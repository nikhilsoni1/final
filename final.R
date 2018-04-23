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
library(earth)
cat("\014")
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
crossValidate<-function(cvtype='kfold',folds=10,dataset,model=NULL,resp,typ=1)
{
  if(typ==1)
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
  else if(typ==2)
  {
    df<-dataset
    df$knum<-sample(1:folds,nrow(df),replace = T)
    dnt<-c('knum',resp)
    l <- vector("list", 2)
    rmse_kfold<-0
    for(i in 1:folds)
    {
      df.train<-df[!df$knum==i,]
      df.test<-df[df$knum==i,]
      model<-bartMachine(df.train[,!names(df.train) %in% dnt],df.train[,resp],verbose=F)
      pred<-predict(model,df.test[,!names(df.test) %in% dnt])
      rmse_kfold<-cbind(rmse_kfold,rmse(df.test[,resp],pred))
      print(i)
    }
    l[[1]]<-rmse_kfold[,-1]
    l[[2]]<-mean(rmse_kfold[,-1])
    return (l)
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
rmse.cv<-data.frame(matrix(ncol=3,nrow=0))
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
png(filename='plots/0.1corrplot.png',width=20,height=20,units='in',res=300)
corrplot.mixed(cor(master.train[,-c(1,3,13,18)],use="complete.obs"),lower="number",
         upper="circle",diag='l',tl.pos="lt")
dev.off()
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

# MARS1----
mars1<-earth(KWHSPC~.,df.train)
rmse <- rbind(rmse,data.frame('Model'="MARS1",'RMSE.IS'=rmse(df.train$KWHSPC,mars1$fitted.values),
                              'RMSE.OS'=rmse(predict(mars1,df.test),df.test$KWHSPC),stringsAsFactors = F))
png(filename = "plots/22.mars1_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=df.train$KWHSPC,y=mars1$fitted.values,xlab="Actual",ylab="Predicted",
     main="Y vs. Y-hat")
dev.off()
png(filename = "plots/23.mars1_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(mars1$residuals,main="Normality Plot for Residuals")
qqline(mars1$residuals)
dev.off()
png(filename = "plots/24.mars1_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=mars1$residuals,x=mars1$fitted.values,xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()

# GAM1----
scope<-gam.scope(df.train,response = which(colnames(df.train)==resp),smoother = 's',
                 arg=c("df=4","df=6","d=8"),form=T)
form<-as.formula(paste0(resp,"~", paste0(colnames(df.train[,!names(df.train) %in% resp]),collapse="+")))
obj<-gam(formula=form,family=gaussian,data=df.train)
gam1<-step.Gam(object=obj,scope = scope,direction = 'both',trace = F)
rm(scope,form,obj)
rmse <- rbind(rmse,data.frame('Model'="GAM1",'RMSE.IS'=rmse(gam1$y,gam1$fitted.values),
                              'RMSE.OS'=rmse(predict(gam1,df.test),df.test$KWHSPC),stringsAsFactors = F))
png(filename = "plots/25.gam1_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=gam1$y,y=gam1$fitted.values,xlab="Actual",ylab="Predicted",
     main="Y vs. Y-hat")
dev.off()
png(filename = "plots/26.gam1_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(gam1$residuals,main="Normality Plot for Residuals")
qqline(gam1$residuals)
dev.off()
png(filename = "plots/27.gam1_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=gam1$residuals,x=gam1$fitted.values,xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()

temp<-df.train
temp$knum<-sample(1:10,nrow(temp),replace = T)
model_RMSE<-0
for(i in 1:10)
{
  temp.train<-temp[!temp$knum==i,]
  temp.test<-temp[temp$knum==i,]
  temp.model<-earth(formula=gam1$formula,data=temp.train)
  temp.pred<-predict(temp.model,temp.test[,-c(1,2,3,4)])
  model_RMSE<-cbind(model_RMSE,rmse(temp.test[,resp],temp.pred))
  print(i)
}
gam1.cv<-vector("list", 2)
gam1.cv[[1]]<-model_RMSE[,-1]
gam1.cv[[2]]<-mean(model_RMSE[,-1])
rm(temp.train,temp.test,temp,temp.pred,temp.model,model_RMSE,i)
rmse.cv<-rbind(rmse.cv,data.frame('Model'='GAM1','RMSE.IS'=gam1.cv[[2]],
                                  'RMSE.OS'=rmse(predict(gam1,df.test),df.test[,resp])))

# BART2----
png(filename = "plots/28.bart2_var_selection_by_permute.png",width=10,height=10,units = 'in',res=300)
bart1.cv.important<-var_selection_by_permute(bart1.cv, 
                                num_reps_for_avg = 10, num_permute_samples = 100, 
                                num_trees_for_permute = 20, alpha = 0.05, 
                                plot = TRUE, num_var_plot = Inf, bottom_margin = 10)
dev.off()
imp<-bart1.cv.important$important_vars_local_col_nums
bart2<-bartMachine(X=df.train[,imp],df.train[,resp],serialize = T)
rmse <- rbind(rmse,data.frame('Model'='BART2','RMSE.IS'=rmse(bart2$y,bart2$y_hat_train),
                              'RMSE.OS'=rmse(predict(bart2,df.test[,imp]),df.test[,resp])))
png(filename = "plots/29.bart2_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=bart2$y,y=bart2$y_hat_train,xlab="Actual",ylab="Predicted",
     main="Y vs. Y-hat")
dev.off()
png(filename = "plots/30.bart2_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(bart2$residuals,main="Normality Plot for Residuals")
qqline(bart2$residuals)
dev.off()
png(filename = "plots/31.bart2_yyhat_credible.png",width=10,height=10,units = 'in',
    res=300)
plot_y_vs_yhat(bart2, credible_intervals = TRUE)
dev.off()
png(filename = "plots/32.bart2_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=bart2$residuals,x=bart2$y,xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()

imp<-append(imp,which(names(df.train)==resp))
bart2.cv<-crossValidate(dataset = df.train[,imp],resp=resp,typ = 2)
rmse.cv<-rbind(rmse.cv,data.frame('Model'='bart2','RMSE.IS'=bart2.cv[[2]],
                                  'RMSE.OS'=rmse(predict(bart2,df.test[,imp[!imp %in% which(names(df.train)==resp)]]),df.test[,resp])))
rm(imp)

# MARS2----
mars1.importance<-evimp(mars1,trim = F)
png(filename = "plots/33.mars1_evimp.png",width=10,height=10,units = 'in',
    res=300)
plot(mars1.importance)
dev.off()
imp<-names(mars1.importance[mars1.importance[,'nsubsets']>0,'nsubsets'])
form<-as.formula(paste0(resp,'~',paste0(imp,collapse = '+')))
mars2<-earth(formula=form,data=df.train)
rmse <- rbind(rmse,data.frame('Model'="MARS2",'RMSE.IS'=rmse(df.train$KWHSPC,mars2$fitted.values),
                              'RMSE.OS'=rmse(predict(mars2,df.test),df.test$KWHSPC),stringsAsFactors = F))
png(filename = "plots/34.mars2_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=df.train$KWHSPC,y=mars2$fitted.values,xlab="Actual",ylab="Predicted",
     main="Y vs. Y-hat")
dev.off()
png(filename = "plots/35.mars2_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(mars2$residuals,main="Normality Plot for Residuals")
qqline(mars2$residuals)
dev.off()
png(filename = "plots/36.mars2_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=mars2$residuals,x=mars2$fitted.values,xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()
png(filename = "plots/37.mars2_resid_v_predictors.png",width=10,height=10,units = 'in',
    res=300)
pairs(mars2$residuals~.,data=df.train[,imp])
dev.off()
mars2.cv<-crossValidate(dataset=df.train,model=mars2,resp=resp)

temp<-df.train
temp$knum<-sample(1:10,nrow(temp),replace = T)
model_RMSE<-0
for(i in 1:10)
{
  temp.train<-temp[!temp$knum==i,]
  temp.test<-temp[temp$knum==i,]
  temp.model<-earth(formula=form,data=temp.train)
  temp.pred<-predict(temp.model,temp.test[,-c(1,2,3,4)])
  model_RMSE<-cbind(model_RMSE,rmse(temp.test[,resp],temp.pred))
  print(i)
}
mars2.cv<-vector("list", 2)
mars2.cv[[1]]<-model_RMSE[,-1]
mars2.cv[[2]]<-mean(model_RMSE[,-1])
rm(temp.train,temp.test,temp,temp.pred,temp.model,model_RMSE,i)
rmse.cv<-rbind(rmse.cv,data.frame('Model'='mars2','RMSE.IS'=mars2.cv[[2]],
                                  'RMSE.OS'=rmse(predict(mars2,df.test),df.test[,resp])))







# Final Model----
imp<-bart1.cv.important$important_vars_local_col_nums
bart2<-bartMachineCV(X=df.train[,imp],y=df.train[,resp],serialize = TRUE,k_folds = 10)
rmse <- rbind(rmse,data.frame('Model'='BART2','RMSE.IS'=rmse(bart2$y,bart2$y_hat_train),
                              'RMSE.OS'=rmse(predict(bart2,df.test[,imp]),df.test[,resp])))
png(filename = "plots/38.bart2_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=bart2$y,y=bart2$y_hat_train,xlab="Actual",ylab="Predicted",main="Y vs. Y-hat")
dev.off()
png(filename = "plots/39.bart2_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(bart2$residuals,main="Normality Plot for Residuals")
qqline(bart2$residuals)
dev.off()
png(filename = "plots/40.bart2_yyhat_credible.png",width=10,height=10,units = 'in',
    res=300)
plot_y_vs_yhat(bart2, credible_intervals = TRUE)
dev.off()
png(filename = "plots/41.bart2_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=bart2$residuals,x=bart2$y,xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()
png(filename = "plots/42.bart2_var_imp.png",width=10,height=10,units = 'in',
    res=300)
investigate_var_importance(bart2, num_replicates_for_avg = 20)
dev.off()
png(filename = "plots/43.bart2_dolecol.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart2, 1, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/44.bart2_kwh.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart2, 2, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/45.bart2_dolelwth.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart2, 3, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/46.bart2_dolelsph.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart2, 4, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()

# bart3----
bart3.include<-c('DOLELCOL','KWH','DOLELWTH','DOLELSPH','DOLLAREL','TOTUCSQFT',
                 'TOTCSQFT','TOTSQFT','AIA_ZONE','HDD30YR')
bart3<-bartMachineCV(X=df.train[,names(df.train) %in% bart3.include],y=df.train[,resp],serialize=T)
rmse <- rbind(rmse,data.frame('Model'='BART3','RMSE.IS'=rmse(bart3$y,bart3$y_hat_train),
                              'RMSE.OS'=rmse(predict(bart3,df.test[,names(df.train) %in% bart3.include]),df.test[,resp])))
png(filename = "plots/47.bart3_y_yhat.png",width=10,height=10,units = 'in',
    res=300)
plot(x=bart3$y,y=bart3$y_hat_train,xlab="Actual",ylab="Predicted",main="Y vs. Y-hat")
dev.off()
png(filename = "plots/48.bart3_resid_normal.png",width=10,height=10,units = 'in',
    res=300)
qqnorm(bart3$residuals,main="Normality Plot for Residuals")
qqline(bart3$residuals)
dev.off()
png(filename = "plots/49.bart3_yyhat_credible.png",width=10,height=10,units = 'in',
    res=300)
plot_y_vs_yhat(bart3, credible_intervals = TRUE)
dev.off()
png(filename = "plots/50.bart3_resid_v_pred.png",width=10,height=10,units = 'in',
    res=300)
plot(y=bart3$residuals,x=bart3$y,xlab='Predicted',ylab='Residuals',main='Residuals vs. Predicted')
abline(h=0)
dev.off()
png(filename = "plots/51.bart3_var_imp.png",width=10,height=10,units = 'in',
    res=300)
investigate_var_importance(bart3, num_replicates_for_avg = 20)
dev.off()
png(filename = "plots/52.bart3_TOTUCSQFT.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart3, 1, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/53.bart3_DOLELSPH.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart3, 2, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/54.bart3_TOTCSQFT.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart3, 3, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/55.bart3_DOLLAREL.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart3, 4, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/56.bart3_DOLELWTH.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart3, 5, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/57.bart3_HDD30YR.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart3, 6, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/58.bart3_KWH.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart3, 7, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/59.bart3_DOLELCOL.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart3, 8, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename ="plots/60.bart3_TOTSQFT.png",width=10,height=10,units = 'in',res=300)
pd_plot(bart3, 9, 
        levs = c(0.05, seq(from = 0.1, to = 0.9, by = 0.1), 0.95), 
        lower_ci = 0.025, upper_ci = 0.975, prop_data = 1)
dev.off()
png(filename = "plots/61.bart3_var_selection_by_permute.png",width=10,height=10,units = 'in',res=300)
bart3.cv.important<-var_selection_by_permute(bart3, 
                                             num_reps_for_avg = 10, num_permute_samples = 100, 
                                             num_trees_for_permute = 20, alpha = 0.05, 
                                             plot = TRUE, num_var_plot = Inf, bottom_margin = 10)
dev.off()

