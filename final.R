# init----
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load(paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/final.Rdata'))

# include----
library(corrplot)
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

# dfcreate----
master<-read.csv("inp/recs2009_public.csv")
master.az<-master[which(master$REPORTABLE_DOMAIN==24),]
ucol<-read.csv('inp/usecol.csv')
ucol<-(ucol$Variable.Name)
ucol<-as.character(ucol)
master.az<-master.az[,ucol]
master<-master.az # Subset with only data for AZ and useful columns
rm(master.az,ucol)
save(list=ls(all=T),file='final.RData')