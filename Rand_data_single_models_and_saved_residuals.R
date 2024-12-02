### Import the libraries
library(fabletools)
library(tidyverse)
library(datasets)
library(zoo)
library(tsibble)
library(tseries)
library(lubridate)
library(fpp3)
library(tsibbledata)
library(feasts)
library(GGally)
library(grid)
library(gridExtra)
library(forecast)
library(Metrics)
library(caTools)
library(TSstudio)
library(dplyr)
library(Rssa)
library("plyr")
rm(list=ls())

###############################Rands data

## read and visualize the data

df1<-read.csv('DEXSFUS.csv')
head(df1)
tail(df1)
df1$DEXSFUS<-as.numeric(df1$DEXSFUS)
df1$DATE<-as.Date(df1$DATE)
###ggplot
df1%>%ggplot(aes(x = DATE, y = DEXSFUS))+
  geom_path(colour = "blue")+
  ggtitle("Spot Exchange Rate")+
  theme(plot.title = element_text(size = 9))+ labs(x = "Date", y = "ZAR/USD")

## Split the data into train and test 
ts1<-ts(df1$DEXSFUS,frequency = 5, start=c(2015,1,1))

split1 <- ts_split(ts1, sample.out = 20)

train1 <- split1$train
test1<-split1$test


###Arima model

fit_ar1<-auto.arima(train1)
summary(fit_ar1)

# Forecast the next 20 days 
for_ar1<-forecast(fit_ar1,h=20)
autoplot(for_ar1)

plot(train1)
lines(fitted(fit_ar1), col="red")


##Exponential Smoothing 

fit_es1<-ets(train1)
summary(fit_es1)

# Forecast the next 20 days
for_es1<-forecast(fit_es1,h=20)
autoplot(for_es1)

plot(train1)
lines(fitted(fit_es1), col="red")


### SSA-LRF model

#Perform SSA on the train data
S1<-ssa(train1,
        kind = c("1d-ssa"),
        circular = FALSE,
        svd.method = c( "svd"),
        force.decompose = TRUE)

summary(S1)

# Correlation between components
w1 <- wcor(S1, groups = 1:30)
plot(w1, grid = c(14,15))

# Signals and noise (51 first components as signal)
r1 <- reconstruct(S1,groups = list(c(1:51)))
plot(r1,plot.method = "xyplot")

plot(train1,type = 'l')
lines(r1$F1,col='red')

# Forecast the next 20 days
rfor1<- rforecast(S1, groups = list(c(1:51)), len = 20, only.new=FALSE)

matplot(data.frame(c(train1, rep(NA, 20)), rfor1), type="l")

forecast::accuracy(test1,rfor1[2068:2087])
### Now we are using the ARIMA and ES models on the time series signal, which will be further used to build SSA-ARIMA-LSTM and SSA-ES-LSTM.

fit_tar1<-auto.arima(r1$F1,seasonal = FALSE)
fit_tes1<-ets(r1$F1)

summary(fit_tar1)
summary(fit_tes1)

#Forecast the next 20 days
for_tar1<-forecast(fit_tar1,h=20)
for_tes1<-forecast(fit_tes1,h=20)

########### Now we will uses all of the previously build models to perform a short term daily forcast
### The forecast is therefore performed based on the test data.

##Arima
d_ar1<-forecast(test1,model = fit_ar1)
summary(d_ar1)

##ES
d_es1<-forecast(test1,model = fit_es1)
summary(d_es1)

##SSA-LRF
# For this case we forecast on day ahead with the different SSA models
rfor1_<-0
for (i in seq_along(test1)) {
  b <- ts1[1:(i+2066)]
  S1<-ssa(b,
          kind = c("1d-ssa"),
          circular = FALSE,
          svd.method = c( "svd"),
          force.decompose = TRUE)
  rfor1a<- rforecast(S1, groups = list(c(1:51)), len = 1, only.new=FALSE)
  rfor1_<-c(rfor1_,tail(rfor1a, n=1))
}

rfor1_<-rfor1_[-1] # We collect the 20 forecast

### We use SSA on the entire time series this will be useful to build the SSA-LSTM-LSTM models  

Sa<-ssa(ts1,
        kind = c("1d-ssa"),
        circular = FALSE,
        svd.method = c( "svd"),
        force.decompose = TRUE)

summary(Sa)

## Reconstruct the signal
ra <- reconstruct(Sa,groups = list(c(1:51)))
plot(ra,plot.method = "xyplot")
plot(ts1,type = 'l')
lines(ra$F1,col='red')

##### Now we are using the ARIMA and ES models on the  signal
#We fit other models because we will need the residuals of the signal forecast for daily forecast 

dfit_tar1<-auto.arima(ra$F1[1:2067],seasonal = FALSE)
dfit_tes1<-ets(ra$F1[1:2067])

summary(dfit_tar1)
summary(dfit_tes1)

# Forecasting 
d_tar1<-forecast(ra$F1[2068:2087],model=dfit_tar1)
d_tes1<-forecast(ra$F1[2068:2087],model=dfit_tes1)


###We save each model and its residuals for further analysis in python

## Creation of the dataframe for the 20 days forcasting models 

#fit models 
fit_arim1<-c(fit_ar1$fitted,for_ar1$mean)
fit_esp1<-c(fit_es1$fitted,for_es1$mean)
fit_lrf1<-c(rfor1)

fit_tarim1<-c(fit_tar1$fitted,for_tar1$mean)
fit_tesp1<-c(fit_tes1$fitted,for_tes1$mean)

#Residuals
res_arim1<-c(fit_ar1$residuals,test1-for_ar1$mean)
res_esp1<-c(fit_es1$residuals,test1-for_es1$mean)
res_lrf1<-c(ts1-rfor1)

res_tarim1<-c(fit_tar1$residuals,test1-for_tar1$mean)
res_tesp1<-c(fit_tes1$residuals,test1-for_tes1$mean)

####
Rand <-data.frame(fit_arim1,fit_esp1,fit_lrf1,res_arim1,res_esp1,res_lrf1,
                  fit_tarim1,fit_tesp1,res_tarim1,res_tesp1)

merge_1<- cbind(df1, Rand)
head(merge_1)

write.csv(merge_1, "data_final_Rand.csv", row.names = FALSE)

# Creation of the DataFrame Considering the Short-Term Daily Forecasting Approach of the 20-Day Forecasting Models

#fit models 
dfit_arim1<-c(fit_ar1$fitted,d_ar1$fitted)
dfit_esp1<-c(fit_es1$fitted,d_es1$fitted)
dfit_lrf1<-c(r1$F1,rfor1_)

dfit_tarim1<-c(dfit_tar1$fitted,d_tar1$fitted)
dfit_tesp1<-c(dfit_tes1$fitted,d_tes1$fitted)


#Residuals
dres_arim1<-c(fit_ar1$residuals,test1-as.vector(d_ar1$fitted))
dres_esp1<-c(fit_es1$residuals,test1-as.vector(d_es1$fitted))
dres_lrf1<-c(residuals(r1), test1-rfor1_)

dres_tarim1<-c(dfit_tar1$residuals,ra$F1[2068:2087]-d_tar1$fitted)
dres_tesp1<-c(dfit_tes1$residuals,ra$F1[2068:2087]-d_tes1$fitted)

signal1<-ra$F1
noice1<-residuals(ra)


##
Rand_ <-data.frame(dfit_arim1,dfit_esp1,dfit_lrf1,dres_arim1,dres_esp1,dres_lrf1,
                  dfit_tarim1,dfit_tesp1,dres_tarim1,dres_tesp1,signal1,noice1)

merge_1<- cbind(df1, Rand_)
head(merge_1)

write.csv(merge_1, "data_final_day_Rand.csv", row.names = FALSE)

