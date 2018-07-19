## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.show='hold'----------------------------------------------------


## ---- fig.show='hold'----------------------------------------------------
# Clopidogrel vs Aspirin trial
cases<-c(939,1021)
person_yrs<-c(17636,17519)
patients<-c(9599,9586)   

## ---- fig.show='hold'----------------------------------------------------
library(metalr)
ici.rr(cases = cases,patients = patients,person_yrs = person_yrs)

## ---- fig.show='hold'----------------------------------------------------
library(metalr)
data("sample_metarr_data")
head(sample_metarr_data)

## ---- fig.show='hold'----------------------------------------------------
metalr_rr(sample_metarr_data,refval=0,num_iter=3000,increm=0.001,method="random")

## ---- fig.show='hold'----------------------------------------------------
library(metalr)
data("statindata")
statindata

## ---- fig.show='hold'----------------------------------------------------
ici.or(idata = statindata[6,2:5])

## ---- fig.show='hold'----------------------------------------------------
metalr_rr(statindata[,2:5],refval=0,num_iter=3000,increm=0.001,method="fixed")

