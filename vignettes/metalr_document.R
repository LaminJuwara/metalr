## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.show='hold'----------------------------------------------------
if (!require("devtools"))install.packages("devtools")
devtools::install_github("laminjuwara/metalr")

## ---- fig.show='hold'----------------------------------------------------
# Clopidogrel vs Aspirin trial dataset
cases<-c(939,1021)
person_yrs<-c(17636,17519)
patients<-c(9599,9586)   

