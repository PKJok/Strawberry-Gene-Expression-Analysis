library(gvlma)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
setwd("C:/Users/ASUS/OneDrive/Desktop")
bio1<- read.csv("gvl.csv", header = TRUE)

head(bio1)
tail(bio1)

# attach
attach(bio1)
detach(pos = 3) 
detach(pos = 4)
bio$rep<- as.factor(rep)
bio$trt<-as.factor(trt)

# test assumption of anova 
require(gvlma)
fit=lm(FPKM~rep+trt)
gvmodelfit= gvlma(fit)
gvmodelfit
boxplot(FPKM)
