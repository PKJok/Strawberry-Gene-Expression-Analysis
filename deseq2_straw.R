library(dplyr)
library(tidyr)
library(DESeq2)
library(tidyverse)
setwd("C:/Users/ASUS/OneDrive/Desktop/Project/straw")

# read CSV file
dat<- read.csv("GSE283248_gene_fpkm.csv", sep = ";")
dim(dat)
head(dat)
