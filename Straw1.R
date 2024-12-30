library(dplyr)
library(tidyr)
library(GEOquery)
library(tidyverse)
setwd("C:/Users/ASUS/OneDrive/Desktop/Project/straw")

# read CSV file
dat<- read.csv("GSE283248_gene_fpkm.csv", sep = ";")
dim(dat)

# get metadata
gse<- getGEO(GEO="GSE283248", GSEMatrix = TRUE)
dim(gse)
head(gse)

#get pheno data

metadata<- pData(phenoData(gse[[1]]))
head(metadata)

metadata.subset<-metadata%>%
  subset(., select=c(1,12))%>%
  rename("treatment"="characteristics_ch1.2")%>%
  mutate(treatment= gsub("treatment:","",treatment))%>%
  mutate(title= gsub("ATCC 8739,","", title))


head(metadata.subset)  

# add additional column for creating our same column in both data set

metadata.subset$sample<- paste(metadata.subset$title)
head(metadata.subset)



metadata_subset<- metadata.subset%>%
  mutate(sample=gsub("control, rep1", "C1A", sample))%>%
  mutate(sample=gsub("control, rep2", "C1B", sample))%>%
  mutate(sample=gsub("control, rep3", "C1C", sample))%>%
  mutate(sample=gsub("thermal, 60°C, rep1", "T2A", sample))%>%
  mutate(sample=gsub("thermal 60°C, rep2", "T2B", sample))%>%
  mutate(sample=gsub("thermal, 60°C, rep3", "T2C", sample))%>%
  mutate(sample=gsub("thermal, 72°C, rep1", "T3A", sample))%>%
  mutate(sample=gsub("thermal 72°C, rep2", "T3B", sample))%>%
  mutate(sample=gsub("thermal, 72°C, rep3", "T3C", sample))%>%
  mutate(sample=gsub("high pressure, 300 MPa, rep1", "H4A", sample))%>%
  mutate(sample=gsub("high pressure, 300 MPa, rep2", "H4B", sample))%>%
  mutate(sample=gsub("high pressure, 300 MPa, rep3", "H4C", sample))%>%
  mutate(sample=gsub("high pressure, 400 MPa, rep1", "H5A", sample))%>%
  mutate(sample=gsub("high pressure, 400 MPa, rep2", "H5B", sample))%>%
  mutate(sample=gsub("high pressure, 400 MPa, rep3", "H5C", sample))%>%
  mutate(sample=gsub("pulsed electric field, 50 kJ/kg, rep1", "P6A", sample))%>%
  mutate(sample=gsub("pulsed electric field, 50 kJ/kg, rep2", "P6B", sample))%>%
  mutate(sample=gsub("pulsed electric field, 50 kJ/kg, rep3", "P6C", sample))%>%
  mutate(sample=gsub("pulsed electric field, 70 kJ/kg, rep1", "P7A", sample))%>%
  mutate(sample=gsub("pulsed electric field, 70 kJ/kg, rep2", "P7B", sample))%>%
  mutate(sample=gsub("pulsed electric field, 70 kJ/kg, rep3", "P7C", sample))
  
colnames(dat)
dim(metadata_subset)  


# extract essentials from dat
dat_subset<- subset(dat, select= c(23,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22))

# change it to long format
dat.long<- dat_subset%>%
  gather(key= "ID", value = "FPKM", -gene_name)

tail(dat.long)
typeof(dat.long$ID)
unique(dat.long$ID)

typeof(metadata_subset$sample)
unique(metadata_subset$sample)
metadata_subset<- metadata_subset%>%
  mutate(sample = str_trim(sample, side = "left"))
unique(metadata_subset$sample)

dat.long<- dat.long%>%
  left_join(., metadata_subset, by = c("ID" = "sample"))
tail(dat.long)
unique(dat.long$treatment)
unique(dat.long$title)
typeof(dat.long$FPKM)
 
#save dat.long
#write.csv(dat.long, "dat_long.csv", row.names = FALSE)

# save this bio file
bio<- dat.long%>%
  subset(., select= c(2,5,3))%>%
  rename("rep"= "ID")%>%
  rename("trt"="treatment")
write.csv(bio, "bio.csv", row.names = FALSE)

#explore the data
# treatment effect in the FPKM
trt_FPKM<- dat.long%>%
  group_by(treatment)%>%
  summarize(mean_FPKM= mean(FPKM), medain_FPKM= median(FPKM))
#plot this in box plot and barplot
ggplot(trt_FPKM, aes(x= treatment, y= mean_FPKM))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_col()+
  labs(title = "Bar Plot of Treatment vs mean_FPKM", x="Treatments")
  
library(scales)
# box plot 
dat.long%>%
  filter(treatment==" pulsed electric field, 50 kJ/kg, 6 kV/cm" )%>%
  mutate(FPKM = as.integer(FPKM))%>%
  ggplot(., aes(x=treatment  , y= FPKM))+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_boxplot()+
  scale_y_continuous()+
  labs(title = "Boxplot of Treatment vs FPKM", x="Treatments")


#NATURAL LOG TRANSFORMATION

dat.long$log_FPKM <- log(dat.long$FPKM)

#log transformed box plot
dat.long%>%
  filter(treatment==" pulsed electric field, 50 kJ/kg, 6 kV/cm" )%>%
  ggplot(., aes(x=treatment  , y= log_FPKM))+
  theme_bw()+
  geom_boxplot()

#log box plots of all treatment
dat.long%>%
  ggplot(., aes(x=treatment  , y= log_FPKM))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_boxplot()

# explore genes now 
dat.long%>%
  filter(gene_name== "yfdY"|gene_name=="yjfY"|gene_name=="yjbJ")%>%
  ggplot(., aes(x= gene_name , y= log_FPKM))+
  theme_bw()+
  geom_boxplot()+
  ggtitle("box plot of 3 random  genes")

# scatter plot of the genes
scat.long<- dat.long%>%
  subset(., select=c(1,2,5,6))

scat.long %>% 
  filter(gene_name == "yfdY" | gene_name == "yjfY") %>% 
  pivot_wider(names_from = gene_name, values_from = log_FPKM) %>% 
  ggplot(., aes(x=yfdY, y= yjfY, color= treatment))+
  geom_point()+
  theme_bw()+
  #geom_smooth(method = "lm", se= FALSE)+
  labs(title = "yjfY vs yfdY")

# heat map
genes.of.interest<- c("nlpl", "hokB","cspB","ispE", "yjfY","yfdY")
scat.long%>%
  filter(gene_name %in% genes.of.interest)%>%
  ggplot(., aes(x= ID, y= gene_name, fill = log_FPKM))+
  geom_tile()+
  scale_fill_gradient(low="white",high = "red")+
  labs(title= "Heat map of genes according to samples")

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
