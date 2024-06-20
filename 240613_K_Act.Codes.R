###Library packages####
library(tidyverse)
library(magrittr)
library(vegan)
library(viridis)
library(cowplot)
library(ggsignif)
library(ggrepel)
library(scales)
library(ggbiplot)
library(corrplot)
library(psych)
library(compositions) #Centered log ratio (clr) transformation


####Alpha-, beta-diversity####
###packages for Alpha-, beta-diversity
library(vegan)


###Alpha-diversity (Hill numbers) #Number of species: Number of ASV
##Data to be used
#1.reads.count.rarefy.df
head(reads.count.rarefy.df[1:5])

#Number of species (Richness)
richness.div <- reads.count.rarefy.df %>%
  dplyr::mutate(richness = apply(reads.count.rarefy.df, 1, function(x) sum(x > 0))) %>%
  dplyr::select(richness)
write.table(richness.div, "./data/richness.div.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)

#Shannon's index
shannon.div0 <- vegan::diversity(reads.count.rarefy.df,index="shannon")
shannon.div <- as.data.frame(shannon.div0) %>%
  dplyr::rename(shannon = shannon.div0)
write.table(shannon.div, "./data/shannon.div.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)

#Simpson's index
simpson.div0 <- vegan::diversity(reads.count.rarefy.df,index="simpson")
simpson.div <- as.data.frame(simpson.div0) %>%
  dplyr::rename(simpson = simpson.div0)
write.table(simpson.div, "./data/simpson.div.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)

##Calculate Hill numbers (D)
#D = number of species (richness.div) when q = 0.
#D = exp(Shannon's index) when q = 1.
#D = 1/(1-(Simpson's index)) when q = 2.
D0.div <- richness.div %>%
  dplyr::rename(D0 = richness) %>%
  dplyr::mutate(SampleID = rownames(.))

D1.div <- exp(shannon.div) %>%
  dplyr::rename(D1 = shannon) %>%
  dplyr::mutate(SampleID = rownames(.))

D2.div <- as.data.frame(1/(1-(simpson.div))) %>%
  dplyr::rename(D2 = simpson) %>%
  dplyr::mutate(SampleID = rownames(.))

Hill.div <- inner_join(D0.div, D1.div, by = "SampleID")
Hill.div <- inner_join(Hill.div, D2.div, by = "SampleID") %>%
  dplyr::select(SampleID, D0, D1, D2)
write.table(Hill.div, "./data/Hill.div.txt", quote = F, sep = "\t", col.names = T, row.names = F, append = F)


###Two-way ANOVA
#using JMP Pro (MA, USA)

##Summarize the result of Two-way ANOVA (Layer*Fraction)
Adiv.LSmeans <- read.table("./data/Hill.div.LSmean.txt",header=T,sep="\t")
Adiv.LSmeans <- Adiv.LSmeans %>% 
  dplyr::mutate(high.SD = Adiv.LSmeans$Least.Sq.Mean + Adiv.LSmeans$Std.Error,
                low.SD = Adiv.LSmeans$Least.Sq.Mean - Adiv.LSmeans$Std.Error) %>% 
  transform(Layer = factor(Layer, levels = c("ML", "Z1%")), 
            Fraction = factor(Fraction, levels = c("SK", "SS", "FL")))
head(Adiv.LSmeans)

Adiv.LSmeans_D0 <- Adiv.LSmeans %>%
  dplyr::filter(D == "D0")
Adiv.LSmeans_D1 <- Adiv.LSmeans %>%
  dplyr::filter(D == "D1")
Adiv.LSmeans_D2 <- Adiv.LSmeans %>%
  dplyr::filter(D == "D2")



###Beta-diversity (NMDS plot, PERMDISP, PERMANOVA)
##Data to be used
#1.reads.hel.df
head(reads.hel.df[1:5])
m_reads.hel.df <- as.matrix(reads.hel.df)

##NMDS
set.seed(123)
nmdsD <- metaMDS(m_reads.hel.df, distance = "bray", trymax =100, autotransform = F)
nmdsD$stress
#0.1734072
plot(nmdsD)
dev.off()

#Output NMDS
data.scoresD <- as.data.frame(scores(nmdsD, "sites")) %>% 
  dplyr::mutate(sample = rownames(.)) %>%
  separate(sample, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_") %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "Sinking")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "Suspended")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "Free-living")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  transform(Layer= factor(Layer, levels = c("ML", "Z1%"))) %>%
  transform(Fraction= factor(Fraction, levels = c("Sinking", "Suspended", "Free-living")))
head(data.scoresD)



##PERMDISP
reads.hel.dist <- metaMDSdist(m_reads.hel.df, distance = "bray", trymax =100, autotransform = F)
reads.hel.dist.factors <- as.data.frame(as.matrix(reads.hel.dist)) %>%
  dplyr::mutate(sample = rownames(.)) %>%
  dplyr::select(sample) %>%
  separate(sample, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_")
head(reads.hel.dist.factors)
Fraction.disper <- betadisper(reads.hel.dist, reads.hel.dist.factors$Fraction)
Fraction.disper.test <- permutest(Fraction.disper, perm = 999)
Fraction.disper.res <- as.data.frame(Fraction.disper.test$tab)
Fraction.disper.dis <- Fraction.disper$distances
Fraction.disper.dis.df <- as.data.frame(Fraction.disper.dis) %>%
  dplyr::rename(dis.center = Fraction.disper.dis) %>%
  dplyr::mutate(sample = rownames(.)) %>%
  separate(sample, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_")
#Fraction.disper.res:The result of PERMDISP (Factor =  Fraction)
write.table(Fraction.disper.res, "./data/Fraction.disper.res.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)

#
Layer.disper <- betadisper(reads.hel.dist, reads.hel.dist.factors$Layer)
Layer.disper.test <- permutest(Layer.disper, perm = 999)
Layer.disper.res <- as.data.frame(Layer.disper.test$tab)
Layer.disper.dis <- Layer.disper$distances
Layer.disper.dis.df <- as.data.frame(Layer.disper.dis) %>%
  dplyr::rename(dis.center = Layer.disper.dis) %>%
  dplyr::mutate(sample = rownames(.)) %>%
  separate(sample, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_")
#Layer.disper.res:The result of PERMDISP (Factor =  Layer)
write.table(Layer.disper.res, "./data/Layer.disper.res.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)

  
#PERMDISP plot
Dispersal_plot <- ggplot(Layer.disper.dis.df) +
  geom_boxplot(aes(x = Layer, y = dis.center)) +
  labs(x = "Fraction", y = "distance to center") +
  scale_fill_manual(values = c("white", "darkgrey"))+
  NULL
graphics.off()
plot(Dispersal_plot)
dev.off()

Dispersal_plot <- ggplot(Fraction.disper.dis.df) +
  geom_boxplot(aes(x = Fraction, y = dis.center)) +
  labs(x = "Fraction", y = "distance to center") +
  scale_fill_manual(values = c("white", "darkgrey"))+
  NULL
graphics.off()
plot(Dispersal_plot)
dev.off()


##PERMANOVA
m_reads.hel.dist <- as.matrix(reads.hel.dist)
reads.hel.dist.factors

PERMANOVA.dist <- adonis2(m_reads.hel.dist ~ Layer+Fraction, data = reads.hel.dist.factors, permutations = 999, method="bray")
PERMANOVA.dist
PERMANOVA.dist <- as.data.frame(PERMANOVA.dist)
write.table(PERMANOVA.dist, "./data/PERMANOVA.dist.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)

#####


###Differential abundance analysis####
###packages for ANCOM-BC (Differential abundance analysis)
library(MASS)
library(nlme)
library(phyloseq)
library(readr)
source("scripts/ancom_bc.R") #scripts for ANCOM-BC function from "https://github.com/FrederickHuangLin/ANCOM-BC"

###1st ANCOM-BC (1st Differential abundance analysis)####
##Data to be used
#1.reads.count.pre.df
head(reads.count.pre.df[1:5])

##Import feature.table (Pre-processed count data)
First.feature.tab <- reads.count.pre.df
head(First.feature.tab)
dim(First.feature.tab)
#1046 ASVs, 48 samples
#Check the data is a data frame of matrix
class(First.feature.tab)

##Make sample type data; row is sample
First.meta.data <- colnames(First.feature.tab)
First.meta.data <- as.data.frame(First.meta.data) %>%
  dplyr::rename(SampleID = First.meta.data) %>%
  separate(SampleID, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_", remove = F) %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "Sinking")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "Suspended")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "Free-living")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  transform(Fraction= factor(Fraction, levels = c("Free-living", "Suspended", "Sinking"))) %>%
  transform( Layer= factor(Layer, levels = c("ML", "Z1%")))
rownames(First.meta.data) = First.meta.data$SampleID
#Check the data is a data frame
class(First.meta.data)

##Pre-processing
#feature.table = First.feature.tab
#meta.data = First.meta.data
#neg.lb = FALSE #Only ASVs that do not occur at all in groups is detected as ASVs with structural zero
lib.cut = 1000 #Cut samples without 1000 read library
zero.cut = 0.90 #Cut ASVs with >90% 0 occurrence
sample.var = "SampleID" #Character. The name of column storing sample IDs.
group.var = "Fraction" #Character. The name of the group indicator.
pre.process.1st = feature_table_pre_process(First.feature.tab, First.meta.data, sample.var, 
                                                group.var, zero.cut, lib.cut, neg.lb = FALSE)
feature.table.1st = pre.process.1st$feature.table
group.name.1st = pre.process.1st$group.name
group.ind.1st = pre.process.1st$group.ind
struc.zero.1st = pre.process.1st$structure.zeros
dim(feature.table.1st)
#841 ASVs 48 samples

##ANCOM-BC main function
#P-value adjust method = "BH"
#tol.EM, max.iterNum, perNum are default values.
#alpha (significance level) = 0.05
adj.method = "BH"
tol.EM = 1e-5
max.iterNum = 100
perNum = 1000
alpha = 0.05
#Record start time and end time.
start_time <- Sys.time()
out.1st = ANCOM_BC(feature.table.1st, group.name.1st, group.ind.1st, struc.zero.1st,
                       adj.method, tol.EM, max.iterNum, perNum, alpha)
end_time <- Sys.time()
end_time - start_time

##Output the results of 1st ANCOM-BC
res.1st = cbind(ASV = rownames(out.1st$feature.table), out.1st$res)
write_csv(res.1st, "./data/ANCOMBCres.1st.csv")
res.ANCOM_BC.1st = data.frame(ASV=rownames(out.1st$feature.table), out.1st$res, 
                                  struc.zero.1st[rownames(out.1st$feature.table), ], 
                                  row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
alpha.adj.1st=0.05/nrow(res.ANCOM_BC.1st)
critic.val.1st=qnorm(1-alpha.adj.1st/2)


##Summarize the results of 1st ANCOM-BC
#Output all data
LogFC.1st_All <- res.ANCOM_BC.1st %>% dplyr::rename(
  "LogFC_SKFL" = "mean.difference (Sinking - Free-living)", 
  "LogFC_SSFL" = "mean.difference (Suspended - Free-living)", 
  "LogFC_SKSS" = "mean.difference (Sinking - Suspended)",
  "se_SKFL" = "se (Sinking - Free-living)",
  "se_SSFL" = "se (Suspended - Free-living)",
  "se_SKSS" = "se (Sinking - Suspended)",
  "structural.zero.FL" = "structural.zero (Free-living)",
  "structural.zero.SS" = "structural.zero (Suspended)",
  "structural.zero.SK" = "structural.zero (Sinking)") %>%
  dplyr::mutate(ci.lo.adj_SKFL = LogFC_SKFL - critic.val.1st*se_SKFL,
                ci.up.adj_SKFL = LogFC_SKFL + critic.val.1st*se_SKFL,
                ci.lo.adj_SSFL = LogFC_SSFL - critic.val.1st*se_SSFL,
                ci.up.adj_SSFL = LogFC_SSFL + critic.val.1st*se_SSFL,
                ci.lo.adj_SKSS = LogFC_SKSS - critic.val.1st*se_SKSS,
                ci.up.adj_SKSS = LogFC_SKSS + critic.val.1st*se_SKSS)
head(LogFC.1st_All)
nrow(LogFC.1st_All) #841ASVs


##Classify ASVs to Category PA, FL, and FL/PA according to the results of 1st ANCOM-BC
#Select ASVs with Non-differential abundance (Category FL/PA)
PA_FL.ASV.1st.list <- LogFC.1st_All %>% dplyr::filter(diff.abn == "FALSE") %>% 
  dplyr::mutate(Category = "FL/PA") 
nrow(PA_FL.ASV.1st.list) #158ASVs

#Select ASVs with differential abundance
DA.ASV.1st.list <- LogFC.1st_All %>% dplyr::filter(diff.abn == "TRUE")
nrow(DA.ASV.1st.list) #683ASVs

#Select ASVs:SK<=FL or SS<=FL
EitherFL.ASV.1st.list <- DA.ASV.1st.list %>% dplyr::filter(LogFC_SKFL <= 0 | LogFC_SSFL <= 0)
nrow(EitherFL.ASV.1st.list) #470

#Select FL-ASVs (Category FL): SK<FL & SS<FL
FL.ASV.1st.list <- EitherFL.ASV.1st.list %>% dplyr::filter(LogFC_SKFL < 0 & LogFC_SSFL < 0) %>% 
  dplyr::mutate(Category = "FL") 
nrow(FL.ASV.1st.list) #360

#onlySK-ASVs:SK>SS & SK>FL (SS or FL have structural zero)
onlySK.ASV.1st.list <- EitherFL.ASV.1st.list %>% dplyr::filter(LogFC_SKFL > 0 & LogFC_SKSS > 0) %>% 
  dplyr::mutate(Category = "PA") 
nrow(onlySK.ASV.1st.list) #43

#onlySS-ASVs:SS>SK & SS>FL (SK or FL have structural zero)
onlySS.ASV.1st.list <- EitherFL.ASV.1st.list %>% dplyr::filter(LogFC_SSFL > 0 & LogFC_SKSS < 0) %>% 
  dplyr::mutate(Category = "PA")
nrow(onlySS.ASV.1st.list) #67 

#Select ASVs:SK>FL and SS>FL
BothPA.ASV.1st.list <- DA.ASV.1st.list %>% dplyr::filter(LogFC_SKFL > 0 & LogFC_SSFL > 0) %>% 
  dplyr::mutate(Category = "PA")
nrow(BothPA.ASV.1st.list) #213

#Gather PA ASVs (Category PA):onlySK, onlySS, and BothPA
PA.ASV.1st.list <- dplyr::bind_rows(BothPA.ASV.1st.list, onlySK.ASV.1st.list, onlySS.ASV.1st.list)
nrow(PA.ASV.1st.list) #323ASVs
write.table(PA.ASV.1st.list, "./data/PA.ASV.1st.list.txt", quote = F, sep = "\t", row.names=F, col.names = T, append = F)
#PA.ASV.1st.list:Category PA-ASVs

#Category PA and FL ASVs
PAandFL.ASV.1st.list <- dplyr::bind_rows(PA.ASV.1st.list, FL.ASV.1st.list)
nrow(PAandFL.ASV.1st.list) #683ASVs
unique(PA_FL.ASV.1st.list$Category)

#Category PA, FL, and FL/PA ASVs list
DA.1st.Class.ASV.1st.list <- dplyr::bind_rows(PAandFL.ASV.1st.list, PA_FL.ASV.1st.list)
DA.1st.Class.ASV.1st.list <- replace_na(DA.1st.Class.ASV.1st.list, replace = list(Category = "FL/PA"))
nrow(DA.1st.Class.ASV.1st.list) #841ASVs
unique(DA.1st.Class.ASV.1st.list$Category)
write.table(DA.1st.Class.ASV.1st.list, "./data/DA.1st.Class.ASV.1st.list.txt", quote = F, sep = "\t", row.names=F, col.names = T, append = F)
#DA.1st.Class.ASV.1st.list:The DA results of 841 ASVs

#FL/PA and FL ASVs
FL.PA_FL.ASV.list <- DA.1st.Class.ASV.1st.list %>%
  dplyr::filter(Category == "FL/PA" | Category == "FL")
nrow(FL.PA_FL.ASV.list) #518ASVs
write.table(FL.PA_FL.ASV.list, "./data/FL.PA_FL.ASV.list.txt", quote = F, sep = "\t", row.names=F, col.names = T, append = F)
#FL.PA_FL.ASV.list:The DA results of Category FL/PA and FL-ASVs


#####


###2nd ANCOM-BC (1st Differential abundance analysis)####
##Data to be used
#1.reads.count.pre.df
#2.PA.ASV.list
head(reads.count.pre.df[1:5])
head(PA.ASV.list)

##10% prevalence filter performed to 32 samples (Filter ASVs that do not occur more than 4 times)
reads.count.2ndpre0 <- reads.count.pre.df %>% 
  dplyr::select(contains("_Sink"), contains("_NonS")) %>%
  dplyr::mutate(ASV = rownames(.))
CategoryPA.ASVs <- PA.ASV.list %>% 
  dplyr::select(ASV)
reads.count.2ndpre <- left_join(CategoryPA.ASVs, reads.count.2ndpre0, by = "ASV")
rownames(reads.count.2ndpre) = reads.count.2ndpre$ASV
dim(reads.count.2ndpre)
#32 samples (Sinking and Suspended particle fractions) + ASV column and PA ASVs
#reads.count.2ndpre:before Prevalence filter

reads.count.2ndpre.check0 <- reads.count.2ndpre %>%
  dplyr::select(-ASV) %>%
  mutate_all(function(x) as.numeric(as.character(x)))
reads.count.2ndpre.check <- reads.count.2ndpre.check0 %>%
  dplyr::mutate(sum.occurrence = apply(reads.count.2ndpre.check0, 1, function(x) sum(x > 0))) %>% #Total occurrence of ASVs
  dplyr::filter(sum.occurrence >= 4) %>% #test
  dplyr::select(-sum.occurrence)
dim(reads.count.2ndpre.check)
#321ASVs, 32samples
#reads.count.2ndpre.check; columns are samples (selected samples) and rows are ASVs occurring at least 5 times.

reads.count.2ndpre.df <- reads.count.2ndpre.check
head(reads.count.2ndpre.df[1:5])
dim(reads.count.2ndpre.df) #1046ASVs, 48samples
#reads.count.2ndpre.df; Prevalence filtered read count data. rows are ASVs and columns are samples
write.table(reads.count.2ndpre.df, "./data/reads.count.2ndpre.df.txt", quote = F, sep = "\t", col.names = NA, row.names = TRUE, append = F)


##Import feature.table (Pre-processed count data)
Second.feature.tab <- reads.count.2ndpre.df
head(Second.feature.tab)
dim(Second.feature.tab)
#321 ASVs, 32 samples
#Check the data is a data frame of matrix
class(Second.feature.tab)


##Make sample type data; row is sample
Second.meta.data <- colnames(Second.feature.tab)
Second.meta.data <- as.data.frame(Second.meta.data) %>%
  dplyr::rename(SampleID = Second.meta.data) %>%
  separate(SampleID, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_", remove = F) %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "Sinking")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "Suspended")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "Free-living")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  transform(Fraction= factor(Fraction, levels = c("Free-living", "Suspended", "Sinking"))) %>%
  transform( Layer= factor(Layer, levels = c("ML", "Z1%")))
rownames(Second.meta.data) = Second.meta.data$SampleID
#Check the data is a data frame
class(Second.meta.data)


##Pre-processing
#feature.table = Second.feature.tab
#meta.data = Second.meta.data
#neg.lb = FALSE #Only ASVs that do not occur at all in groups is detected as ASVs with structural zero
lib.cut = 1000 #Cut samples without 1000 read library
zero.cut = 1.0 #Cut ASVs with >90% 0 occurrence
sample.var = "SampleID" #Character. The name of column storing sample IDs.
group.var = "Fraction" #Character. The name of the group indicator.
pre.process.2nd = feature_table_pre_process(Second.feature.tab, Second.meta.data, sample.var, 
                                            group.var, zero.cut, lib.cut, neg.lb = FALSE)
feature.table.2nd = pre.process.2nd$feature.table
group.name.2nd = pre.process.2nd$group.name
group.ind.2nd = pre.process.2nd$group.ind
struc.zero.2nd = pre.process.2nd$structure.zeros
dim(feature.table.2nd)
#321 ASVs 32 samples


##ANCOM-BC main function
#P-value adjust method = "BH"
#tol.EM, max.iterNum, perNum are default values.
#alpha (significance level) = 0.05
adj.method = "BH"
tol.EM = 1e-5
max.iterNum = 100
perNum = 1000
alpha = 0.05
#Record start time and end time.
start_time <- Sys.time()
out.2nd = ANCOM_BC(feature.table.2nd, group.name.2nd, group.ind.2nd, struc.zero.2nd,
                   adj.method, tol.EM, max.iterNum, perNum, alpha)
end_time <- Sys.time()
end_time - start_time


##Output the results of 2nd ANCOM-BC
res.2nd = cbind(ASV = rownames(out.2nd$feature.table), out.2nd$res)
write_csv(res.2nd, "./data/ANCOMBCres.2nd.csv")
res.ANCOM_BC.2nd = data.frame(ASV=rownames(out.2nd$feature.table), out.2nd$res, 
                              struc.zero.2nd[rownames(out.2nd$feature.table), ], 
                              row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
alpha.adj.2nd=0.05/nrow(res.ANCOM_BC.2nd)
critic.val.2nd=qnorm(1-alpha.adj.2nd/2)


##Summarize the results of 2nd ANCOM-BC
#Output all data
LogFC.2nd_All <- res.ANCOM_BC.2nd %>% dplyr::rename(
  "LogFC_SKSS" = "mean.difference (Sinking - Suspended)",
  "se_SKSS" = "se (Sinking - Suspended)",
  "structural.zero.SS" = "structural.zero (Suspended)",
  "structural.zero.SK" = "structural.zero (Sinking)") %>%
  dplyr::mutate(ci.lo.adj_SKSS = LogFC_SKSS - critic.val.2nd*se_SKSS,
                ci.up.adj_SKSS = LogFC_SKSS + critic.val.2nd*se_SKSS)
head(LogFC.2nd_All)
nrow(LogFC.2nd_All) #321ASVs


##Classify ASVs to Category PA to Category SK, SS, and SK/SS according to the results of 2nd ANCOM-BC
#ASVs with Non-differential abundance
SK_SS.ASV.2nd.list <- LogFC.2nd_All %>% dplyr::filter(diff.abn == "FALSE") %>% 
  dplyr::mutate(Category = "SK/SS")
nrow(SK_SS.ASV.2nd.list) #165

#ASVs with differential abundance
DA.ASV.2nd.list <- LogFC.2nd_All %>% dplyr::filter(diff.abn == "TRUE")
names(DA.ASV.2nd.list)
nrow(DA.ASV.2nd.list) #156

#Select SK-ASVs (Category SK): SK>SS
SK.ASV.list <- DA.ASV.2nd.list %>% dplyr::filter(LogFC_SKSS > 0) %>% 
  dplyr::mutate(Category = "SK")
names(SK.ASV.list)
nrow(SK.ASV.list) #103

#Select SS-ASVs (Category SS): SK<SS
SS.ASV.list <- DA.ASV.2nd.list %>% dplyr::filter(LogFC_SKSS < 0) %>% 
  dplyr::mutate(Category = "SS")
names(SS.ASV.list)
nrow(SS.ASV.list) #53

#Marge SKSS-ASVs
marged_SKSS.ASV.list <- dplyr::bind_rows(SK.ASV.list, SS.ASV.list)
nrow(marged_SKSS.ASV.list) #156

#Category SK, SS, and SK/SS ASVs list
DA.2nd.Class.ASV.list <- dplyr::bind_rows(marged_SKSS.ASV.list, SK_SS.ASV.2nd.list)
nrow(DA.2nd.Class.ASV.list) #321ASVs
unique(DA.2nd.Class.ASV.list$Category)
write.table(DA.2nd.Class.ASV.list, "./data/DA.2nd.Class.ASV.list.txt", quote = F, sep = "\t", row.names=F, col.names = T, append = F)
#DA.2nd.Class.ASV.list:The DA results of 321 ASVs 


###Summarize data of classified ASV (Category SK, SS, SK/SS, FL/PA, and FL)####
##Data to be used
head(DA.2nd.Class.ASV.list)
head(FL.PA_FL.ASV.list)
Category.SK.SS.SK_SS <- DA.2nd.Class.ASV.list %>%
  dplyr::select(ASV, Category)
Category.FL.PA_FL <- FL.PA_FL.ASV.list %>%
  dplyr::select(ASV, Category)
DA.res0 <- rbind(Category.SK.SS.SK_SS, Category.FL.PA_FL)
DA.res1 <- left_join(DA.res0, FL.PA_FL.ASV.list, by = "ASV") %>%
  dplyr::select(ASV, Category.x, LogFC_SSFL, LogFC_SKFL, q.val) %>%
  dplyr::rename(q.val.1st = q.val, Category = Category.x) %>%
  dplyr::mutate(star.1st = ifelse(q.val.1st<.001, "***", 
                                  ifelse(q.val.1st<.01, "**",
                                         ifelse(q.val.1st<.05, "*", ""))))
DA.res2 <- left_join(DA.res1, DA.2nd.Class.ASV.list, by = "ASV") %>%
  dplyr::select(ASV, Category.x, LogFC_SSFL, LogFC_SKFL, q.val.1st, star.1st, LogFC_SKSS, q.val) %>%
  dplyr::rename(q.val.2nd = q.val, Category = Category.x) %>%
  dplyr::mutate(star.2nd = ifelse(q.val.2nd<.001, "***", 
                                  ifelse(q.val.2nd<.01, "**",
                                         ifelse(q.val.2nd<.05, "*", ""))))
DA.res.df <- DA.res2
head(DA.res.df)
dim(DA.res.df)
#839ASVs;1046ASVs decreased to 841ASVs at 1st ANCOM-BC Pre-processing. 2ASVs decreased at 2nd Prevalence filter.
#DA.res.df:ASV and Category
#####


####Redundancy analysis####
###packages for Zero permutation and CLR transformation
library(compositions)
library(zCompositions)
library(vegan)

##Data to be used
#1.reads.count.2ndpre.df:SK, SS, and SK/SS ASVs reads count data
#2.DA.2nd.Class.ASV.list:The results of 2nd ANCOM-BC
#3.Env.data
head(reads.count.2ndpre.df[1:5])
dim(reads.count.2ndpre.df)
#321ASVs and 32samples
head(DA.2nd.Class.ASV.list)
head(Env.data)

###Pre-processing (Zeros permutation and CLR transformation)####
##Convert SK, SS, and SK/SS ASVs reads count data to 100 closed data
SKSS.close100 <- as.data.frame(t(reads.count.2ndpre.df)) %>%
  sweep(1, rowSums(.), "/") # RA for relative abundance
SKSS.close100[1:5, 1:5]
SKSS.close100.df <- SKSS.close100*100
rowSums(SKSS.close100.df)
dim(SKSS.close100.df)
#SKSS.close100.df:columns are sampple and rows are ASVs

##Check occurrence pattern of ASVs
head(SKSS.close100.df[1:5])
sum(SKSS.close100.df == 0) #There are 5890 cells with 0
sum(apply(SKSS.close100.df, 2, sum) == 0) #There isn't ASV that never occurring.

##Check occurrence pattern using "zPatterns" function in package "ZCompositions" 
SKSS.reads.pattern.ID <- zPatterns(SKSS.close100.df, label=0, bar.colors = c("red", "blue"), bar.labels = TRUE, cell.colors = c("green", "white"), cell.labels = c("Nondetected", "Observed"), cex.axis=0.8)

##Define thresholds
#dl=Thresholds (Min value in the data)
SKSS.close100.df_rmZero <- SKSS.close100.df
SKSS.close100.df_rmZero[SKSS.close100.df_rmZero == 0] <- NA #remave 0
#Min value table
SKSS.close100.df_min <- min(SKSS.close100.df_rmZero, na.rm = T)
dlMin = rep(SKSS.close100.df_min, 321) #321 ASVs

##Zero permutation using function “multRepl” in package "ZCompositions" 
#delta = 0.65 (default)
#z.warning:90％ zeros
SKSS.close100.df_multReplDif <- multRepl(SKSS.close100.df, label = 0, dl = dlMin,  z.warning=0.90)
dim(SKSS.close100.df_multReplDif)
#Add row names
perZero.SKSS.close100.df <- SKSS.close100.df_multReplDif
rownames(perZero.SKSS.close100.df) = rownames(SKSS.close100.df)
#Check min value
min(perZero.SKSS.close100.df) #min value = dlMin*delta(0.65)
#Check the value of data
check100 <- apply(perZero.SKSS.close100.df, 1, sum)
summary(check100)
write.table(perZero.SKSS.close100.df, "./data/perZero.SKSS.close100.df.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)
#perZero.SKSS.close100.df:Zero perutated SKSS.close100.df

##clr transformation with compositions's function "clr" in package "compositions"
SKSS.CLR <- clr(perZero.SKSS.close100.df)
SKSS.CLR.df <- as.data.frame(SKSS.CLR)
write.table(SKSS.CLR.df, "./data/SKSS.CLR.df.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)
#SKSS.CLR.df: CLR transformed data


###Make data sets####
##Make environmental data for RDA analysis
#select sinking and suspended particle fractions for RDA
head(Env.data)
Env.data.RDA0 <- Env.data %>% 
  dplyr::filter(str_detect(SampleID, "_Sink") | str_detect(SampleID, "_NonS")) %>%
  dplyr::select(SampleID, CSP_TEP, TEP_POC, CSP_POC, Depth, Temperature, Nitrate.con.) %>%
  dplyr::mutate(NAdata = apply(., 1, function(x) sum(is.na(x))))

#Check frequency distributions
Env.data.RDA0.remNA.sel <- Env.data.RDA0 %>% 
  dplyr::filter(NAdata == 0) %>%
  dplyr::select(-SampleID, -NAdata)
pairs(Env.data.RDA0.remNA.sel, lower.panel = NULL, pch=1, col="blue")
ggplot(gather(Env.data.RDA0.remNA.sel), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')

#Check Collinearity (correlations between factors)
cor(Env.data.RDA0.remNA.sel)
Env.RDA.remNA.cor <- as.matrix(cor(Env.data.RDA0.remNA.sel))
Env.RDA.remNA.cor
Env.RDA.remNA.cor.test <- psych::corr.test(Env.data.RDA0.remNA.sel, y = NULL, use = "pairwise",method="pearson",adjust="BH", 
                                     alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
Env.RDA.remNA.cor.test$r
Env.RDA.remNA.cor.test$p
Env.RDA.remNA.cor.test$p.adj
cor_Env.data.RDA0.remNA.sel <- as.data.frame(cor(Env.data.RDA0.remNA.sel))
write.table(cor_Env.data.RDA0.remNA.sel, "./data/cor_Env.data.RDA0.remNA.sel.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)
write.table(as.data.frame(Env.RDA.remNA.cor.test$p), "./data/cor.test_Env.data.RDA0.remNA.sel.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)
write.table(as.data.frame(Env.RDA.remNA.cor.test$p.adj), "./data/cor.test.adj_Env.data.RDA0.remNA.sel.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)
heatmap(abs(cor(Env.data.RDA0.remNA.sel)), 
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

#Select explanatory values and remove sample with NA
Env.data.RDA.NACheck <- Env.data.RDA0 %>% 
  dplyr::select(CSP_TEP, TEP_POC, Depth, Temperature, Nitrate.con.) %>%
  dplyr::mutate(NAdata = apply(., 1, function(x) sum(is.na(x))))
Env.data.RDA <- Env.data.RDA.NACheck %>% 
  dplyr::filter(NAdata == 0) %>%
  dplyr::select(-NAdata)
head(Env.data.RDA)
dim(Env.data.RDA)
#32samples 
summary(Env.data.RDA)
write.table(Env.data.RDA, "./data/Env.data.RDA.txt", quote = F, sep = "\t", col.names = T, append = F)
#Env.data.RDA:The environmental data (Particle condition & Surrounding seawater environment

#Standardize (Scaling and centering) 
Env.data.RDA.sc <- vegan::decostand(Env.data.RDA, method = "standardize")
#Variables are now centered around a mean of 0 and scaled to have a standard deviation of 1
round(apply(Env.data.RDA.sc, 2, mean), 1)
apply(Env.data.RDA.sc, 2, sd)
pairs(Env.data.RDA.sc, lower.panel = NULL, pch=1, col="blue")
#Env.data.RDA.sc:Standardized environmental data (Particle condition & Surrounding seawater environment

##Remove samples with NA from CLR transformed data
SKSS.CLR.RDA0 <- SKSS.CLR.df %>%
  dplyr::mutate(SampleID = rownames(.)) #Add SampleID
RDA.env.samples <- rownames(Env.data.RDA.sc)
RDA.env.samples <- as.data.frame(RDA.env.samples) %>% 
  dplyr::rename(SampleID = RDA.env.samples)
SKSS.CLR.RDA <- left_join(RDA.env.samples, SKSS.CLR.RDA0)
rownames(SKSS.CLR.RDA) = SKSS.CLR.RDA$SampleID
SKSS.CLR.RDA %<>% dplyr::select(-SampleID)
head(SKSS.CLR.RDA)
dim(SKSS.CLR.RDA)
#32samples and 321ASVs
#SKSS.CLR.RDA:CLR transformed data. row (samples) are sorted same as Env.data.RDA.sc


###Redundancy analysis (RDA)####
###RDA using package "vagan"
##First RDA model (explanatory variables;CSP_TEP, TEP_POC, Depth, and Temperature, Nitrate.con.)
RDA.model.first <- vegan::rda(SKSS.CLR.RDA ~ ., data = Env.data.RDA.sc)
summary(RDA.model.first)
RDA.model.first
#Constrained Proportion: variance of Y explained by 39.1(%)
#Unconstrained Proportion: unexplained variance in 60.9(%)
RsquareAdj(RDA.model.first)
#the adjusted R2 (corrected for the number of explanatory variables): 0.3140191


##Select explanatory variables by forward selection of variables:
RDA.model.fwd.sel <- vegan::ordiR2step(vegan::rda(SKSS.CLR.RDA ~ 1, data = Env.data.RDA.sc),
                                    scope = formula(RDA.model.first), 
                                    direction = "forward",
                                    R2scope = TRUE,
                                    pstep = 1000,
                                    trace = FALSE) 
#Check the new model with forward-selected variables
RDA.model.fwd.sel$call
#Depth, Temperature, TEP_POC, and Nitrate.con. were selected.


##Significant RDA model
RDA.model.signif <- rda(SKSS.CLR.RDA ~  Depth + Temperature + TEP_POC + 
                                Nitrate.con., data = Env.data.RDA.sc) 
RDA.model.signif
#Constrained Proportion: variance of Y explained by 36.0(%)
#Unconstrained Proportion: unexplained variance in 64.0(%)
RsquareAdj(RDA.model.signif)
#the adjusted R2 (corrected for the number of explanatory variables): 0.264591


##Save the summary of RDA model
RDA.model.tab <- vegan::anova.cca(RDA.model.signif, step = 1000)
RDA.ExpVar.tab <- vegan::anova.cca(RDA.model.signif, step = 1000, by = "term")
RDA.axes.tab <- vegan::anova.cca(RDA.model.signif, step = 1000, by = "axis")
write.table(RDA.model.tab, "./data/RDA.model.tab.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)
write.table(RDA.ExpVar.tab, "./data/RDA.ExpVar.tab.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)
write.table(RDA.axes.tab, "./data/RDA.axes.tab.txt", quote = F, sep = "\t", col.names = NA, row.names = T, append = F)

##Save the correlation of explanatory variables and scores of sites(samples) and species(ASVs).
axes.perc <- round(100*(summary(RDA.model.signif)$cont$importance[2, 1:2]), 1) #Contributions of axes (%)
species.score <- as.data.frame(scores(RDA.model.signif, display="species", choices=c(1,2), scaling=1)) #Species distribution
sites.score <- as.data.frame(scores(RDA.model.signif, display="sites", choices=c(1,2), scaling=1)) #Sites distribution
ExpVar.cor <- as.data.frame(scores(RDA.model.signif, display="bp", choices=c(1,2), scaling=2)) #Angles of explanatory variables reflect their correlation

##Marge meta data
#ASV data
species.score0 <- species.score %>%
  dplyr::mutate(ASV = rownames(.))
species.score1 <- inner_join(species.score0, DA.res.df, by = "ASV")
species.score2 <- inner_join(species.score1, ASV_16S.Phylum.RE, by = "ASV")
species.score3 <- inner_join(species.score2, ASV_TaxID, by = "ASV")
species.score.df <- species.score3 %>%
  dplyr::mutate(Taxa = gsub(Taxa, pattern="_",replacement = " ")) %>%
  transform(Category = factor(Category, levels = c("SK", "SS", "SK/SS"))) %>%
  transform(Taxa= factor(Taxa, levels = Taxonomy.level)) %>%
  dplyr::select(ASV, everything())
#species.score.df;ASVs, RDA1, RDA2, Category, and Taxonomy
write.table(species.score.df, "./table/species.score.df.txt", quote = F, sep = "\t", col.names = T, row.names = F, append = F)

species.score.Alpha <- species.score.df %>% dplyr::filter(Taxa == "Alphaproteobacteria")
species.score.Gamma <- species.score.df %>% dplyr::filter(Taxa == "Gammaproteobacteria")
species.score.Bac <- species.score.df %>% dplyr::filter(Taxa == "Bacteroidia")
species.score.Plan <- species.score.df %>% dplyr::filter(Taxa == "Planctomycetota")
species.score.Bdell <- species.score.df %>% dplyr::filter(Taxa == "Bdellovibrionota")
species.score.Desul <- species.score.df %>% dplyr::filter(Taxa == "Desulfobacterota")



#Sample data
sites.score.df <- sites.score %>%
  dplyr::mutate(SampleID = rownames(.)) %>%
  separate(SampleID, c("CruiseID", "Station", "Layer", "Fraction"), sep = "_", remove = F) %>%
  dplyr::mutate(Fraction = gsub(Fraction, pattern="Sink",replacement = "SK fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="NonS",replacement = "SS fraction")) %>% 
  dplyr::mutate(Fraction = gsub(Fraction, pattern="FL",replacement = "FL fraction")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L1",replacement = "ML")) %>%
  dplyr::mutate(Layer = gsub(Layer, pattern="L2",replacement = "Z1%")) %>%
  transform(Layer= factor(Layer, levels = c("ML", "Z1%"))) %>%
  transform(Fraction= factor(Fraction, levels = c("SK fraction", "SS fraction", "FL fraction"))) %>%
  unite(Fraction_Layer, c("Fraction", "Layer"), sep = "-", remove = F)
head(sites.score.df)
#sites.score.df:SampleID, RDA1, RDA2, Station, Layer, Fraction, and Fraction_Layer

#Explanatory variables
ExpVar.cor.df <- ExpVar.cor %>%
  dplyr::mutate(Variables = rownames(.)) %>%
  dplyr::mutate(Variables = gsub(Variables, pattern="TEP_POC",replacement = "TEP : POC ratio")) %>%
  dplyr::mutate(Variables = gsub(Variables, pattern="Nitrate.con.",replacement = "Nitrate conc.")) %>%
  dplyr::mutate(text.posi.RDA1 = c(-0.4, -0.8, -1.0, 0.25)) %>%
  dplyr::mutate(text.posi.RDA2 = c(0.95, -0.5, 0.2, 0.8))
#ExpVar.cor.df:Explanatory variables, RDA1, RDA2, text.position

#Ranges of RDA1 and RDA
scores.df <- dplyr::bind_rows(species.score, sites.score, ExpVar.cor)
range(scores.df$RDA1)
range(scores.df$RDA2)


#####



