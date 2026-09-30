
#############################################################################################################
# Estimate genomic heritability and Principal Component Analyses of NIRS and UAV data (secondary traits and latent variables from autoencoder)
# 
#############################################################################################################

rm(list=ls())

library(data.table)
library(rrBLUP)
library(sommer)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(MM4LMM)
library(matrixcalc)

setwd("D:/sauvegarde/Projets/2022_SPS_PhenoMaize/Manuscrit_drones/")

source("Script_Publication/BuildMMEstInputs.R") # Script to run MM4LMM

#############################################################################################################
# Load data
#############################################################################################################

# genotyping data
geno=fread("Data_Publications/GenoArticle1Alizarine_maize.csv")

indgeno=as.character(geno$V1)
geno = as.matrix(geno[,-1])
rownames(geno) = indgeno

geno=geno/2

# phenotyping data
pheno2 = read.csv("Data_Publications/Phenotypes_maize.csv",header=T)
pheno2 = pheno2[,-1]

flint=unique(pheno2$CodeParentC) # flint parents
dente=unique(pheno2$CodeParentD) # dent parents

# Pretreated Spectra from the nursery (Moulon 2024)
nurs_leaf_24 = readRDS("Data_Publications/spectra_transformed_maize_leaf.Rds")
nurs_grain_24 = readRDS("Data_Publications/spectra_transformed_maize_Grain.Rds")

# UAV secondary traits

drone = read.csv("Data_Publications/adjmeans_secondarytraits_maize.csv")
rownames(drone) = drone[,1]
drone = as.matrix(drone[,-1])

herit_drone = readRDS("Data_Publications/Heritabilities_secondarytraits_maize.Rds")
herit_drone = unlist(herit_drone)

# UAV latent variable from autoencoder

vec_s=read.csv("Data_Publications/predictions_AE_maize.csv")

rownames(vec_s)=vec_s[,1]
vec_s = vec_s[,-1]

herit_s=read.csv("Data_Publications/Heritabilities_AE_LatentVariables_maize.csv") # heritability of the latent variables from the autoencoder

#############################################################################################################
# Compute kinships and covariance matrices
#############################################################################################################

# list of phenomic data
spectres = list(nurs_leaf=nurs_leaf_24$norm_der1$NIRS,nurs_grain=nurs_grain_24$norm_der1$NIRS,UAV_traits=t(drone),VJEPA=t(vec_s))

for (l in 1:length(spectres)) {colnames(spectres[[l]])[match("PSMC1287 (F9024)",colnames(spectres[[l]]))] <- "PSMC1287"}

commun = names(which(table(c(colnames(spectres[[1]]),colnames(spectres[[2]]),colnames(spectres[[3]]),colnames(spectres[[4]])))==4))

spectres = lapply(spectres,function(a) a[,match(commun,colnames(a))])

ind=colnames(spectres[[1]])

indD = ind[which(substr(ind,1,4)=="PSMD")] # Dent parental lines

indF = ind[which(substr(ind,1,4)=="PSMC")] # Flint parental lines

# Phenomic Covariance matrices
############################################################

spectres_cs = lapply(spectres,function(a) scale(t(a))) # Center and scale each wavelength or trait

H_mat= lapply(spectres_cs,function(a) a%*%t(a)/ncol(a)) # Compute covariance matrices between parental lines

# Split into a Dent and Flint covariance matrix
H_matD = lapply(H_mat,function(a) a[match(indD,rownames(a)),match(indD,colnames(a))]) # Covariance matrix between dent lines
H_matF = lapply(H_mat,function(a) a[match(indF,rownames(a)),match(indF,colnames(a))]) # Covariance matrix between flint lines

spectres_Dent = lapply(spectres,function(a) {a = a[,match(indD,colnames(a))] } )
spectres_Flint = lapply(spectres,function(a) {a = a[,match(indF,colnames(a))] } )

spectres_Dent_cs = lapply(spectres_Dent,function(a) scale(t(a))) 
spectres_Flint_cs = lapply(spectres_Flint,function(a) scale(t(a))) 

H_matD= lapply(spectres_Dent_cs,function(a) a%*%t(a)/ncol(a)) # Compute covariance matrices
H_matF= lapply(spectres_Flint_cs,function(a) a%*%t(a)/ncol(a)) # Compute covariance matrices

# Genomic kinship : Dent
############################################################

whichdent  = match(indD,rownames(geno))
genoD = geno[whichdent,]

Freq=apply(genoD,2,function(a) mean(a,na.rm=T))

Zdose=genoD
Freq=apply(Zdose,2,function(a) mean(a))
Denom=Freq*(1-Freq)

for (t in 1:ncol(Zdose)) {
  Zdose[,t]=Zdose[,t]-Freq[t]
}

kinD=Zdose%*%t(Zdose)/(sum(Denom))

# Genomic kinship : Flint
##############################

whichflint  = match(indF,rownames(geno))
genoF = geno[whichflint,]

Freq=apply(genoF,2,function(a) mean(a,na.rm=T))

Zdose=genoF
Freq=apply(Zdose,2,function(a) mean(a))
Denom=Freq*(1-Freq)

for (t in 1:ncol(Zdose)) {
  Zdose[,t]=Zdose[,t]-Freq[t]
}

kinF=Zdose%*%t(Zdose)/(sum(Denom))

# Scale the kinship and covariance matrices
##############################

kinF <- kinF/(matrix.trace(kinF)/ncol(kinF))
kinD <- kinD/(matrix.trace(kinD)/ncol(kinD))

H_matD = lapply(H_matD,function(a) a/(matrix.trace(a)/ncol(a)))
H_matF = lapply(H_matF,function(a) a/(matrix.trace(a)/ncol(a)))

phenomicDent = lapply(spectres,function(a) a[,match(indD,colnames(a))])
phenomicFlint = lapply(spectres,function(a) a[,match(indF,colnames(a))])

#############################################################################################################
# Estimation of genomic heritability for each phenomic predictor
#############################################################################################################

# Dent

phenomic = phenomicDent
kinm =kinD

herit=list()
  
for (i in 1:length(phenomic)) {
  print(i)
  herit_tmp = rep(NA,nrow(phenomic[[i]]))
  for (t in 1:nrow(phenomic[[i]])) {
  
    mm = mixed.solve(y=phenomic[[i]][t,],K=kinm)
    herit_tmp[t] = mm$Vu/(mm$Vu+mm$Ve)
  }
  herit[[i]] = herit_tmp
}

heritD = herit
names(heritD) =c("NIRS_leaf","NIRS_grain","UAV (traits)","UAV (LV)")


# Flint

phenomic = phenomicFlint
kinm =kinF

herit=list()

for (i in 1:length(phenomic)) {
  print(i)
  herit_tmp = rep(NA,nrow(phenomic[[i]]))
  for (t in 1:nrow(phenomic[[i]])) {
    
    mm = mixed.solve(y=phenomic[[i]][t,],K=kinm)
    herit_tmp[t] = mm$Vu/(mm$Vu+mm$Ve)
  }
  herit[[i]] = herit_tmp
}

heritF = herit
names(heritF) =c("NIRS_leaf","NIRS_grain","UAV (traits)","UAV (LV)")


###########################################################
# Number of PCA axes with h²>0.1, for each kind of phenomic predictor
###########################################################

# PCA on the Dent parental lines

resuACP_D = rep(NA,length(phenomicDent))
herit_ax_D_tmp = matrix(NA,68,length(phenomicDent))
herit_ax_D = rep(NA,length(phenomicDent))

for (i in 1:length(phenomicDent)) {  # Loop on the different kinds of phenomic predictors
  
res <- prcomp(t(phenomicDent[[i]]), center = TRUE, scale = TRUE)
var_explained <- res$sdev^2
prop_var <- var_explained / sum(var_explained)
test = max(which(cumsum(prop_var)<=0.9))+1 # test the number of PCA axes to reach 90% of explained variance
if (length(test)==0) {test=1}

resuACP_D[i]=test

for (tr in 1:ncol(res$rotation)) {
  
  mm = mixed.solve(y=res$x[,tr],K=kinD)
  herit_ax_D_tmp[tr,i] = mm$Vu/(mm$Vu+mm$Ve)
  
}

herit_ax_D[i] = sum(herit_ax_D_tmp[,i]>=0.1)

}

herit_ax_D

# PCA on the Flint parental lines

resuACP_F = rep(NA,length(phenomicFlint))
herit_ax_F_tmp = matrix(NA,78,length(phenomicFlint))
herit_ax_F = rep(NA,length(phenomicFlint))

for (i in 1:length(phenomicFlint)) {  # Loop on the different kinds of phenomic predictors
  
  res <- prcomp(t(phenomicFlint[[i]]), center = TRUE, scale = TRUE)
  var_explained <- res$sdev^2
  prop_var <- var_explained / sum(var_explained)
  test = max(which(cumsum(prop_var)<=0.9))+1 # test the number of PCA axes to reach 90% of explained variance
  resuACP_F[i]=test
  
  for (tr in 1:ncol(res$rotation)) {
    
    mm = mixed.solve(y=res$x[,tr],K=kinF)
    herit_ax_F_tmp[tr,i] = mm$Vu/(mm$Vu+mm$Ve)
    
  }
  
  herit_ax_F[i] = sum(herit_ax_F_tmp[,i]>=0.1)
}

herit_ax_F

###########################################################
# Synthesis of the results
###########################################################

# summary stat on genomic heritability
resD <- sapply(heritD, summary)
resF <- sapply(heritF, summary)

resultat = cbind(resD,resF)

# total nb of varialbes
NvarD = sapply(spectres_Dent_cs,function(a) dim(a)[2])
NvarF = sapply(spectres_Flint_cs,function(a) dim(a)[2])
Nvar=c(NvarD,NvarF)

resultat = rbind(Nvar,resultat)

# Nb variables avec h2>0.4
highh2_D <- sapply(heritD, function(a) sum(a>=0.4))
highh2_F <- sapply(heritF, function(a) sum(a>=0.4))

resultat = rbind(resultat,c(highh2_D,highh2_F))

# Nb of axes PCA with h2>=0.1

resuACP = c(herit_ax_D,herit_ax_F)

resultat = rbind(resultat,resuACP)

resultat = resultat[-match(c("Min.","Mean","Max."),rownames(resultat)),]

rownames(resultat) = c("Nb variables","1st Qu.","Median","3rd Qu.","Nb h²>0.4","Nb PCA Axes")


 # write.csv(resultat,file="Figure_Publication/SummaryStatistics.csv")



# end