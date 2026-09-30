
#############################################################################################################
# Genomic and phenomic predictions of the hybrid values of nursery-grown parental lines
# Scenario CV_newFamilies
#############################################################################################################

rm(list=ls())

# Load packages
library(data.table)
library(rrBLUP)
library(sommer)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(matrixcalc)

setwd("D:/sauvegarde/Projets/2022_SPS_PhenoMaize/Manuscrit_drones/")

#############################################################################################################
#############################################################################################################
# Load and prepare data
#############################################################################################################
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

# trait heritabilities (from Lorenzi et al. 2022)

herit = c(0.89,0.93,0.93,0.86)
names(herit) = c("DMY","DMC","Silking","MFU")

# Coding of the Flint and Dent parents of each hybrid

flint=unique(pheno2$CodeParentC) # flint parents
dente=unique(pheno2$CodeParentD) # dent parents

# Pretreated Spectra from the nursery (Moulon 2024)
nurs_leaf_24 = readRDS("Data_Publications/spectra_transformed_maize_leaf.Rds")
nurs_grain_24 = readRDS("Data_Publications/spectra_transformed_maize_Grain.Rds")

# Secondary traits from drone images 

drone = read.csv("Data_Publications/adjmeans_secondarytraits_maize.csv")
rownames(drone) = drone[,1]
drone = as.matrix(drone[,-1])

herit_drone = readRDS("Data_Publications/Heritabilities_secondarytraits_maize.Rds")
herit_drone = unlist(herit_drone)
drone = drone[,herit_drone>=0.4]

# Drone images (latent variables VJEPA)

vec_s=read.csv("Data_Publications/predictions_AE_maize.csv")

rownames(vec_s)=vec_s[,1]
vec_s = vec_s[,-1]

herit_s=read.csv("Data_Publications/Heritabilities_AE_LatentVariables_maize.csv") # heritability of the latent variables from VJEPA
vec_s = vec_s[,herit_s[,2]>=0.4]

#############################################################################################################
#############################################################################################################
# Compute kinships and phenomic covariance matrices of the Dent and Flint parental lines
#############################################################################################################
#############################################################################################################

# list of pre-treated phenomic data
phenomic = list(nurs_leaf=nurs_leaf_24$norm_der1$NIRS,nurs_grain=nurs_grain_24$norm_der1$NIRS,drone=t(drone),AE_tmax=t(vec_s))

for (l in 1:length(phenomic)) {colnames(phenomic[[l]])[match("PSMC1287 (F9024)",colnames(phenomic[[l]]))] <- "PSMC1287"}

# keep only the genotypes having all kinds of phenomic variables
commun = names(which(table(c(colnames(phenomic[[1]]),colnames(phenomic[[2]]),colnames(phenomic[[3]]),colnames(phenomic[[4]])))==4))

phenomic = lapply(phenomic,function(a) a[,match(commun,colnames(a))])

# Hyperspectral Covariance matrices
############################################################

# name of the genotypes
ind=colnames(phenomic[[1]])

# Center and scale each wavelength, drone trait or latent variable
phenomic = lapply(phenomic,function(a) scale(t(a))) 

# Compute covariance matrices
H_mat= lapply(phenomic,function(a) a%*%t(a)/ncol(a))

# Dent genotypes
indD = colnames(H_mat[[1]])[which(substr(colnames(H_mat[[1]]),1,4)=="PSMD")]

# Flint genotypes
indF = colnames(H_mat[[1]])[which(substr(colnames(H_mat[[1]]),1,4)=="PSMC")]

# Extract Dent and Flint phenomic covariance matrices
H_matD = lapply(H_mat,function(a) a[match(indD,rownames(a)),match(indD,colnames(a))]) # Covariance matrix between dent lines
H_matF = lapply(H_mat,function(a) a[match(indF,rownames(a)),match(indF,colnames(a))]) # Covariance matrix between flint lines


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

# Dent kinship:
kinD=Zdose%*%t(Zdose)/(sum(Denom)) # VanRaden estimate

# Genomic kinship : Flint
##############################
whichflint  = match(indF,rownames(geno))
genoF = geno[whichflint,]

Freq=apply(genoF,2,function(a) mean(a,na.rm=T))

sum(is.na(genoF)) # No missing value

Zdose=genoF
Freq=apply(Zdose,2,function(a) mean(a))
Denom=Freq*(1-Freq)

for (t in 1:ncol(Zdose)) {
  Zdose[,t]=Zdose[,t]-Freq[t]
}

# Flint kinship:
kinF=Zdose%*%t(Zdose)/(sum(Denom)) # VanRaden estimate

# Scale the kinships and covariance matrices
##############################

kinF <- kinF/(matrix.trace(kinF)/ncol(kinF))
kinD <- kinD/(matrix.trace(kinD)/ncol(kinD))

H_matD = lapply(H_matD,function(a) a/(matrix.trace(a)/ncol(a)))
H_matF = lapply(H_matF,function(a) a/(matrix.trace(a)/ncol(a)))

# Keep the phenotypes of the hybrids whose both parents are NIRSed
pheno2=pheno2[which(pheno2$CodeParentC%in%indF & pheno2$CodeParentD%in%indD),] 

H_matD= lapply(H_matD,function(a) a[match(pheno2$CodeParentD,rownames(a)),match(pheno2$CodeParentD,colnames(a))])
H_matF= lapply(H_matF,function(a) a[match(pheno2$CodeParentC,rownames(a)),match(pheno2$CodeParentC,colnames(a))])

kinG_D = kinD[match(pheno2$CodeParentD,rownames(kinD)),match(pheno2$CodeParentD,colnames(kinD))]
kinG_F = kinF[match(pheno2$CodeParentC,rownames(kinF)),match(pheno2$CodeParentC,colnames(kinF))]


#############################################################################################################
#############################################################################################################
# Predictions 
#############################################################################################################
#############################################################################################################


#######################################################################################################
# Scenario 2 : Predict completely new individuals: the parents of the predicted hybrids are from different
#         families than the parents of the hybrids in the training set (CV_newFamilies)
#######################################################################################################

trait=c("DMY","DMC","Silking","MFU")
Ntrait = length(trait)

# Names of your prediction models
modname=c("GdGf","HdHf_nursL","HdHf_nursG","HdHf_nursGL","HdHf_drone","LV","Trait_LV")


famillesF=substr(pheno2$CodeParentC,5,5)
famillesD=substr(pheno2$CodeParentD,5,5)
famillesDF = paste(famillesD,famillesF,sep="_")
famillesDF_unik = unique(famillesDF)

Nfam=length(famillesDF_unik)

Nvalid <- rep(NA,Nfam)  # Nb of hybrids in each new family
for (rep in 1:Nfam) {
  valid = which(famillesDF==famillesDF_unik[rep])
  Nvalid[rep]= length(valid)
}

# determine the combinations of families for which the number of training genotypes is at least 100 
# and the number of validation genotypes is at least equal to 30
Fam1 = combn(famillesDF_unik,1)
Fam2 = combn(famillesDF_unik,2)
Fam3 = combn(famillesDF_unik,3)
Fam4 = combn(famillesDF_unik,4)
Fam5 = combn(famillesDF_unik,5)

Fam = list(Fam1,Fam2,Fam3,Fam4,Fam5)

CombOK = list()

for (t in 1:length(Fam)) { # Loop on the nb of families combined in the validation set
  print(t)
  Ncomb = ncol(Fam[[t]])
  
  Size = matrix(NA,2,ncol(Fam[[t]]))
  rownames(Size) = c("valid","training")
  
  for (l in 1:ncol(Fam[[t]])) { # Loop on each combination
    
  Size[1,l] =   sum(Nvalid[match(Fam[[t]][,l],famillesDF_unik)]) # Nb of hybrids in the validation set
  FamD_valid = substr(Fam[[t]][,l],1,1)
  FamF_valid = substr(Fam[[t]][,l],3,3)    
  
  Fam_training = famillesDF_unik[which(substr(famillesDF_unik,1,1)%notin%FamD_valid & substr(famillesDF_unik,3,3)%notin%FamF_valid) ] # Families of the training set
  
  Size[2,l] =   sum(Nvalid[match(Fam_training,famillesDF_unik)]) # Nb of hybrids in the training set
  }

  CombOK[[t]] = Fam[[t]][,which(Size[1,]>=20 & Size[2,]>=120)] # Combination of families with the necessary conditions
}

CombOK[[1]] = as.matrix(CombOK[[1]])
CombOK[[4]] = as.matrix(CombOK[[4]])

# Predictions
##############################################

CombOK = CombOK[which(unlist(lapply(CombOK,ncol))!=0)]

NcombOK = sum(unlist(lapply(CombOK,ncol)))
NcombOK

acc <-matrix(NA,NcombOK,Ntrait)

accu <- list(acc,acc,acc,acc,acc,acc,acc)

Nvalid=rep(NA,NcombOK)
Ncalib=rep(NA,NcombOK)

cpt=0
for (cb in 1:length(CombOK)) { # loop on the different number of families combined in the validation set
  print(paste("Nb: ",cb,sep=""))

  for (l in 1:ncol(CombOK[[cb]])) {   # Loop on the different combinations of families in the validation set
  cpt =cpt+1  
  valid = which(famillesDF%in%CombOK[[cb]][,l]) # Predicted genotypes
  hideD = which(famillesD%in%substr(CombOK[[cb]][,l],1,1))
  hideF = which(famillesF%in%substr(CombOK[[cb]][,l],3,3))
  hide = unique(c(hideD,hideF))
  
  calib = seq(1,nrow(pheno2))[-hide]
  
  Nvalid[cpt]= length(valid)
  Ncalib[cpt]= length(calib)
  
  for (tr in 1:Ntrait) { # loop on the traits
    print(tr)
    
    trait2pred=pheno2[,match(trait[tr],colnames(pheno2))]
    traitNA=trait2pred
    traitNA[hide] =NA
    
    # data preparation
    
    traitcal=data.frame(genoF=pheno2$CodeParentC,genoD=pheno2$CodeParentD,pheno=traitNA) 
    traitcal[hide,3]=NA   # Hide validation set
    traitcal$genoF=as.factor(traitcal$genoF)
    traitcal$genoF2=as.factor(traitcal$genoF)
    traitcal$genoF3=as.factor(traitcal$genoF)
    traitcal$genoD=as.factor(traitcal$genoD)
    traitcal$genoD2=as.factor(traitcal$genoD)
    traitcal$genoD3=as.factor(traitcal$genoD)
    
    
    # REML Model 1 (genomic)
    mixmod_sommer <- mmer(pheno ~ 1,
                          random=~vs(genoD,Gu=kinG_D)+   vs(genoF,Gu=kinG_F),
                          rcov=~units,
                          data=traitcal,tolParInv = 1e-02,verbose=FALSE)
    
    # Predictions
    pred1=mixmod_sommer$U$`u:genoD`$pheno[match(traitcal$genoD,names(mixmod_sommer$U$`u:genoD`$pheno))] + 
      mixmod_sommer$U$`u:genoF`$pheno[match(traitcal$genoF,names(mixmod_sommer$U$`u:genoF`$pheno))] 
    accu[[1]][cpt,tr] = cor(pred1[valid],trait2pred[valid])  # Predictive ability
    
    
    # REML model 2 (Phenomic - NIRS - nursery leaves)
    mixmod_sommer <- mmer(pheno ~ 1,
                          random=~vs(genoD,Gu=H_matD[[1]])+   vs(genoF,Gu=H_matF[[1]]),
                          rcov=~units,
                          data=traitcal,tolParInv = 1e-02,verbose=FALSE)
    
    # Predictions
    pred1=mixmod_sommer$U$`u:genoD`$pheno[match(traitcal$genoD,names(mixmod_sommer$U$`u:genoD`$pheno))] + 
      mixmod_sommer$U$`u:genoF`$pheno[match(traitcal$genoF,names(mixmod_sommer$U$`u:genoF`$pheno))] 
    accu[[2]][cpt,tr] = cor(pred1[valid],trait2pred[valid])  # Predictive ability
    
    
    # REML model 3 (Phenomic - NIRS - nursery grains)
    mixmod_sommer <- mmer(pheno ~ 1,
                          random=~vs(genoD,Gu=H_matD[[2]])+   vs(genoF,Gu=H_matF[[2]]),
                          rcov=~units,
                          data=traitcal,tolParInv = 1e-02,verbose=FALSE)
    
    # Predictions
    pred1=mixmod_sommer$U$`u:genoD`$pheno[match(traitcal$genoD,names(mixmod_sommer$U$`u:genoD`$pheno))] + 
      mixmod_sommer$U$`u:genoF`$pheno[match(traitcal$genoF,names(mixmod_sommer$U$`u:genoF`$pheno))] 
    accu[[3]][cpt,tr] = cor(pred1[valid],trait2pred[valid])  # Predictive ability
    
    # REML model 4 (Phenomic - NIRS - nursery leaves + grains)
    mixmod_sommer <- mmer(pheno ~ 1,
                          random=~vs(genoD,Gu=H_matD[[1]])+ vs(genoD2,Gu=H_matD[[2]])+ vs(genoF,Gu=H_matF[[1]]) +  vs(genoF2,Gu=H_matF[[2]]),
                          rcov=~units,
                          data=traitcal,tolParInv = 1e-02,verbose=FALSE)
    
    # Predictions
    pred1=mixmod_sommer$U$`u:genoD`$pheno[match(traitcal$genoD,names(mixmod_sommer$U$`u:genoD`$pheno))] + 
      mixmod_sommer$U$`u:genoD2`$pheno[match(traitcal$genoD2,names(mixmod_sommer$U$`u:genoD2`$pheno))] + 
      mixmod_sommer$U$`u:genoF`$pheno[match(traitcal$genoF,names(mixmod_sommer$U$`u:genoF`$pheno))] +
      mixmod_sommer$U$`u:genoF2`$pheno[match(traitcal$genoF2,names(mixmod_sommer$U$`u:genoF2`$pheno))]
    accu[[4]][cpt,tr] = cor(pred1[valid],trait2pred[valid])  # Predictive ability
    
    # REML Model 5 (Phenomic - Drone - secondary traits)
    mixmod_sommer <- mmer(pheno ~ 1,
                          random=~vs(genoD,Gu=H_matD[[3]])+   vs(genoF,Gu=H_matF[[3]]),
                          rcov=~units,
                          data=traitcal,tolParInv = 1e-02,verbose=FALSE)
    
    # Predictions
    pred1=mixmod_sommer$U$`u:genoD`$pheno[match(traitcal$genoD,names(mixmod_sommer$U$`u:genoD`$pheno))] + 
      mixmod_sommer$U$`u:genoF`$pheno[match(traitcal$genoF,names(mixmod_sommer$U$`u:genoF`$pheno))] 
    accu[[5]][cpt,tr] = cor(pred1[valid],trait2pred[valid])  # Predictive ability
    
    # REML Model 6 (Phenomic - Drone - VJEPA)
    mixmod_sommer <- mmer(pheno ~ 1,
                          random=~vs(genoD,Gu=H_matD[[4]])+   vs(genoF,Gu=H_matF[[4]]),
                          rcov=~units,
                          data=traitcal,tolParInv = 1e-02,verbose=FALSE)
    
    # Predictions
    pred1=mixmod_sommer$U$`u:genoD`$pheno[match(traitcal$genoD,names(mixmod_sommer$U$`u:genoD`$pheno))] + 
      mixmod_sommer$U$`u:genoF`$pheno[match(traitcal$genoF,names(mixmod_sommer$U$`u:genoF`$pheno))] 
    accu[[6]][cpt,tr] = cor(pred1[valid],trait2pred[valid])  # Predictive ability
    
    
    # REML Model 7 (Phenomic - Drone - secondary traits + VL from VJEPA)
    mixmod_sommer <- mmer(pheno ~ 1,
                          random=~vs(genoD,Gu=H_matD[[3]])+ vs(genoD2,Gu=H_matD[[4]])+ vs(genoF,Gu=H_matF[[3]]) +  vs(genoF2,Gu=H_matF[[4]]),
                          rcov=~units,
                          data=traitcal,tolParInv = 1e-02,verbose=FALSE)
    
    # Predictions
    pred1=mixmod_sommer$U$`u:genoD`$pheno[match(traitcal$genoD,names(mixmod_sommer$U$`u:genoD`$pheno))] + 
      mixmod_sommer$U$`u:genoD2`$pheno[match(traitcal$genoD2,names(mixmod_sommer$U$`u:genoD2`$pheno))] + 
      mixmod_sommer$U$`u:genoF`$pheno[match(traitcal$genoF,names(mixmod_sommer$U$`u:genoF`$pheno))] +
      mixmod_sommer$U$`u:genoF2`$pheno[match(traitcal$genoF2,names(mixmod_sommer$U$`u:genoF2`$pheno))]
    accu[[7]][cpt,tr] = cor(pred1[valid],trait2pred[valid])  # Predictive ability
    
  }
  }
}


# Compute accuracy by dividing predictive ability by square root of heritability
for (tr in 1:length(accu)) {
  accu[[tr]][,1] = accu[[tr]][,1]/sqrt(herit[1]) 
  accu[[tr]][,2] = accu[[tr]][,2]/sqrt(herit[2]) 
  accu[[tr]][,3] = accu[[tr]][,3]/sqrt(herit[3]) 
  accu[[tr]][,4] = accu[[tr]][,4]/sqrt(herit[4]) 
}

# Summarize results

table=sapply(accu,function(a) apply(a,2,mean))
colnames(table)=modname
rownames(table)=trait
round(table,2)

# write.table(table,file="Figure_Publication/Accuracy_CV_newFamilies.csv")

# Produce figure
##############################################

accu2 = accu
modname <- c("GBLUP","NIRS (leaf)","NIRS (grain)","NIRS (leaf+grain)","UAV (traits)","UAV (LV)","UAV (traits+LV)")
names(accu2) = modname
for (i in 1:length(accu2)) {colnames(accu2[[i]]) = trait}

df_long <- lapply(seq_along(accu2), function(i) {
  as.data.frame(accu2[[i]]) %>%
    mutate(rep = 1:n()) %>%
    pivot_longer(
      cols = -rep,
      names_to = "trait",
      values_to = "PA"
    ) %>%
    mutate(vol = modname[i])
}) %>%
  bind_rows()

colnames(df_long)=c("rep","trait","PA","model") # PA = Prediction accuracy

df_long$model = factor(df_long$model,levels=modname)

df_long$trait[which(df_long$trait=="FLOF")]="Silking"
df_long$trait[which(df_long$trait=="UFL")]="MFU"

df_long$trait = factor(df_long$trait,levels=c("DMY","DMC","Silking","MFU"))

# Split GS and PS
df_genomic <- df_long[df_long$model == "GBLUP", ]
df_phenomic <- df_long[df_long$model != "GBLUP", ]

# plot
p <- ggplot() +
  
  # Bloc 1 : Genomic
  geom_boxplot(
    data = df_genomic,
    aes(x = model, y = PA, fill = model)
  ) +
  scale_fill_manual(
    name = "Genomic prediction",
    values = c("GBLUP" = "#009E73")
  ) +
  
  # 👉 reset du fill
  new_scale_fill() +
  
  # Bloc 2 : Phenomic
  geom_boxplot(
    data = df_phenomic,
    aes(x = model, y = PA, fill = model)
  ) +
  scale_fill_manual(
    name = "Phenomic prediction",
    values = c(
      "NIRS (leaf)" = "#F0E442",
      "NIRS (grain)" = "#E69F00",
      "NIRS (leaf+grain)" = "#D55E00",
      "UAV (traits)" = "#56B4E9",
      "UAV (LV)" = "#0072B2",
      "UAV (traits+LV)" = "#5C4FA3"
    )
  ) +
  
  facet_wrap(~trait, ncol = 2) +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL, y = "Prediction accuracy")+
  
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box.background = element_rect(colour = "black", linewidth = 0.5),
    legend.key.size = unit(1, "lines"),
    axis.text.x = element_blank(),
    panel.spacing = unit(2, "lines"),
    strip.text = element_text(size = 16) ,
    legend.text = element_text(size = 16),   # Taille du texte des éléments
    legend.title = element_text(size = 16),
    axis.title.y = element_text(
      size = 14,
      face = "bold")
  )


p <- p +
  labs(
    title = "Scenario CV_newFamilies",
    y = "Prediction accuracy"
  ) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5,
      colour = "grey30"
    )
  )

print(p)

# ggsave("Figure_Publication/Figure_Accuracy_CV_newFamilies.jpg",
# plot = p,
# width = 12,
# height = 8,
# units = "in",
# dpi = 300,
# quality = 100)


# end