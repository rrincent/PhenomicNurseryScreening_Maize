
#############################################################################################################
# Phenomic predictions with different UAV acquisition dates, nb of dates, using UAV secondary traits
#############################################################################################################

rm(list=ls())

library(data.table)
library(rrBLUP)
library(sommer)
library(tidyverse)
library(gtools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(matrixcalc)

setwd("D:/sauvegarde/Projets/2022_SPS_PhenoMaize/Manuscrit_drones/")

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

# trait heritabilities (from Lorenzi et al. 2022)

herit = c(0.89,0.93,0.93,0.86)
names(herit) = c("DMY","DMC","Silking","MFU")

# UAV secondary traits

drone = readRDS("Data_Publications/adjmeans_drone_per_flight_maize.Rds")

for (i in 1:length(drone)) {rownames(drone[[i]])[which(rownames(drone[[i]])=="psmc1287_f9024")] = "psmc1287"}

rnam = rownames(drone[[1]])
drone_all = matrix(unlist(drone),nrow=length(rnam))
rownames(drone_all) = toupper(rnam)


#############################################################################################################
# Compute kinships and covariance matrices
#############################################################################################################


alldata = lapply(drone,function(a) t(a))
alldata[[7]]=t(drone_all) # combination of all dates

alldata = lapply(alldata,function(a) {colnames(a)=toupper(colnames(a))
return(a)})

for (l in 1:length(alldata)) {colnames(alldata[[l]])[match("PSMC1287 (F9024)",colnames(alldata[[l]]))] <- "PSMC1287"}
for (l in 1:length(alldata)) {colnames(alldata[[l]])[match("PSMC1287_F9024",colnames(alldata[[l]]))] <- "PSMC1287"}

commun = read.csv("Data_Publications/Geno_Alldata_maize.csv",sep=";") # names of the parental lines with phenomic data
commun = commun[,2]

alldata = lapply(alldata,function(a) a[,match(commun,colnames(a))])

# Phenomic covariance matrices
############################################################

ind=colnames(alldata[[1]])

alldata_cs = lapply(alldata,function(a) scale(t(a))) # Center and scale each wavelength or trait

H_mat= lapply(alldata_cs,function(a) a%*%t(a)/ncol(a)) # Compute covariance matrices

indD = colnames(H_mat[[1]])[which(substr(colnames(H_mat[[1]]),1,4)=="PSMD")]

indF = colnames(H_mat[[1]])[which(substr(colnames(H_mat[[1]]),1,4)=="PSMC")]


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

#############################################################################################################
# Predictions
#############################################################################################################

trait=c("DMY","DMC","Silking","MFU")
Ntrait = length(trait)

pheno2=pheno2[which(pheno2$CodeParentC%in%indF & pheno2$CodeParentD%in%indD),] # Keep the phenotypes of the hybrids whose both parents have phenomic data

H_matD= lapply(H_matD,function(a) a[match(pheno2$CodeParentD,rownames(a)),match(pheno2$CodeParentD,colnames(a))])
H_matF= lapply(H_matF,function(a) a[match(pheno2$CodeParentC,rownames(a)),match(pheno2$CodeParentC,colnames(a))])

kinG_D = kinD[match(pheno2$CodeParentD,rownames(kinD)),match(pheno2$CodeParentD,colnames(kinD))]
kinG_F = kinF[match(pheno2$CodeParentC,rownames(kinF)),match(pheno2$CodeParentC,colnames(kinF))]

#######################################################################################################
# Effect of the number of flights : Scenario CV_related
#######################################################################################################

Nrep=200 # Nb of repetitions

acc <-matrix(NA,Nrep,Ntrait)
accu <- list(acc,acc,acc,acc,acc,acc,acc)
accu2 <- accu


for (rep in 1:Nrep) {
  print(paste("Repetion Nb: ",rep,sep=""))
  
  calib = sample(nrow(pheno2),round(4*nrow(pheno2)/5))  # Calibration set
  valid = seq(1,nrow(pheno2))[-calib]                   # Validation set
  
  for (tr in 1:Ntrait) {
    print(tr)
    
    trait2pred=pheno2[,match(trait[tr],colnames(pheno2))]
    traitNA=trait2pred
    traitNA[valid] =NA
    
    # Data preparation
    traitcal=data.frame(genoF=pheno2$CodeParentC,genoD=pheno2$CodeParentD,pheno=traitNA) 
    traitcal[valid,3]=NA   # we hide the phenotypes of the validation set
    traitcal$genoF=as.factor(traitcal$genoF)
    traitcal$genoF2=as.factor(traitcal$genoF)
    traitcal$genoF3=as.factor(traitcal$genoF)
    traitcal$genoD=as.factor(traitcal$genoD)
    traitcal$genoD2=as.factor(traitcal$genoD)
    traitcal$genoD3=as.factor(traitcal$genoD)
    
    for (fly in 1:7) { 
      
      if (fly<=6) { spl = sample(6,fly) } # Random samples among the UAV acquisition dates
      if (fly==7) { spl = c(2,4,6) } 
      
      H_matD_sel = H_matD[spl]
      H_matD_sel <- Reduce("+", H_matD_sel) / length(H_matD_sel) # Average the corresponding Dent covariance matrices
      H_matF_sel = H_matF[spl]
      H_matF_sel <- Reduce("+", H_matF_sel) / length(H_matF_sel) # Average the corresponding Flint covariance matrices 
      
      # Fit prediction model
      mixmod_sommer <- mmer(pheno ~ 1,
                            random=~vs(genoD,Gu=H_matD_sel)+   vs(genoF,Gu=H_matF_sel),
                            rcov=~units,
                            data=traitcal,tolParInv = 1e-02,verbose = FALSE)
      
      # Phenomic predictions
      
      pred1=mixmod_sommer$U$`u:genoD`$pheno[match(traitcal$genoD,names(mixmod_sommer$U$`u:genoD`$pheno))] + 
        mixmod_sommer$U$`u:genoF`$pheno[match(traitcal$genoF,names(mixmod_sommer$U$`u:genoF`$pheno))] 
      
      # Compute predictive ability
      accu2[[fly]][rep,tr] = cor(pred1[valid],trait2pred[valid])
      
    } # end fly
  }
}

# Compute prediction accuracy by dividing predictive ability by square root of trait heritability
for (fly in 1:7) { 
  for (tr in 1:4) {
  accu2[[fly]][,tr] = accu2[[fly]][,tr]/sqrt(herit[tr])
  }
}

# Average across repetitions
table=sapply(accu2,function(a) apply(a,2,mean))
colnames(table)=seq(1,7)
rownames(table)=trait
round(table,2)

# saveRDS(accu2,file="Figure_Publication/Accuracy_NbFlights.Rds")

accu2 = lapply(accu2, function(a) {
  colnames(a)=trait
  return(a)})


df <- lapply(seq_along(accu2), function(i) {
  as.data.frame(accu2[[i]]) %>%
    setNames(trait) %>%
    mutate(vol = i)
}) %>%
  bind_rows() %>%
  pivot_longer(cols = all_of(trait),
               names_to = "trait",
               values_to = "value")

df$vol[df$vol=="7"]="3opt"

df$vol <- factor(df$vol, levels = c("1","2","3","4","5","6","3opt"))

# Load accuracy of each acquisition date
FbyF = readRDS("Figure_Publication/Accu_FlightByFlight.Rds")
FbyF=sapply(FbyF,function(a) apply(a,2,mean))
FbyF = FbyF[,-7] # Remove the accuracy obtained with all flights together
rownames(FbyF)=trait
round(FbyF,2)

hline_df <- data.frame(  # Horizontal lines corresponding to the 6 individual flights
  trait = rep(c("DMY", "DMC", "Silking", "MFU"), each = 6),
  hline = as.numeric(t(FbyF[,1:6]))
)

df$trait <- factor(df$trait, levels = c("DMY","DMC","Silking","MFU"))

df <- df %>%
  group_by(trait) %>%
  mutate(is_last = ifelse(vol == "3opt", "last", "other")) %>%
  ungroup()

hline_df$trait <- factor(hline_df$trait,
                         levels = levels(df$trait))


# Generate Figure
#########################################

p <- ggplot(df, aes(x = vol, y = value)) +
  geom_boxplot(aes(fill = is_last), color = "black") + 
  facet_wrap(~trait, scales = "free_y") +
  geom_hline(
    data = hline_df,
    aes(yintercept = hline),
    linetype = "dashed",
    color = "red",
  ) +
  scale_fill_manual(values = c("other" = "lightblue", "last" = "indianred")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,face ="bold",size=16,colour = "grey40"),
    legend.position = "none" 
  ) +
  labs(title = "Scenario CV_related",x = "Number of acquisition dates", y = "Prediction accuracy")
p


# ggsave("Figure_Publication/Number_of_Flights_CV_related.jpg",
#        plot = p,
#        width = 12,
#        height = 8,
#        units = "in",
#        dpi = 300,
#        quality = 100)

# #######################################################################################################
# # Effect of the number of flights : Scenario CV_newFamilies
# #######################################################################################################
# 
# famillesF=substr(pheno2$CodeParentC,5,5)
# famillesD=substr(pheno2$CodeParentD,5,5)
# famillesDF = paste(famillesD,famillesF,sep="_")
# famillesDF_unik = unique(famillesDF)
# 
# Nfam=length(famillesDF_unik)
# 
# Nvalid <- rep(NA,Nfam)  # Nb of hybrids in each new family
# for (rep in 1:Nfam) {
#   valid = which(famillesDF==famillesDF_unik[rep])
#   Nvalid[rep]= length(valid)
# }
# 
# # determine the combinations of families for which the number of training genotypes is at least 100 
# # and the number of validation genotypes is at least equal to 30
# Fam1 = combn(famillesDF_unik,1)
# Fam2 = combn(famillesDF_unik,2)
# Fam3 = combn(famillesDF_unik,3)
# Fam4 = combn(famillesDF_unik,4)
# Fam5 = combn(famillesDF_unik,5)
# 
# Fam = list(Fam1,Fam2,Fam3,Fam4,Fam5)
# 
# CombOK = list()
# 
# for (t in 1:length(Fam)) { # Loop on the nb of families combined in the validation set
#   print(t)
#   Ncomb = ncol(Fam[[t]])
#   
#   Size = matrix(NA,2,ncol(Fam[[t]]))
#   rownames(Size) = c("valid","training")
#   
#   for (l in 1:ncol(Fam[[t]])) { # Loop on each combination
#     
#     Size[1,l] =   sum(Nvalid[match(Fam[[t]][,l],famillesDF_unik)]) # Nb of hybrids in the validation set
#     FamD_valid = substr(Fam[[t]][,l],1,1)
#     FamF_valid = substr(Fam[[t]][,l],3,3)    
#     
#     Fam_training = famillesDF_unik[which(substr(famillesDF_unik,1,1)%notin%FamD_valid & substr(famillesDF_unik,3,3)%notin%FamF_valid) ] # Families of the training set
#     
#     Size[2,l] =   sum(Nvalid[match(Fam_training,famillesDF_unik)]) # Nb of hybrids in the training set
#   }
#   
#   CombOK[[t]] = Fam[[t]][,which(Size[1,]>=30 & Size[2,]>=100)] # Combination of families with the necessary conditions
# }
# 
# CombOK[[1]] = as.matrix(CombOK[[1]])
# CombOK[[4]] = as.matrix(CombOK[[4]])
# 
# CombOK = CombOK[which(unlist(lapply(CombOK,ncol))!=0)]
# 
# NcombOK = sum(unlist(lapply(CombOK,ncol)))
# 
# acc <-matrix(NA,NcombOK,Ntrait)
# 
# accu <- list(acc,acc,acc,acc,acc,acc,acc)
# 
# Nvalid=rep(NA,NcombOK)
# Ncalib=rep(NA,NcombOK)
# 
# 
# # Generate the different combinations of UAV acquisition dates
# appariements = list()
# 
# for (nFly in 1:6) {
#   appariements[[nFly]]<- combinations(6, nFly, v = 1:6) 
# }
# 
# appariements[[7]]<- t(as.matrix(c(2,4,6)))
# 
# acc1 <-matrix(NA,NcombOK*nrow(appariements[[1]]),Ntrait)
# acc2 <-matrix(NA,NcombOK*nrow(appariements[[2]]),Ntrait)
# acc3 <-matrix(NA,NcombOK*nrow(appariements[[3]]),Ntrait)
# acc4 <-matrix(NA,NcombOK*nrow(appariements[[4]]),Ntrait)
# acc5 <-matrix(NA,NcombOK*nrow(appariements[[5]]),Ntrait)
# acc6 <-matrix(NA,NcombOK*nrow(appariements[[6]]),Ntrait)
# acc7 <-matrix(NA,NcombOK*nrow(appariements[[7]]),Ntrait)
# 
# accu <- list(acc1,acc2,acc3,acc4,acc5,acc6,acc7)
# 
# 
# for (fly in 1:7) { # Loop on the number of UAV acquisition dates
#   if (fly<=6) {print(paste("Nb of dates: ", fly))}
#   if (fly==7) {print("3 optimized dates")}
#   
#   appariements_fl <- appariements[[fly]] # The different possibilities to combine the x acquisition dates
#   Nb_appariements = nrow(appariements_fl)
#   
#   cpt=0
#   for (app in 1:Nb_appariements) {
#     
#     spl = appariements_fl[app,]
#     H_matD_sel = H_matD[spl]
#     H_matD_sel <- Reduce("+", H_matD_sel) / length(H_matD_sel) # average the corresponding Dent covariance matrices
#     H_matF_sel = H_matF[spl]
#     H_matF_sel <- Reduce("+", H_matF_sel) / length(H_matF_sel)  # average the corresponding Flint covariance matrices
#     
# 
#     for (cb in 1:length(CombOK)) { # loop on the different number of families combined in the validation set
#       print(paste("Nb: ",cb,sep=""))
#       
#       for (l in 1:ncol(CombOK[[cb]])) {   # Loop on the different combinations of families in the validation set
#         cpt =cpt+1  
#         valid = which(famillesDF%in%CombOK[[cb]][,l]) # Predicted genotypes
#         hideD = which(famillesD%in%substr(CombOK[[cb]][,l],1,1))
#         hideF = which(famillesF%in%substr(CombOK[[cb]][,l],3,3))
#         hide = unique(c(hideD,hideF))
#         
#         calib = seq(1,nrow(pheno2))[-hide]
#         
#         Nvalid[cpt]= length(valid)
#         Ncalib[cpt]= length(calib)
#       
#       for (tr in 1:Ntrait) { # Loop on the traits
#         
#         trait2pred=pheno2[,match(trait[tr],colnames(pheno2))]
#         traitNA=trait2pred
#         traitNA[hide] =NA  # hide the phenotypes of the validation set
#         
#         traitcal=data.frame(genoF=pheno2$CodeParentC,genoD=pheno2$CodeParentD,pheno=traitNA) 
#         traitcal$genoF=as.factor(traitcal$genoF)
#         traitcal$genoD=as.factor(traitcal$genoD)
#         
#         # Fit prediction model
#         mixmod_sommer <- mmer(pheno ~ 1,
#                               random=~vs(genoD,Gu=H_matD_sel)+   vs(genoF,Gu=H_matF_sel),
#                               rcov=~units,
#                               data=traitcal,tolParInv = 1e-02,verbose=FALSE)
#         
#         
#         # Phenomic predictions
#         pred1=mixmod_sommer$U$`u:genoD`$pheno[match(traitcal$genoD,names(mixmod_sommer$U$`u:genoD`$pheno))] + 
#           mixmod_sommer$U$`u:genoF`$pheno[match(traitcal$genoF,names(mixmod_sommer$U$`u:genoF`$pheno))] 
#         
#         # Compute predictive ability
#         accu[[fly]][cpt,tr] = cor(pred1[valid],trait2pred[valid])
#       }
#       } 
#     }
#   }
# }
# 
# 
# # Compute prediction accuracy by dividing predictive ability by square root of trait heritability
# for (fly in 1:7) { 
#   for (tr in 1:4) {
#     accu[[fly]][,tr] = accu[[fly]][,tr]/sqrt(herit[tr])
#   }
# }
# 
# table=sapply(accu,function(a) apply(a,2,mean))
# rownames(table)=trait
# round(table,2)
# 
# # write.table(table,file="Figure_Publication/tableAccuracy_newFamDF_NbFlights.csv")
# 
# 
# accu = lapply(accu, function(a) {
#   colnames(a)=trait
#   return(a)})
# 
# df <- lapply(seq_along(accu), function(i) {
#   as.data.frame(accu[[i]]) %>%
#     setNames(trait) %>%
#     mutate(vol = i)
# }) %>%
#   bind_rows() %>%
#   pivot_longer(cols = all_of(trait),
#                names_to = "trait",
#                values_to = "value")
# 
# df$vol[df$vol=="7"]="3opt"
# 
# df$vol <- factor(df$vol, levels = c("1","2","3","4","5","6","3opt"))
# 
# # Prediction accuracies by acquisition dates
# FbyF = read.csv("Figure_Publication/tableAccuracy_newFamDF_ByFlights.csv",sep=" ",check.names=F)
# FbyF = FbyF[,-7] # remove prediction accuracy obtained when all acquisition dates are combined
# 
# hline_df <- data.frame(  # Horizontal lines corresponding to the 6 individual flights
#   trait = rep(trait, each = 6),
#   hline = as.numeric(t(FbyF[,1:6]))
# )
# 
# df$trait = factor(df$trait,levels = c("DMY","DMC","Silking","MFU"))
# 
# df <- df %>%
#   group_by(trait) %>%
#   mutate(is_last = ifelse(vol == "3opt", "last", "other")) %>%
#   ungroup()
# 
# hline_df$trait <- factor(hline_df$trait,
#                          levels = levels(df$trait))
# 
# 
# 
# 
# p <- ggplot(df, aes(x = vol, y = value)) +
#   geom_boxplot(aes(fill = is_last), color = "black") +
#   facet_wrap(~trait, scales = "free_y") +
#   geom_hline(data = hline_df,
#              aes(yintercept = hline),
#              linetype = "dashed",
#              color = "red") +
#   scale_fill_manual(values = c("other" = "lightblue", "last" = "indianred")) +
#   coord_cartesian(ylim = c(-0.5, 1)) +
#   theme_bw(base_size = 14) +
#   theme(
#     strip.background = element_rect(fill = "white"),
#     panel.grid.major = element_line(color = "grey90"),
#     panel.grid.minor = element_blank(),
#     plot.title = element_text(hjust = 0.5,face ="bold",size=16,colour = "grey40"),
#     legend.position = "none" 
#   ) +
#   labs(title = "Scenario CV_newFamilies",x = "Number of acquisition dates", y = "Prediction accuracy")
# 
# p
# 
# # ggsave("Figure_Publication/Number_of_Flights_newFamilies.jpg",
# #        plot = p,
# #        width = 12,
# #        height = 8,
# #        units = "in",
# #        dpi = 300,
# #        quality = 100)
# 



# end