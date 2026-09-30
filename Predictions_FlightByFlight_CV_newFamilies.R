
#############################################################################################################
#############################################################################################################
# Phenomic predictions of the different drone flights, taken individually
#############################################################################################################
#############################################################################################################

rm(list=ls())

library(data.table)
library(rrBLUP)
library(sommer)
library(tidyverse)
library(gtools)
library(matrixcalc)
library(patchwork)

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

# Secondary traits from drone images 

drone = readRDS("Data_Publications/adjmeans_drone_per_flight_maize.Rds")

for (i in 1:length(drone)) {rownames(drone[[i]])[which(rownames(drone[[i]])=="psmc1287_f9024")] = "psmc1287"}

rnam = rownames(drone[[1]])
drone_all = matrix(unlist(drone),nrow=length(rnam))
rownames(drone_all) = toupper(rnam)

herit_drone = readRDS("Data_Publications/Heritabilities_secondarytraits_maize.Rds")
herit_drone2 = unlist(herit_drone)
drone_all = drone_all[,herit_drone2>=0.4]

for (y in 1:length(drone)) {drone[[y]] = drone[[y]][,which(herit_drone[[y]]>=0.4)]}

drone[[7]] = drone_all

for (y in 1:length(drone)) {rownames(drone[[y]]) = toupper(rownames(drone[[y]])) 
drone[[y]] = t(drone[[y]])
}

alldata = drone

rm(drone_all)
rm(drone)

commun = read.csv("Data_Publications/Geno_Alldata_maize.csv",sep=";") # names of the parental lines with phenomic data
commun = commun[,2]

alldata = lapply(alldata,function(a) a[,match(commun,colnames(a))])

#############################################################################################################
#############################################################################################################
# Compute kinships and phenomic covariance matrices of the Dent and Flint parental lines
#############################################################################################################
#############################################################################################################


# Phenomic covariance matrices
############################################################

alldata_cs = lapply(alldata,function(a) scale(t(a))) # Center and scale each secondary traits (from the drone images)

# Covariance matrices (one per flight, and one for all flights together)
H_mat= lapply(alldata_cs,function(a) a%*%t(a)/ncol(a)) # Compute covariance matrices

# Dent parental lines
indD = colnames(H_mat[[1]])[which(substr(colnames(H_mat[[1]]),1,4)=="PSMD")]

# Flint parental lines
indF = colnames(H_mat[[1]])[which(substr(colnames(H_mat[[1]]),1,4)=="PSMC")]

# Covariance matrix between the Dent parents
H_matD = lapply(H_mat,function(a) a[match(indD,rownames(a)),match(indD,colnames(a))]) # Covariance matrix between dent lines
# Covariance matrix between the Flint parents
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


#############################################################################################################
#############################################################################################################
# Predictions 
#############################################################################################################
#############################################################################################################

trait=c("DMY","DMC","Silking","MFU")
Ntrait = length(trait)

pheno2=pheno2[which(pheno2$CodeParentC%in%indF & pheno2$CodeParentD%in%indD),] # Keep the phenotypes of the hybrids whose both parents have phenomic data

H_matD= lapply(H_matD,function(a) a[match(pheno2$CodeParentD,rownames(a)),match(pheno2$CodeParentD,colnames(a))])
H_matF= lapply(H_matF,function(a) a[match(pheno2$CodeParentC,rownames(a)),match(pheno2$CodeParentC,colnames(a))])

kinG_D = kinD[match(pheno2$CodeParentD,rownames(kinD)),match(pheno2$CodeParentD,colnames(kinD))]
kinG_F = kinF[match(pheno2$CodeParentC,rownames(kinF)),match(pheno2$CodeParentC,colnames(kinF))]

# Names of the dates
modname=c("06/28","07/11","07/24","08/13","09/05","10/02","All Flights")

#######################################################################################################
# Predictions in Scenario 2 : Predict completely new individuals: the parents of the predicted hybrids are from different
#         families than the parents of the hybrids in the training set (CV_newFamilies)
#######################################################################################################

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

CombOK = CombOK[which(unlist(lapply(CombOK,ncol))!=0)]

NcombOK = sum(unlist(lapply(CombOK,ncol)))

acc <-matrix(NA,NcombOK,Ntrait)

accu <- list(acc,acc,acc,acc,acc,acc,acc)

Nvalid=rep(NA,NcombOK)
Ncalib=rep(NA,NcombOK)


for (fly in 1:7) { # Loop on the number of flights (+1 for all flights together)
  print(paste("Fly Nb: ", fly))
  
  H_matD_sel = H_matD[[fly]]
  H_matF_sel = H_matF[[fly]]
  
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
      
      for (tr in 1:Ntrait) {
        
        trait2pred=pheno2[,match(trait[tr],colnames(pheno2))]
        traitNA=trait2pred
        traitNA[hide] =NA
        
        traitcal=data.frame(genoF=pheno2$CodeParentC,genoD=pheno2$CodeParentD,pheno=traitNA) 
        traitcal[hide,3]=NA   # Hide validation set
        traitcal$genoF=as.factor(traitcal$genoF)
        traitcal$genoD=as.factor(traitcal$genoD)
        
        # REML model 1
        mixmod_sommer <- mmer(pheno ~ 1,
                              random=~vs(genoD,Gu=H_matD_sel)+   vs(genoF,Gu=H_matF_sel),
                              rcov=~units,
                              data=traitcal,tolParInv = 1e-02,verbose=FALSE)
        
        
        # Les predictions
        pred1=mixmod_sommer$U$`u:genoD`$pheno[match(traitcal$genoD,names(mixmod_sommer$U$`u:genoD`$pheno))] + 
          mixmod_sommer$U$`u:genoF`$pheno[match(traitcal$genoF,names(mixmod_sommer$U$`u:genoF`$pheno))] 
        accu[[fly]][cpt,tr] = cor(pred1[valid],trait2pred[valid])
        if (var(pred1[valid])==0) {accu[[fly]][cpt,tr] = 0}
        
      } 
    }
  }
}


table=sapply(accu,function(a) apply(a,2,mean))
colnames(table)=modname
rownames(table)=trait
round(table,2)

# Compute accuracy by dividing predictive ability by square root of trait heritability

for (trt in 1:Ntrait) { 
  table[1,] = table[1,]/sqrt(herit[1])
  table[2,] = table[2,]/sqrt(herit[2])
  table[3,] = table[3,]/sqrt(herit[3])
  table[4,] = table[4,]/sqrt(herit[4])
}

#write.table(table,file="Figure_Publication/tableAccuracy_newFamilies_ByFlights.csv")

df <- as.data.frame(table[,1:6])
df$trait <- rownames(df)

df_long <- df %>%
  pivot_longer(
    cols = -trait,
    names_to = "date",
    values_to = "valeur"
  )

df_long$trait <- factor(df_long$trait,
                        levels = c("DMY", "DMC", "Silking", "MFU"))

last_date <- "10/02"

temp <- df_long[df_long$date==last_date,]
df_last <- mutate(temp, numero = 1:n())

df_last$trait <- factor(df_last$trait,
                        levels = c("DMY", "DMC", "Silking", "MFU"))

df_long$date <- as.Date(df_long$date, format = "%m/%d")
df_last$date <- as.Date(df_last$date, format = "%m/%d")
df_last$valeur = table[,7]


## Generate Figure
###########################

cols <- c(
  "DMY"      = "#0072B2",
  "DMC"      = "#F28E2B",
  "Silking"  = "#009E73",
  "MFU"      = "#CC79A7"
)

p <-
  ggplot(df_long,
         aes(date,
             valeur,
             colour = trait,
             group = trait)) +
  
  ## horizontal line
  geom_hline(
    yintercept = 0,
    colour = "grey55",
    linewidth = 0.45,
    linetype = "22"
  ) +
  
  ## curves
  geom_line(
    linewidth = 1.2,
    lineend = "round"
  ) +
  
  ##  Final points (accuracy with all flights taken together)
  geom_point(
    data = df_last,
    aes(fill = trait),
    shape = 21,
    size = 4.5,
    stroke = 0.7,
    colour = "black",
    show.legend = FALSE
  ) +
  
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks",
    expand = expansion(mult = c(0.02,0.02))
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.12,0.05))
  ) +
  
  labs(
    x = "UAV acquisition dates",
    y = "Prediction accuracy"
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    axis.text = element_text(
      size = 14,
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    legend.position = "top",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 14
    ),
    
    panel.border = element_blank(),
    
    plot.margin = margin(
      10,
      10,
      15,
      10
    )
  ) +
  geom_rect(
    inherit.aes = FALSE,
    xmin = as.Date("07/24",format="%m/%d"),
    xmax = as.Date("07/31",format="%m/%d"),
    ymin = -Inf,
    ymax = Inf,
    fill = "#E69F00",
    alpha = 0.01
  ) +
  
  geom_rect(
    inherit.aes = FALSE,
    xmin = as.Date("08/11",format="%m/%d"),
    xmax = as.Date("08/27",format="%m/%d"),
    ymin = -Inf,
    ymax = Inf,
    fill = "#D55E00",
    alpha = 0.01
  )

p


# saveRDS(p,file="Figure_Publication/AccuperFlight_newFamilies.Rds")

# Combine the two figures :
################################################################################

p1 = readRDS("Figure_Publication/Accuracy_per_flight2.Rds")
p2 = readRDS("Figure_Publication/AccuperFlight_newFamilies.Rds")

p1 <- p1 +
  labs(title = "Scenario CV_related") +
  theme(
    plot.title = element_text(
      colour = "grey40",
      size = 16,
      face = "bold",
      hjust = 0.5
    )
  )

p2 <- p2 +
  labs(title = "Scenario CV_newFamilies") +
  theme(
    plot.title = element_text(
      colour = "grey40",
      size = 16,
      face = "bold",
      hjust = 0.5
    )
  )

p2 <- p2 + labs(y = NULL)+ theme(legend.position = "none")

p3 = (p1 | p2) +
  plot_annotation(
    tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20,face="bold")
  )

p3

# ggsave("Figure_Publication/Accuracy_perFlight_combined.jpg",
#        plot = p3,
#        width = 12,
#        height = 8,
#        units = "in",
#        dpi = 300,
#        quality = 100)

# end