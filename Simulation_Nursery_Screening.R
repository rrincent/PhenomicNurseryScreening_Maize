##################################################################################
# Script to estimate the risk of losing best genotypes when screening nursery 
# with genomic or phenomic prediction
##################################################################################

n_cores <- 50

options(bitmapType = "cairo")

library(MASS)
library(parallel)
library(ggplot2)
library(scico)
library(dplyr)
library(patchwork)
library(grid)


# The different level of prediction accuracy (= correlation between TBV and predictions)
corl_seq = seq(0.1,0.9,by=0.02)

# Proportion of genotypes eliminated based on predictions
elim_rate = seq(0.1,0.9,by=0.02)

# Size of the nursery
Nind = 2000

# Number of repetitions of the process
Nrep=500000


Nbest_1 = 1 # Best genotype
Nbest_2 = 3 # Three best genotypes
Nbest_3 = 30 # 50 best genotypes


ComputeRisk <- function(i, Nbest_1,nbest_2,N_best3,corl_seq,elim_rate,Nrep) {

  saver_1 = matrix(0,length(corl_seq ),length(elim_rate))
  saver_2 = matrix(0,length(corl_seq ),length(elim_rate))
  propbest = matrix(0,length(corl_seq ),length(elim_rate))
  
for (crl in 1:length(corl_seq)) {
  print(crl)
   
    sigma = matrix(1,2,2)
    sigma[1,2] <- sigma[2,1] <-corl_seq[crl]
    # Generate predictions and TBV based on a given level of accuracy
    test = mvrnorm(n=Nind,mu=c(0,0),Sigma = sigma)
    colnames(test) =c("Pred","TBV")
  
  for (elim in 1:length(elim_rate)) {
    
    eliminate = elim_rate[elim]

      # Apply screening
      remove = order(test[,1],decreasing = FALSE)[1:(eliminate*Nind)]
      
      # Test is the Nbest_1 genotypes are eliminated or not
      best = order(test[,2],decreasing = TRUE)[1:Nbest_1]
      if (sum(best%in%remove) == Nbest_1) {saver_1[crl,elim] = saver_1[crl,elim] + 1}
      
      # Test if at least one of the Nbest_2 genotypes is eliminated or not
      best = order(test[,2],decreasing = TRUE)[1:Nbest_2]
      if (sum(best%in%remove) >=1) {saver_2[crl,elim] = saver_2[crl,elim] + 1}
      
      # Compute the number of Nbest_3 genotypes eliminated
      best = order(test[,2],decreasing = TRUE)[1:Nbest_3]
      propbest[crl,elim] = propbest[crl,elim] + sum(best%in%remove)
      
    }}

gc()
list(saver_1 = saver_1,saver_2 = saver_2,propbest = propbest)
}


results <- mclapply(seq(1,Nrep),
                       ComputeRisk,
                    Nbest_1=Nbest_1,nbest_2=Nbest_2,N_best3=Nbest_3,
                    corl_seq=corl_seq,elim_rate=elim_rate,Nrep=Nrep,
                    mc.cores = n_cores)

adjmeans1 <- matrix(0, length(corl_seq), length(elim_rate))
rownames(adjmeans1) <- corl_seq
colnames(adjmeans1) <- elim_rate
adjmeans3 <- adjmeans2 <- adjmeans1 

for (i in 1:Nrep) {
  adjmeans1 <-   adjmeans1 + results[[i]]$saver_1
  adjmeans2 <- adjmeans2 + results[[i]]$saver_2
  adjmeans3 <- adjmeans3 + results[[i]]$propbest
}

adjmeans1 <-   adjmeans1/Nrep
adjmeans2 <- adjmeans2/Nrep
adjmeans3 <- adjmeans3/Nrep

saveRDS(adjmeans1,  "Risk_Best1.Rds")
saveRDS(adjmeans2,  "Risk_Best2.Rds")
saveRDS(adjmeans3,  "Risk_PropBest.Rds")


# #############################################################################
# # Produce Figure from the output
# #############################################################################
# #############################################################################
# 
# Risk_Best1 = readRDS("Risk_Best1.Rds")
# Risk_Best2 = readRDS("Risk_Best2.Rds")
# Risk_PropBest = readRDS("Risk_PropBest.Rds")
# 
# Risk_Best1 = Risk_Best1*100
# Risk_Best2 = Risk_Best2*100
# Risk_PropBest = Risk_PropBest*100
# 
# # The different level of prediction accuracy (= correlation between TBV and predictions)
# corl_seq = seq(0.1,0.9,by=0.01)
# 
# # Proportion of genotypes eliminated based on predictions
# elim_rate = seq(0.1,0.9,by=0.01)
# 
# 
# # Draw Figure Risk 1
# #################################################################################
# 
# df = data.frame(x = rep(corl_seq,length(elim_rate)),y = rep(elim_rate,each=length(corl_seq)), z=as.numeric(Risk_Best1))
# 
# df$z_class <- cut(
#   df$z,
#   breaks = c(0,2,5,10,25,50,100),
#   include.lowest = TRUE,
#   labels = c(
#     "0–2",
#     "2–5",
#     "5–10",
#     "10–25",
#     "25–50",
#     "50–100"
#   )
# )
# 
# p1 = ggplot(df, aes(x,y)) +
#   geom_raster(
#     aes(fill=z_class),
#     interpolate=TRUE
#   ) +
#   scale_fill_manual(
#     values = rev(scico(6, palette="batlow")),
#     name = "Probability (%)"
#   ) +
#   geom_contour(
#     aes(z=z),
#     breaks = c(0,2,5,10,25,50,100),
#     colour="white",
#     linewidth=1
#   ) +
#   
#   scale_x_continuous(
#     breaks = seq(min(df$x), max(df$x), by = 0.1)
#   ) +
#   scale_y_continuous(
#     breaks = seq(min(df$y), max(df$y), by = 0.1)
#   ) +
#   coord_fixed() +
#   theme_classic() +
#   labs(
#     x="",
#     y="Elimination Rate",
#     title = "Risk of losing the best genotype",
#     fill = "Probability (%)"
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     legend.position = "none",
#     legend.title = element_text(size=10),
#     plot.title = element_text(
#       face = "bold",
#       size = 12,
#       hjust = 0.5,
#       colour = "grey30"
#     )
#   )
# 
# p1
# 
# # Draw Figure Risk 2
# #################################################################################
# df = data.frame(x = rep(corl_seq,length(elim_rate)),y = rep(elim_rate,each=length(corl_seq)), z=as.numeric(Risk_Best2))
# 
# 
# # Draw Figure Risk 3
# #################################################################################
# df = data.frame(x = rep(corl_seq,length(elim_rate)),y = rep(elim_rate,each=length(corl_seq)), z=as.numeric(Risk_PropBest))
# 
# df$z_class <- cut(
#   df$z,
#   breaks = c(0,2,5,10,25,50,100),
#   include.lowest = TRUE,
#   labels = c(
#     "0–2 %",
#     "2–5 %",
#     "5–10 %",
#     "10–25 %",
#     "25–50 %",
#     "50–100 %"
#   )
# )
# 
# p3 = ggplot(df, aes(x,y)) +
#   geom_raster(
#     aes(fill=z_class),
#     interpolate=TRUE
#   ) +
#   scale_fill_manual(
#     values = rev(scico(6, palette="batlow")),
#     name = "Probability (%)"
#   ) +
#   geom_contour(
#     aes(z=z),
#     breaks = c(0,2,5,10,25,50,100),
#     colour="white",
#     linewidth=1
#   ) +
#   
#   scale_x_continuous(
#     breaks = seq(min(df$x), max(df$x), by = 0.1)
#   ) +
#   scale_y_continuous(
#     breaks = seq(min(df$y), max(df$y), by = 0.1)
#   ) +
#   coord_fixed() +
#   theme_classic() +
#   labs(
#     x="",
#     y="",
#     title = "% discarded from the top 30",
#     fill = "Probability (%)"
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     legend.position = "right",
#     legend.title = element_text(size=10),
#     plot.title = element_text(
#       face = "bold",
#       size = 12,
#       hjust = 0.5,
#       colour = "grey30"
#     )
#   )
# 
# p3
# 
# # Assemble the figures
# 
# library(patchwork)
# 
# p1 <- p1 +
#   labs(x = NULL)
# 
# p3 <- p3 +
#   labs(x = NULL,y=NULL)
# 
# # Ajouter uniquement A/B
# p4 <- (p1 | p3) +
#   plot_annotation(
#     tag_levels = "A"
#   ) &
#   theme(
#     plot.tag = element_text(
#       face = "bold",
#       size = 16
#     ),
#     plot.tag.position = c(0.02, 1)
#   )
# 
# p5 <- p4 +
#   plot_annotation(
#     caption = "Prediction Accuracy"
#   ) &
#   theme(
#     plot.caption = element_text(
#       size = 14,
#       hjust = 0.45
#     )
#   )
# 
# p5
# 
# 
# 
# ggsave("FigureRisk.jpg",
#        plot = p5,
#        width = 12,
#        height = 8,
#        units = "in",
#        dpi = 300,
#        quality = 100)
# 
# 
