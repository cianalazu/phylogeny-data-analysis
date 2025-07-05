
library(tidyr)
library(ggplot2)
library(dplyr)
library(MCMCglmm)

endcoms <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/endcoms.csv")


# Assuming 'endcoms' and 'divergence_times' both have a 'comdist' column
endcoms <- merge(endcoms, divergence_times, by = "comdist", all = TRUE)
endcoms$divergence_times_ma[is.na(endcoms$divergence_times_ma)] <- 0



sortcommeans$disttoname <- c("Lemna", "Wolffia", "Spirodela", "Elodea", "Ceratophyllum", "Nymphaea", "Riccia")[sortcommeans$comdist + 1]


#________________________________________________________________________________________________________

#____________________creation of endcoms_minlemna____________________________________________________________________

# Filter out rows where percentother is 0 or NA, and remove rows with missing values in relevant columns
endcoms_minlemna <- subset(endcoms, percentother != 0 & !is.na(percentother))

# Keep only complete cases for the variables of interest
vars_to_keep <- c("log_divergence_time_sq", "RGR", "log_divergence_time", "percentother", "improvement", "comdist", "divergence_times_ma")
endcoms_minlemna <- endcoms_minlemna[complete.cases(endcoms_minlemna[, vars_to_keep]), ]

#save endcoms_minlemna to csv
write.csv(endcoms_minlemna, "C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/endcoms_minlemna.csv", row.names = FALSE)

# _______________PLOT IM LOOKING FOR!_______________________________
#calculates the average Relative Growth Rate (RGR) for each unique group in comdist, and stores the result in cdmeans
cdmeans <- tapply(endcoms$RGR, endcoms$comdist, mean)
cdses <- tapply(endcoms$RGR, endcoms$comdist, se)


png("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/improvement_meansEndComsRGR.png", width = 7, height = 3, units = "in", res = 1200, bg = "white")
plot(improvement_means~c(0:6), pch=1,cex=2, xlab="",ylab="", ylim=range(c(improvement_means+improvement_ses,improvement_means-improvement_ses)))
arrows(c(0:6),y0=improvement_means-improvement_ses,y1=improvement_means+improvement_ses,length=0,lwd=2)
dev.off()

png("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/comdistEndComsRGR.png", width = 7, height = 3, units = "in", res = 1200, bg = "white")
plot(cdmeans~c(0:6), pch=1,cex=2, xlab="",ylab="", ylim=range(c(cdmeans+cdses,cdmeans-cdses)))
arrows(c(0:6),y0=cdmeans-cdses,y1=cdmeans+cdses,length=0,lwd=2)
dev.off()

modcomdistfig <- MCMCglmm(cdmeans ~ comdist -1 , data = endcoms, nitt = 10000, burnin = 5000, thin = 10,verbose=F)


MCMCglmm(Distance ~ Group , data = plot_within_df, nitt = 10000, burnin = 5000, thin = 10,verbose=F)
summary(modcomdistfig)




#model of divergence times w/ end coms INCLUDES LEMNA
modelmc <- MCMCglmm(RGR ~  divergence_times_ma + divergence_time_sq, data = endcoms, nitt = 500000, burnin = 5000, thin = 50,verbose=F)
summary(modelmc)

solmc <- modelmc$Sol

# Create a sequence of values for comdist and its square for predictions
divlin_seq <- seq(min(endcoms$divergence_times_ma), max(endcoms$divergence_times_ma), length.out = 1000)
divsq_seq <- divlin_seq^2

mean_pred <- sapply(1:length(divlin_seq), function(z) mean(solmc[,1] + solmc[,2] * divlin_seq[z] + solmc[,3] *divsq_seq[z]) )
hpdi_pred <- sapply(1:length(divlin_seq), function(z) HPDinterval(solmc[,1] + solmc[,2] * divlin_seq[z] + solmc[,3] *divsq_seq[z], 0.95) )

# Print the summary of the model
summary(modelmc)

png("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/comdistEndComsRGR_pred.png", width = 7, height = 3, units = "in", res = 1200, bg = "white")
plot(cdmeans~sort(unique(endcoms$divergence_times_ma)), pch=1,cex=2, xlab="",ylab="", ylim=range(c(cdmeans+cdses,cdmeans-cdses)))
arrows(sort(unique(endcoms$divergence_times_ma)),y0=cdmeans-cdses,y1=cdmeans+cdses,length=0,lwd=2)
lines(divlin_seq, mean_pred, col="red", lwd=2)
dev.off()





#________________________________________________________________________________________
#_______model with derived communities - RGR NO LEMNA_____________________________________


rgrpercentotherend <-  MCMCglmm(RGR ~ percentother, data = endcoms_minlemna, verbose=F, nitt = 100000, burnin = 5000, thin = 10)

summary(rgrpercentotherend)

rgrlogdivergenceend <- MCMCglmm(RGR ~ log_divergence_time, data = endcoms_minlemna, verbose=F, nitt = 100000, burnin = 5000, thin = 10)
summary(rgrlogdivergenceend)

interaction_model <- MCMCglmm(RGR ~ percentother + log_divergence_time + log_divergence_time_sq + percentother*log_divergence_time_sq, data = endcoms_minlemna, verbose=F, nitt = 100000, burnin = 5000, thin = 10)
summary(interaction_model)


#________________________________________________________________________________________

modelmclog <- MCMCglmm(RGR ~  log_divergence_time + log_divergence_time_sq, data = endcoms_minlemna,  nitt = 100000, burnin = 5000, thin = 10,verbose=F)
solmclog <- modelmclog$Sol
# Print the summary of the model
summary(modelmclog)

ldivlin_seq <- seq(min(endcoms_minlemna$log_divergence_time), max(endcoms_minlemna$log_divergence_time), length.out = 1000)
ldivsq_seq <- ldivlin_seq^2

lmean_pred <- sapply(1:length(ldivlin_seq), function(z) mean(solmclog[,1] + solmclog[,2] * ldivlin_seq[z] + solmclog[,3] *ldivsq_seq[z]) )

#plot i was looking for but with a line fitted to the log divergence time
png("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/endcoms_log_divergence_time_nolemna.png", width = 7, height = 5, units = "in", res = 1200, bg = "white")
plot(cdmeans[-1]~sort(unique(endcoms_minlemna$log_divergence_time)), pch=1,cex=2, xlab="log Divergence Time",ylab="Average RGR of Derived Communities", ylim=range(c(cdmeans+cdses,cdmeans-cdses)))
arrows(sort(unique(endcoms_minlemna$log_divergence_time)),y0=cdmeans[-1]-cdses[-1],y1=cdmeans[-1]+cdses[-1],length=0,lwd=2)
lines(ldivlin_seq, lmean_pred, col="red", lwd=2)
dev.off()


#________________________________________________________________________________________
#_________________________IMPROVEMENT OF DERIVED COMMUNITIES FIT TO MODEL_______________________

improvement_means <- tapply(endcoms_minlemna$improvement, endcoms_minlemna$comdist, mean)
improvement_ses <- tapply(endcoms_minlemna$improvement, endcoms_minlemna$comdist, se)

modelmclogImpr <- MCMCglmm(improvement ~  log_divergence_time + log_divergence_time_sq, data = endcoms_minlemna,  nitt = 100000, burnin = 5000, thin = 10,verbose=F)
solmclogImpr <- modelmclogImpr$Sol

# Print the summary of the model
summary(modelmclogImpr)

ldivlin_seq <- seq(min(endcoms_minlemna$log_divergence_time), max(endcoms_minlemna$log_divergence_time), length.out = 1000)
ldivsq_seq <- ldivlin_seq^2

lmean_pred_Impr <- sapply(1:length(ldivlin_seq),
                     function(z) mean(solmclogImpr[,1] + solmclogImpr[,2] * ldivlin_seq[z] + solmclogImpr[,3] *ldivsq_seq[z]) )

#save plot
png("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/improvement_meansEndComs_divergence.png", width = 7, height = 3, units = "in", res = 1200, bg = "white")
plot(improvement_means~sort(unique(endcoms_minlemna$log_divergence_time)), 
     pch=1,cex=2, xlab="log Divergence Time",ylab="Average Improvement of Derived Communities", 
     ylim=range(c(improvement_means+improvement_ses,improvement_means-improvement_ses)))
arrows(sort(unique(endcoms_minlemna$log_divergence_time)),
       y0=improvement_means-improvement_ses,y1=improvement_means+improvement_ses,length=0,lwd=2)
lines(ldivlin_seq, lmean_pred_Impr, col="red", lwd=2)


