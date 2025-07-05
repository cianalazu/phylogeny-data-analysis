#Load in libraries
library(MCMCglmm)
library(ggplot2)

#Load in data
enddat <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/enddat.csv")
endcoms <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/endcoms.csv")
endcoms_minlemna <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/endcoms_minlemna.csv")
enddat_no_lemna <- enddat[enddat$comdist > 0, ]
exp1_data <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/exp1_data_v2")

#filter exp1_data to only be data associated with round_num = 10
exp1_data_round10 <- exp1_data[exp1_data$round_num == 10, ]


#sums of squares function
ssbyvar <- function(response,category.vec){ 
  means <- tapply(response,category.vec,mean,na.rm=T)
  ssresid <- sum(sapply(sort(unique(category.vec)), function(z) sum( (response[category.vec==z] - means[names(means)==z])^2,na.rm=T )))
  sstot <- sum((response-mean(response,na.rm=T))^2,na.rm=T)
  sst <- (sstot-ssresid)
  return(sst/sstot)
}

#________________________________________________________________________________
#sum of squares - experiment 1 data


#RGR ~ round_num
print(ssbyvar(exp1_data$RGR, exp1_data$round_num))*100
round_variances <- tapply(exp1_data$RGR, exp1_data$round_num, var, na.rm = TRUE)
print(round_variances)*100

plot(as.numeric(names(round_variances)), round_variances,
     type = "b", pch = 19, col = "blue",
     xlab = "Round Number", ylab = "Variance in RGR",
     main = "RGR Variance Across Rounds")

summary(MCMCglmm(RGR ~ round_num, data = exp1_data, nitt = 10000, burnin = 500, thin = 10, verbose=F))
tapply(exp1_data$RGR, exp1_data$round_num, mean, na.rm = TRUE)

plot(round_num_means <- tapply(exp1_data$RGR, exp1_data$round_num, mean, na.rm = TRUE),
     type = "b", pch = 19, col = "blue", xlab = "Round", ylab = "Mean RGR",
     main = "RGR Across Rounds")
abline(lm(RGR ~ round_num, data = exp1_data), col = "red", lty = 2)

#_______________________________________________________________________________
#final_experiment data - sum of squares

#enddat is only data from the last day
print(ssbyvar(enddat_no_lemna$RGR, enddat_no_lemna$stend))
print(ssbyvar(enddat_no_lemna$RGR, enddat_no_lemna$comdist))
print(ssbyvar(enddat_no_lemna$RGR, enddat_no_lemna$percentother))
print(ssbyvar(enddat_no_lemna$RGR, enddat_no_lemna$treatment))
print(ssbyvar(enddat_no_lemna$RGR, enddat_no_lemna$comnum))
print(ssbyvar(enddat_no_lemna$RGR, paste(enddat_no_lemna$comnum,enddat_no_lemna$stend)))



# sum of squares - only derived communities NO LEMNA
print(ssbyvar(endcoms_minlemna$improvement, endcoms_minlemna$comnum))*100
print(ssbyvar(endcoms_minlemna$improvement, endcoms_minlemna$comdist))*100
#print(ssbyvar(endcoms_minlemna$improvement, endcoms_minlemna$treatment))*100
print(ssbyvar(endcoms_minlemna$improvement, endcoms_minlemna$percentother))*100

print(ssbyvar(endcoms_minlemna$improvement, paste(endcoms_minlemna$percentother,endcoms_minlemna$comdist )) )*100

#only derived communities LEMNA INCLUDED
print(ssbyvar(endcoms$improvement, endcoms$percentother))*100
print(ssbyvar(endcoms$improvement, endcoms$log_divergence_time))*100
print(ssbyvar(endcoms$improvement, endcoms$log_divergence_time_sq))*100
print(ssbyvar(endcoms$improvement, endcoms$comdist))*100
print(ssbyvar(endcoms$improvement, paste(endcoms$comdist,endcoms$percentother)))*100

print(ssbyvar(endcoms$improvement, endcoms$comnum))*100

#___________________________________________________________________________________

# plots

trtonly <- endcoms[,c("comnum","treatment")] #[]#subset to treatment columns
datonly <- endcoms[,c("RGR","improvement")]#subset to measured responses
#for(i in 1000){
  randorder <- sample(1:nrow(endcoms),nrow(endcoms),replace=F)
  randdat <- datonly[randorder,]
  randComExpl <- ssbyvar(datonly$improvement, trtonly$comnum)
#}

final_experiment_data <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/final_experiment_data.csv")
enddat_no_lemna <- merge(enddat_no_lemna, divergence_times, by = "comdist")





#___________________________________________________________________________________
#________________________ANOVAs___________________________________________________
#more work with anovas, can look at evolution paper
anova(lm(improvement~treatment, data=endcoms_minlemna))

anova(lm(improvement~comnum, data=endcoms_minlemna))

anova(lm(RGR~treatment, data=enddat_no_lemna,na.action = na.omit))

anova(lm(RGR~comnum, data=enddat_no_lemna,na.action = na.omit))

anova(lm(RGR~stend, data=enddat_no_lemna,na.action = na.omit))


anova(lm(improvement~percentother, data=endcoms_minlemna))

anova(lm(RGR~comdist, data=enddat_no_lemna,na.action = na.omit))

anova(lm(RGR~divergence_times_ma, data=enddat_no_lemna,na.action = na.omit))


