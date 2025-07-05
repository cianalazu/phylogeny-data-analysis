
library(ggplot2)
library(dplyr)
library(gridExtra)
library(broom)
library(MCMCglmm)
library(cowplot)

setwd("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1")

# Replace the path with your actual file path
# data <- read.csv("path/to/your/dataset.csv")
exp1_data <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/exp1_data.csv")
final_experiment_data <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/final_experiment_data.csv")
# have to load manually: syn_com_placement_details <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/syn_com_placement_details.csv", header= "True")
exp1_avg_of_reps_per_round <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/exp1_avg_of_reps_per_round.csv")
divergence_times <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/divergence_times.csv")

# summary( MCMCglmm(RGR~other_host_source, data=exp1_data, verbose=T, nitt= 10000, burnin = 500, thin = 10))
# 
# summary( MCMCglmm(RGR~round_num, data=exp1_data, verbose=T, nitt= 10000, burnin = 500, thin = 10))

#_________________________________________________________________________________________

se <- function(x) {sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)])) }

final_experiment_data$platenum <- gsub("[^[:digit:]]", "",final_experiment_data$treatment)
final_experiment_data$stend <- gsub("[[:digit:]]", "",final_experiment_data$treatment)

final_experiment_data$comnum <- final_experiment_data$platenum
final_experiment_data$comnum[final_experiment_data$stend=="End"] <- sapply(
  final_experiment_data$platenum[final_experiment_data$stend=="End"], 
  function(z) syn_com_placement_details$Community.num[syn_com_placement_details$Cycle.10.Plate.num==z])

comdist<- tapply(exp1_data$distance_from_lemna,exp1_data$syn_com,mean)
perother <- tapply(exp1_data$percent_other,exp1_data$syn_com,mean)

final_experiment_data$comdist <- comdist[final_experiment_data$comnum]
final_experiment_data$percentother <- perother[final_experiment_data$comnum]

#add in divergence_times column into final_experiment_data based on shared column comdist in divergence_times df

# Assuming both dataframes have a column named 'comdist'
final_experiment_data <- final_experiment_data %>%
  left_join(divergence_times, by = "comdist")






#startmeans is the average pixel area of each ancestral community (avg of day 1 and day 10)
#endmeans is the average pixel area of each derived community (avg of day 1 and day 10)
startmeans <- tapply(final_experiment_data$area[final_experiment_data$stend=="Start"], final_experiment_data$comnum[final_experiment_data$stend=="Start"],mean) 
endmeans <- tapply(final_experiment_data$area[final_experiment_data$stend=="End"], paste(final_experiment_data$comnum,final_experiment_data$platenum)[final_experiment_data$stend=="End"],mean) 


enddatall <- final_experiment_data[final_experiment_data$day==10,]
enddatall$startarea <- sapply(1:nrow(enddatall), function(z) final_experiment_data$area[
      which(final_experiment_data$day==1 & 
              final_experiment_data$well==enddatall$well[z] &
              final_experiment_data$plate==enddatall$plate[z]
            )])
enddat<- enddatall[enddatall$startarea>0,]

#RGR per day calculation
#((exp1_data$end_area - exp1_data$start_area) / exp1_data$start_area)/exp1_data$daysinround
enddat$RGR <- ((enddat$area-enddat$startarea) / enddat$startarea)/10


# summary( MCMCglmm(RGR~stend, data=enddat,verbose=T, nitt = 100000, burnin = 5000, thin = 10 ))


#save enddat as csv
write.csv(enddat, "C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/enddat.csv", row.names = FALSE)


#startRGR is the RGR of the ancestral communities
 endcoms <- enddat[enddat$stend=="End",]
 
 endcoms$startRGR <- sapply(1:nrow(endcoms), function(z)  mean(enddat$RGR[which(enddat$comnum == endcoms$comnum[z] & enddat$stend=="Start" )])  )
 
 startcoms <- enddat[enddat$stend=="Start",]
 allmeans <- c(tapply(endcoms$RGR,endcoms$treatment,mean,na.rm=T),tapply(startcoms$RGR,startcoms$comnum,mean,na.rm=T))
 allses <- c(tapply(endcoms$RGR,endcoms$treatment,se),tapply(startcoms$RGR,startcoms$comnum,se))
 startendv <- c(rep("E",times=75),rep("S",times=25))
 comdistmn <- c(tapply(endcoms$comdist,endcoms$treatment,mean,na.rm=T),tapply(startcoms$comdist,startcoms$comnum,mean,na.rm=T))
 perothmn <- c(tapply(endcoms$percentother,endcoms$treatment,mean,na.rm=T),tapply(startcoms$percentother,startcoms$comnum,mean,na.rm=T))
 sortcommeans <- data.frame(allmeans = allmeans,allses=allses,stend= startendv,comdist=comdistmn,percentother=perothmn)[order(allmeans),]
 sortcommeans$disttoname <- c("Lemna","Wolffia","Spirodela","Elodea","Ceratophyllum","Nymphaea","Riccia")[sortcommeans$comdist+1]
 
 
 # Create the log-transformed divergence_time
 endcoms$log_divergence_time <- log(endcoms$divergence_times_ma+1)#add 1 so as to not take log of 0
 endcoms$log_divergence_time_sq <- log(endcoms$divergence_times_ma+1)^2 #add 1 so as to not take log of 0
 endcoms$divergence_time_sq <- endcoms$divergence_times_ma^2

  #making column with improvement score
 endcoms$improvement <- (endcoms$RGR - endcoms$startRGR) / abs(endcoms$startRGR)
 #save endcoms
 write.csv(endcoms, "C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/endcoms.csv", row.names = FALSE)
 

 
 #average and standard error of improvement score for each percent_other
 
avgimprovement_df <- endcoms |> 
   filter(percentother %in% c(0, 50, 75)) |> 
   group_by(percentother) |> 
   summarise(
     mean_improvement = mean(improvement, na.rm = TRUE),
     se_improvement = sd(improvement, na.rm = TRUE) / sqrt(n())
   )
 
 # Step 2: Plot with SE error bars
improvementscorefinal <- ggplot(avgimprovement_df, aes(x = factor(percentother), y = mean_improvement)) +
   geom_point(size = 4, shape = 21, fill = "black", color = "black") +  # circle as mean
   geom_errorbar(aes(ymin = mean_improvement - se_improvement,
                     ymax = mean_improvement + se_improvement),
                 width = 0) +  # line through the circle = error bar
   labs(
     x = "% Initial microbiome from other",
     y = "Derived communities improvement",
     title = ""
   ) +
   theme_cowplot()
 #save figure
ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/avg_improvement_by_percent_other.png", width = 5, height = 4, dpi = 300, bg = "white")
 
#t-test between percentother and improvement

 #________________________________________________________________________________________
 #red and black plot showing average RGR in community, with error bars
 png("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/AverageRGRinComm.png",units="in",height=4,width=15, res=1200)
 
 par(mar = c(7, 4, 0, 1))
 
 plot(sortcommeans$allmeans ~ c(1:100),
      pch = c(17, 18, 16)[as.numeric(as.factor(sortcommeans$percentother))],
      col = as.numeric(as.factor(sortcommeans$stend)),
      ylim = range(c(sortcommeans$allmeans - sortcommeans$allses, 
                     sortcommeans$allmeans + sortcommeans$allses)),
      xaxt = "n", ylab = "Average RGR in Community", xlab = ""
 )
 
 arrows(c(1:100), sortcommeans$allmeans - sortcommeans$allses, 
        y1 = sortcommeans$allmeans + sortcommeans$allses, 
        length = 0, lwd = 1, 
        col = as.numeric(as.factor(sortcommeans$stend)))
 
 axis(at = c(1:100), labels = sortcommeans$disttoname, 
      side = 1, las = 2, cex.axis = 0.5)
 
 legend(-1, 1.75, c("Ancestral", "Derived"), fill = c("red", "black"), bty = "n")
 legend(-1, 1, c("0%", "50%", "75%"), pch = c(17, 18, 16), bty = "n")
 text(-1, 1, "% Initial microbiome from other host", pos = 4)
 
dev.off()
 



se <- function(x) {sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)])) }

cdmeans <- tapply(endcoms$RGR, endcoms$comdist, mean)
cdses <- tapply(endcoms$RGR, endcoms$comdist, se)
png("~/comdistEndComsRGR.png",units="in",height=2,width=4, res=1200)
plot(cdmeans~c(0:6), pch=1,cex=2, xlab="",ylab="", ylim=range(c(cdmeans+cdses,cdmeans-cdses)))
arrows(c(0:6),y0=cdmeans-cdses,y1=cdmeans+cdses,length=0,lwd=2)
dev.off()

pomeans <- tapply(endcoms$RGR, endcoms$percentother, mean)
poses <- tapply(endcoms$RGR, endcoms$percentother, se)
png("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/percentotherEndComsRGR.png",units="in",height=3.5,width=4.25, res=1200)#can set bg to transparent
par(mar=c(3,3,0,1))
plot(pomeans~c(0:2), pch=16,cex=2, xlab="",ylab="", xlim=c(-0.5,2.5), xaxt="n", ylim=range(c(1.1*(pomeans+poses),0.9*(pomeans-poses))))
arrows(c(0:2),y0=pomeans-poses,y1=pomeans+poses,length=0,lwd=2)
  axis(side=1, at=c(0,1,2),labels=c("0","50","75"))
  mtext("%Initial microbiome from other host", line=2, side=1)
  mtext("Derived communities RGR", line=2, side=2)
dev.off()




hist(endcoms$improvement,breaks=60)#c(-25,-16,-8,-4,-2,-1,-0.5,0.25,0,0.25,0.5,1,2,4,8,16,25)) #,xlim=c(-10,10),xlab="Improvement in RGR",main="Distribution of Improvement Scores")
hist(endcoms$RGR,breaks=60)#c(-25,-16,-8,-4,-2,-1,-0.5,0.25,0,0.25,0.5,1,2,4,8,16,25)) #,xlim=c(-10,10),xlab="Improvement in RGR",main="Distribution of Improvement Scores")
plot(endcoms$RGR ~ endcoms$startRGR)
plot(endcoms$RGR ~ endcoms$percentother)
plot(endcoms$RGR ~ endcoms$comdist)

summary( MCMCglmm(improvement~ comdist + percentother + comdist:percentother, data=endcoms,verbose=F,  nitt = 100000, burnin = 5000, thin = 10))

summary( MCMCglmm(improvement~ log_divergence_time + percentother + log_divergence_time:percentother, data=endcoms,verbose=F,  nitt = 100000, burnin = 5000, thin = 10))

summary (MCMCglmm(improvement~ percentother, data=endcoms,verbose=F, nitt=100000, burnin = 5000, thin = 10))
tapply(endcoms$improvement, endcoms$percentother, mean)

summary (MCMCglmm(RGR~ percentother, data=endcoms,verbose=F, nitt=100000, burnin = 5000, thin = 10))
tapply(endcoms$RGR, endcoms$percentother, mean)


summary( MCMCglmm(RGR~startRGR + log_divergence_time:startRGR + percentother:startRGR, data=endcoms,verbose=F, nitt = 10000, burnin = 5000, thin = 10))


#_______________________________________________________________________________________
# ANOVA of improvement by percentother
anova_model <- aov(improvement ~ factor(percentother), data = endcoms)
summary(anova_model)


#doing all the checks youre supposed to do before running an ANOVA
# Step 1: Check normality of improvement
shapiro.test(endcoms$improvement) #violates, so shouldnt do anova

library(dunn.test)
#Kruskal-Wallis test instead?
kruskal.test(improvement ~ factor(percentother), data = endcoms)


dunn_result <- dunn.test(endcoms$improvement, endcoms$percentother, method = "bonferroni")

# View the result
summary(dunn_result)


#_______________________________________________________________________________________

png("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/start_vs_endmeans.png",
    units = "in", height = 7, width = 8, res = 1200)

# Set up the plot with no points and suppress the default axis
plot(c(startmeans, endmeans) ~ c(rep(1, length(startmeans)), rep(2, length(endmeans))),
     pch = NA, xaxt = "n", xlab = "", ylab = "Mean")

# Add only the two desired x-axis ticks
axis(1, at = c(1, 2), labels = c("", ""))

# Draw the arrows
arrows(1, rep(startmeans, each = 3), 2, endmeans, length = 0)

dev.off()




growthslodat <- data.frame(syncom=rep(1:25,each=3),
                           replicate = rep(c("A","B","C"), times=25), 
                           growthslope= rep(NA, times= 75))
growthslodat$otherhost <- sapply(1:nrow(growthslodat), function(z) 
            exp1_avg_of_reps_per_round$other_host_source[exp1_avg_of_reps_per_round$syn_com==growthslodat$syncom[z]][1])
growthslodat$percentother <- sapply(1:nrow(growthslodat), function(z) 
  exp1_avg_of_reps_per_round$percentother[exp1_avg_of_reps_per_round$syn_com==growthslodat$syncom[z]][1])
growthslodat$otherdistance <- sapply(1:nrow(growthslodat), function(z) 
  exp1_avg_of_reps_per_round$distance_from_lemna[exp1_avg_of_reps_per_round$syn_com==growthslodat$syncom[z]][1])

for(s in 1:25){
  for(r in c("A","B","C")){
          dat <- exp1_avg_of_reps_per_round[which(exp1_avg_of_reps_per_round$syn_com == s & exp1_avg_of_reps_per_round$rep == r),]
          growthslodat$growthslope[which(growthslodat$syncom==s & growthslodat$replicate==r)] <- lm(dat$avg_of_rep~dat$round_num)$coef[2]
  }
}

 summary(MCMCglmm(growthslope~percentother*otherdistance,data=growthslodat,verbose=F,nitt=100000,thin=100,burnin=1000))
 summary(lm(growthslope~percentother*otherdistance,data=growthslodat,verbose=F))


# new figure, showing two panels one 50% and one 75% other host. Non-filled dots for control

png("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/all_lemna_vs_50_75.png", width = 7, height = 3, units = "in", res = 1200, bg = "transparent")
layout(matrix(1:2, ncol = 2))  # Two side-by-side plots

par(mar = c(4, 0, 1, 0))       # Inner margins: bottom, left, top, right
par(oma = c(0, 4, 0, 1))       # Outer margins

for (p in c(50, 75)) {
  plot(growthslope ~ otherdistance, 
       data = growthslodat[growthslodat$percentother == p, ],
       yaxt = "n", ylim = c(-200, 700), pch = 16, ylab = "", xlab = "", xlim = c(0, 6))
  
  if (p == 50) {
    axis(side = 2)
    mtext("Mean growth rate in pixels", side = 2, line = 2)
  }
  
  if (p == 75) {
    mtext("Level of phylogenetic distance", side = 1, line = 2, at = -0.25)
    
    # Add legend INSIDE the plot in top right
    legend(
      "topright",
      legend = c("50% non-host associated", "75% non-host associated", "Control (0%)"),
      pch = c(16, 16, 1),          # filled circles for 50% and 75%, open circle for control
      col = c("black", "black", "black"),
      pt.cex = 1.2,
      bty = "n",
      cex = 0.8,
      horiz = FALSE,
      inset = c(0, 0)
    )
  }
  
  # Add control group points
  points(growthslodat$growthslope[growthslodat$percentother == 0] ~ 
           growthslodat$otherdistance[growthslodat$percentother == 0], 
         pch = 1)
}

dev.off()




