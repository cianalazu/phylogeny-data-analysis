#for the figures i had made, turn into rgr, split into phylogenetic levels (7 plots)
#make them all share a y-axis, phylogenetic title label,
#lines should be colored by percent other

library(tidyr)
library(ggplot2)
library(dplyr)
library(MCMCglmm)

#set working directory
setwd("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1")

exp1_data <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/exp1_data.csv")

print(unique(exp1_data$percent_other))

# Clean up the 'rep' column: Remove leading/trailing spaces and convert to factor
exp1_data$rep <- trimws(exp1_data$rep)
exp1_data$rep <- factor(exp1_data$rep)

print(unique(exp1_data$percent_other))

# Check the levels to confirm that the space is removed
print(levels(exp1_data$rep))

exp1_data <- exp1_data |>  
  filter(exp1_data$start_area > 0)

print(unique(exp1_data$percent_other))

# Filter out rows where 'start_area' is less than or equal to zero
exp1_data <- exp1_data[exp1_data$start_area > 0, ]
exp1_data$daysinround <- sapply(exp1_data$round_num, function(z) ifelse(z==4, 2,5))
exp1_data$RGR <- ((exp1_data$end_area - exp1_data$start_area) / exp1_data$start_area)/exp1_data$daysinround

print(unique(exp1_data$percent_other))

#_______________________________________________________________________________
#_____________RGR over rounds plot______________________________________________

exp1_data_avgRGR <- exp1_data |>
  group_by(round_num, distance_from_lemna, rep) |>
  summarise(RGR = mean(RGR))


rgr_trajectories_plot <- exp1_data_avgRGR |>
  ggplot(aes(x = round_num, y = RGR, color = rep)) +
  geom_line() +
  ggtitle("RGR by Distance from Lemna") +
  xlab("Round Number") +
  ylab("RGR") + 
  scale_x_continuous(breaks = 1:10) +  # Ensure x-axis shows round number
  scale_y_continuous(limits = c(-0.02,.57)) + # Set y-axis limits
  facet_wrap(~ distance_from_lemna, scales = "free_y", ncol =1) +
  theme_minimal()

#save plot
ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/rgr_trajectories_plot.png", plot = rgr_trajectories_plot, width = 10, height = 15, dpi = 300, bg ="white")

#_______________________________________________________________________________





#_______________________________________________________________________________


exp1_data$log_divergence_time <- log(exp1_data$divergence_times_ma + 1)
#square log_divergence_time
exp1_data$log_divergence_time_sq <- exp1_data$log_divergence_time^2
#improvement score column

print(unique(exp1_data$percent_other))

# Step 1: Aggregate RGR to ensure no duplicates
summarized_data <- exp1_data %>%
  filter(round_num %in% c(1, 10)) %>%
  group_by(syn_com, rep, round_num) %>%
  summarise(RGR = mean(RGR, na.rm = TRUE), .groups = "drop")

# Step 2: Pivot to wide format
improvement_data <- summarized_data %>%
  pivot_wider(names_from = round_num, values_from = RGR, names_prefix = "RGR_day") %>%
  mutate(improvement = (RGR_day10 - RGR_day1) / abs( RGR_day1))

# Step 3: Join back to original exp1_data
exp1_data <- exp1_data %>%
  left_join(improvement_data %>% select(syn_com, rep, improvement),
            by = c("syn_com", "rep"))

#no lemna in dataframe
exp1_data_minlemna <- exp1_data %>% 
  filter(distance_from_lemna != 0)


#save exp1_data to csv
write.csv(exp1_data, "C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/exp1_data_v2", row.names = FALSE)


#model with exp1 improvement scores LEMNA INCLUDED
improvement_percentother <- MCMCglmm(improvement ~ percent_other, data = exp1_data, verbose=T, nitt = 10000, burnin = 500, thin = 10)
summary(improvement_percentother)
tapply(exp1_data$improvement, exp1_data$percent_other, mean)


exp1avgimprovement_df <- exp1_data |> 
  filter(percent_other %in% c(0, 50, 75)) |> 
  group_by(percent_other) |> 
  summarise(
    mean_improvement = mean(improvement, na.rm = TRUE),
    se_improvement = sd(improvement, na.rm = TRUE) / sqrt(n())
  )

# Step 2: Plot with SE error bars
ggplot(exp1avgimprovement_df, aes(x = factor(percent_other), y = mean_improvement)) +
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
#save plot
ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/exp1_avg_improvement_plot.png", width = 6, height = 4, dpi = 300, bg ="white")



#models with exp1 improvement scores MIN_LEMNA
# improvement_percentother_minlemna <- MCMCglmm(improvement ~ percent_other, data = exp1_data_minlemna, verbose=T, nitt = 10000, burnin = 500, thin = 10)
# summary(improvement_percentother_minlemna)
# tapply(exp1_data_minlemna$improvement, exp1_data_minlemna$percent_other, mean)

improvement_log_divergence_min_lemna <- MCMCglmm(improvement ~ log_divergence_time, data = exp1_data_minlemna, verbose=T, nitt = 10000, burnin = 500, thin = 10)
summary(improvement_log_divergence_min_lemna)
tapply(exp1_data_minlemna$improvement, exp1_data_minlemna$log_divergence_time, mean)

improvement_log_divergence_min_lemnasquared <- MCMCglmm(improvement ~ log_divergence_time + log_divergence_time_sq, data = exp1_data_minlemna, verbose=T, nitt = 10000, burnin = 500, thin = 10)
summary(improvement_log_divergence_min_lemnasquared)
tapply(exp1_data_minlemna$improvement, exp1_data_minlemna$log_divergence_time_sq, mean)
tapply(exp1_data_minlemna$improvement, exp1_data_minlemna$log_divergence_time, mean)




#_______________________________________________________________________________
#_____________Improvement score by divergence fit to model__________________________________
# Function to calculate standard error
se <- function(x) {sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)])) }

exp1_improvement_means <- tapply(exp1_data_minlemna$improvement, exp1_data_minlemna$distance_from_lemna, mean)
exp1_improvement_ses <- tapply(exp1_data_minlemna$improvement, exp1_data_minlemna$distance_from_lemna, se)

exp1_modelmclog <- MCMCglmm(improvement ~ log_divergence_time + log_divergence_time_sq, data = exp1_data_minlemna,  nitt = 10000, burnin = 500, thin = 10,verbose=T)
exp1_solmclog <- exp1_modelmclog$Sol

# Print the summary of the model
summary(exp1_modelmclog)

ldivlin_seq <- seq(min(exp1_data_minlemna$log_divergence_time), max(exp1_data_minlemna$log_divergence_time), length.out = 1000)
ldivsq_seq <- ldivlin_seq^2

exp1_lmean_pred <- sapply(1:length(ldivlin_seq),
                          function(z) mean(exp1_solmclog[,1] + exp1_solmclog[,2] * ldivlin_seq[z] + exp1_solmclog[,3] *ldivsq_seq[z]) )



#save plot
png("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/exp1_improvement_by_divergence.png", width = 800, height = 600, res = 100, bg = "white")
plot(exp1_improvement_means~sort(unique(exp1_data_minlemna$log_divergence_time)), 
     pch=1,cex=2, xlab="log Divergence Time",ylab="Average Improvement of Derived Communities", 
     ylim=range(c(exp1_improvement_means+exp1_improvement_ses,exp1_improvement_means-exp1_improvement_ses)))
arrows(sort(unique(exp1_data_minlemna$log_divergence_time)),
       y0=exp1_improvement_means-exp1_improvement_ses,y1=exp1_improvement_means+exp1_improvement_ses,length=0,lwd=2)
lines(ldivlin_seq, exp1_lmean_pred, col="red", lwd=2)
dev.off()








#________PLOTS______________________________________________

par(mfrow=c(1,7))
par(oma=c(0,5,0,0))
par(mar=c(4,0,1,0))
for(d in 0:6){
  datdist <- exp1_data_avg[exp1_data_avg$distance_from_lemna==d,]
  plot(1~c(1),pch=NA,ylim=c(0,1),xlim=c(1,10),ylab="",xlab="",yaxt="n")
  for(s in unique(datdist$syn_com)){
    for(o in unique(datdist$percent_other)){
      lines(datdist$RGR[datdist$syn_com==s & datdist$percent_other==o ]~datdist$round_num[datdist$syn_com==s & datdist$percent_other==o ],
            col=rgb(o/100,0,0))
    }}
}


par(mfrow=c(1,7))
par(oma=c(0,5,0,0))
par(mar=c(4,0,1,0))
for(d in 0:6){
  datdist <- exp1_data_avg[exp1_data_avg$distance_from_lemna==d,]
  plot(1~c(1),pch=NA,ylim=c(0,0.5),xlim=c(1,10),ylab="",xlab="",yaxt="n")
  for(o in unique(datdist$percent_other)){
    lines(datdist$RGR[datdist$percent_other==o ]~datdist$round_num[datdist$percent_other==o ],
          col=rgb(o/100,0,0))
  }
  if(d==0){
    axis(2)
    mtext("RGR",side=2,line=2)
  }
  if(d==3){mtext("Round", side=1,line=2)}
}

#______________________________________________________________________________________





