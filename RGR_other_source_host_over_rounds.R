# R script to analyze and visualize RGR trajectories across different host sources
# lost the plot of RGR for the seven host types (stacked plots)
# x axis should be round_num 1-10
# y-axis is the RGR
# each panel is other_host_source

library(dplyr)
library(ggplot2)
library(MCMCglmm)

# 1. Clean data: drop rows with NA in relevant columns
samples_df_clean <- samples_df %>%
  filter(!is.na(rep), !is.na(RGR), !is.na(round_num), !is.na(other_host_source))

# 2. Calculate mean RGR per rep, round_num, other_host_source
rep_means_df <- samples_df_clean %>%
  group_by(other_host_source, rep, round_num) %>%
  summarise(mean_RGR = mean(RGR, na.rm = TRUE), .groups = "drop")

# 3. Calculate overall mean and SE per round_num and other_host_source (from raw data)
mean_se_df <- samples_df_clean %>%
  group_by(other_host_source, round_num) %>%
  summarise(
    mean_RGR = mean(RGR, na.rm = TRUE),
    se_RGR = sd(RGR, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )



# Ensure rep is a factor
rep_means_df$rep <- as.factor(rep_means_df$rep)

# Create a new column in rep_means_df for the legend grouping
rep_means_df <- rep_means_df %>%
  mutate(group = "Replicates")

# Create a new column in mean_se_df for the legend grouping
mean_se_df$group <- "Mean"

ggplot() +
  # Replicate lines (pink), with group for legend
  geom_line(data = rep_means_df,
            aes(x = round_num, y = mean_RGR, group = rep, color = group),
            alpha = 0.6) +
  geom_point(data = rep_means_df,
             aes(x = round_num, y = mean_RGR, group = rep, color = group),
             alpha = 0.6) +
  
  # Mean line (medium blue) with error bars
  geom_line(data = mean_se_df,
            aes(x = round_num, y = mean_RGR, color = group),
            size = 1.2) +
  geom_errorbar(data = mean_se_df,
                aes(x = round_num, ymin = mean_RGR - se_RGR, ymax = mean_RGR + se_RGR, color = group),
                width = 0.2) +
  
  # Stack panels vertically with one column
  facet_wrap(~ other_host_source, ncol = 1) +
  
  scale_x_continuous(breaks = 1:10) +
  
  scale_color_manual(
    name = "Legend",
    values = c("Replicates" = "#FF1493", "Mean" = "#0073C2")
  ) +
  
  labs(
    x = "Round Number",
    y = "RGR",
    title = "Replicate RGR Trajectories with Mean ± SE by Host Source"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(), # cleaner panel labels
    strip.text = element_text(size = 12)
  )

# Save the plot
ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/outputs/RGR_trajectories_by_host_source.png", width = 20, height = 16, bg = "white")



lm_model <- lm(RGR ~ other_host_source + round_num, data = samples_df_clean)
summary(lm_model)


prior <- list(
  R = list(V = 1, nu = 0.002),        # residual variance
  G = list(G1 = list(V = 1, nu = 0.002))  # for random effect rep
)

model_mcmc <- MCMCglmm(
  fixed = RGR ~ other_host_source + round_num,
  random = ~ rep,
  data = samples_df_clean,
  prior = prior,
  nitt = 13000, burnin = 3000, thin = 10,
  verbose = FALSE
)
summary(model_mcmc)

