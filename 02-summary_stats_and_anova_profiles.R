#-----------------------------------------------------------------------------#
# The R script performs statistical analysis on scuba diving profile data.    #
# It calculates mean temperature, standard deviation of temperature, and      #
# the slope of temperature vs. depth for each filename, State, and Sector.    #
# The script then conducts ANOVA tests to compare these variables across      #
# different Sectors after a specific state and across different States        #
# within the "South" sector. Post-hoc analysis using Tukey's HSD test         #
# follows each ANOVA to identify significant differences between group means. #
#-----------------------------------------------------------------------------#

# Load necessary libraries
library(tidyverse)
library(broom)
library(reshape2)
library(ggpubr)

# Read in the dataset, please mind the path on local machine

scuba_profiles <- read_csv("profile_data.csv")

# Create a summary table with mean temperature, standard deviation of temperature,
# and the slope of temperature vs. depth for each combination of dive_id, State, and Sector.
summary_table <- scuba_profiles %>%
     group_by(dive_id, state, sector) %>%
     nest() %>%
     mutate(
          # Apply a generalized linear model for each group to calculate the slope
          model = map(data, ~ glm(temperature ~ depth, data = .x)),
          # Extract the coefficients from each model
          tidied = map(model, tidy),
          # Calculate mean temperature for each group
          mean_temperature = map_dbl(data, ~ mean(.$temperature, na.rm = TRUE)),
          # Calculate standard deviation of temperature for each group
          sd_temperature = map_dbl(data, ~ sd(.$temperature, na.rm = TRUE))
     ) %>%
     # Clean the data frame by selecting relevant columns
     select(-data, -model) %>%
     # Expand the tidied data
     unnest(tidied) %>%
     # Filter rows where the term is 'depth'
     filter(term == "depth") %>%
     # Remove unnecessary columns and rename 'estimate' to 'slope'
     select(-term, -statistic, -p.value, -std.error) %>%
     rename(slope = estimate) %>%
     # Remove the grouping structure
     ungroup()

# Display the summary table
summary_table



## Statistical comparisons

# Between Sectors of After state ------------------------------------------

# Filter data for comparisons between sectors after a specific state
sectors <-
     summary_table %>%
     filter(state == "After") %>% 
     mutate(sector = as.factor(sector), state = as.factor(state))

# Reshape data to long format for ANOVA
long_format_data <-
     melt(sectors, id.vars = c("dive_id", "state", "sector"))

# Perform ANOVA for mean temperature by sector
anova_mean_temp <-
     aov(value ~ sector , data = long_format_data[long_format_data$variable == "mean_temperature", ])
# Display summary of ANOVA
summary(anova_mean_temp)

# Perform ANOVA for standard deviation of temperature by sector
anova_sd_temp <-
     aov(value ~ sector, 
         data = long_format_data[long_format_data$variable == "sd_temperature", ])
# Display summary of ANOVA
summary(anova_sd_temp)

# Perform ANOVA for slope by sector
anova_slope <-
     aov(value ~ sector, data = long_format_data[long_format_data$variable == "slope", ])
# Display summary of ANOVA
summary(anova_slope)


# Between States of South Sector ------------------------------------------


# Filter data for comparisons between states in the South sector
states <-
     summary_table %>% 
     filter(sector == "South") %>% 
     mutate(sector = as.factor(sector), state = as.factor(state))

# Reshape data to long format for ANOVA
long_format_data <-
     melt(states, id.vars = c("dive_id", "state", "sector"))

# Perform ANOVA for mean temperature by state
anova_mean_temp <-
     aov(value ~ state , data = long_format_data[long_format_data$variable == "mean_temperature", ])
# Display summary of ANOVA
summary(anova_mean_temp)

# Perform ANOVA for standard deviation of temperature by state
anova_sd_temp <-
     aov(value ~ state, data = long_format_data[long_format_data$variable == "sd_temperature", ])
# Display summary of ANOVA
summary(anova_sd_temp)

# Perform ANOVA for slope by state
anova_slope <-
     aov(value ~ state, data = long_format_data[long_format_data$variable == "slope", ])
# Display summary of ANOVA
summary(anova_slope)

# Perform Tukey's post-hoc test for mean temperature
TukeyHSD(anova_mean_temp)
# Perform Tukey's post-hoc test for standard deviation of temperature
TukeyHSD(anova_sd_temp)
# Perform Tukey's post-hoc test for slope
TukeyHSD(anova_slope)

summary(anova_slope)
TukeyHSD(anova_slope)
report::report(anova_slope)


p1 <- states %>% 
     mutate(state = factor(state, levels = c("Before", "During", "After"))) %>% 
     ggplot(aes(x = state, y = mean_temperature, fill = state)) +
     geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), width = 0.7) +
     geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(0.8), width = 0.25) +
     coord_cartesian(ylim = c(28, 30)) +  # Setting y-axis limits
     scale_fill_manual(values = c("royalblue", "firebrick", "orange")) +
     theme_bw() +
     labs(x = "", y = "Mean Temperature °C") +
     theme(legend.position = "")
p1
p2 <- states %>% 
     mutate(state = factor(state, levels = c("Before", "During", "After"))) %>% 
     ggplot(aes(x = state, y = sd_temperature, fill = state)) +
     geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), width = 0.7) +
     geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(0.8), width = 0.25) +
     #coord_cartesian(ylim = c(28, 30)) +  # Setting y-axis limits
     scale_fill_manual(values = c("royalblue", "firebrick", "orange")) +
     theme_bw() +
     labs(x = "", y = "Thermal variation °C") +
     theme(legend.position = "")
p2
p3 <- states %>% 
     mutate(state = factor(state, levels = c("Before", "During", "After"))) %>% 
     ggplot(aes(x = state, y = slope, fill = state)) +
     geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), width = 0.7) +
     geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(0.8), width = 0.25) +
     #coord_cartesian(ylim = c(28, 30)) +  # Setting y-axis limits
     scale_fill_manual(values = c("royalblue", "firebrick", "orange")) +
     theme_bw() +
     labs(x = "", y = "Thermal slope (°C/m)") +
     theme(legend.position = "")
p3

library(patchwork)
p1+p2+p3 + plot_annotation(tag_levels = "A")





p1 <- sectors %>% 
     ggplot(aes(x = sector, y = mean_temperature, fill = sector)) +
     geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), width = 0.7) +
     geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(0.8), width = 0.25) +
     coord_cartesian(ylim = c(25, 30)) +  # Setting y-axis limits
     scale_fill_manual(values = c("gray50", "gray90")) +
     theme_bw() +
     labs(x = "", y = "Mean Temperature °C") +
     theme(legend.position = "")
p1
p2 <- sectors %>% 
     ggplot(aes(x = sector, y = sd_temperature, fill = sector)) +
     geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), width = 0.7) +
     geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(0.8), width = 0.25) +
     #coord_cartesian(ylim = c(28, 30)) +  # Setting y-axis limits
     scale_fill_manual(values = c("gray50", "gray90")) +
     theme_bw() +
     labs(x = "", y = "Thermal variation °C") +
     theme(legend.position = "")
p2
p3 <-sectors %>% 
     ggplot(aes(x = sector, y = slope, fill = sector)) +
     geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.8), width = 0.7) +
     geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(0.8), width = 0.25) +
     #coord_cartesian(ylim = c(28, 30)) +  # Setting y-axis limits
     scale_fill_manual(values = c("gray50", "gray90")) +
     theme_bw() +
     labs(x = "", y = "Thermal slope (°C/m)") +
     theme(legend.position = "")
p3

library(patchwork)
p1+p2+p3 + plot_annotation(tag_levels = "A")





### END OF SCRIPT //

