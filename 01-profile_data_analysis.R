#-----------------------------------------------------------------------------#
# The R script performs statistical analysis on scuba diving profile data.    # 
# The script performs data wrangling, and visualizes                          #
# temperature profiles and anomalies. Additionally, it conducts linear and    #
# non-linear modeling to understand the effects of various environmental      #
# factors on temperature.                                                     #
#-----------------------------------------------------------------------------#

# Loading necessary libraries
library(tidyverse)
library(mgcv)
library(patchwork)
library(broom)
library(reshape2)

# Reading in profile data
profile_data <- read_csv("profile_data.csv")

# Temperature Profile Analysis in the South Sector
# -------------------------------------------------

# Calculating mean and standard deviation of temperature
p_ba <- profile_data %>% 
     filter(sector == "South") %>% 
     mutate(state = factor(state, levels = c("Before", "During", "After")),
            depth2 = round(depth)) %>% 
     group_by(depth2, state) %>% 
     summarise(mean_temperature = mean(temperature, na.rm = TRUE), 
               sd_temperature = sd(temperature, na.rm = TRUE)) %>% 
     ggplot(aes(x = mean_temperature, y = -depth2, col = state, fill = state)) +
     geom_ribbon(aes(xmin = mean_temperature - sd_temperature, 
                     xmax = mean_temperature + sd_temperature), 
                 alpha = 0.2, col = NA) +
     geom_line(orientation = "y") +
     scale_color_manual(values = c("royalblue", "firebrick", "orange")) +
     scale_fill_manual(values = c("royalblue", "firebrick", "orange")) +
     labs(x = "Temperature (°C)", y = "Depth (m)") +
     theme_bw() +
     theme(legend.position = "top", strip.background = element_blank())
p_ba

# Analyzing Temperature Anomalies
# --------------------------------

# Preparing data for temperature anomaly analysis
ad <- profile_data %>% 
     filter(sector == "South") %>% 
     mutate(state = factor(state, levels = c("Before", "During", "After")),
            depth2 = round(depth)) %>% 
     group_by(depth2, state) %>% 
     summarise(temperature = mean(temperature, na.rm = TRUE)) %>% 
     pivot_wider(names_from = "state", values_from = "temperature") %>% 
     mutate(During = During - Before, After = After - Before) %>% 
     select(depth2, During, After) %>% 
     pivot_longer(During:After, names_to = "state", values_to = "anomaly") %>% 
     mutate(state = factor(state, levels = c("During", "After"))) %>% 
     ggplot(aes(x = anomaly, y = -depth2)) +
     geom_line(aes(col = state, fill = state), orientation = "y") +
     scale_color_manual(values = c("firebrick", "orange")) +
     labs(x = "Temperature anomalies (°C)", y = "Depth (m)") +
     scale_x_continuous(breaks = seq(0, 5.5, by = 1)) +
     theme_bw() +
     theme(legend.position = "", strip.background = element_blank())
ad

# Exponential Model Fitting
# --------------------------

# Preparing data for modeling
to_model <- profile_data %>% 
     filter(sector == "South") %>% 
     mutate(state = factor(state, levels = c("Before", "During", "After")),
            depth2 = round(depth)) %>% 
     group_by(depth2, state) %>% 
     summarise(temperature = mean(temperature, na.rm = TRUE)) %>% 
     pivot_wider(names_from = "state", values_from = "temperature") %>% 
     mutate(During = During - Before, After = After - Before) %>% 
     select(depth2, During, After) %>% 
     pivot_longer(During:After, names_to = "state", values_to = "anomaly") %>% 
     mutate(state = factor(state, levels = c("During", "After"))) 

# Fitting a linear model and using it as a starting point for non-linear fitting
lm_model <- lm(log(anomaly) ~ depth2, data = to_model)
start_a <- exp(coef(lm_model)[1])  # Exponential of the intercept
start_b <- coef(lm_model)[2]       # Slope

# Fitting an exponential model
model <- nls(anomaly ~ a * exp(b * depth2), data = to_model, start = list(a = start_a, b = start_b))

# Predicting values using the model
depth_range <- seq(min(to_model$depth2, na.rm = TRUE), max(to_model$depth2, na.rm = TRUE), length.out = 100)
predicted_anomaly <- predict(model, newdata = data.frame(depth2 = depth_range))

# Extracting the model equation
coefficients <- coef(model)
equation <- paste("y =", round(coefficients["a"], 2), "* e^(", round(coefficients["b"], 2), "* x)")

# Calculating R-squared value
predicted_values <- predict(model, newdata = to_model)
rss <- sum((to_model$anomaly - predicted_values)^2, na.rm = TRUE)
tss <- sum((to_model$anomaly - mean(to_model$anomaly, na.rm = TRUE))^2, na.rm = TRUE)
r_squared <- round(1 - (rss / tss), 2)
r_squared

# Plotting the curve with annotations
p2 <- ggplot(to_model) +
     geom_line(aes(y = anomaly, x = -depth2, col = state, fill = state)) +
     geom_line(data = data.frame(depth2 = depth_range, anomaly = predicted_anomaly), 
               aes(y = anomaly, x = -depth2), col = "black") +
     scale_color_manual(values = c("firebrick", "orange")) +
     labs(y = "Temperature anomalies (°C)", x = "Depth (m)") +
     theme_bw() +
     theme(legend.position = "", strip.background = element_blank(),
           plot.background = element_blank()) +
     coord_flip()
p2


# Analysis of Coral Data
# -----------------------

# Reading coral data and creating a plot
# this part depends on the script 03-reef_community_analysis.R which 
# creates the model coefficients, here is added just to reproduce Fig. 2
# please refer to that for further details. 


pcorals <- read_rds("outputs/model_coefficients.RDS") %>% 
     mutate(Sector = factor(Sector, levels = c("South", "North")),
            Depth = factor(Depth, levels = c("Deep", "Shallow"))) %>% 
     ggplot(aes(x = EffectSize, y = Depth, color = Type, 
                shape = Sector, group = interaction(Sector, Depth))) +
     geom_vline(xintercept = 0) +
     geom_errorbarh(aes(xmin = EffectSize - Error, xmax = EffectSize + Error), 
                    position = position_dodge(width = 0.25),
                    height = 0.1, linewidth = 1) +
     geom_point(size = 3, position = position_dodge(width = 0.25), stroke = 1) +
     scale_color_manual(values = c("royalblue", "firebrick")) +
     scale_shape_manual(values = c(16, 17), guide = "none") +
     xlim(-1.5, 1.8) +
     labs(y = "", x = "Effect Size", color = "") +
     theme_bw() +
     theme(strip.background = element_blank(), 
           plot.background = element_blank(), 
           legend.position = "top")

# Combining Plots
# ----------------

# Combining previous plots with annotations
combined_plot <- p_ba + p2 + pcorals + plot_annotation(tag_levels = "A")
combined_plot

# Analysis for Northern and Southern Sectors
# ------------------------------------------

# Temperature Profile Analysis after 'After' State
p_ns <- profile_data %>% 
     group_by(dive_id) %>% 
     mutate(anomaly = temperature - seas) %>% 
     filter(state == "After") %>% 
     mutate(sector = factor(sector, levels = c("North", "South"))) %>% 
     mutate(depth2 = round(depth)) %>% 
     group_by(depth2, sector) %>% 
     summarise(mean_temperature = mean(temperature, na.rm = TRUE), 
               sd_temperature = sd(temperature, na.rm = TRUE)) %>% 
     ggplot(aes(x = mean_temperature, y = -depth2, shape = sector, col = sector)) +
     geom_ribbon(aes(xmin = mean_temperature - sd_temperature, 
                     xmax = mean_temperature + sd_temperature), 
                 alpha = 0.2, col = NA) +
     geom_line(orientation = "y") +
     geom_point(size = 2, alpha = 0.5) +
     scale_shape_manual(values = c(24, 21)) +
     scale_color_manual(values = c("black", "gray60")) +
     labs(x = "Temperature (°C)", y = "Depth (m)") +
     theme_bw() +
     theme(legend.position = "top", strip.background = element_blank())
p_ns

# Linear and Non-linear Modelling
# --------------------------------

# Modelling for South Sector
tomod <- profile_data %>% 
     filter(sector == "South")

combined_model_south <- lm(temperature ~ state * depth + state * I(depth^2), data = tomod)
summary(combined_model_south)
report::report(combined_model_south)

# Modelling for 'After' State across Sectors
tomod <- profile_data %>% 
     filter(state == "After")

combined_model_after <- lm(temperature ~ sector * depth + sector * I(depth^2), data = tomod)
summary(combined_model_after)
report::report(combined_model_after)

