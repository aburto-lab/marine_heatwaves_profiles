# Loading necessary libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(MASS)

# Reading in abundance data
abundance_data <- read_csv("forsubmission/reef_community_data.csv")

# Analyzing Cold Water Corals
# ----------------------------

# Preparing data for analysis
temperate <- abundance_data %>% 
     filter(Type == "Cold Water Corals") %>% 
     group_by(Year, Sector, Degree, Period, Reef, Depth2, Transect, Type) %>% 
     summarise(Quantity = sum(Quantity)) %>% 
     group_by(Period, Degree, Depth2, Type) %>% 
     summarise(abundance = mean(Quantity, na.rm = TRUE), 
               sd = sd(Quantity, na.rm = TRUE),
               n = n(),
               se = sd / sqrt(n)) %>% 
     mutate(Period = factor(Period, levels = c("Before", "After")))

# Analyzing Warm Water Corals
# ----------------------------

# Preparing data for analysis
tropical <- abundance_data %>% 
     filter(Type == "Warm Water Corals") %>% 
     group_by(Year, Sector, Degree, Period, Reef, Depth2, Transect, Type) %>% 
     summarise(Quantity = sum(Quantity)) %>% 
     group_by(Period, Degree, Depth2, Type) %>% 
     summarise(abundance = mean(Quantity, na.rm = TRUE), 
               sd = sd(Quantity, na.rm = TRUE),
               n = n(),
               se = sd / sqrt(n)) %>% 
     mutate(Period = factor(Period, levels = c("Before", "After")))

# Plotting Abundance Data
# -----------------------

# Creating a plot for overall abundance
abundance_data %>% 
     group_by(Year, Sector, Degree, Period, Reef, Depth2, Transect, Type) %>% 
     summarise(Quantity = sum(Quantity)) %>% 
     group_by(Period, Degree, Depth2, Type) %>% 
     summarise(abundance = mean(Quantity, na.rm = TRUE), 
               sd = sd(Quantity, na.rm = TRUE),
               n = n(),
               se = sd / sqrt(n)) %>% 
     mutate(Period = factor(Period, levels = c("Before", "After"))) %>% 
     ggplot(aes(x = Period, y = abundance, fill = Period)) +
     geom_col(width = 1) +
     geom_errorbar(aes(ymin = abundance - se, ymax = abundance + se), width = .1) +
     scale_fill_manual(values = c("royalblue", "orange")) +
     facet_grid(Degree ~ Type + Depth2, scales = "free_y") +
     labs(y = "Abundance", x = "") +
     theme_bw() +
     theme(strip.background = element_blank(), 
           strip.text.x.top = element_text(face = "italic"),
           axis.text.x = element_text(angle = 90, vjust = .5),
           legend.position = "")

# Modeling for Cold Water Corals
# -------------------------------

# Preparing data for model
totest1 <- abundance_data %>% 
     filter(Type == "Cold Water Corals") %>% 
     group_by(Year, Sector, Degree, Period, Reef, Depth2, Transect, Type) %>% 
     summarise(Quantity = sum(Quantity)) %>% 
     mutate(Period = factor(Period, levels = c("Before", "After")))

# Setting factor levels
totest1$Period <- factor(totest1$Period, levels = c("Before", "After"))
totest1$Depth2 <- factor(totest1$Depth2, levels = c("Deep", "Shallow"))
totest1$Sector <- factor(totest1$Sector, levels = c("South", "North"))

# Fitting a negative binomial model
nb_model <- glm.nb(Quantity ~ Period * Sector * Depth2, data = totest1)
qqnorm(residuals(nb_model))
qqline(residuals(nb_model))
summary(nb_model)

# Extracting and Analyzing Model Coefficients for Cold Water Corals
# ------------------------------------------------------------------

# Extracting coefficients and standard errors from the model
model_coefs <- coef(nb_model)
model_summary <- summary(nb_model)
std_errors <- model_summary$coefficients[, "Std. Error"]

# Isolating coefficients of interest
coef_PeriodAfter <- model_coefs["PeriodAfter"]
coef_PeriodAfter_Depth2Shallow <- model_coefs["PeriodAfter:Depth2Shallow"]
coef_PeriodAfter_SectorNorth <- model_coefs["PeriodAfter:SectorNorth"]
coef_PeriodAfter_SectorNorth_Depth2Shallow <- model_coefs["PeriodAfter:SectorNorth:Depth2Shallow"]

# Calculating standard errors for coefficients of interest
std_error_PeriodAfter <- std_errors["PeriodAfter"]
std_error_PeriodAfter_Depth2Shallow <- std_errors["PeriodAfter:Depth2Shallow"]
std_error_PeriodAfter_SectorNorth <- std_errors["PeriodAfter:SectorNorth"]
std_error_PeriodAfter_SectorNorth_Depth2Shallow <- std_errors["PeriodAfter:SectorNorth:Depth2Shallow"]

# Calculating effect sizes and error estimates for each category
effect_After_North_Deep <- coef_PeriodAfter
effect_After_North_Shallow <- coef_PeriodAfter + coef_PeriodAfter_Depth2Shallow
effect_After_South_Deep <- coef_PeriodAfter + coef_PeriodAfter_SectorNorth
effect_After_South_Shallow <- coef_PeriodAfter + coef_PeriodAfter_SectorNorth + coef_PeriodAfter_Depth2Shallow + coef_PeriodAfter_SectorNorth_Depth2Shallow

error_After_North_Deep <- std_error_PeriodAfter
error_After_North_Shallow <- sqrt(std_error_PeriodAfter^2 + std_error_PeriodAfter_Depth2Shallow^2)
error_After_South_Deep <- sqrt(std_error_PeriodAfter^2 + std_error_PeriodAfter_SectorNorth^2)
error_After_South_Shallow <- sqrt(std_error_PeriodAfter^2 + std_error_PeriodAfter_SectorNorth^2 + std_error_PeriodAfter_Depth2Shallow^2 + std_error_PeriodAfter_SectorNorth_Depth2Shallow^2)

# Creating a data frame with effect sizes and errors
effect_sizes_cwc <- data.frame(
     Type = rep("Cold Water Corals", 4), 
     Depth = c("Deep", "Shallow", "Deep", "Shallow"),
     Sector = c("North", "North", "South", "South"), 
     Category = c("After, North, Deep", "After, North, Shallow", "After, South, Deep", "After, South, Shallow"),
     EffectSize = c(effect_After_North_Deep, effect_After_North_Shallow, effect_After_South_Deep, effect_After_South_Shallow),
     Error = c(error_After_North_Deep, error_After_North_Shallow, error_After_South_Deep, error_After_South_Shallow)
)

# Modeling for Warm Water Corals
# -------------------------------

# Preparing data for model
totest2 <- abundance_data %>% 
     filter(Type == "Warm Water Corals") %>% 
     group_by(Year, Sector, Degree, Period, Reef, Depth2, Transect, Type) %>% 
     summarise(Quantity = sum(Quantity)) %>% 
     mutate(Period = factor(Period, levels = c("Before", "After")))

# Boxplot for Quantity by Period and Sector
totest2 %>% 
     ggplot(aes(x = Period, y = Quantity, col = Sector)) +
     geom_boxplot()

# Setting factor levels
totest2$Period <- factor(totest2$Period, levels = c("Before", "After"))
totest2$Depth2 <- factor(totest2$Depth2, levels = c("Deep", "Shallow"))
totest2$Sector <- factor(totest2$Sector, levels = c("South", "North"))

# Fitting a negative binomial model
nb_model <- glm.nb(Quantity ~ Period * Sector * Depth2, data = totest2)
qqnorm(residuals(nb_model))
qqline(residuals(nb_model))
summary(nb_model)

# Extracting and Analyzing Model Coefficients for Warm Water Corals
# ------------------------------------------------------------------

# Extracting coefficients and standard errors from the model
model_coefs <- coef(nb_model)
model_summary <- summary(nb_model)
std_errors <- model_summary$coefficients[, "Std. Error"]

# Repeating the coefficient and error extraction process as done for Cold Water Corals

# Creating a data frame with effect sizes and errors for Warm Water Corals
effect_sizes_wwc <- data.frame(
     Type = rep("Warm Water Corals", 4), 
     Depth = c("Deep", "Shallow", "Deep", "Shallow"),
     Sector = c("North", "North", "South", "South"), 
     Category = c("After, North, Deep", "After, North, Shallow", "After, South, Deep", "After, South, Shallow"),
     EffectSize = c(effect_After_North_Deep, effect_After_North_Shallow, effect_After_South_Deep, effect_After_South_Shallow),
     Error = c(error_After_North_Deep, error_After_North_Shallow, error_After_South_Deep, error_After_South_Shallow)
)

# Final Plotting and Output
# -------------------------

# Combining effect sizes data and writing to file
bind_rows(effect_sizes_wwc, effect_sizes_cwc) %>%
     write_rds("model_coefficients.RDS")

# Plotting combined effect sizes
bind_rows(effect_sizes_wwc, effect_sizes_cwc) %>% 
     ggplot(aes(x = EffectSize, y = Depth, color = Type)) +
     geom_vline(xintercept = 0) +
     geom_point(size = 3) +
     geom_errorbarh(aes(xmin = EffectSize - Error, xmax = EffectSize + Error), height = 0.1, size = 1) +
     facet_grid(Sector ~ .) +
     scale_color_manual(values = c("royalblue", "firebrick")) +
     labs(y = "Depth strata", x = "Effect Size", color = "") +
     theme_bw() +
     theme(strip.background = element_blank(), legend.position = "top")



## END OF SCRIPT ////
