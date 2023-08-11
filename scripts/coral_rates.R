## Coral diversification

# set directory
setwd("C:/Users/lars-/Desktop/Uni/Master/2. Semester/Phylogenetics/Phylogenetics project")

devtools::install_github("https://github.com/revbayes/RevGadgets")
# install.packages(RevGadgets)
library(RevGadgets)
library(ape)
library(ggplot2)
library(tibble)

tree <- read.nexus("data/Scleractinia_species_level.tre")


# path files
speciation_time_file <- paste0("output/corals_EBD_Corr_speciation_times.log")
speciation_rate_file <- paste0("output/corals_EBD_Corr_speciation_rates.log")
extinction_time_file <- paste0("output/corals_EBD_Corr_extinction_times.log")
extinction_rate_file <- paste0("output/corals_EBD_Corr_extinction_rates.log")

rates <- processDivRates(speciation_time_log = speciation_time_file,
                         speciation_rate_log = speciation_rate_file,
                         extinction_time_log = extinction_time_file,
                         extinction_rate_log = extinction_rate_file,
                         burnin = 0.25,
                         summary = "median")

ex_df <- read.table(extinction_time_file, sep ="\t", header = TRUE)
#pdf("output/EBD_Corr.pdf")
#par(mfrow=c(2,2))
#plotDivRates(rev_out, predictor.ages=co2_age, predictor.var=co2, use.geoscale=TRUE)
#dev.off()

# the T values as a reference in our plot
T <- c(21.37600, 15.89934, 17.31874, 16.43609, 20.36971, 20.65145, 17.78779, 16.63133, 18.84282, 16.15000, 18.62500, 17.05000, 27.56098, 18.85000, 29.78249, 31.82526, 35.50000, 27.45749, 25.60000, 25.30000, 24.06004, 26.95000, 30.55000, 24.25000, 21.32500, 22.45000, 24.25000, 21.10000, 21.04708, 23.80000, 21.10000, 20.87500, 24.25000, 14.80000, 18.85000, 18.85000, 22.84327, 21.55000,
21.00411, 23.80000, 19.30000, 17.95000, 17.59783, 20.20000, 17.05000, 15.70000, 22.50625, 20.65000, 13.90000, 18.72000, 16.80857, 17.50000, 25.47000, 20.05000, 20.25429, 18.69483, 19.89531,
20.20000, 20.62614, 24.91147, 21.24230, 22.84717, 23.58299, 20.15921, 19.98539, 16.23663, 15.29296, 23.90600, 23.30763, 27.21080, 30.68000, 22.79000, 26.36300, 27.85192, 29.20420)

MAX_VAR_AGE <- 250
NUM_INTERVALS <- length(T)
T_age <- MAX_VAR_AGE * (1:NUM_INTERVALS) / NUM_INTERVALS
predictor.ages <- T_age
predictor.var <- T



T_df <- tibble("T" = T, "T_age" = T_age)


p1 <- plotDivRates(rates)
p1
## We can add the T measurements on top of the ggplot object
coeff1 <- max(T_df$T)*2 ## Change the relative plotting scale between rate units and co2 units (ppm)
p2 <- p1 +
  geom_line(data=T_df, mapping = aes(x = T_age, y = T/coeff1), inherit.aes = FALSE) +
  scale_y_continuous(
    name = "Rate", ## First y-axis (left)
    sec.axis = sec_axis(~.*coeff1, name="T (°C)") ## Second y-axis (right)
  )
p2
ggsave(paste0("output/Temperature.pdf"), plot = p2, width = 150, height = 120, units = "mm")

summary(p2)

#################
# CO2

# path files
speciation_time_file <- paste0("output/CO2_corals_EBD_Corr_speciation_times.log")
speciation_rate_file <- paste0("output/CO2_corals_EBD_Corr_speciation_rates.log")
extinction_time_file <- paste0("output/CO2_corals_EBD_Corr_extinction_times.log")
extinction_rate_file <- paste0("output/CO2_corals_EBD_Corr_extinction_rates.log")

rates <- processDivRates(speciation_time_log = speciation_time_file,
                         speciation_rate_log = speciation_rate_file,
                         extinction_time_log = extinction_time_file,
                         extinction_rate_log = extinction_rate_file,
                         burnin = 0.25,
                         summary = "median")

ex_df <- read.table(extinction_time_file, sep ="\t", header = TRUE)
#pdf("output/EBD_Corr.pdf")
#par(mfrow=c(2,2))
#plotDivRates(rev_out, predictor.ages=co2_age, predictor.var=co2, use.geoscale=TRUE)
#dev.off()

# the CO2 values as a reference in our plot
CO2 <- c(424.8073,  794.7754,  850.7455,  761.6740,  775.4860,  968.6500, 1246.0200, 1390.9620, 1439.3850, 1742.0267, 2196.3980,
              2490.7600, 2163.2750, 1358.7550,  779.9031,  557.5220,  569.8060,  583.7780,  526.7770,  561.9113,  894.3447, 1782.6520,
              2014.1950, 2037.1800, 2255.3680)
       

MAX_VAR_AGE <- 250
NUM_INTERVALS <- length(CO2)
CO2_age <- MAX_VAR_AGE * (1:NUM_INTERVALS) / NUM_INTERVALS
predictor.ages <- CO2_age
predictor.var <- CO2



CO2_df <- tibble("CO2" = CO2, "CO2_age" = CO2_age)


p3 <- plotDivRates(rates)
p3
## We can add the T measurements on top of the ggplot object
coeff2 <- max(CO2_df$CO2)*2 ## Change the relative plotting scale between rate units and co2 units (ppm)
p4 <- p3 +
  geom_line(data=CO2_df, mapping = aes(x = CO2_age, y = CO2/coeff2), inherit.aes = FALSE) +
  scale_y_continuous(
    name = "Rate", ## First y-axis (left)
    sec.axis = sec_axis(~.*coeff2, name="CO2 (ppm)") ## Second y-axis (right)
  )
p4
ggsave(paste0("output/CO2.pdf"), plot = p4, width = 150, height = 120, units = "mm")

#######################
## Sr isotopes

# path files
speciation_time_file <- paste0("output/Sr_corals_EBD_Corr_speciation_times.log")
speciation_rate_file <- paste0("output/Sr_corals_EBD_Corr_speciation_rates.log")
extinction_time_file <- paste0("output/Sr_corals_EBD_Corr_extinction_times.log")
extinction_rate_file <- paste0("output/Sr_corals_EBD_Corr_extinction_rates.log")

rates <- processDivRates(speciation_time_log = speciation_time_file,
                         speciation_rate_log = speciation_rate_file,
                         extinction_time_log = extinction_time_file,
                         extinction_rate_log = extinction_rate_file,
                         burnin = 0.25,
                         summary = "median")

ex_df <- read.table(extinction_time_file, sep ="\t", header = TRUE)
#pdf("output/EBD_Corr.pdf")
#par(mfrow=c(2,2))
#plotDivRates(rev_out, predictor.ages=co2_age, predictor.var=co2, use.geoscale=TRUE)
#dev.off()

# the Sr values as a reference in our plot
Sr <- c(0.7090184, 0.7087190, 0.7081701, 0.7078171, 0.7077420, 0.7077501, 0.7078148, 0.7076480, 0.7074318, 0.7073783, 0.7074178,
        0.7072761, 0.7073982, 0.7073793, 0.7071852, 0.7069474, 0.7069454, 0.7072605, 0.7071754, 0.7075384, 0.7077311, 0.7078536,
        0.7077222, 0.7076284, 0.7078606)


MAX_VAR_AGE <- 250
NUM_INTERVALS <- length(Sr)
Sr_age <- MAX_VAR_AGE * (1:NUM_INTERVALS) / NUM_INTERVALS
predictor.ages <- Sr_age
predictor.var <- Sr



Sr_df <- tibble("Sr" = Sr, "Sr_age" = Sr_age)


p5 <- plotDivRates(rates)
p5
## We can add the Sr measurements on top of the ggplot object ## second axis adden or in powerpoint, values around 0.700
Sr_df$scaled <- Sr_df$Sr*100 - 70.55 ## Change the relative plotting scale between rate units and co2 units (ppm)
coeff3 <- max(Sr_df$Sr)*100 - 70.55
p6 <- p5 +
  geom_line(data=Sr_df, mapping = aes(x = Sr_age, y = scaled), inherit.aes = FALSE) +
  scale_y_continuous(
    name = "Rate", ## First y-axis (left)
   # sec.axis = sec_axis(trans = ~.) ## Second y-axis (right)
  )
p6

## Sr plot
plot(Sr_age, Sr, type = "l", lty = "dashed", ylab = "87Sr/86Sr", xlab = "Age (Ma)", xlim = c(max(Sr_age), min(Sr_age)))

