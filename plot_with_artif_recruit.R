load("parameters_GLV_plussishek_adapteddiscfunc_200.RData")

##### implement artificial recruitement #####

species_names <- colnames(
  parameters_GLV_plussishek_adapteddiscfunc_200[[1]]$interaction_matrix)

source("adapted_calculate_projections.R")
source("adapted_ode_solve.R")

projections_200_nonscaled_upperbounds_recruitment <-
  adapted_calculate_projections(parameters_GLV_plussishek_adapteddiscfunc_200,
                                initial_condition = c(6/20000000,
                                                      5000/20000000, #1
                                                      90000/20000000, #2
                                                      1500/20000000, #3
                                                      10000/20000000, #4
                                                      90000/20000000, #5
                                                      90000000/20000000, #6
                                                      90000000/20000000, #7
                                                      900000/20000000), #8
                                t_window = c(0, 10),
                                derivative = EEMtoolbox::derivative_func,
                                scaled = FALSE,
                                species_names = species_names)

save(projections_200_nonscaled_upperbounds_recruitment,
     file = "projections_200_nonscaled_upperbounds_recruitment.RData")

##### plot with own way #####

library(ggplot2)
library(dplyr)
abundance <- group_by(projections_200_nonscaled_upperbounds_recruitment,
                      time,
                      species)
abundance <- summarise(abundance,
                       median_pop = median(pop),
                       upper = quantile(pop, 0.975),
                       lower = quantile(pop, 0.025))

#basic plot
ggplot(abundance) +
  guides(fill = guide_legend(title = "Species"),
         color = guide_legend(title = "Species")) +
  xlab("Time") +
  ylab("Abundance") +
  facet_wrap( ~ species) +
  geom_ribbon(aes(x = time,
                  ymin = lower,
                  ymax = upper,
                  color = species,
                  fill = species),
              alpha = 0.2) +
  geom_line(aes(x = time, y = median_pop, color = species),
            linewidth = 1.5) +
  theme_bw()
