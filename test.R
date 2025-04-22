## =============================================================================
## 1. Build the equilibrium ecosystem ####
## =============================================================================
species_eq <- c("seabirds",
                "terrestrial crabs",
                "carnivorous crabs",
                "cane spiders",
                "geckos",
                "cockroaches",
                "terrestrial arthropods",
                "native trees")

interaction_matrix_eq <- matrix(c(-1, 0,-1, 0, 0, 0, 0, 1,
                                   0,-1,-1, 0, 0, 0, 0, 1,
                                   1, 1,-1, 0, 0, 1, 1, 0,
                                   0, 0, 0,-1, 1, 1, 1, 0,
                                   0, 0, 0,-1,-1, 0, 1, 0,
                                   0, 0,-1,-1, 0,-1, 0, 1,
                                   0, 0,-1,-1,-1, 0,-1, 1,
                                   1,-1,-1, 0, 0, 0, 1,-1),
                                ncol = 8, nrow = 8,
                                dimnames = list(species_eq, species_eq),
                                byrow = TRUE)

upper_gr_eq <- c(1.1, 1.5, 1.5, 0.39, 0.49, 3.0, 3.0, 3.0)

divider <- 100000

lower_100m2 <- c(30, #"seabirds",
                 200, #"terrestrial crabs",
                 20, #"carnivorous crabs",
                 100, #"cane spiders",
                 300, #"geckos",
                 900, #"cockroaches",
                 50000, #"terrestrial arthropods",
                 30)/divider #"native trees"

upper_100m2 <- c(70, #"seabirds",
                 400, #"terrestrial crabs",
                 40, #"carnivorous crabs",
                 400, #"cane spiders",
                 500, #"geckos",
                 1500, #"cockroaches",
                 150000, #"terrestrial arthropods",
                 50)/divider #"native trees" ---> all per 100 m^2

EEM_eq <- EEM(interaction_matrix_eq,
              upper_bounds_growth_rate = upper_gr_eq,
              algorithm = "SMC-EEM",
              disc_func = function(data) {
                adapted_discrepancy_continuous_sum(
                  data,
                  target_lower = c(30, #"seabirds",
                                   200, #"terrestrial crabs",
                                   20, #"carnivorous crabs",
                                   100, #"cane spiders",
                                   300, #"geckos",
                                   900, #"cockroaches",
                                   50000, #"terrestrial arthropods",
                                   30)/100000, #"native trees",
                  target_upper = c(70, #"seabirds",
                                   400, #"terrestrial crabs",
                                   40, #"carnivorous crabs",
                                   400, #"cane spiders",
                                   500, #"geckos",
                                   1500, #"cockroaches",
                                   150000, #"terrestrial arthropods",
                                   50)/100000) #"native trees"
              },
              n_ensemble = 5000,
              n_cores = 3)

save(EEM_eq, file = "EEM_eq.RData")

source("add_species_names.R")

EEM_eq <- add_species_names(EEM_eq,
                            c("seabirds",
                              "terrestrial crabs",
                              "carnivorous crabs",
                              "cane spiders",
                              "geckos",
                              "cockroaches",
                              "terrestrial arthropods",
                              "native trees"))

save(EEM_eq, file = "test_EEM_eq.RData")

## =============================================================================
## 2. check if equilibrium values are between upper and lower bounds  ####
## =============================================================================

source("select_EEM_outputs.R")

EEM_outputs_eq <- select_EEM_outputs(ensemble = EEM_eq,
                                     target_lower = lower_100m2,
                                     target_upper = upper_100m2)

#it works

## =============================================================================
## 3. plot projections native system, scaled  ####
## =============================================================================
EEM_eq_short <- EEM_eq[1:5]

source("adapted_calculate_projections.R")
source("adapted_ode_solve.R")

mean_eq <- (upper_100m2 + lower_100m2) / 2

eq_projections_scaled <-
  adapted_calculate_projections(EEM_eq_short,
                                initial_condition = rep(1, length(species_eq)),
                                t_window = c(0, 20),
                                scaled = TRUE,
                                species_names = species_eq,
                                multiplier = divider)

source("adapted_plot_projections.R")

plot_eq_scaled <- adapted_plot_projections(
  data = eq_projections_scaled,
  title = " Projections native system, scaled")

#it works

## =============================================================================
## 4. plot projections native system, not scaled  ####
## =============================================================================

mean_eq <- (upper_100m2 + lower_100m2) / 2

eq_projections <- adapted_calculate_projections(EEM_eq_short,
                                                initial_condition = mean_eq,
                                                t_window = c(0, 20),
                                                scaled = FALSE,
                                                species_names = species_eq,
                                                multiplier = divider)

source("adapted_plot_projections.R")

plot_eq <- adapted_plot_projections(data = eq_projections,
                                    title = "Projections native system",
                                    multiplier = divider)

#it works

## =============================================================================
## 5. add sihek  ####
## =============================================================================

source("add_introduced_species.R")

EEM_sihek <-
  add_introduced_species(native_parameters = EEM_eq,
                         introduced_lower_bound_growth_rate = 0,
                         introduced_upper_bound_growth_rate = 1.1,
                         introduced_self_sign = -1,
                         introduced_row_signs =
                           c(1, 1, 1, 1, 1, 1, 1, 1),
                         introduced_col_signs =
                           c(-1, -1, -1, -1, -1, -1, -1, 0),
                         introduced_k = 0.002/divider) #40/200ha = 0.002/100m^2

## =============================================================================
## 6. check if equilibrium values are between upper and lower bounds  ####
## =============================================================================

EEM_outputs_sihek <- select_EEM_outputs(ensemble = EEM_sihek,
                                        target_lower = c(0.001/divider,
                                                         lower_100m2),
                                        target_upper = c(0.003/divider,
                                                         upper_100m2),
                                        mode = "disturbed")

r <- EEM_sihek[[1]]$growthrates
A <- EEM_sihek[[1]]$interaction_matrix
solve(A, -r)

## =============================================================================
## 7. plot projections with sihek but without artificial recruitment
## =============================================================================

mean_sihek <- c(0.00045/divider, mean_eq) #9/200ha = 0.00045/100m^2

species_sihek <- c("sihek", species_eq)

source("adapted_calculate_projections.R")
source("adapted_ode_solve.R")
sihek_projections <-
  adapted_calculate_projections(EEM_sihek[EEM_outputs_sihek$outputs_selected],
                                initial_condition = mean_sihek,
                                t_window = c(0, 20),
                                scaled = FALSE,
                                species_names = species_sihek,
                                multiplier = divider)

library(ggplot2)
library(dplyr)

abundance_intro <- group_by(sihek_projections,
                            time,
                            species)
abundance_intro <- summarise(abundance_intro,
                             median_pop = median(pop),
                             upper = quantile(pop, 0.975),
                             lower = quantile(pop, 0.025))

#basic plot
plot_sihek <- ggplot(abundance_intro) +
  guides(fill = guide_legend(title = "Species"),
         color = guide_legend(title = "Species")) +
  xlab("Time") +
  ylab("Abundance") +
  facet_wrap( ~ species, scales = "free") +
  geom_ribbon(aes(x = time,
                  ymin = lower,
                  ymax = upper,
                  color = species,
                  fill = species),
              alpha = 0.1,
              linewidth = 0.2) +
  geom_line(aes(x = time, y = median_pop, color = species),
            linewidth = 0.8) +
  ggtitle("Projections with sihek") +
  theme_bw()

## =============================================================================
## 8. plot projections with sihek with artificial recruitment
## =============================================================================

sihek_projections_artrec <-
  adapted_calculate_projections(EEM_sihek[EEM_outputs_sihek$outputs_selected],
                                initial_condition = mean_sihek,
                                t_window = c(0, 20),
                                scaled = FALSE,
                                species_names = species_sihek,
                                init_intervention_amount = 0.00045/divider, #9/200ha = 0.00045/100m^2
                                init_intervention_timepoints = 1,
                                sustain_intervention_amount = 0.00025/divider, #5/200ha = 0.00025/100m^2
                                sustain_intervention_timepoints = c(3, 5, 7,
                                                                    9, 11, 13,
                                                                    15, 17, 19),
                                sustain_intervention_threshold = 0.001/divider, #20/200ha = 0.001/100m^2
                                introduced_species_index = 1,
                                time_step_len = 0.01,
                                multiplier = divider)

library(ggplot2)
library(dplyr)

abundance_sihek_artrec <- group_by(sihek_projections_artrec,
                                   time,
                                   species)
abundance_sihek_artrec <- summarise(abundance_sihek_artrec,
                                    median_pop = median(pop),
                                    upper = quantile(pop, 0.975),
                                    lower = quantile(pop, 0.025))

#add red intercept line
threshholds <- data.frame(yintercepts = c(0.001, 0.002),
                          species = "sihek")

#basic plot
plot_sihek_artrec <- ggplot(abundance_sihek_artrec) +
  guides(fill = guide_legend(title = "Species"),
         color = guide_legend(title = "Species")) +
  xlab("Time") +
  ylab("Abundance") +
  facet_wrap( ~ species, scales = "free") +
  geom_ribbon(aes(x = time,
                  ymin = lower,
                  ymax = upper,
                  color = species,
                  fill = species),
              alpha = 0.1,
              linewidth = 0.2) +
  geom_line(aes(x = time, y = median_pop, color = species),
            linewidth = 0.8) +
  geom_hline(data = threshholds, aes(yintercept = yintercepts),
             color = "red", linewidth = 0.2) +
  ggtitle("Projections with sihek & artificial recruitment") +
  theme_bw()

## =============================================================================
## 9. compare the three plots
## =============================================================================

library(gridExtra)

grid.arrange(plot_eq, plot_sihek, plot_sihek_artrec, ncol = 3)

## =============================================================================
## 10. add coconut palms to the equilibrium system
## =============================================================================

source("add_introduced_species.R")

EEM_palm <-
  add_introduced_species(native_parameters = EEM_eq,
                         introduced_lower_bound_growth_rate = 0,
                         introduced_upper_bound_growth_rate = 3,
                         introduced_self_sign = -1,
                         introduced_row_signs = c(1, -1, -1, 1, 1, 1, 1, -1),
                         introduced_col_signs = rep(-1,8),
                         introduced_k = 106.5/divider) #2130000/200ha = 106.5/100m^2

## =============================================================================
## 11. check if equilibrium values are between upper and lower bounds  ####
## =============================================================================

source("select_EEM_outputs.R")

EEM_outputs_palms <- select_EEM_outputs(ensemble = EEM_palm,
                                        target_lower = c(103.5/divider,
                                                         lower_100m2),
                                        target_upper = c(110.5/divider,
                                                         upper_100m2),
                                        mode = "disturbed")

r <- EEM_palm[[5]]$growthrates
A <- EEM_palm[[5]]$interaction_matrix
solve(A, -r)

## =============================================================================
## 9. plot projections with palm trees and with artificial control
## =============================================================================

initial_palm <- c(0, mean_eq)

species_names_palm <- c("palm trees", species_eq)

source("adapted_calculate_projections.R")
source("adapted_ode_solve.R")

palm_projections <-
  adapted_calculate_projections(EEM_palm[EEM_outputs_palms$outputs_selected],
                                initial_condition = initial_palm,
                                t_window = c(0, 20),
                                scaled = FALSE,
                                species_names = species_names_palm,
                                mode = "removal",
                                init_intervention_amount = 5/divider, #1000/200ha = 50/100m^2
                                init_intervention_timepoints = 1,
                                sustain_intervention_amount = -100/divider, #-2000000/200ha = -100/100m^2
                                sustain_intervention_timepoints = 15:16,
                                sustain_intervention_threshold = 0,
                                introduced_species_index = 1,
                                time_step_len = 0.01,
                                multiplier = divider)

save(test_palm_projections, file = "test_palm_projections.RData")

abundance_palm <- group_by(palm_projections,
                           time,
                           species)
abundance_palm <- summarise(abundance_palm,
                            median_pop = median(pop),
                            upper = quantile(pop, 0.975),
                            lower = quantile(pop, 0.025))

#basic plot
plot_palm <- ggplot(abundance_palm) +
  guides(fill = guide_legend(title = "Species"),
         color = guide_legend(title = "Species")) +
  xlab("Time") +
  ylab("Abundance") +
  facet_wrap( ~ species, scales = "free") +
  geom_ribbon(aes(x = time,
                  ymin = lower,
                  ymax = upper,
                  color = species,
                  fill = species),
              alpha = 0.1,
              linewidth = 0.2) +
  geom_line(aes(x = time, y = median_pop, color = species),
            linewidth = 0.8) +
  ggtitle("Projections with palm trees") +
  theme_bw()


