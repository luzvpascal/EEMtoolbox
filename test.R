## =============================================================================
## 1. Build the equilibrium ecosystem
## =============================================================================
interaction_matrix_nosihek <- matrix(c(-1, 0, -1, 0, 0, 0, 0, 1,
                                       0, -1, -1, 0, 0, 0, 0, 1,
                                       1, 1, -1, 0, 0, 1, 1, 0,
                                       0, 0, 0, -1, 1, 1, 1, 0,
                                       0, 0, 0, -1, -1, 0, 1, 0,
                                       0, 0, -1, -1, 0, -1, 0, 1,
                                       0, 0, -1, -1, -1, 0, -1, 1,
                                       1, -1, -1, 0, 0, 0, 1, -1),
                                     ncol = 8, nrow = 8, byrow = TRUE)

upper_bounds_growth_rate_nosihek = c(1.1, 1.5, 1.5, 0.39, 0.49, 3.0, 3.0, 3.0)

test_EEM_eq <- EEM(interaction_matrix_nosihek,
                   upper_bounds_growth_rate = upper_bounds_growth_rate_nosihek,
                   algorithm = "SMC-EEM",
                   disc_func = function(data) {
                     adapted_discrepancy_continuous_sum(
                       data,
                       target_lower = c(500/200000000, #1
                                        10000/200000000, #2
                                        500/200000000, #3
                                        1000/200000000, #4
                                        10000/200000000, #5
                                        1000000/200000000, #6
                                        1000000/200000000, #7
                                        100000/200000000), #8
                       target_upper = c(5000/200000000, #1
                                        90000/200000000, #2
                                        1500/200000000, #3
                                        10000/200000000, #4
                                        90000/200000000, #5
                                        90000000/200000000, #6
                                        90000000/200000000, #7
                                        900000/200000000)) #8
                   },
                   n_ensemble = 5000,
                   n_cores = 3)

save(test_EEM_eq, file = "test_EEM_eq.RData")

source("add_species_names.R")

test_EEM_eq <- add_species_names(test_EEM_eq,
                                 c("seabirds",
                                   "terrestrial crabs",
                                   "carnivorous crabs",
                                   "cane spiders",
                                   "geckos",
                                   "cockroaches",
                                   "terrestrial arthropods",
                                   "native trees"))

save(test_EEM_eq, file = "test_EEM_eq.RData")

## =============================================================================
## 2. check whether the equilibrium values are between the upper and lower bounds
## =============================================================================

test_arg <- args_function(interaction_matrix = interaction_matrix_nosihek,
                          upper_bounds_growth_rate =
                            upper_bounds_growth_rate_nosihek,
                          model = "GLV")

source("select_EEM_outputs.R")

test_target_lower <- c(500/200000000, #1
                       10000/200000000, #2
                       500/200000000, #3
                       1000/200000000, #4
                       10000/200000000, #5
                       1000000/200000000, #6
                       1000000/200000000, #7
                       100000/200000000) #8

test_target_upper <- c(5000/200000000, #1
                       90000/200000000, #2
                       1500/200000000, #3
                       10000/200000000, #4
                       90000/200000000, #5
                       90000000/200000000, #6
                       90000000/200000000, #7
                       900000/200000000) #8

test_EEM_outputs <- select_EEM_outputs(ensemble = test_EEM_eq,
                                       target_lower = test_target_lower, #8
                                       target_upper = test_target_upper,
                                       sim_args = test_arg)

#it works

## =============================================================================
## 3. plot projections without sihek
## =============================================================================
test_EEM_eq_short <- test_EEM_eq[1:5]

source("adapted_calculate_projections.R")
source("adapted_ode_solve.R")

mean_eq <- (test_target_upper + test_target_lower) / 2

species_names_eq <- colnames(test_EEM_eq_short[[1]]$interaction_matrix)

test_eq_projections <-
  adapted_calculate_projections(test_EEM_eq_short,
                                initial_condition = mean_eq,
                                t_window = c(0, 20),
                                scaled = FALSE,
                                species_names = species_names_eq)

library(ggplot2)
library(dplyr)

test_abundance_eq <- group_by(test_eq_projections,
                              time,
                              species)
test_abundance_eq <- summarise(test_abundance_eq,
                               median_pop = median(pop),
                               upper = quantile(pop, 0.975),
                               lower = quantile(pop, 0.025))

#basic plot
test_plot_eq <- ggplot(test_abundance_eq) +
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
  ggtitle("Projections without sihek") +
  ylim(c(0, 1)) +
  theme_bw()

## =============================================================================
## 4. add sihek
## =============================================================================

source("add_introduced_species.R")

test_EEM_intro <-
  add_introduced_species(native_parameters = test_EEM_eq_short,
                         introduced_lower_bound_growth_rate = 0,
                         introduced_upper_bound_growth_rate = 1.1,
                         introduced_self_sign = -1,
                         introduced_row_signs = c(1, 1, 1, 1, 1, 1, 1, 1),
                         introduced_col_signs = c(-1, -1, -1, -1, -1, -1, -1, 0),
                         introduced_k = 40/200000000)

## =============================================================================
## 5. plot projections with sihek but without artificial recruitment
## =============================================================================

mean_intro <- c(9/200000000, mean_eq)

species_names_intro <- c("sihek",
                         colnames(test_EEM_eq_short[[1]]$interaction_matrix))

test_intro_projections <-
  adapted_calculate_projections(test_EEM_intro,
                                initial_condition = mean_intro,
                                t_window = c(0, 20),
                                scaled = FALSE,
                                species_names = species_names)

library(ggplot2)
library(dplyr)

test_abundance_intro <- group_by(test_intro_projections,
                                 time,
                                 species)
test_abundance_intro <- summarise(test_abundance_intro,
                                  median_pop = median(pop),
                                  upper = quantile(pop, 0.975),
                                  lower = quantile(pop, 0.025))

#basic plot
test_plot_intro <- ggplot(test_abundance_intro) +
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
## 6. plot projections with sihek with artificial recruitment
## =============================================================================

test_intro_projections_artrec <-
  adapted_calculate_projections(test_EEM_intro,
                                initial_condition = mean_intro,
                                t_window = c(0, 20),
                                scaled = FALSE,
                                species_names = species_names_intro,
                                init_release_amount = 9/200000000,
                                init_release_timepoints = 1,
                                sustain_release_amount = 5/200000000,
                                sustain_release_timepoints = c(3, 5, 7, 9, 11,
                                                               13, 15, 17, 19),
                                sustain_release_threshold = 20/200000000,
                                introduced_species_index = 1,
                                time_step_len = 0.01)

library(ggplot2)
library(dplyr)

test_abundance_intro_artrec <- group_by(test_intro_projections_artrec,
                                        time,
                                        species)
test_abundance_intro_artrec <- summarise(test_abundance_intro_artrec,
                                         median_pop = median(pop),
                                         upper = quantile(pop, 0.975),
                                         lower = quantile(pop, 0.025))

#add red intercept line
threshholds <- data.frame(yintercepts = c(20/200000000, 40/200000000),
                          species = "sihek")

#basic plot
test_plot_intro_artrec <- ggplot(test_abundance_intro_artrec) +
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
## 7. compare the three plots
## =============================================================================

library(gridExtra)

grid.arrange(test_plot_eq, test_plot_intro, test_plot_intro_artrec, ncol = 3)

## =============================================================================
## 8.1. change scale of sihek to see impact better
## =============================================================================

test_EEM_intro_high <-
  add_introduced_species(native_parameters = test_EEM_eq_short,
                         introduced_lower_bound_growth_rate = 0,
                         introduced_upper_bound_growth_rate = 1.1,
                         introduced_self_sign = -1,
                         introduced_row_signs = c(1, 1, 1, 1, 1, 1, 1, 1),
                         introduced_col_signs = c(-1, -1, -1, -1, -1, -1, -1, 0),
                         introduced_k = 4)

test_intro_projections_high <-
  adapted_calculate_projections(test_EEM_intro_high,
                                initial_condition = c(1, mean_eq),
                                t_window = c(0, 20),
                                scaled = FALSE,
                                species_names = species_names)

library(ggplot2)
library(dplyr)

test_abundance_intro_high <- group_by(test_intro_projections_high,
                                      time,
                                      species)
test_abundance_intro_high <- summarise(test_abundance_intro_high,
                                       median_pop = median(pop),
                                       upper = quantile(pop, 0.975),
                                       lower = quantile(pop, 0.025))

#basic plot
test_plot_intro_high <- ggplot(test_abundance_intro_high) +
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
  ggtitle("Projections with sihek, high K") +
  theme_bw()

## =============================================================================
## 8.2. change scale of sihek to see impact better + artificial recruitment
## =============================================================================

test_intro_projections_artrec_high <-
  adapted_calculate_projections(test_EEM_intro_high,
                                initial_condition = c(1, mean_eq),
                                t_window = c(0, 20),
                                scaled = FALSE,
                                species_names = species_names,
                                init_release_amount = 1,
                                init_release_timepoints = 1,
                                sustain_release_amount = 1.5,
                                sustain_release_timepoints = c(3, 5, 7, 9, 11,
                                                               13, 15, 17, 19),
                                sustain_release_threshold = 3,
                                introduced_species_index = 1,
                                time_step_len = 0.01)

library(ggplot2)
library(dplyr)

test_abundance_intro_artrec_high <- group_by(test_intro_projections_artrec_high,
                                             time,
                                             species)
test_abundance_intro_artrec_high <- summarise(test_abundance_intro_artrec_high,
                                              median_pop = median(pop),
                                              upper = quantile(pop, 0.975),
                                              lower = quantile(pop, 0.025))

#add red intercept line
threshholds_high <- data.frame(yintercepts = c(3, 4),
                               species = "sihek")

#basic plot
test_plot_intro_artrec_high <- ggplot(test_abundance_intro_artrec_high) +
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
  ggtitle("Projections with sihek & artificial recruitment, high K") +
  theme_bw()

## =============================================================================
## 8.3. compare the two plots
## =============================================================================

library(gridExtra)

grid.arrange(test_plot_intro_high, test_plot_intro_artrec_high, ncol = 2)
