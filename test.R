## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 1. Build the equilibrium ecosystem ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
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

load("EEM_eq.RData")

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 2. check if equilibrium values are between upper and lower bounds  ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

source("select_EEM_outputs.R")

EEM_outputs_eq <- select_EEM_outputs(ensemble = EEM_eq,
                                     target_lower = lower_100m2,
                                     target_upper = upper_100m2)

#it works

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 3. plot projections native system, scaled  ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
EEM_eq_short <- EEM_eq[1:5]

source("adapted_calculate_projections.R")
source("adapted_ode_solve.R")

mean_eq <- (upper_100m2 + lower_100m2) / 2

eq_projections_scaled <-
  adapted_calculate_projections(EEM_eq_short,
                                initial_condition = rep(1, length(species_eq)),
                                t_window = c(0, 30),
                                scaled = TRUE,
                                species_names = species_eq,
                                multiplier = divider)

source("adapted_plot_projections.R")

plot_eq_scaled <- adapted_plot_projections(projections = eq_projections_scaled,
                                           title = " Projections native
                                           system, scaled")

#it works

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 4. plot projections native system, not scaled  ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

mean_eq <- (upper_100m2 + lower_100m2) / 2

eq_projections <- adapted_calculate_projections(EEM_eq_short,
                                                initial_condition = mean_eq,
                                                t_window = c(0, 30),
                                                scaled = FALSE,
                                                species_names = species_eq,
                                                multiplier = divider)

source("adapted_plot_projections.R")

plot_eq <- adapted_plot_projections(projections = eq_projections,
                                    title = "Projections native system")

#it works

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 5. add sihek  ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

source("add_introduced_species.R")

EEM_sihek <-
  add_introduced_species(native_parameters = EEM_eq_short,
                         introduced_lower_bound_growth_rate = 0,
                         introduced_upper_bound_growth_rate = 1.1,
                         introduced_self_sign = -1,
                         introduced_row_signs =
                           c(1, 1, 1, 1, 1, 1, 1, 1),
                         introduced_col_signs =
                           c(-1, -1, -1, -1, -1, -1, -1, 0),
                         introduced_k = 0.002/divider) #40/200ha = 0.002/100m^2

source("add_species_names.R")

EEM_sihek <- add_species_names(EEM_sihek,
                               species_names = c("sihek",
                                                 "seabirds",
                                                 "terrestrial crabs",
                                                 "carnivorous crabs",
                                                 "cane spiders",
                                                 "geckos",
                                                 "cockroaches",
                                                 "terrestrial arthropods",
                                                 "native trees"))

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 6. check if equilibrium values are between upper and lower bounds ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

EEM_outputs_sihek <- select_EEM_outputs(ensemble = EEM_sihek,
                                        target_lower = c(0.0015/divider,
                                                         lower_100m2),
                                        target_upper = c(0.0025/divider,
                                                         upper_100m2),
                                        mode = "disturbed")

#it works

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 7. plot projections with sihek but without artificial recruitment ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

mean_sihek <- c(0.00045/divider, mean_eq) #9/200ha = 0.00045/100m^2

species_sihek <- c("sihek", species_eq)

source("adapted_calculate_projections.R")
source("adapted_ode_solve.R")
sihek_projections <-
  adapted_calculate_projections(EEM_sihek,
                                initial_condition = mean_sihek,
                                t_window = c(0, 30),
                                scaled = FALSE,
                                species_names = species_sihek,
                                multiplier = divider)

plot_sihek <- adapted_plot_projections(projections = sihek_projections,
                                      title = "Projections with sihek")

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 8. plot projections with sihek with artificial recruitment ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

mean_sihek_artrec <- c(0, mean_eq)

sihek_projections_artrec <-
  adapted_calculate_projections(EEM_sihek,
                                initial_condition = mean_sihek_artrec,
                                t_window = c(0, 30),
                                scaled = FALSE,
                                species_names = species_sihek,
                                init_intervention_amount = 0.00045/divider, #9/200ha = 0.00045/100m^2
                                init_intervention_timepoints = c(2,3),
                                sustain_intervention_amount = 0.00025/divider, #5/200ha = 0.00025/100m^2
                                sustain_intervention_timepoints = c(3, 5, 7,
                                                                    9, 11, 13,
                                                                    15, 17, 19),
                                sustain_intervention_threshold = 0.001/divider, #20/200ha = 0.001/100m^2
                                time_step_len = 0.01,
                                multiplier = divider)


plot_sihek_artrec <- adapted_plot_projections(
  projections = sihek_projections_artrec,
  title = "Projections with sihek & artificial recruitment")

#add red intercept line
threshholds <- data.frame(yintercepts = c(0.001, 0.002),
                          species = "sihek")

library(ggplot2)
plot_sihek_artrec <- plot_sihek_artrec +
  geom_hline(data = threshholds, aes(yintercept = yintercepts),
             color = "red", linewidth = 0.2)

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 9. compare the three plots ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

library(gridExtra)

grid.arrange(plot_eq, plot_sihek, plot_sihek_artrec, ncol = 3)

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 10. add coconut palms to the equilibrium system ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

source("add_introduced_species.R")

EEM_palm <-
  add_introduced_species(native_parameters = EEM_eq_short,
                         introduced_lower_bound_growth_rate = 0,
                         introduced_upper_bound_growth_rate = 3,
                         introduced_self_sign = -1,
                         introduced_row_signs = c(1, -1, -1, 1, 1, 1, 1, -1),
                         introduced_col_signs = rep(-1,8),
                         introduced_k = 106.5/divider) #2130000/200ha = 106.5/100m^2

EEM_palm <- add_species_names(EEM_palm,
                              species_names = c("palm trees",
                                                "seabirds",
                                                "terrestrial crabs",
                                                "carnivorous crabs",
                                                "cane spiders",
                                                "geckos",
                                                "cockroaches",
                                                "terrestrial arthropods",
                                                "native trees"))

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 11. check if equilibrium values are between upper and lower bounds  ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

source("select_EEM_outputs.R")

selected_EEM_palms <- select_EEM_outputs(ensemble = EEM_palm,
                                        target_lower = c(106.4/divider,
                                                         lower_100m2),
                                        target_upper = c(106.6/divider,
                                                         upper_100m2),
                                        mode = "disturbed")

#it works

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 12. plot projections with palm trees and with artificial control ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

mean_palm <- c(50/divider, mean_eq) #1000000/200ha = 50/100m^2

species_names_palm <- c("palm trees", species_eq)

source("adapted_calculate_projections.R")
source("adapted_ode_solve.R")

palm_projections <-
  adapted_calculate_projections(EEM_palm,
                                initial_condition = mean_palm,
                                t_window = c(0, 30),
                                scaled = FALSE,
                                species_names = species_names_palm,
                                mode = "removal",
                                init_intervention_amount = -75/divider, #-1500000/200ha = -75/100m^2
                                init_intervention_timepoints = c(20,21),
                                sustain_intervention_amount = -10/divider, #-100000/200ha = -5/100m^2
                                sustain_intervention_timepoints = c(21.5, 22,
                                                                    22.5, 23,
                                                                    23.5, 24,
                                                                    24.5, 25,
                                                                    25.5, 26,
                                                                    26.5, 27,
                                                                    27.5, 28,
                                                                    28.5, 29),
                                sustain_intervention_threshold = 0,
                                intro_species_index = 1,
                                time_step_len = 0.01,
                                multiplier = divider)

plot_palm <- adapted_plot_projections(
  projections = palm_projections,
  title = "Projections with palm trees & artificial control")

#add red intercept line
threshholds <- data.frame(yintercepts = c(0, 106.5),
                          species = "palm trees")

library(ggplot2)
plot_palm <- plot_palm +
  geom_hline(data = threshholds, aes(yintercept = yintercepts),
             color = "red", linewidth = 0.2)

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 13. Merge palm and sihek introductions ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

source("merge_introductions.R")
source("adapted_calculate_projections.R")
source("adapted_ode_solve.R")

signs_sihek_palm <- matrix(c(NA, -1,
                             -1, NA),
                           nrow = 2,
                           ncol = 2,
                           byrow = TRUE)

EEM_merged <- merge_introductions(EEM_intros = pairlist(EEM_sihek,
                                                         EEM_palm),
                                   sign_interaction_intros = signs_sihek_palm,
                                   mode = "recycled")

species_names_merged <- c("sihek", species_names_palm)

merged_projections <-
  adapted_calculate_projections(EEM_merged,
                                initial_condition = c(0, mean_palm),
                                t_window = c(0, 30),
                                scaled = FALSE,
                                species_names = species_names_merged,
                                mode = c(
                                  "recruitment",
                                  "removal"),
                                init_intervention_amount = c(
                                  0.00045/divider, #9/200ha = 0.00045/100m^2
                                  -75/divider), #-1500000/200ha = -75/100m^2
                                init_intervention_timepoints = list(
                                  c(2,3),
                                  c(20,21)),
                                sustain_intervention_amount = c(
                                  0.00025/divider,
                                  -10/divider), #-100000/200ha = -5/100m^2
                                sustain_intervention_timepoints = list(
                                  c(3, 5, 7,
                                    9, 11, 13,
                                    15, 17, 19),
                                  c(21.5, 22,
                                    22.5, 23,
                                    23.5, 24,
                                    24.5, 25,
                                    25.5, 26,
                                    26.5, 27,
                                    27.5, 28,
                                    28.5, 29)),
                                sustain_intervention_threshold = c(
                                  0.001/divider, #20/200ha = 0.001/100m^2
                                  0),
                                intro_species_index = c(
                                  1,
                                  2),
                                time_step_len = 0.01,
                                multiplier = divider)

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 14. check if equilibrium values are between upper and lower bounds ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

source("select_EEM_outputs.R")

selected_EEM_merged <- select_EEM_outputs(ensemble = EEM_merged,
                                          target_lower = c(0.0015/divider,
                                                           106.4/divider,
                                                           lower_100m2),
                                          target_upper = c(0.0025/divider,
                                                           106.6/divider,
                                                           upper_100m2),
                                          mode = "disturbed",
                                          n_intro = 2)

#it works

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 15. plot projections with palm trees & sihek and with artificial control ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

plot_merged <- adapted_plot_projections(
  projections = merged_projections,
  title = "Projections with sihek & palm trees & artificial control")


threshholds <- data.frame(yintercepts = c(0.001, 106.5, #k
                                          0.002, 0), #threshold
                          species = c("sihek", "palm trees"))

library(ggplot2)
plot_merged <- plot_merged +
  geom_hline(data = threshholds, aes(yintercept = yintercepts),
             color = "red", linewidth = 0.2)

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 16. Normalise abundances of native species, merged system ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
merged_norm_projections <- merged_projections
merged_norm_projections[merged_norm_projections$species %in% species_eq,]$pop <-
  merged_norm_projections[
    merged_norm_projections$species %in% species_eq,]$pop / eq_projections$pop

plot_merged_norm <- adapted_plot_projections(
  projections = merged_norm_projections,
  title = "Projections with sihek & palm trees & artificial control, normalised")

threshholds <- data.frame(yintercepts = c(0.001, 106.5, #k
                                          0.002, 0), #threshold
                          species = c("sihek", "palm trees"))

library(ggplot2)
plot_merged_norm <- plot_merged_norm +
  geom_hline(data = threshholds, aes(yintercept = yintercepts),
             color = "red", linewidth = 0.2)

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 17. Normalise abundances of native species, merged system vs palm system ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

palm_norm_projections <- merged_projections
palm_norm_projections[palm_norm_projections$species %in% species_eq,]$pop <-
  palm_norm_projections[palm_norm_projections$species %in% species_eq,]$pop /
  palm_projections[palm_projections$species %in% species_eq,]$pop

plot_palm_norm <- adapted_plot_projections(
  projections = palm_norm_projections,
  title = "Projections with sihek & palm trees & artificial control, normalised")

threshholds <- data.frame(yintercepts = c(0, 106.5),
                          species = "palm trees")

library(ggplot2)
plot_palm_norm <- plot_palm_norm +
  geom_hline(data = threshholds, aes(yintercept = yintercepts),
             color = "red", linewidth = 0.2)

## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## 18. plot only the native species from the normalised merged system ####
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

plot_eq_merged_norm_palm <- adapted_plot_projections(
  projections = palm_norm_projections[palm_norm_projections$species %in% species_eq,],
  title = "merged vs palm system, native species",
  scaled = FALSE)
