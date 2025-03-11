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

parameters_GLV_nosihek_adapteddiscfunc_200 <- EEM(
  interaction_matrix_nosihek,
  upper_bounds_growth_rate = upper_bounds_growth_rate_nosihek,
  algorithm = "SMC-EEM",
  disc_func = function(data) {
    adapted_discrepancy_continuous_sum(data,
                                       target_lower = c(500/20000000, #1
                                                        10000/20000000, #2
                                                        500/20000000, #3
                                                        1000/20000000, #4
                                                        10000/20000000, #5
                                                        1000000/20000000, #6
                                                        1000000/20000000, #7
                                                        100000/20000000), #8
                                       target_upper = c(5000/20000000, #1
                                                        90000/20000000, #2
                                                        1500/20000000, #3
                                                        10000/20000000, #4
                                                        90000/20000000, #5
                                                        90000000/20000000, #6
                                                        90000000/20000000, #7
                                                        900000/20000000)) #8
  },
  n_cores = 3)

 source("add_species_names.R")

parameters_GLV_nosihek_adapteddiscfunc_200 <-
  add_species_names(parameters_GLV_nosihek_adapteddiscfunc_200,
                    c("seabirds",
                      "terrestrial crabs",
                      "carnivorous crabs",
                      "cane spiders",
                      "geckos",
                      "cockroaches",
                      "terrestrial arthropods",
                      "native trees"))

save(parameters_GLV_nosihek_adapteddiscfunc_200,
     file = "parameters_GLV_nosihek_adapteddiscfunc_200.RData")

##### check if equilibrium abundances match

arg <- args_function(interaction_matrix = interaction_matrix_nosihek,
                     upper_bounds_growth_rate =
                       upper_bounds_growth_rate_nosihek,
                     model = "GLV")

# check the equilibrium for one parameter set
reconstructed <-
  reconstruct_matrix_growthrates(parameters_GLV_nosihek_adapteddiscfunc_200[[1]],
                                 arg)

r <- reconstructed$growthrates$growthrates
A <- reconstructed$growthrates$interaction_matrix
equilibrium <- solve(A, -r)

source("select_EEM_outputs.R")
#implement function
select_EEM_outputs(ensemble = parameters_GLV_nosihek_adapteddiscfunc_200,
                   target_lower = c(500/20000000, #1
                                    10000/20000000, #2
                                    500/20000000, #3
                                    1000/20000000, #4
                                    10000/20000000, #5
                                    1000000/20000000, #6
                                    1000000/20000000, #7
                                    100000/20000000), #8
                   target_upper = c(5000/20000000, #1
                                    90000/20000000, #2
                                    1500/20000000, #3
                                    10000/20000000, #4
                                    90000/20000000, #5
                                    90000000/20000000, #6
                                    90000000/20000000, #7
                                    900000/20000000),
                   sim_args = arg)
# it works¨
