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

parameters_GLV_nosihek_adapteddiscfunc <- EEM(interaction_matrix_nosihek,
       upper_bounds_growth_rate = upper_bounds_growth_rate_nosihek,
                                              algorithm = "SMC-EEM",
                                              disc_func = function(data) {
       adapted_discrepancy_continuous_sum(data,
                                          target_lower = c(5000/100000,
                                                           100000/1000000,
                                                           1000/1000000,
                                                           5000/10000,
                                                           100000/10000000,
                                                           10000000/1000000000,
                                                           100000000/10000000000,
                                                           100000/10000000),
                                          target_upper = c(90000/100000,
                                                           9000000/1000000,
                                                           90000/1000000,
                                                           90000/10000,
                                                           9000000/10000000,
                                                           900000000/1000000000,
                                                           9000000000/10000000000,
                                                           9000000/10000000))
                                              },
                                              n_cores = 3)

source("add_species_names.R")

add_species_names(parameters_GLV_nosihek_adapteddiscfunc,
                  c("seabirds",
                    "terrestrial crabs",
                    "carnivorous crabs",
                    "cane spiders",
                    "geckos",
                    "cockroaches",
                    "terrestrial arthropods",
                    "native trees"))

save(parameters_GLV_nosihek_adapteddiscfunc,
     file = "parameters_GLV_nosihek_adapteddiscfunc.RData")

##### check if equilibrium abundances match

arg <- args_function(interaction_matrix = interaction_matrix_nosihek,
                     upper_bounds_growth_rate =
                       upper_bounds_growth_rate_nosihek,
                     model = "GLV")

# check the equilibrium for one parameter set
reconstructed <-
  reconstruct_matrix_growthrates(parameters_GLV_nosihek_adapteddiscfunc[[1]],
                                 arg)

r <- reconstructed$growthrates$growthrates
A <- reconstructed$growthrates$interaction_matrix
equilibrium <- solve(A, -r)

source("select_EEM_outputs.R")
#implement function
select_EEM_outputs(ensemble = parameters_GLV_nosihek_adapteddiscfunc,
                   target_lower = c(5000/100000,
                                    100000/1000000,
                                    1000/1000000,
                                    5000/10000,
                                    100000/10000000,
                                    10000000/1000000000,
                                    100000000/10000000000,
                                    100000/10000000),
                   target_upper = c(90000/100000,
                                    9000000/1000000,
                                    90000/1000000,
                                    90000/10000,
                                    9000000/10000000,
                                    900000000/1000000000,
                                    9000000000/10000000000,
                                    9000000/10000000),
                   sim_args = arg)
# function works but gives nothing if tolerance is that low!
