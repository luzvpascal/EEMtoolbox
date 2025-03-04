source("add_introduced_species.R")

###### add the introduced species to parameter sets produced by EEM #####
parameters_GLV_plussishek_adapteddiscfunc_200 <-
  add_introduced_species(native_parameters =
                           parameters_GLV_nosihek_adapteddiscfunc_200,
                         introduced_lower_bound_growth_rate = 0,
                         introduced_upper_bound_growth_rate = 1.1,
                         introduced_self_sign = -1,
                         introduced_row_signs = c(1, 1, 1, 1, 1, 1, 1, 1),
                         introduced_col_signs = c(-1, -1, -1, -1, -1, -1, -1, 0),
                         interaction_strength_upper_bounds = c(0.5, 0.5, 1, 1,
                                                               0.5, 0.1, 0.1,
                                                               1, 1),
                         interaction_strength_lower_bounds = c(0, 0, 0, 0, 0,
                                                               0, 0, 0, 0))

#extract species names
species_names <- c(colnames(parameters_GLV_plussishek_adapteddiscfunc[[1]]$interaction_matrix)[1:8], "sihek")

##### check the abundance of sihek #####
interaction_matrix_with_sihek <- matrix(c(-1, 0, -1, 0, 0, 0, 0, 1, -1,
                                          0, -1, -1, 0, 0, 0, 0, 1, -1,
                                          1, 1, -1, 0, 0, 1, 1, 0, -1,
                                          0, 0, 0, -1, 1, 1, 1, 0, -1,
                                          0, 0, 0, -1, -1, 0, 1, 0, -1,
                                          0, 0, -1, -1, 0, -1, 0, 1, -1,
                                          0, 0, -1, -1, -1, 0, -1, 1, -1,
                                          1, -1, -1, 0, 0, 0, 1, -1, 0,
                                          1, 1, 1, 1, 1, 1, 1, 1, -1),
                                        ncol = 9, nrow = 9, byrow = TRUE)

arg_withsihek <- args_function(
  interaction_matrix = interaction_matrix_with_sihek,
  upper_bounds_growth_rate = c(upper_bounds_growth_rate_nosihek, 1.1),
  model = "GLV")

# check the equilibrium for one parameter set
reconstructed_withsihek <-
  reconstruct_matrix_growthrates(parameters_GLV_plussishek_adapteddiscfunc_200[[1]],
                                 arg_withsihek)

r <- reconstructed_withsihek$growthrates$growthrates
A <- reconstructed_withsihek$growthrates$interaction_matrix
equilibrium <- solve(A, -r)

##### ------> find a way to: select the values where all abundances are at least 0 ? set something where if abundance is negative, set to 0?? is it needed? can we adapt that in the plot_projection function by adapting the derivative_func?
##### ------> find a way to set abundance of sihek to 9 when all other abundances are set to 1 from the nosihek case
##### ------> run without sihek first, save abundance of all species, take them as initial conditions for the rest?


##### plot projections with sihek in the model #####
plot_projections(parameters = parameters_GLV_plussishek_adapteddiscfunc_200,
                 initial_condition = c(1, 1, 1, 1, 1, 1, 1, 1,
                                       6/40),
                 t_window = c(0,20),
                 scaled = TRUE,
                 species_names = species_names)

##### plot projections first without sihek and then sihek added -> plot_projections with initial conditions = last number form previous plot, scaled = FALSE
