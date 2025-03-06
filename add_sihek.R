source("add_introduced_species.R")
load("parameters_GLV_nosihek_adapteddiscfunc_200.RData")

###### add the introduced species to parameter sets produced by EEM #####
parameters_GLV_plussishek_adapteddiscfunc_200 <-
  add_introduced_species(native_parameters =
                           parameters_GLV_nosihek_adapteddiscfunc_200,
                         introduced_lower_bound_growth_rate = 0,
                         introduced_upper_bound_growth_rate = 1.1,
                         introduced_self_sign = -1,
                         introduced_row_signs = c(1, 1, 1, 1, 1, 1, 1, 1),
                         introduced_col_signs = c(-1, -1, -1, -1, -1, -1, -1, 0),
                         introduced_k = 40/20000000)

#extract species names
species_names <- c(colnames(parameters_GLV_plussishek_adapteddiscfunc[[1]]$interaction_matrix)[1:8], "sihek")

##### check the equilibrium abundance of sihek --> not the exact number but the order of size --> eg same number of decimals#####
eq_which <- which(sapply(parameters_GLV_plussishek_adapteddiscfunc_200, function(x) x$eq_intro_abundance) > 0)
for (i in 1:5000) {
  r <- parameters_GLV_plussishek_adapteddiscfunc_200[[i]]$growthrates
  A <- parameters_GLV_plussishek_adapteddiscfunc_200[[i]]$interaction_matrix
  equilibrium <- solve(A, -r)
  if ((40/20000000)/equilibrium[1] > 0.1
      && (40/20000000)/equilibrium[1] < 10) {
    eq[[i]]$eq_param_number <- i
    eq[[i]]$eq_intro_abundance <- equilibrium[1]
  }
}

eq <- c()
non_eq <- c()
for (i in 1:5000) {
  r <- parameters_GLV_plussishek_adapteddiscfunc_200[[i]]$growthrates
  A <- parameters_GLV_plussishek_adapteddiscfunc_200[[i]]$interaction_matrix
  equilibrium <- solve(A, -r)
  if ((40/20000000)/equilibrium[1] > 0.1
      && (40/20000000)/equilibrium[1] < 10) {
    eq[i] <- equilibrium[1]
  } else {
    non_eq[i] <- equilibrium[1]
  }
}
eq <- eq[!sapply(eq, is.na)]
non_eq <- non_eq[!sapply(non_eq, is.na)]
eq
non_eq

##### ------> find a way to: select the values where all abundances are at least 0 ? set something where if abundance is negative, set to 0?? is it needed? can we adapt that in the plot_projection function by adapting the derivative_func?
##### ------> find a way to set abundance of sihek to 9 when all other abundances are set to 1 from the nosihek case --> this is done --> set eq abundance to K and then percentage of this
##### ------> run without sihek first, save abundance of all species, take them as initial conditions for the rest? --> not needed since taking the eq abundance without sihek? noo! the eq abundances are not the same since sihek is in there. the one should be the 1 of the first EEM(). can we extract that? Yes! because we set them


##### plot projections with sihek in the model #####
plot_projections(parameters = parameters_GLV_plussishek_adapteddiscfunc_200,
                 initial_condition = c(6/40, 1, 1, 1, 1, 1, 1, 1, 1),
                 t_window = c(0,10),
                 scaled = TRUE,
                 species_names = species_names)

##### plot projections first without sihek and then sihek added -> plot_projections with initial conditions = last number form previous plot, scaled = FALSE
