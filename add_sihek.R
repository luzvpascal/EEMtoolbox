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

save(parameters_GLV_plussishek_adapteddiscfunc_200,
     file = "parameters_GLV_plussishek_adapteddiscfunc_200.RData")

##### check the equilibrium abundance of sihek --> not the exact number but the order of size --> eg same number of decimals#####
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
eq
eq[!is.na(eq)]
which(is.na(eq))
non_eq[!is.na(non_eq)]
parameters_GLV_plussishek_adapteddiscfunc_200[which(is.na(eq))]
length(eq[!is.na(eq)])*100/5000

##### ------> find a way to: select the values where all abundances are at least 0 ? set something where if abundance is negative, set to 0?? is it needed? can we adapt that in the plot_projection function by adapting the derivative_func?
##### ------> find a way to set abundance of sihek to 9 when all other abundances are set to 1 from the nosihek case --> this is done --> set eq abundance to K and then percentage of this
##### ------> run without sihek first, save abundance of all species, take them as initial conditions for the rest? --> not needed since taking the eq abundance without sihek? noo! the eq abundances are not the same since sihek is in there. the one should be the 1 of the first EEM(). can we extract that? Yes! because we set them --> just set the lower or upper pool of the ones we set!


##### plot projections with sihek in the model #####

species_names <- c("sihek", colnames(parameters_GLV_plussishek_adapteddiscfunc_200[[1]]$interaction_matrix)[2:9])

plot_projections(parameters = parameters_GLV_plussishek_adapteddiscfunc_200,
                 initial_condition = c(6/40, 1, 1, 1, 1, 1, 1, 1, 1), #--> doesn't make sense! the 1 is the resulting eq!! we need the 1 from the previous state without sihek (EEM param set)
                 t_window = c(0,10),
                 scaled = TRUE,
                 species_names = species_names)

plot_projections(parameters = parameters_GLV_plussishek_adapteddiscfunc_200,
                 initial_condition = c((6/40)/20000000,
                                       5000/20000000, #1
                                       90000/20000000, #2
                                       1500/20000000, #3
                                       10000/20000000, #4
                                       90000/20000000, #5
                                       90000000/20000000, #6
                                       90000000/20000000, #7
                                       900000/20000000), #8
                 t_window = c(0,10),
                 scaled = FALSE,
                 species_names = species_names)

projections_200_nonscaled_upperbounds <- calculate_projections(
  parameters = parameters_GLV_plussishek_adapteddiscfunc_200,
  initial_condition = c(6/20000000,
                        5000/20000000, #1
                        90000/20000000, #2
                        1500/20000000, #3
                        10000/20000000, #4
                        90000/20000000, #5
                        90000000/20000000, #6
                        90000000/20000000, #7
                        900000/20000000), #8
  t_window = c(0,10),
  scaled = FALSE,
  species_names = species_names)

save(projections_200_nonscaled_upperbounds,
     file = "projections_200_nonscaled_upperbounds.RData")

##### plot with own way #####

library(ggplot2)
library(dplyr)
abundance <- group_by(projections_200_nonscaled_upperbounds, time, species)
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

#new plot
ggplot(abundance) +
  guides(fill = guide_legend(title = "Species"),
         color = guide_legend(title = "Species")) +
  xlab("Time") +
  ylab("Abundance") +
  geom_ribbon(aes(x = time,
                  ymin = lower,
                  ymax = upper,
                  color = species,
                  fill = species),
              alpha = 0.2) +
  geom_line(aes(x = time, y = median_pop, color = species),
            linewidth = 1.5) +
  theme_bw()

