select_EEM_outputs <- function(ensemble,
                               target_lower,
                               target_upper,
                               mode = "native",
                               n_intro = 1) {
  # ensemble: a list of parameter sets from EEM.
  # target_equilibrium: the vector of desired species abundances (e.g., c(10000, 100, ...)).
  # tolerance: a threshold on the difference (relative error) between the normalized computed equilibrium and target.
  # mode: the mode of the simulation, either "native" or "disturbed".

  # Prepare an empty list to collect selected parameter sets
  selected <- list()
  unselected <- list()
  outputs <- list(outputs_selected = c(),
                  outputs_unselected = c())

  for (i in seq_along(ensemble)) {
    # Reconstruct the parameters: growth rates and interaction matrix.
    r <- ensemble[[i]]$growthrates
    A <- ensemble[[i]]$interaction_matrix

    # Compute the equilibrium: A %*% N + r = 0  =>  N = solve(A, -r)
    current_eq <- tryCatch(solve(A, -r),
                           error = function(e) rep(NA, length(r)))

    # Check for NA values (if the system is unsolvable, skip this set)
    if (any(is.na(current_eq))) {
      next
    }
    if (mode == "native") {
    # If the equilibrium is within the bounds, select this parameter set.
    if (length(which(current_eq > target_lower)) == length(r) &&
        length(which(current_eq < target_upper)) == length(r)) {
      selected[[length(selected) + 1]] <-
        list(ensemble[[i]], data.frame("equilibrium" = as.vector(current_eq, mode = "numeric"),
                                   "species" = 1:length(r)))
      outputs$outputs_selected <-
        c(outputs$outputs_selected, i)
    } else {
      unselected[[length(unselected) + 1]] <-
        list(ensemble[[i]], data.frame("equilibrium" = as.vector(current_eq, mode = "numeric"),
                                   "species" = 1:length(r)))
      outputs$outputs_unselected <-
        c(outputs$outputs_unselected, i)
    }
    }
  if (mode == "disturbed") {
    # If the equilibrium is within the bounds, select this parameter set.
    if (current_eq[1] > target_lower[1] &&
        current_eq[1] < target_upper[1]) {
      selected[[length(selected) + 1]] <-
        list(ensemble[[i]], data.frame("equilibrium" = as.vector(current_eq, mode = "numeric"),
                                       "species" = 1:length(r)))
      outputs$outputs_selected <-
        c(outputs$outputs_selected, i)
    } else {
      unselected[[length(unselected) + 1]] <-
        list(ensemble[[i]], data.frame("equilibrium" = as.vector(current_eq, mode = "numeric"),
                                       "species" = 1:length(r)))
      outputs$outputs_unselected <-
        c(outputs$outputs_unselected, i)
    }
  }
}
  if (length(selected) == 0) {
    cat("No parameter sets found within the bounds.\n")
    df <- data.frame(equilibrium = numeric(), species = numeric())
    for (i in c(1:length(unselected))) {
      df <- data.frame(equilibrium = c(df$equilibrium,
                                       unselected[[i]][[2]]$equilibrium),
                       species = c(df$species,
                                   unselected[[i]][[2]]$species))
    }
  } else {
    cat("Selected", length(selected), "parameter sets.\n")
  df <- data.frame(equilibrium = numeric(), species = numeric())
  for (i in c(1:length(selected))) {
    df <- data.frame(equilibrium = c(df$equilibrium,
                                     selected[[i]][[2]]$equilibrium),
                     species = c(df$species,
                                 selected[[i]][[2]]$species))
  }
  }

  if (mode == "disturbed") {
    for (i in 1:length(df$species)) {
      if (df$species[i] %in% seq_len(n_intro)) {
        df$species[i] <- paste("intro",
                               df$species[i], sep = "_")
      }
    }
  }

  uf <- data.frame(upper = target_upper, species = unique(df$species))
  lf <- data.frame(lower = target_lower, species = unique(df$species))

  a <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df, ggplot2::aes(x = species, y = equilibrium,
                                                color = "equilibrium")) +
    ggplot2::geom_point(data = uf, ggplot2::aes(x = species, y = upper,
                                                color = "upper")) +
    ggplot2::geom_point(data = lf, ggplot2::aes(x = species, y = lower,
                                                color = "lower")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45, vjust = 0.7)) +
    ggplot2::labs(x = "Species", y = "Computed equilibrium",
                  colors = "Legends")

  print(a)
  return(outputs)
}
