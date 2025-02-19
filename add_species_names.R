add_species_names <- function(parameter,
                              species_names) {
  lapply(parameter, function(x) {
    colnames(x$interaction_matrix) <- c(species_names)
    rownames(x$interaction_matrix) <- c(species_names)
    return(x)
  })
}
