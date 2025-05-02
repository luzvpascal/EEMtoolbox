#' @title Merge multiple single-species introductions into one system
#' @param EEM_intros a pairlist of objects obtained with add_introduced_species()
#' @param sign_interaction_intros one number or a matrix of signs for the interactions between introduced species and native species
#' @param mode "updated" or "recycled". If "updated", the interaction signs are updated for each projection. If "recycled", the interaction signs are recycled in all projections
#' @return a list of EEM objects




merge_introductions <- function(EEM_intros,
                                sign_interaction_intros = 0,
                                mode = "updated") {

  if (is.vector(EEM_intros)) {
    stop("parameter EEM_intros should be a pairlist of EEM objects obtained with pairlist(EEM_object1, EEM_object2, ...)")
  }

  #number of native species
  n_native <- ncol(EEM_intros[[1]][[1]]$interaction_matrix) - 1

  #number of introduced species
  n_intro <- length(EEM_intros)

  tot_species <- n_native + n_intro

  n_proj <- length(EEM_intros[[1]])
  for (i in seq_len(n_intro)) {
    if (length(EEM_intros[[i]]) != n_proj) {
      stop("All objects in EEM_intros should have the same number of projections")
    }
  }

  A_full <- vector(mode = "pairlist", length = n_proj)

  species_names_intro <- c()
  species_names_native <- colnames(
    EEM_intros[[1]][[1]]$interaction_matrix)[2:(n_native + 1)]
  for (i in seq_len(n_intro)) {
    if (!is.null(colnames(EEM_intros[[i]][[1]]$interaction_matrix))) {
      species_names_intro <- c(
        species_names_intro,
        colnames(EEM_intros[[i]][[1]]$interaction_matrix)[1])
    }
  }

  # Create a new interaction matrix for the merged system
  A_new <- matrix(NA,
                  nrow = tot_species,
                  ncol = tot_species,
                  dimnames = list(c(species_names_intro, species_names_native),
                                  c(species_names_intro, species_names_native)))

  r_full <- rep(NA, tot_species)

  for (i in seq_len(n_proj)) {
    A_full[[i]] <- list(growthrates = r_full,
                        interaction_matrix = A_new)
  }
    for (i in seq_len(n_proj)) {
      A_full[[i]]$interaction_matrix[(n_intro + 1):tot_species,
                                     (n_intro + 1):tot_species] <-
      EEM_intros[[1]][[i]]$interaction_matrix[2:(n_native + 1),
                                              2:(n_native + 1)]
      A_full[[i]]$growthrates[(n_intro + 1):tot_species] <-
  EEM_intros[[1]][[i]]$growthrates[2:(n_native + 1)]

    for (j in seq_len(n_intro)) {
      A_full[[i]]$interaction_matrix[j, (n_intro + 1):tot_species] <-
        EEM_intros[[j]][[i]]$interaction_matrix[1, 2:(n_native + 1)]
      A_full[[i]]$interaction_matrix[(n_intro + 1):tot_species, j] <-
        EEM_intros[[j]][[i]]$interaction_matrix[2:(n_native + 1), 1]
      A_full[[i]]$interaction_matrix[j, j] <- EEM_intros[[j]][[i]]$interaction_matrix[1, 1]

      A_full[[i]]$growthrates[j] <- EEM_intros[[j]][[i]]$growthrates[1]
    }
    }
  if (mode == "recycled") {
    interaction_term <-
      runif(sign_interaction_intros)*sign_interaction_intros
  }

  if (length(sign_interaction_intros) != 1) {
    for (i in seq_len(n_proj)) {
      for (j in seq_len(n_intro)) {
        if (mode == "updated") {
          interaction_term <-
            runif(sign_interaction_intros)*sign_interaction_intros
        }
        A_full[[i]]$interaction_matrix[
          j, which(is.na(A_full[[i]]$interaction_matrix[j,]))] <-
          interaction_term[j,which(!is.na(interaction_term[j,]))]
      }
    }
  } else if (length(sign_interaction_intros) == 1 &&
             sign_interaction_intros != 0) {
    for (i in seq_len(n_proj)) {
      A_full[[i]]$interaction_matrix[is.na(A_full[[i]]$interaction_matrix)] <-
        sign_interaction_intros
    }
  }

  return(A_full)
         }
