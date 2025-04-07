add_introduced_species <- function(native_parameters,
                                   introduced_lower_bound_growth_rate = 1,
                                   introduced_upper_bound_growth_rate = 5,
                                   introduced_self_sign,
                                   introduced_row_signs,
                                   introduced_col_signs,
                                   introduced_k) {
  #get the number of native species
  n_native <- length(native_parameters[[1]]$growthrates)
  extended_parameters <-
    lapply(seq_along(native_parameters),
           function(i) {
             x <- native_parameters[[i]]
             #sample the introduced species growth rates
             introduced_growth_rate <-
               runif(1,
                     min = introduced_lower_bound_growth_rate,
                     max = introduced_upper_bound_growth_rate)
             #sample interactions between natives and introduced species
             introduced_col <-
               sapply(seq_len(n_native),
                      function(j) {
                        if (introduced_col_signs[j] != 0) {
                          introduced_col_signs[j]*runif(1)
                        } else {
                          0
                        }
                      })
             introduced_row <-
               sapply(seq_len(n_native),
                      function(j) {
                        if (introduced_row_signs[j] != 0) {
                          introduced_row_signs[j]*runif(1)
                        } else{
                          0
                        }
                      })
             #---> until here it's the same

             # Let r_native and A be the native growth rates and interaction matrix
             r_native <- x$growthrates #new
             A <- x$interaction_matrix #new
             # Compute the full native equilibrium when the introduced species is at its target:
             # Solve: A * N_native + introduced_col * introduced_k + r_native = 0 -> net per capita growth rate of each species is 0 if abundance of introduced species is k
             eq_with_intro <- tryCatch(-solve(A, introduced_col * introduced_k + r_native),
                                       error = function(e) rep(NA, n_native))
             if (any(is.na(eq_with_intro))) {
               stop("Parameter set", i, "Failed to compute full native equilibrium. No unique solution to Ax = b.")
             }

             # Now compute the self-interaction for the introduced species.
             # The equilibrium condition for the introduced species (full system) is:
             #   r_int + sum(introduced_row * eq_with_intro) + a_int,int * introduced_k = 0.
             # Solve for a_int,int:
             introduced_self <- (-introduced_growth_rate -
                                   sum(introduced_row *
                                         eq_with_intro)) / introduced_k

             # Optionally check that the sign is as desired:
             if (sign(introduced_self_sign) != sign(introduced_self)) {
               print(paste("Parameter set", i, "has inconsistent introduced self interaction sign. Removing this parameter set."))
               return(NULL)
             }
             extended_growthrates <- c(x$growthrates, introduced_growth_rate)
             extended_interaction_matrix <-
               matrix(0, nrow = n_native + 1, ncol = n_native + 1)
             extended_interaction_matrix[2:(n_native + 1), 2:(n_native + 1)] <-
               A
             extended_interaction_matrix[2:(n_native + 1), 1] <- introduced_col
             extended_interaction_matrix[1, 2:(n_native + 1)] <- introduced_row
             extended_interaction_matrix[1, 1] <- introduced_self
             if (!is.null(colnames(A))) {
               colnames(extended_interaction_matrix) <-
                 c("Introduced species", colnames(A))
             } else {
               colnames(extended_interaction_matrix) <-
                 c("Introduced species", seq_len(ncol(A)))
             }
             if (!is.null(rownames(A))) {
               rownames(extended_interaction_matrix) <-
                 c("Introduced species", rownames(A))
             } else {
               rownames(extended_interaction_matrix) <-
                 c("Introduced species", seq_len(nrow(A)))
             }
             return(list(growthrates = extended_growthrates,
                         interaction_matrix = extended_interaction_matrix))
           })
  # Filter out parameter sets that returned NULL
  extended_parameters <- extended_parameters[!sapply(extended_parameters, is.null)]
  return(extended_parameters)
}
