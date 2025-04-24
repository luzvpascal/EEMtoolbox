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
             r_native <- x$growthrates #new
             A <- x$interaction_matrix #new
             #sample the introduced species growth rates
             introduced_growth_rate <-
               runif(1,
                     min = introduced_lower_bound_growth_rate,
                     max = introduced_upper_bound_growth_rate)
             #sample interactions between natives and introduced species
             introduced_col <-
               sapply(seq_len(n_native),
                      function(j) {
                        ifelse(introduced_col_signs[j] != 0,
                               introduced_col_signs[j] * runif(1),
                               0)
                      })
             introduced_row <-
               sapply(seq_len(n_native),
                      function(j) {
                        ifelse(introduced_row_signs[j] != 0,
                               introduced_row_signs[j] * runif(1),
                               0)
                      })

             # Compute self‐interaction for full system equilibrium = k
             invA <- solve(A)   # A^{-1}
             vec <- (r_native + introduced_col * introduced_k)
             numerator <- as.numeric(introduced_row %*% (invA %*% vec)) - introduced_growth_rate
             introduced_self <- numerator / introduced_k

             # Check that the sign is as desired:
             if (sign(introduced_self) != sign(introduced_self_sign)) {
               cat("Introduced self-interaction sign of parameter set", i, "does not match the desired sign.\n")
             }
             extended_growthrates <- c(introduced_growth_rate, r_native)
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
