add_introduced_species <- function(native_parameters,
                                   introduced_lower_bound_growth_rate = 1,
                                   introduced_upper_bound_growth_rate = 5,
                                   introduced_self_sign,
                                   introduced_row_signs,
                                   introduced_col_signs,
                                   # interaction_strength_bounds is a matrix with the lower and upper bounds for the interaction strengths of the introduced species with the native species. The columns are the species. The first row is the lower bound, and the second row is the upper bound. The minimum interaction strength is 0 and the maximum is 1.
                                   interaction_strength_lower_bounds,
                                   interaction_strength_upper_bounds) {
  n_native <- length(native_parameters[[1]]$growthrates)
  extended_parameters <-
    lapply(native_parameters,
           function(x) {
             introduced_growth_rate <-
               runif(1,
                     min = introduced_lower_bound_growth_rate,
                     max = introduced_upper_bound_growth_rate)
             introduced_col <-
               sapply(seq_len(n_native),
                      function(i) {
                        if (introduced_col_signs[i] != 0) {
                          introduced_col_signs[i]*runif(1,
                                                        interaction_strength_lower_bounds[i],
                                                        interaction_strength_upper_bounds[i])
                        } else {
                          0
                        }
                      })
             introduced_row <-
               sapply(seq_len(n_native),
                      function(i) {
                        if (introduced_row_signs[i] != 0) {
                          introduced_row_signs[i]*runif(1,
                                                        interaction_strength_lower_bounds[i],
                                                        interaction_strength_upper_bounds[i])
                        } else{
                          0
                        }
                      })
             introduced_self <-
               if (introduced_self_sign != 0) {
                 introduced_self_sign*runif(1,
                                            interaction_strength_lower_bounds[1],
                                            interaction_strength_upper_bounds[1])
               } else {
                 0
               }
             extended_growthrates <- c(x$growthrates, introduced_growth_rate)
             extended_interaction_matrix <-
               matrix(0, nrow = n_native + 1, ncol = n_native + 1)
             extended_interaction_matrix[1:n_native, 1:n_native] <-
               x$interaction_matrix
             extended_interaction_matrix[1:n_native, n_native + 1] <-
               introduced_col
             extended_interaction_matrix[n_native + 1, 1:n_native] <-
               introduced_row
             extended_interaction_matrix[n_native + 1, n_native + 1] <-
               introduced_self
             if (is.null(colnames(x$interaction_matrix)) == FALSE) {
               colnames(extended_interaction_matrix) <-
                 c(colnames(x$interaction_matrix),
                   "Introduced species")
             } else {
               colnames(extended_interaction_matrix) <-
                 c(c(1:ncol(x$interaction_matrix)),
                   "Introduced species")
             }
             if (is.null(rownames(x$interaction_matrix)) == FALSE) {
               rownames(extended_interaction_matrix) <-
                 c(rownames(x$interaction_matrix),
                   "Introduced species")
             } else {
               rownames(extended_interaction_matrix) <-
                 c(c(1:ncol(x$interaction_matrix)),
                   "Introduced species")
             }
             return(list(growthrates = extended_growthrates,
                         interaction_matrix = extended_interaction_matrix))
           })
  return(extended_parameters)
}
