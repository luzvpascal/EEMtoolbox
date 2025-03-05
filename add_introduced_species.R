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
    lapply(native_parameters,
           function(x) {
             #sample the introduced species growth rates
             introduced_growth_rate <-
               runif(1,
                     min = introduced_lower_bound_growth_rate,
                     max = introduced_upper_bound_growth_rate)
             #sample interactions between natives and introduced species
             introduced_col <-
               sapply(seq_len(n_native),
                      function(i) {
                        if (introduced_col_signs[i] != 0) {
                          introduced_col_signs[i]*runif(1)
                        } else {
                          0
                        }
                      })
             introduced_row <-
               sapply(seq_len(n_native),
                      function(i) {
                        if (introduced_row_signs[i] != 0) {
                          introduced_row_signs[i]*runif(1)
                        } else{
                          0
                        }
                      })
             #compute the self interaction of the introduced species
             ## -> we need to compute the required self-interaction to ensure
             ##that the introduced species's equilibrium value is K
             ## --> we need introduced_growth_rate + native_effect + introduced_self * introduced_k = 0
             ## --> introduced_self = -(introduced _growth_rate + native_effect)/introduced_k
             ## the native_effect is the sum of the effects from natives on the
             ## introduced species, calculated as
             ## sum_{j=1}^{n_native} (introduced_row[j]*native_eq[j]))
             ## native_eq is the native equilibrium values of the system
             native_eq <- solve(x$interaction_matrix, -x$growthrates)
             native_effect <- sum(introduced_row*native_eq)
             introduced_self <- -(introduced_growth_rate + native_effect)/introduced_k
             if (sign(introduced_self_sign) != sign(introduced_self)){
               print("The introduced self interaction sign is not consistent with
                    the required equilibrium value")
             }
             extended_growthrates <- c(x$growthrates, introduced_growth_rate)
             extended_interaction_matrix <-
               matrix(0, nrow = n_native + 1, ncol = n_native + 1)
             extended_interaction_matrix[2:(n_native+1), 2:(n_native+1)] <-
               x$interaction_matrix
             extended_interaction_matrix[2:(n_native+1), 1] <-
               introduced_col
             extended_interaction_matrix[1, 2:(n_native+1)] <-
               introduced_row
             extended_interaction_matrix[1, 1] <-
               introduced_self
             if (is.null(colnames(x$interaction_matrix)) == FALSE) {
               colnames(extended_interaction_matrix) <-
                 c("Introduced species", colnames(x$interaction_matrix))
             } else {
               colnames(extended_interaction_matrix) <-
                 c("Introduced species", c(1:ncol(x$interaction_matrix)))
             }
             if (is.null(rownames(x$interaction_matrix)) == FALSE) {
               rownames(extended_interaction_matrix) <-
                 c("Introduced species", rownames(x$interaction_matrix))
             } else {
               rownames(extended_interaction_matrix) <-
                 c("Introduced species", c(1:ncol(x$interaction_matrix)))
             }
             return(list(growthrates = extended_growthrates,
                         interaction_matrix = extended_interaction_matrix))
           })
  return(extended_parameters)
}
