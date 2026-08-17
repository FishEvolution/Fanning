

import.will <- function(data) {
  #rearrange data so each row is a taxa, each column is one of (4 * n_harmonics) coefficients
  rows <- seq_len(nrow(data)) #extract string with rows indices
  coefmat <- matrix(NA, ncol = nrow(data) * 2, nrow = ncol(data)/2)
  for(i in 1:nrow(coefmat)) {

    A <- data[rows[rows%%2==1],(i*2)-1]
    B <- data[rows[rows%%2==1],i*2]
    C <- data[1 + rows[rows%%2==1],(i*2)-1]
    D <- data[1 + rows[rows%%2==1],i*2]

    coefmat[i,] <- c(A, B, C, D)
  }

  #name columns (coefficient:n_harmonic)
  colnames(coefmat) <- c(paste0("A", 1:(ncol(coefmat)/4)),
                         paste0("B", 1:(ncol(coefmat)/4)),
                         paste0("C", 1:(ncol(coefmat)/4)),
                         paste0("D", 1:(ncol(coefmat)/4)))

  return(coefmat)
}








