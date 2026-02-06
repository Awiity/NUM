simpson <- function(x_base, coefs, a, b, m = 100) {
  h <- (b - a) / m
  x_vals <- seq(a, b, by = h)
  y_vals <- sapply(x_vals, function(v) poly_val(x_base, coefs, v))
  
  s <- y_vals[1] + y_vals[m + 1] + 
    4 * sum(y_vals[seq(2, m, by = 2)]) + 
    2 * sum(y_vals[seq(3, m - 1, by = 2)])
  return(s * h / 3)
}
simpson_hilbert <- function(k, m = 100) {
  a <- 0; b <- 1
  h <- (b - a) / m
  t <- seq(a, b, by = h)
  y <- t^k
  s <- y[1] + y[m + 1] + 
    4 * sum(y[seq(2, m, by = 2)]) + 
    2 * sum(y[seq(3, m - 1, by = 2)])
  return(s * h / 3)
}