midPointRule <- \(f, a, b, n=1){
  h <- (b-a)/n
  return(h*sum(f(a+h*(1:n)-0.5*h)))
}
Romberg <- function(f, a, b, n = 2){
  m <- 2^n
  h <- (b-a)/m
  x <- a + h*(1:(m-1))
  y <- f(x)
  res <- numeric(n)
  
  odd <- seq(1, m, 2)
  pocetHodnot <- m
  for(i in n:1){
    h <- 2*h
    pocetHodnot <- pocetHodnot/2
    sekvence <- odd[1:pocetHodnot]
    res[i] <- h * sum(y[sekvence])
    y <- y[-sekvence]
  }
  A <- matrix(0, nrow=n, ncol=n)
  A[, 1] <- res
  nasobic <- 1
  for(i in 2:n){
    nasobic <- nasobic*4
    A[i:n, i] <- (nasobic*A[i:n, i-1]-A[(i-1):(n-1), i-1])/(nasobic-1)
  }
  return(A)
  # 1, 3, 5, 7, 9, 11, 13, 15
  # 2, 6, 10, 14
  # 4, 12
  # 8
}
RombergJ <- function(f, a, b, n = 2){
  m <- 2^n
  h <- (b-a)/m
  x <- a + h*(1:(m-1))
  y <- f(x)
  res <- numeric(n)
  
  sekvence <- seq(1, m, 2)
  pocetHodnot <- m/2
  for(i in n:1){
    h <- 2*h
    # print(sekvence)
    res[i] <- h * sum(y[sekvence])
    pocetHodnot <- pocetHodnot/2
    sekvence <- 2*sekvence[1:pocetHodnot]
  }
  A <- matrix(0, nrow=n, ncol=n)
  A[, 1] <- res
  nasobic <- 1
  for(i in 2:n){
    nasobic <- nasobic*4
    A[i:n, i] <- (nasobic*A[i:n, i-1]-A[(i-1):(n-1), i-1])/(nasobic-1)
  }
  return(A)
}
print(Romberg(sin, 0, pi, n = 5))
print(RombergJ(sin, 0, pi, n = 5))