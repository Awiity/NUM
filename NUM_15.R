Horner <- function(x, a, N){
  y <- a[N]
  for(i in (N-1):1) y <- y*x+a[i]
  return(y)
}
MNCforPolynomials <- function(x, y, n){
  A <- matrix(0, n, n)
  b <- numeric(n)
  
  for(r in 1:n){
    b[r] <- sum(y*x^{r-1})
    for(s in 1:n){
      A[r, s] <- sum(x^{r+s-2})
    }
  }
  return(solve(A, b))
}
MNCforPolynomialsBetter <- function(x, y, n){
  m <- length(x)
  X <- matrix(1, n, m)
  for(r in 2:n) X[r, ] <- X[r-1, ]*x
  A <- X%*%t(X)
  b <- X%*%y
  return(solve(A, b))
}

n <- 6
a <- c(0.3,6,-0.2,-32,0.1,32)
x <- seq(-1, 1, 0.000001)
m <- length(x)
y <- Horner(x, a, n) #*runif(m, 0.99, 1.01)
# plot(x, y)

n <- 5
tb <- Sys.time()
coef <- MNCforPolynomials(x, y, n)
te <- Sys.time()
print(coef)
print(te-tb)
# lines(x, Horner(x, coef, n), col='blue')

tb <- Sys.time()
coef <- MNCforPolynomialsBetter(x, y, n)
te <- Sys.time()
print(coef)
print(te-tb)
# lines(x, Horner(x, coef, n), col='green')