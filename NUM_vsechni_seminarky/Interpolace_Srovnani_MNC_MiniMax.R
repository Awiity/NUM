Horner <- function(x, a, N){
  y <- a[N]
  for(i in (N-1):1) y <- y*x+a[i]
  return(y)
}

# interpolace polynomomem \sum_{i} a_i x^i
# matice soustavy: Vandermondova

N <- 14
x <- runif(N, -2,2)
y <- sin(x)
t <- seq(min(x), max(x), length.out=100)

plot(x, y, col='red')

A <- matrix(1, N, N)
for(i in 2:N) A[,i] <- A[,i-1]*x
a <- solve(A, y)
lines(t, Horner(t, a, N), col='blue')

## pro N=2
# A <- matrix(c(1, 1, x), 2, 2)
# a <- solve(A, y)
# lines(t, a[1] + a[2] * t, col='blue')

## pro N=3
# A <- matrix(c(1, 1, 1, x, x*x), 3, 3)
# a <- solve(A, y)
# lines(t, a[1] + a[2] * t + a[3] * t*t, col='blue')

# Newton (viz minule)
NewtonPolynomialCoefs <- function(x, y){
  N <- length(x)
  a <- y
  if(N == 1) return(a)
  a[2] <- (y[2]-y[1])/(x[2]-x[1])
  if(N == 2) return(a)
  for(i in 3:N){
    suma <- y[i]
    soucin <- 1
    for(j in 1:(i-1)){
      suma <- suma - a[j]*soucin
      soucin <- soucin*(x[i]-x[j])
    }
    a[i] <- suma/soucin
  }
  return(a)
}
NewtonPolynomialValues <- function(a, x, t){
  N <- length(a)
  y <- rep(a[N], length(t))
  if(N > 1){
    for(i in (N-1):1) y <- y*(t-x[i])+a[i]
  }
  return(y)
}
lines(t, NewtonPolynomialValues(
  NewtonPolynomialCoefs(x, y), x, t), col='green')
## pro N=2
# a <- c(y[1], (y[2]-y[1])/(x[2]-x[1]))
# lines(t, a[1] + a[2]*(t-x[1]), col='green')

## pro N=3
# a <- c(a, (y[3]-(a[1]+a[2]*(x[3]-x[1])))/
#          ((x[3]-x[1])*(x[3]-x[2])))
# lines(t, a[1] +
#         a[2]*(t-x[1]) +
#         a[3]*(t-x[1])*(t-x[2]), col='green')

# Langrange (viz 13. Ĺ™Ă­jna)
Lagrange <- function(xk, yk, x){
  # https://en.wikipedia.org/wiki/Lagrange_polynomial
  k <- length(xk)
  L <- 0
  for(j in 1:k){
    lj <- 1
    for(m in 1:k){
      if(m != j) lj <- lj*(x-xk[m])/(xk[j]-xk[m])
    }
    L <- L + yk[j]*lj
  }
  return(L)
}

lines(t, Lagrange(x, y, t), col='cyan')

## pro N=2
# w1 <- (x[2]-t)/(x[2]-x[1])
# w2 <- (t-x[1])/(x[2]-x[1])
# lines(t, w1*y[1] + w2*y[2], col='violet')

## pro N=3
# w1 <- (t-x[2])*(t-x[3])/((x[1]-x[2])*(x[1]-x[3]))
# w2 <- (t-x[1])*(t-x[3])/((x[2]-x[1])*(x[2]-x[3]))
# w3 <- (t-x[1])*(t-x[2])/((x[3]-x[1])*(x[3]-x[2]))
# lines(t, w1*y[1] + w2*y[2] + w3*y[3], col='violet')

# lineĂˇrnĂ­ aproximace
## LSA
A <- matrix(c(N, sum(x), sum(x), sum(x*x)), 2)
a <- solve(A, c(sum(y), sum(x*y)))
lines(t, a[1]+a[2]*t, col='darkgreen', lwd = 2)
abline(a, col='yellow')

## Monte Carlo MiniMax
aerror <- abs((a[1]+a[2]*x)-y)
maxaerror <- max(aerror)

maxerror <- maxaerror
bestb <- a
for(i in 1:100000){
  b <- runif(2, -0.5, 1.5)*a
  berror <- abs((b[1]+b[2]*x)-y)
  maxberror <- max(berror)
  if(maxberror < maxerror){
    print(i)
    maxerror <- maxberror
    bestb <- b
  }
}
abline(bestb, col='brown')
print(maxaerror)
print(maxerror)