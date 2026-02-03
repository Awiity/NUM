
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

# 3*x^2 - 5*x - 7 --> c(-7, -5, 3)
p <- function(a, x){
  N <- length(a)
  res <- 0
  for(i in 1:N) res <- res + a[i]*x^{i-1}
  return(res)
}
Horner <- function(a, x){
  N <- length(a)
  # if(N == 1) return(rep(a, length(x)))
  if(N == 1) a
  res <- a[N]
  for(i in (N-1):1) res <- res*x + a[i]
  return(res)
}
NewtonHorner <- function(a, x, eps=0.00000001){
  N <- length(a)
  repeat{
    y <- 0
    yd <- 0
    for(i in N:2){
      y <- y*x + a[i]
      yd <- yd*x + y
    }
    y <- y*x + a[1]
    xnew <- x-y/yd
    if(abs(xnew-x) < eps) return(x)
    x <- xnew
  }
}


plot(sin, xlim = c(0, 2*pi), col='red')
xk <- runif(10, 0, 2*pi)
yk <- sin(xk)
points(xk, yk, col='darkgreen')
x <- seq(0, 2*pi, 0.01)
lines(x, Lagrange(xk, yk, x), col='darkgreen', lwd=2)

a <- c(-7, -5, 3)
x <- seq(-3, 3, 0.001)
y <- p(a, x)
plot(x, y, type='l', col='red')
y <- Horner(a, x)
lines(x, y, col='blue')
abline(h=0)
xroot <- NewtonHorner(a, -2)
points(xroot, 0, col='red')
xroot <- NewtonHorner(a, 200)
points(xroot, 0, col='red')

a <- c(1, 2)
y <- p(a, x)
plot(x, y, type='l', col='red')
y <- Horner(a, x)
lines(x, y, col='blue')
abline(h=0)

a <- c(1)
y <- p(a, x)
plot(x, y, type='l', col='red')
y <- Vectorize(Horner, 'x')(a, x)
lines(x, y, col='blue')
abline(h=0)