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
NewtonPolynomialCoef <- function(x, y, a, i){
  N <- length(x)
  if(i == 1) return(y[1])
  if(i == 2) return((y[2]-y[1])/(x[2]-x[1]))
  suma <- y[i]
  soucin <- 1
  for(j in 1:(i-1)){
    suma <- suma - a[j]*soucin
    soucin <- soucin*(x[i]-x[j])
  }
  return(suma/soucin)
}
NewtonPolynomialValues <- function(a, x, t){
  N <- length(a)
  y <- rep(a[N], length(t))
  if(N > 1){
    for(i in (N-1):1) y <- y*(t-x[i])+a[i]
  }
  return(y)
}

N <- 4
x <- seq(0, 2*pi, length.out = N)
y <- sin(x)
# y <- exp(x)
plot(x, y)
t <- seq(0, 2*pi, 0.01)
a <- NewtonPolynomialCoefs(x, y)
z <- NewtonPolynomialValues(a, x, t)
lines(t, z, col='red')

N <- N + 1
x <- c(x, 3)
y <- sin(x)
y[N] <- y[N]+0.2
points(x[N], y[N], col='blue', pch=17)

a <- c(a, NewtonPolynomialCoef(x, y, a, N))
z <- NewtonPolynomialValues(a, x, t)
lines(t, z, col='violet')

# DalĹˇĂ­ test
xod <- 1
xdo <- 2

N <- 5
x <- sample(seq(xod, xdo, length.out = N))
y <- exp(x)
plot(x, y)

t <- seq(xod, xdo, 0.01)
a <- NULL
for(i in 1:N){
  a <- c(a, NewtonPolynomialCoef(x, y, a, i))
  z <- NewtonPolynomialValues(a, x, t)
  lines(t, z, col=i+1)
}

# Gauss-Seidel
N <- 10
x <- runif(N, -1, 1)
A <- matrix(runif(N*N), N)
diag(A) <- 100*diag(A)
b <- A %*% x

print(x)
print(solve(A, b))

Jacobi <- function(A, b){
  N <- length(b)
  Drec <- 1/diag(A)
  diag(A) <- 0
  x <- rep(0, N)
  for(i in 1:10) x <- (b-A%*%x)*Drec
  return(x)
}
print(Jacobi)

print(Jacobi(A, b) - solve(A, b))