NewtonPolAdd <- function(x, y, a, i){
  if(i == 1) return(y[i])
  if(i == 2) return((y[2]-y[1])/(x[2]-x[1]))
  nasobitko <- 1
  suma <- 0
  for(j in 1:(i-1)){
    suma <- suma + a[j]*nasobitko
    nasobitko <- nasobitko * (x[i]-x[j])
  }
  return((y[i] - suma)/nasobitko)
}
NewtonPolynomial <- function(x, y){
  N <- length(x)
  a <- y
  for(i in 2:N) a[i] <- NewtonPolAdd(x, y, a, i)
  return(a)
}

N <- 3
x <- 1:N
y <- sin(x)
plot(x, y)
a <- NewtonPolynomial(x, y)
print(a)
t <- seq(1, N, 0.01)
z <- a[1] + a[2]*(t - x[1]) + a[3]*(t - x[1])*(t - x[2])
lines(t, z)


Jacobi <- function(A, b, maxiter = 1000, eps = 0.0000000001){
  N <- length(b)
  xold <- xnew <- rep(0, N)
  for(iter in 1:maxiter){
    for(i in 1:N){
      xnew[i] <- b[i]
      for(j in 1:N){
        if(j != i) xnew[i] <- xnew[i] - A[i,j]*xold[j]
      }
      xnew[i] <- xnew[i]/A[i,i]
    }
    error <- sum((xnew - xold)^2)
    if(is.infinite(error)) return(NA)
    if(error < eps){
      cat("Number of iterations (Jacobi): ", iter, "\n")
      return(xnew)
    }
    xold <- xnew
  }
  return(NA)
}
JacobiVer2 <- function(A, b, maxiter = 1000, eps = 0.0000000001){
  N <- length(b)
  D <- diag(A)
  diag(A) <- 0
  xold <- xnew <- rep(0, N)
  for(iter in 1:maxiter){
    xnew <- (b - A%*%xold)/D
    error <- sum((xnew - xold)^2)
    if(is.infinite(error)) return(NA)
    if(error < eps){
      cat("Number of iterations (Jacobi): ", iter, "\n")
      return(xnew)
    }
    xold <- xnew
  }
  return(NA)
}
GS <- function(A, b, maxiter = 1000, eps = 0.0000000001){
  N <- length(b)
  xold <- xact <- rep(0, N)
  for(iter in 1:maxiter){
    for(i in 1:N){
      xact[i] <- b[i]
      for(j in 1:N){
        if(j != i) xact[i] <- xact[i] - A[i,j]*xact[j]
      }
      xact[i] <- xact[i]/A[i,i]
    }
    error <- sum((xact - xold)^2)
    if(is.infinite(error)) return(NA)
    if(error < eps){
      cat("Number of iterations (Gauss-Seidel): ", iter, "\n")
      return(xact)
    }
    xold <- xact
  }
  return(NA)
}
GaussSeidel <- function(A, b, maxiter = 1000, eps = 0.0000000001){
  N <- length(b)
  xold <- xact <- rep(0, N)
  for(iter in 1:maxiter){
    for(i in 1:N) xact[i] <- (b[i] - sum(A[i,-i]*xact[-i]))/A[i,i]
    error <- sum((xact - xold)^2)
    if(is.infinite(error)) return(NA)
    if(error < eps){
      cat("Number of iterations (Gauss-Seidel): ", iter, "\n")
      return(xact)
    }
    xold <- xact
  }
  return(NA)
}

N <- 1000
A <- matrix(runif(N*N, -1, 1), N)
diag(A) <- diag(A)*1.01^N
b <- A %*% 1:N

# print(solve(A, b))

tb <- Sys.time()
x <- Jacobi(A, b, eps = 0.000001)
te <- Sys.time()
print(te-tb)
print(x)

tb <- Sys.time()
x <- JacobiVer2(A, b, eps = 0.000001)
te <- Sys.time()
print(te-tb)
print(x)

tb <- Sys.time()
x <- GS(A, b, eps = 0.000001)
te <- Sys.time()
print(te-tb)
# print(x)

tb <- Sys.time()
x <- GaussSeidel(A, b, eps = 0.000001)
te <- Sys.time()
print(te-tb)
# print(x)