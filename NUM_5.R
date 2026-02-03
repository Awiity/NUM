Chebyshev <- function(N){
  a0 <- rep(0, N)
  if(N==1) return(c(1))
  a1 <- a0
  if(N==2) return(c(1,0))
  a <- a0
  a0[1] <- 1
  a1[2] <- 1
  for(i in 1:(N-2)){
    a <- 2*c(0,a1[1:(N-1)])-a0
    a0 <- a1
    a1 <- a
  }
  return(rev(a))
}
Horner <- function(a, x){
  N <- length(a)
  y <- a[1]
  if(N > 1){
    for(i in 2:N){
      y <- y*x+a[i]
    }
  }
  return(y)
}
HornerD <- function(a, x){
  N <- length(a)
  y <- a[1]
  yd <- a[1]
  if(N > 2){
    for(i in 2:(N-1)){
      y <- y*x+a[i]
      yd <- yd*x+y
    }
    y <- y*x+a[N]
  }else if(N == 2){
    return(c(y*x+a[2], yd))
  }else if(N == 1){
    return(c(y, 0))
  }
  return(c(y, yd))
}
NewtonHorner <- function(a, x, eps = 0.0000001){
  N <- length(a)
  repeat{
    y <- a[1]
    yd <- a[1]
    if(N > 2){
      for(i in 2:(N-1)){
        y <- y*x+a[i]
        yd <- yd*x+y
      }
      y <- y*x+a[N]
    }
    xnew <- x-y/yd
    if(abs(xnew-x) < eps) return(xnew)
    x <- xnew
  }
}
Lagrange <- function(x, y, alpha){
  N <- length(x)
  suma <- 0
  for(i in 1:N){
    pitko <- 1
    for(j in 1:N){
      if(j != i) pitko <- pitko * (alpha - x[j]) /(x[i]-x[j])
    }
    suma <- suma + y[i]*pitko
  }
  return(suma)
}

coef <- Chebyshev(11)
x <- seq(-1, 1, 0.01)
y <- Horner(coef, x)
plot(x, y, type = 'l', col='violet', lwd=4)
abline(h = 0)
x1 <- NewtonHorner(coef, 0.5)
points(x1,0)
x2 <- NewtonHorner(coef, -0.5)
points(x2,0)
x3 <- NewtonHorner(coef, -0.75)
points(x3,0)
x4 <- NewtonHorner(coef, -0.99)
points(x4,0)
x5 <- NewtonHorner(coef, 0.99)
points(x5,0)
x6 <- NewtonHorner(coef, 0.2)
points(x6,0)

z <- c(x1, x2, x3, x4, x5, x6)+runif(6, -0.1, 0.1)
y <- Horner(coef, z)
points(z, y, col='blue')
lines(x, Lagrange(z, y, x), col='blue', lwd=2)

m <- length(x)
z <- x[sample(1:m, 11)]
y <- Horner(coef, z)
points(z, y, col='green')
lines(x, Lagrange(z, y, x), col='green', lwd=2)


print(HornerD(c(11, -4, 3, 2), 3))
print(HornerD(c(11, -4), 3))
print(HornerD(c(11), 3))