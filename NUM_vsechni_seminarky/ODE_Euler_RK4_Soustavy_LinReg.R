rm(list=ls())
EulerI <- function(x, y, h, f){return(y+h*f(x,y))}
EulerIIuprostred <- function(x, y, h, f){
  hpul <- 0.5*h
  return(y+h*f(x+hpul,y+hpul*f(x,y)))
}
EulerIIprumerna<- function(x, y, h, f){
  return(y+0.5*h*(f(x,y)+f(x+h, y+h*f(x,y))))
}
EulerIIprumernaOpakovane<- function(x, y, h, f){
  hpul <- 0.5*h
  k1 <- f(x,y)
  xNew <- x+h
  yNew <- y+h*k1
  yNew <- y+hpul*(k1+f(xNew, yNew))
  yNew <- y+hpul*(k1+f(xNew, yNew))
  # yNew <- y+hpul*(k1+f(xNew, yNew))
  # yNew <- y+hpul*(k1+f(xNew, yNew))
  # yNew <- y+hpul*(k1+f(xNew, yNew))
  # yNew <- y+hpul*(k1+f(xNew, yNew))
  return(yNew)
}
RK4 <- function(x, y, h, f){
  k1 <- f(x, y)
  hpul <- 0.5*h
  xnapul <- x + hpul
  k2 <- f(xnapul, y+hpul*k1)
  k3 <- f(xnapul, y+hpul*k2)
  k4 <- f(x+h, y+h*k3)
  return(y+h*(k1+2*(k2+k3)+k4)/6)
}
f <- function(t, y){return(0.3*y)}
F <- function(t, Y){return(c(Y[2], -0.3*Y[1]))}

tLast <- 7
dt <- 1
t <- seq(0, tLast, dt)
Nsteps <- length(t)
yI <- numeric(Nsteps)
yIIu <- numeric(Nsteps)
yIIp <- numeric(Nsteps)
yIV <- numeric(Nsteps)

y0 = 20
yI[1] <- y0
yIIu[1] <- y0
yIIp[1] <- y0
yIV[1] <- y0
for(i in 2:Nsteps){
  yI[i] <- EulerI(t[i-1], yI[i-1], dt, f)
  yIIu[i] <- EulerIIuprostred(t[i-1], yIIu[i-1], dt, f)
  yIIp[i] <- EulerIIprumernaOpakovane(t[i-1], yIIp[i-1], dt, f)
  yIV[i] <- RK4(t[i-1], yIV[i-1], dt, f)
}
plot(t, 20*exp(0.3*t), type = 'l', lwd=2)
lines(t, yI, col='blue')
lines(t, yIIu, col='green')
lines(t, yIIp, col='cyan')
lines(t, yIV, col='violet')


tLast <- 7
dt <- 0.5
t <- seq(0, tLast, dt)
Nsteps <- length(t)
Y <- matrix(0, nrow=Nsteps, ncol=2)
Y0 = c(-10, 0)
Y[1, ] <- Y0
for(i in 2:Nsteps){
  Y[i, ] <- EulerIIprumernaOpakovane(t[i-1], Y[i-1, ], dt, F)
}
w <- sqrt(0.3)
plot(t, -10*cos(t*w), type = 'l')
lines(t, Y[, 1], col="blue")
lines(t, Y[, 2], col="red")
# lines(t, 5*sin(t*w), col='violet')


GaussElimination <- function(A, b){
  n <- length(b)
  Ab <- cbind(A, b)
  # primy chod
  for(k in 1:(n-1)){
    for(i in (k+1):n){
      j <- (k+1):(n+1)
      kc <- Ab[i,k]/Ab[k,k]
      Ab[i, j] <- Ab[i, j] -kc*Ab[k, j]
    }
  }
  # zpetny chod
  x <- Ab[, n+1]
  x[n] <- x[n]/Ab[n,n]
  for(i in (n-1):1){
    j <- (i+1):n
    x[i] <- (x[i] - sum(Ab[i,j]*x[j]))/Ab[i,i]
  }
  return(x)
}
GaussEliminationPivoting <- function(A, b){
  n <- length(b)
  Ab <- cbind(A, b)
  # primy chod
  for(k in 1:(n-1)){
    kmax <- which.max(abs(Ab[k:n, k])) + k - 1
    if(kmax != k){
      j <- k:(n+1)
      pamatuj <- Ab[k, j]
      Ab[k, j] <- Ab[kmax, j]
      Ab[kmax, j] <- pamatuj
    }
    for(i in (k+1):n){
      j <- (k+1):(n+1)
      kc <- Ab[i,k]/Ab[k,k]
      Ab[i, j] <- Ab[i, j] -kc*Ab[k, j]
    }
  }
  # zpetny chod
  x <- Ab[, n+1]
  x[n] <- x[n]/Ab[n,n]
  for(i in (n-1):1){
    j <- (i+1):n
    x[i] <- (x[i] - sum(Ab[i,j]*x[j]))/Ab[i,i]
  }
  return(x)
}

Jacobi <- function(A, b, epsMin = 0.000000000000000001, epsMax = 10000000){
  n <- length(b)
  x <- numeric(n)
  LpU <- A
  diag(LpU) <- 0
  recD <- 1/diag(A)
  repeat{
    xNew <- (b-LpU%*%x)*recD
    eps <- sum((xNew-x)^2)
    if(eps < epsMin){return(c(xNew))}
    else if(eps > epsMax){stop("Nekonverguje.")}
    x <- xNew
  }
}


n <- 8
A <- matrix(runif(n*n), n, n)
# A[2,1:2] <- A[1,1:2]
diag(A) <- diag(A)*1000
b <- runif(n)
x <- numeric(n)
GaussEliminationPivoting(A, b)
solve(A, b)
Jacobi(A, b)

n <- 10
x <- 1:n
y <- x*x
plot(x, y)
sumx <- sum(x)
sumy <- sum(y)
sumxx <- sum(x*x)
sumxy <- sum(x*y)

detA <- n*sumxx-sumx*sumx
detq <- sumy*sumxx-sumx*sumxy
detk <- n*sumxy-sumx*sumy

q <- detq/detA
k <- detk/detA
abline(a=q, b=k)