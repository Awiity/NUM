rm(list=ls())

f <- \(t, n){return(-lambda*n)}
F <- \(t, Y){
  # return(c(Y[2], -Y[1]))
  return(c(Y[2], -Y[1]-0.01*Y[2]))
}
EulerI <- function(f, x, y, h){
  return(y+h*f(x,y))
}
EulerIIa <- function(f, x, y, h){
  hpul <- 0.5*h
  return(y+h*f(x+hpul,y+hpul*f(x,y)))
}
EulerIIb <- function(f, x, y, h){
  return(y+0.5*h*(f(x,y)+f(x+h,y+h*f(x,y))))
}
RK4 <- function(f, x, y, h){
  hpul <- 0.5*h
  xpul <- x+hpul
  k1 <- f(x,y)
  k2 <- f(xpul,y+hpul*k1)
  k3 <- f(xpul,y+hpul*k2)
  k4 <- f(x+h, y+h*k3)
  return(y+h*(k1+2*(k2+k3)+k4)/6)
}

lambda <- 2
n0 <- 1000
tmin <- 0
tmax <- 2
dt <- 0.05
t <- seq(tmin, tmax, dt)
steps <- length(t)
n <- numeric(steps)
n[1] <- n0
nIIa <- numeric(steps)
nIIa[1] <- n0
nIIb <- numeric(steps)
nIIb[1] <- n0
nRK4 <- numeric(steps)
nRK4[1] <- n0
for(step in 2:steps){
  n[step] <- EulerI(f, t[step-1], n[step-1], dt)
  nIIa[step] <- EulerIIa(f, t[step-1], nIIa[step-1], dt)
  nIIb[step] <- EulerIIb(f, t[step-1], nIIb[step-1], dt)
  nRK4[step] <- RK4(f, t[step-1], nRK4[step-1], dt)
}

plot(t, n, type='l', col='red')
lines(t, nIIa, col='cyan')
lines(t, nIIb, col='blue')
lines(t, nRK4, col='green')
lines(t, n0*exp(-lambda*t), col='black')



Y0 <- c(-1, 0)
tmin <- 0
tmax <- 100
dt <- 0.05
t <- seq(tmin, tmax, dt)
steps <- length(t)
Y <- matrix(0, nrow=steps, ncol=2)
Y[1, ] <- Y0
for(step in 2:steps){
  Y[step, ] <- RK4(F, t[step-1], Y[step-1, ], dt)
}

plot(t, Y[,1], type='l', col='red')
lines(t, Y[,2], type='l', col='green')
A <- Y0[1]
w <- 1
# lines(t, A*cos(w*t))