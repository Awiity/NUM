rm(list=ls())

f <- \(t, P){return(r*P*(1-P/K))}
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

r <- 0.025
K <- r/0.002
P0 <- 12
A <- (K-P0)/P0
tmin <- 0
tmax <- 2000
dt <- 0.05
t <- seq(tmin, tmax, dt)
steps <- length(t)
P <- numeric(steps)
P[1] <- P0

for(step in 2:steps){
  P[step] <- RK4(f, t[step-1], P[step-1], dt)
}

plot(t, P, type='l', col='red')
lines(t, K/(1+A*exp(-r*t)), col='green')