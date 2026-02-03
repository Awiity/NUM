# Interpolace pĹ™Ă­mkou
x <- c(1, 2)
y <- c(4,11)

t <- seq(x[1], x[2], 0.01)

plot(x, y, col='red')

## Vandermonde
A <- cbind(c(1,1), x)
print(A)
a <- solve(A, y)
print(a)

lines(t, a[1]+a[2]*t, col='blue')

## Newton
a <- c(y[1], (y[2]-y[1])/(x[2]-x[1]))
lines(t, a[1]+a[2]*(t-x[1]), col='green')


## Lagrange
w1 <- (x[2]-t)/(x[2]-x[1])
w2 <- (t-x[1])/(x[2]-x[1])
lines(t, w1*y[1]+w2*y[2], col='violet')

# LSA: linear
x <- 1:5
y <- c(3.2, 1.4, 2.7, -0.8, 11.5)
A <- matrix(c(length(x), sum(x), sum(x), sum(x*x)), 2, 2)
a <- solve(A, c(sum(y), sum(x*y)))
plot(x,y)
t <- seq(min(x), max(x), 0.01)
lines(t, a[1]+a[2]*t, col='blue')

# Interpolace parabolou
N <- 3
x <- 1:N
y <- sample(-N:N, N)

t <- seq(min(x), max(x), 0.01)

plot(x, y, col='red')

## Vandermonde
Horner <- function(a, x, N){
  y <- a[N]
  for(i in (N-1):1) y <- y*x+a[i]
  return(y)
}

N <- length(x)
A <- matrix(1, N, N)
for(i in 2:N) A[, i] <- A[, i-1]*x
a <- solve(A, y)
lines(t, Horner(a, t, N), col='blue')



## Newton
a <- c(y[1], (y[2]-y[1])/(x[2]-x[1]))
a <- c(a, (y[3]-(a[1]+a[2]*(x[3]-x[1])))/((x[3]-x[1])*(x[3]-x[2])))
lines(t, a[1]+
        a[2]*(t-x[1])+
        a[3]*(t-x[1])*(t-x[2]), col='green')


## Lagrange
w1 <- ((t-x[2])*(t-x[3]))/((x[1]-x[2])*(x[1]-x[3]))
w2 <- ((t-x[1])*(t-x[3]))/((x[2]-x[1])*(x[2]-x[3]))
w3 <- ((t-x[1])*(t-x[2]))/((x[3]-x[1])*(x[3]-x[2]))
lines(t, w1*y[1]+w2*y[2]+w3*y[3], col='violet')