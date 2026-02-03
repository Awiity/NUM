midPointRule <- function(f, a, b, n = 1){
  h <- (b-a)/n
  return(h*sum(f(h*(1:n)+a-0.5*h)))
}
puleniIntervalu <- function(f, a, b){
  fa <- f(a)
  fb <- f(b)
  if(fa*fb < 0){
    repeat{
      c <- (a+b)/2
      if(c == a || c == b) return(c)
      fc <- f(c)
      if(fa*fc < 0){
        b <- c
        fb <- fc
      }else{
        a <- c
        fa <- fc
      }
    }
  }else{
    stop("Chybny pocatecni interval.")
  }
}
res <- puleniIntervalu(\(x) x*x-2, 1, 2)
print(res)
print(sqrt(2))

res <- midPointRule(\(x) exp(-x*x), 0, 1.175, 100000)
print(res)

res <- puleniIntervalu(\(p) midPointRule(\(x) exp(-x*x), 0, p, 100000) - 0.8, 1, 2)
print(res)

# ------------------------------------
Horner <- function(a, x){
  n <- length(a)
  res <- a[n]
  if(n > 1) for(i in (n-1):1) res <- res*x+a[i]
  return(res)
}
Vandermonde <- function(x, y, n){
  A <- matrix(1, n, n)
  for(i in 2:n) A[, i] <- A[, i-1] * x
  return(solve(A, y))
}

n <- 3
x <- runif(n, 0, 2*pi)
y <- sin(x)
plot(x, y, col = 'red')
coef <- Vandermonde(x, y, n)
t <- seq(0, 2*pi, length.out = 100)
lines(t, Horner(coef, t), col='green')

a <- numeric(3)
a[1] <- y[1]
a[2] <- (y[2]-y[1])/(x[2]-x[1])
a[3] <- (y[3]-a[1]-a[2]*(x[3]-x[1]))/((x[3]-x[1])*(x[3]-x[2]))
lines(t, (a[3]*(t-x[2])+a[2])*(t-x[1])+a[1], col='blue')


z <- y[1]*(t-x[2])*(t-x[3])/((x[1]-x[2])*(x[1]-x[3]))+
     y[2]*(t-x[1])*(t-x[3])/((x[2]-x[1])*(x[2]-x[3]))+
     y[3]*(t-x[1])*(t-x[2])/((x[3]-x[1])*(x[3]-x[2]))
lines(t, z, col='cyan')