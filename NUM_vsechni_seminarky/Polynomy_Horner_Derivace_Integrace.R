MidPointRule <- function(f, a, b, n){
  h <- (b-a)/n
  return(h*sum(f(a+h*(1:n)-0.5*h)))
}
TrapezoidalRule <- function(f, a, b, n){
  h <- (b-a)/n
  suma <- 0.5*(f(a)+f(b))
  if(N > 1) suma <- suma+sum(f(a+h*(1:(n-1))))
  return(h*suma)
}
SimpsonsRule <- function(f, a, b, n){
  h <- (b-a)/n
  suma <- f(a)+f(b)+4*sum(f(a+h*(1:n)-0.5*h))
  if(N > 1) suma <- suma+2*sum(f(a+h*(1:(n-1))))
  return(h*suma/6)
}
Horner <- function(x, a, N){
  y <- 0
  yd <- 0
  ydd <- 0
  for(i in N:3){
    y <- y*x + a[i]
    yd <- yd*x + y
    ydd <- ydd*x + yd
  }
  y <- y*x + a[2]
  yd <- yd*x + y
  y <- y*x + a[1]
  return(list(y=y, yd=yd, ydd=ydd*2))
}
DerivativeOneSide <- \(f, x, h) (f(x+h)-f(x))/h
DerivativeBothSide <- \(f, x, h) 0.5*(f(x+h)-f(x-h))/h
SecondDerivative <- \(f, x, h) (f(x+h)-2*f(x)+f(x-h))/(h*h)

a <- c(2,-1,3,-4,5,-6)
# a <- c(2,-1,3)
N <- length(a)
x <- 1
h <- 100.001
print("------------------------------")
res <- DerivativeOneSide(\(x) Horner(x, a, N)$y, x, h = h)
print(res)
res <- DerivativeBothSide(\(x) Horner(x, a, N)$y, x, h = h)
print(res)
print(Horner(x, a, N)$yd)
print("------------------------------")
res <- SecondDerivative(\(x) Horner(x, a, N)$y, x, h = h)
print(res)
print(Horner(x, a, N)$ydd)
print("------------------------------")

# print(sin(1))
# print(cos(1))
# print(DerivativeOneSide(sin, 1, h=0.0000001))
# print(-sin(1))
# print(SecondDerivative(sin, 1, h=0.0001))
# print("------------------------------")

od <- 1
do <- 2
res <- Horner(do, a, N)$y - Horner(od, a, N)$y
print(res)
n <- 10
res <- MidPointRule(\(x) Horner(x, a, N)$yd, od, do, n=n)
print(res)
res <- TrapezoidalRule(\(x) Horner(x, a, N)$yd, od, do, n=n)
print(res)
res <- SimpsonsRule(\(x) Horner(x, a, N)$yd, od, do, n=n)
print(res)