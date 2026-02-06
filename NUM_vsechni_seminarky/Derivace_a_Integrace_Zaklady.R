derivativeI <- \(f, x, h)(f(x+h)-f(x))/h
derivativeII <- \(f, x, h) 0.5*(f(x+h)-f(x-h))/h
secondDerivative <- \(f, x, h) (f(x+h)-2*f(x)+f(x-h))/(h*h)
midPointRule <- \(f, a, b, n=1){
  h <- (b-a)/n
  return(h*sum(f(a+h*(1:n)-0.5*h)))
}
trapezoidalRule <- \(f, a, b, n=1){
  h <- (b-a)/n
  suma <- 0.5*(f(a)+f(b))
  if(n > 1) suma <- suma+sum(f(a+h*(1:(n-1))))
  return(h*suma)
}
SimpsonRule <- \(f, a, b, n=1){
  h <- (b-a)/n
  suma <- f(a)+f(b)+4*(sum(f(a+h*(1:n)-h/2)))
  if(n > 1) suma <- suma+2*sum(f(a+h*(1:(n-1))))
  return(h*suma/6)
}

x <- seq(0, 2*pi, length.out = 200)

# prvni testovani
plot(cos, xlim=c(0, 2*pi))
lines(x, derivativeI(sin, x, h=0.00000001), col='red')
h <- 0.5
lines(x, derivativeI(sin, x, h), col='green')
lines(x, derivativeII(sin, x, h), col='blue')
lines(x, -secondDerivative(cos, x, h), col='cyan')

# druhe testovani
x <- 0
h <- (1:100)*0.001
plot(h, derivativeI(sin, x, h))
abline(h=1)

x <- 0
h <- 2*2^{-(1:20)}
plot(h, derivativeI(sin, x, h), log='x')
abline(h=1)

x <- 0
h <- 0.0001*2^{-(1:20)}
plot(h, secondDerivative(cos, x, h), log='x')
abline(h=-1)

# vypocet integralu
a <- 0
b <- pi/2
plot(cos, xlim=c(a, b))
print(sin(b)-sin(a))
n <- 100
print(midPointRule(cos, a, b, n))
print(trapezoidalRule(cos, a, b, n))
print(SimpsonRule(cos, a, b, n))