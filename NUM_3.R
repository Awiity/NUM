Newton <- function(f, fd, x0, iter=1000){
  for(i in 1:iter){
    x <- x0-f(x0)/fd(x0)
    x0 <- x
  }
  return(x)
}
NewtonII <- function(f, fd, x0, eps=0.00000001){
  repeat{
    x <- x0-f(x0)/fd(x0)
    if(abs(x-x0) < eps) break
    x0 <- x
  }
  return(x)
}
NewtonIII <- function(f, fd, x0, maxiter=1000, eps=0.00000001){
  for(i in 1:maxiter){
    x <- x0-f(x0)/fd(x0)
    if(abs(x-x0) < eps) break
    x0 <- x
  }
  return(list(x=x, iter=i))
}

Secny <- function(f, x0, h=0.0001, maxiter=1000, eps=0.00000001){
  for(i in 1:maxiter){
    x <- x0-2*h*f(x0)/(f(x0+h)-f(x0-h))
    if(abs(x-x0) < eps) break
    x0 <- x
  }
  return(list(x=x, iter=i))
}

Steffensen <- function(f, x0, maxiter=1000, eps=0.00000001){
  for(i in 1:maxiter){
    h <- f(x0)
    x <- x0-h*f(x0)/(f(x0+h)-f(x0))
    if(abs(x-x0) < eps) break
    x0 <- x
  }
  return(list(x=x, iter=i))
}

Halley <- function(f, fd, fdd, x0, maxiter=1000, eps=0.00000001){
  for(i in 1:maxiter){
    y <- f(x0)
    yd <- fd(x0)
    x <- x0-y*yd/(yd*yd-0.5*y*fdd(x0))
    if(abs(x-x0) < eps) break
    x0 <- x
  }
  return(list(x=x, iter=i))
}
Bisection <- function(f, a, b){
  fa <- f(a)
  fb <- f(b)
  if(fa*fb < 0){
    repeat{
      c <- (a+b)/2
      if(c==a || c==b) break
      fc <- f(c)
      if(fa*fc < 0){
        b <- c
        fb <- fc
      }else{
        a <- c
        fa <- fc
      }
    }
    return(c)
  }else{
    stop("Neplati: fa*fb < 0!")
  }
}

RegulaFalsi <- function(f, a, b){
  fa <- f(a)
  fb <- f(b)
  if(fa*fb < 0){
    repeat{
      c <- (fa*b-fb*a)/(fa-fb)
      if(c==a || c==b) break
      fc <- f(c)
      if(fa*fc < 0){
        b <- c
        fb <- fc
      }else{
        a <- c
        fa <- fc
      }
    }
    return(c)
  }else{
    stop("Neplati: fa*fb < 0!")
  }
}

FixedPointIteration <- function(f, x0, maxiter=1000, eps=0.00000001){
  for(i in 1:maxiter){
    x <- f(x0)
    if(abs(x-x0) < eps) break
    x0 <- x
  }
  return(list(x=x, iter=i))
}

print(Newton(\(x) x*x-2, \(x) 2*x, 2))
print(NewtonII(\(x) x*x-2, \(x) 2*x, 2))
print(NewtonIII(\(x) x*x-2, \(x) 2*x, 2))
print(Newton(function(x) x*x-3, function(x) 2*x, 2))

res <- NewtonIII(\(x) cos(x)-x, \(x) -sin(x)-1, 2)
print(res$iter)
print(res$x)

res <- Secny(\(x) cos(x)-x, 2)
print(res$iter)
print(res$x)

res <- Steffensen(\(x) cos(x)-x, 2)
print(res$iter)
print(res$x)

res <- Halley(\(x) cos(x)-x, \(x) -sin(x)-1, \(x) -cos(x), 2)
print("Halley")
print(res$iter)
print(res$x)

res <- FixedPointIteration(\(x) cos(x), 2)
print("FixedPointIteration")
print(res$iter)
print(res$x)

res <- Bisection(\(x) cos(x)-x, 0.5, 1)
print(res)

res <- RegulaFalsi(\(x) cos(x)-x, 0.5, 1)
print(res)

res <- RegulaFalsi(\(x) x*x-3, 1, 2)
print(res)