Bisection <- function(f, xleft, xright){
  yleft <- f(xleft)
  yright <- f(xright)
  if(yleft * yright < 0){
    # counter <- 0
    repeat{
      xcenter <- (xleft + xright)/2
      # counter <- counter + 1
      if(xcenter == xleft || xcenter == xright) break
      ycenter <- f(xcenter)
      if(yleft * ycenter < 0){
        xright <- xcenter
        yright <- ycenter
      }else{
        xleft <- xcenter
        yleft <- ycenter
      }
    }
  }else{
    stop("Chyba na vstupu.")
  }
  # print(counter)
  return(xcenter)
}

RegulaFalsi <- function(f, xleft, xright){
  yleft <- f(xleft)
  yright <- f(xright)
  if(yleft * yright < 0){
    # counter <- 0
    repeat{
      xnew <- xleft - yleft * (xright - xleft)/(yright - yleft)
      # counter <- counter + 1
      if(xnew == xleft || xnew == xright) break
      ynew <- f(xnew)
      if(yleft * ynew < 0){
        xright <- xnew
        yright <- ynew
      }else{
        xleft <- xnew
        yleft <- ynew
      }
    }
  }else{
    stop("Chyba na vstupu.")
  }
  # print(counter)
  return(xnew)
}

FixedPointIteration <- function(f, xstart, maxiter = 1000, eps = 0.000000001){
  plot(f, xlim=c(-pi/2, pi/2), lwd=2, col="blue")
  abline(0, 1, lwd=2, col="black")
  segments(xstart, 0, xstart, f(xstart), col = "red")
  xold <- xstart
  for(i in 1:maxiter){
    xnew <- f(xold)
    segments(xold, xnew, xnew, xnew, col = "red")
    segments(xnew, xnew, xnew, f(xnew), col = "red")
    if(abs(xnew - xold) < eps) break
    xold <- xnew
  }
  return(list(x = xnew, iter = i, conv = i < maxiter))
}

NewtonMethod <- function(f, fd, xstart, maxiter = 1000, eps = 0.000000001){
  xold <- xstart
  for(i in 1:maxiter){
    xnew <- xold - f(xold)/fd(xold)
    if(abs(xnew - xold) < eps) break
    xold <- xnew
  }
  return(list(x = xnew, iter = i, conv = i < maxiter))
}

NewtonMethodB <- function(f, xstart, h = 0.001, maxiter = 1000, eps = 0.000000001){
  xold <- xstart
  for(i in 1:maxiter){
    xnew <- xold - 2*h*f(xold)/(f(xold+h)-f(xold-h))
    if(abs(xnew - xold) < eps) break
    xold <- xnew
  }
  return(list(x = xnew, iter = i, conv = i < maxiter))
}

print(sqrt(2))

print("Bisekce")
res <- Bisection(\(x) x*x-2, 1, 2)
print(res)

print("RegulaFalsi")
res <- RegulaFalsi(\(x) x*x-2, 1, 2)
print(res)

print("Bisekce")
res <- Bisection(\(x) cos(x) - x, 0, 1)
print(res)

print("FixedPointIteration")
res <- FixedPointIteration(\(x) cos(x), -1)
print(res)

print("NewtonMethod")
res <- NewtonMethod(\(x) cos(x) - x, \(x) -sin(x) - 1, -1)
print(res)

print("NewtonMethodB")
res <- NewtonMethodB(\(x) cos(x) - x, -1)
print(res)