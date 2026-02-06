GausovaEliminaceV1 <- function(A, b, N){
  Ab <- cbind(A, b)
  # primy chod
  for(p in 1:(N-1)){
    for(i in (p+1):N){
      for(j in (p+1):(N+1)){
        nasobek <- Ab[i, p]/Ab[p, p]
        Ab[i, j] <- Ab[i, j] - nasobek * Ab[p, j]
      }
    }
  }
  # zpetny chod
  x <- numeric(N)
  A <- Ab[, -(N+1)]
  b <- Ab[, N+1]
  x[N] <- b[N]/A[N, N]
  for(i in (N-1):1) x[i] <- (b[i] - sum(A[i, (i+1):N]*x[(i+1):N]))/A[i, i]
  return(x)
}
GausovaEliminaceV2 <- function(A, b, N){
  Ab <- cbind(A, b)
  # primy chod
  for(p in 1:(N-1)){
    prec <- 1/Ab[p, p]
    for(i in (p+1):N){
      j <- (p+1):(N+1)
      Ab[i, j] <- Ab[i, j] - Ab[i, p] * prec * Ab[p, j]
    }
  }
  # zpetny chod
  x <- numeric(N)
  A <- Ab[, -(N+1)]
  b <- Ab[, N+1]
  x[N] <- b[N]/A[N, N]
  for(i in (N-1):1){
    j <- (i+1):N
    x[i] <- (b[i] - sum(A[i, j]*x[j]))/A[i, i]
  }
  return(x)
}
GausovaEliminaceSPivotaci <- function(A, b, N){
  Ab <- cbind(A, b)
  # primy chod
  for(p in 1:(N-1)){
    pmax <- which.max(abs(Ab[p:N, p])) + p - 1
    if(pmax != p){
      j <- p:(N+1)
      aux <- Ab[p, j]
      Ab[p, j] <- Ab[pmax, j]
      Ab[pmax, j] <- aux
    }
    prec <- 1/Ab[p, p]
    for(i in (p+1):N){
      j <- (p+1):(N+1)
      Ab[i, j] <- Ab[i, j] - Ab[i, p] * prec * Ab[p, j]
    }
  }
  # zpetny chod
  x <- numeric(N)
  A <- Ab[, -(N+1)]
  b <- Ab[, N+1]
  x[N] <- b[N]/A[N, N]
  for(i in (N-1):1){
    j <- (i+1):N
    x[i] <- (b[i] - sum(A[i, j]*x[j]))/A[i, i]
  }
  return(x)
}
LUrozklad <- function(A, N){
  for(p in 1:(N-1)){
    prec <- 1/A[p, p]
    for(i in (p+1):N){
      j <- (p+1):N
      nasobek <- A[i, p] * prec
      A[i, j] <- A[i, j] - nasobek * A[p, j]
      A[i, p] <- nasobek
    }
  }
  return(A)
}
LUzpetnyChod <- function(LU, b, N){
  # L y = b
  y <- numeric(N)
  y[1] <- b[1]
  for(i in 2:N){
    j <- 1:(i-1)
    y[i] <- b[i] - sum(LU[i, j]*y[j])
  }
  # U x = y
  x <- numeric(N)
  x[N] <- y[N]/LU[N, N]
  for(i in (N-1):1){
    j <- (i+1):N
    x[i] <- (y[i] - sum(LU[i, j]*x[j]))/LU[i, i]
  }
  return(x)
}

N <- 4
A <- matrix(runif(N*N, -1, 1), N)
# A[2, 1:2] <- A[1, 1:2]
b <- runif(N, -1, 1)

print(GausovaEliminaceSPivotaci(A, b, N))
print(solve(A, b))

LU <- LUrozklad(A, N)
print(LUzpetnyChod(LU, b, N))
# print(LUrozklad(A, N))
# U <- LU
# U[lower.tri(U)] <- 0
# L <- LU
# L[upper.tri(L)] <- 0
# diag(L) <- 1
# print("----------------")
# L %*% U
# A
# print("----------------")
# print(LUzpetnyChod(LU, b, N))