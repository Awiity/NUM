GaussElimination <- function(A, b){
  N <- length(b)
  Ab <- cbind(A, b)
  # pĹ™Ă­mĂ˝ chod
  for(p in 1:(N-1)){
    u <- - 1/Ab[p, p]
    for(r in (p+1):N){
      s <- (p+1):(N+1)
      Ab[r, s] <- Ab[r, s] + Ab[r, p] * u * Ab[p, s]
    }
  }
  # zpÄ›tnĂ˝ chod
  x <- b
  x[N] <- Ab[N, N+1]/Ab[N, N]
  for(r in (N-1):1){
    s <- (r+1):N
    x[r] <- (Ab[r, N+1] - sum(Ab[r, s] * x[s]))/Ab[r, r]
  }
  return(x)
}
GaussEliminationPivoting <- function(A, b){
  N <- length(b)
  Ab <- cbind(A, b)
  # pĹ™Ă­mĂ˝ chod
  for(p in 1:(N-1)){
    imax <- which.max(abs(Ab[p:N, p])) + p - 1
    if(imax != p){
      s <- p:(N+1)
      aux <- Ab[imax, s]
      Ab[imax, s] <- Ab[p, s]
      Ab[p, s] <- aux
    }
    u <- - 1/Ab[p, p]
    for(r in (p+1):N){
      s <- (p+1):(N+1)
      Ab[r, s] <- Ab[r, s] + Ab[r, p] * u * Ab[p, s]
    }
  }
  # zpÄ›tnĂ˝ chod
  x <- b
  x[N] <- Ab[N, N+1]/Ab[N, N]
  for(r in (N-1):1){
    s <- (r+1):N
    x[r] <- (Ab[r, N+1] - sum(Ab[r, s] * x[s]))/Ab[r, r]
  }
  return(x)
}
LUrozklad <- function(A){
  N <- dim(A)[1]
  LU <- A
  for(p in 1:(N-1)){
    w <- 1/LU[p, p]
    for(r in (p+1):N){
      s <- (p+1):N
      u <- LU[r, p] * w
      LU[r, s] <- LU[r, s] - u * LU[p, s]
      LU[r, p] <- u
    }
  }
  return(LU)
}

N <- 5
A <- matrix(runif(N*N), N, N)
A[1, 1] <- 0
b <- runif(N)

N <- 3
A <- matrix(c(2, -3, 11, 7, 2, -8, -5, 1, 0), N, N, byrow = TRUE)
b <- c(43, -5, -20)


print(GaussEliminationPivoting(A, b))
print(solve(A, b))
LU <- LUrozklad(A)
print(LU)