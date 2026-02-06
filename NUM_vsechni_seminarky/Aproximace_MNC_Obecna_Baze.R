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
Phi <- function(x){
  n <- 10
  res <- numeric(n)
  res[1] <- 1
  for(i in 2:n) res[i] <- res[i-1]*x
  return(res)
}
Theta <- function(x){
  return(sapply(1:100, function(k) sin(k * x)))
}
MNC <- function(x, y, Phi){
  A <- sapply(x, Phi)
  return(GaussEliminationPivoting(A%*%t(A), A%*%y))
}
n <- 1000
x <- sort(runif(n, 0, 2*pi))
y <- sin(x) + 2*sin(5*x)
plot(x, y, col='red')
a <- MNC(x, y , Phi)
lines(x, sapply(x, function(x) sum(a * Phi(x))), col='blue')

a <- MNC(x, y , Theta)
# plot(a, type='l')
lines(x, sapply(x, function(x) sum(a * Theta(x))), col='cyan')