# log likelihood function
logL <- function(alpha, x, y) {
  # -alpha_0 - alpha_1*x_i is used twice, so compute it separately
  alpha_sum <- -alpha[1] - alpha[2]*x
  
  # the (1-y_i)(-alpha_0 - alpha_1*x_i) term
  term1 <- (1-y)*alpha_sum
  
  # the log(1 + exp(-alpha_0 - alpha_1*x_i)) term
  term2 <- log(1 + exp(alpha_sum))
  
  total <- sum(term1 - term2)
  return (total)
}

# gradient of log likelihood
grad_logL <- function(alpha, x, y) {
  # exp(-alpha_0 - alpha_1*x_1) is used several times
  exp_alpha <- exp(-alpha[1] - alpha[2]*x)
  
  # partial with respect to alpha_0
  # -(1-y_i) + exp(-alpha_0 - alpha_1*x_1) /(1+exp(-alpha_0 - alpha_1*x_1))
  p0 <- -(1-y) + exp_alpha/(1 + exp_alpha)
  
  #partial with respect to alpha_1
  # this works out to be x_i times the alpha_0 partial
  p1 <- x*p0
  
  # sum across i
  p0 <- sum(p0)
  p1 <- sum(p1)
  
  return (c(p0, p1))
}

hessian_logL <- function(alpha, x, y)
{
  # we actually don't need y
  a1 <- alpha[1]; a2 <- alpha[2]

  h <- matrix(0, ncol=2, nrow=2)
  nx <- length(x)
  
  for (i in 1:nx) {
    v <- c(1, x[i])
    numerator <- exp(-a1-a2*x[i])
    denom <- (1 + exp(-a1-a2*x[i]))^2
    h <- h - numerator/denom * (v %*% t(v))
  }
  
  return (h)
}


norm <- function(x)
{
  sqrt(sum(x^2))
}


Newton <- function(start, tol=10^-10, max.iter=10^4)
{
  d <- read.table("o_ring_data.txt", header=T)
  x <- d$Temp; y <- d$Failure
  
  alpha <- start
  
  for (i in 1:10) {
    grad <- grad_logL(alpha, x, y) 
    hessian <- hessian_logL(alpha, x, y)
    alpha <- alpha - solve(hessian, grad)
  
    cat("iteration ", i, alpha, "\n")
  }
  
  return (alpha)
}

### code
d <- read.table("o_ring_data.txt", header=T)
x <- d$Temp; y <- d$Failure
















