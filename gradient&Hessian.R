U <- function(p)
{
  # get the molecule positions
  m1 <- p[1:2]
  m2 <- p[3:4]
  
  # compute r-squared and then r
  r2 <- sum((m2 - m1)^2)
  r <- sqrt(r2)
  
  Ur <- 1/r^12 - 2/r^6
  
  return (Ur)
}

gradU <- function(p)
{
  # get the molecule positions
  m1 <- p[1:2]
  m2 <- p[3:4]
  
  # get the (x,y) of each molecule
  x1 <- m1[1]
  x2 <- m2[1]
  y1 <- m1[2]
  y2 <- m2[2]
  
  # compute r-squared and then r
  r2 <- sum((m2 - m1)^2)
  r <- sqrt(r2)
  
  # gradient
  factor <- 12*(-1/r^14 + 1/r^8)
  vec <- c(x1 - x2, y1 - y2, x2 - x1, y2 - y1)
  
  g <- factor*vec
  return (g)
}

HU <- function(p)
{
  # get the molecule positions
  m1 <- p[1:2]
  m2 <- p[3:4]
  
  # get the (x,y) of each molecule
  x1 <- m1[1]
  x2 <- m2[1]
  y1 <- m1[2]
  y2 <- m2[2]
  
  # compute r-squared and then r
  r2 <- sum((m2 - m1)^2)
  r <- sqrt(r2)
  
  # summand 1 in equation (5) 
  factor1 <- 12*14/r^16 - 12*8/r^10
  vec <- c(x1 - x2, y1 - y2, x2 - x1, y2 - y1)
  summand1 <- factor1*vec %*% t(vec)
  
  # summand 2
  factor2 <- 12*(-1/r^14 + 1/r^8)
  mat <- matrix(c(1,0,-1,0,
                  0,1,0,-1,
                  -1,0,1,0,
                  0,-1,0,1), nrow=4, ncol=4)
  summand2 <- factor2*mat
  
  
  H <- summand1 + summand2
  return (H)
}

###################################
Energy <- function(p)
{
  d <- length(p)
  num_mol <- d/2
  
  e <- 0
  for (i in 2:num_mol) {
    for (j in 1:(i-1)) {
      coor <- pair_coordinates(i,j)
      p_ij <- p[coor]
      e <- e + U(p_ij)
    }
  }
  
  # add p1 term
  e <- e + p[1]^2 + p[2]^2
  
  # add y term
  e <- e + p[4]^2
  
  return (e)
}

pair_coordinates <- function(i,j)
{
  i_coordinates <- c(2*(i-1)+1, 2*(i-1)+2)
  j_coordinates <- c(2*(j-1)+1, 2*(j-1)+2)
  return (c(i_coordinates, j_coordinates))
}

gradE <- function(p)
{
  # find the dimension
  d <- length(p)
  num_mol <- d/2
  
  g <- rep(0, d)
  for (i in 2:num_mol) {
    for (j in 1:(i-1)) {
      coor <- pair_coordinates(i,j)
      p_ij <- p[coor]
      g[coor] <- g[coor] + gradU(p_ij)
    }
  }
  
  # the ||p_1||^2 term 
  x1 <- p[1]
  y1 <- p[2]
  p1_g <- c(2*x1, 2*y1, rep(0, d-2))
  g <- g + p1_g
  
  # the y_2^2 term
  y2 <- p[4]
  y2_g <- c(0,0,0,2*y2,rep(0, d-4))
  g <- g + y2_g
  
  return (g)
}

HE <- function(p) 
{
  d <- length(p)
  num_mol <- d/2
  
  H <- matrix(0, nrow=d, ncol=d)
  
  for (i in 2:num_mol) {
    for (j in 1:(i-1)) {
      coor <- pair_coordinates(i,j)
      p_ij <- p[coor]
      H[coor,coor] <- H[coor,coor] + HU(p_ij)
    }
  }
  
  # the ||p_1||^2 term 
  p1_m <- matrix(c(0,2,0,2), nrow=2, ncol=2)
  H[1:2,1:2] <- H[1:2,1:2] + p1_m
  
  # the y_2^2 term
  H[4,4] <- H[4,4] + 2
  
  return (H)
}

#######################################################
norm <- function(x)
{
  return (sqrt(sum(x^2)))
}

Newton <- function(p, track_iteration) {
  
  g <- gradE(p)
  i <- 1
  I <- diag(length(g))
  s <- 0
  
  while (norm(g) > 1E-6) {
    h <- HE(p)
    g <- gradE(p)
    
    # start with no Hessian modification
    # In class we discussed more sophisticated methods
    # that remember hm from iteration to iteration, but
    # here I'm taking the simplest approach
    hm <- 0
    
    # make sure we can invert h before we try,
    # modify the Hessian if needed
    while (kappa(h + hm*I) > 1E15)
      hm <- 2*(hm + .1)
    h <- h + hm*I
    
    # ok, now make sure we are going in a direction of
    # descent
    d <- -solve(h,g)
    while(sum(d*g) > 0) {
      hm <- 2*hm
      h <- h + hm*I
      d <- -solve(h,g)
    }
    
    # always useful to track your target function value
    # to check for descent
    cE <- Energy(p)
    if (track_iteration)
      cat("iteration", i, "Energy =", cE, "\n")
    
    s <- 1
    while(cE < Energy(p + s*d))
      s <- s/2
    
    p <- p + s*d
    i <- i + 1
  }
  
  return (p)
}

