nn <- read.table("nn.txt", header = T)

attach(nn)
plot(x1, x2, col = (2 - y))

get_beta <- function(alpha, i) {
  return(alpha[(2*i-1):(2*i)])
}

get_gamma <- function(alpha, i, m) {
  return(alpha[((i-1)*m+1):(i*m)])
}

omega <- function(w) {
  return(1/(1+exp(-w)))
}

NN <- function(x, a, m) {
  alpha_beta <- a[1:(2*m)]
  alpha_beta_0 <- a[(2*m+1):(3*m)]
  alpha_gamma <- a[(3*m+1):(5*m)]
  alpha_gamma_0 <- a[(5*m+1):(5*m+2)]
  z <- matrix(0, ncol = m, nrow = 1947)
  for (i in 1:m) {
    z[,i] <- omega(alpha_beta_0[i] + get_beta(alpha_beta, i)[1] * x[,1] + 
                     get_beta(alpha_beta, i)[2] * x[,2])
  }
  t <- matrix(0, ncol = 2, nrow = 1947)
  for (i in 1:2) {
    for (j in 1:m) {
      t[,i] <- t[,i] + get_gamma(alpha_gamma, i, m)[j] * z[,j]
    }
    t[,i] <- omega(alpha_gamma_0[i] + t[,i])
  }
  p <- NULL
  for (i in 1:1947) {
    p[i] <- exp(t[i,1]) / (exp(t[i,1]) + exp(t[i,2]))
  }
  return(p)
}

nn <- as.matrix(nn)

x <- nn[,1:2]

m <- 4

set.seed(123)
a <- runif((5*m+2),0,1)

logL <- function(x, a, m) {
  sum <- 0
  p <- NN(x, a, m)
  for (i in 1:1947) {
    sum <- sum + (1 - nn[i,3]) * log(p[i]) + nn[i,3] * log(1-p[i])
  }
  return(sum)
}

logL(x, a, m)

grad_logL <- function(x, a, m) {
  result <- rep(0, (5*m+2))
  h <- 0.0001
  a_new <- a
  for (i in 1:(5*m+2)) {
    a_new[i] <- a[i] + h
    result[i] <- (logL(x, a_new, m) - logL(x, a, m)) / h
    a_new <- a
  }
  return(result)
}

grad_logL(x, a, m)

norm <- function(x)
{
  sqrt(sum(x^2))
}

steepest_ascent_with_backtrack <- function(x, a, m, 
                                           epsilon=10^-3, max_iter=5000) {
  iter <- 0
  
  while(norm(grad_logL(x, a, m)) > epsilon & iter < max_iter) {
    iter <- iter + 1
    gf <- grad_logL(x, a, m)
    
    # direction, see above for definition of norm
    d <- gf/norm(gf)
    
    # step size with backtracking
    s <- 1
    while(logL(x, a+s*d, m) < logL(x, a, m))
      s <- s/2
    
    # update
    a <- a + s*d
    
  }
  
  return(list(a = a, iter = iter))
  
}

result <- steepest_ascent_with_backtrack(x, a, m)
a_best <- result$a

F_value <- NULL
p <- seq(0.2, 0.8, 0.1)
p_final <- NN(x, a_best, m)
for (i in 1:7) {
  F_value <- ifelse(p_final > p[i], 1, 0)
  plot(x[,1], x[,2], col = (2 - F_value))
}

library(igraph)

A <- as.matrix(read.table("karate club.txt", header = F))

g = graph.adjacency(A, mode = "undirected", weighted = NULL, diag = FALSE)

k <- rowSums(A)
m <- sum(k)/2
B <- A - 1/(2*m) * k %*% t(k)

n <- nrow(B)
V <- diag(n)

for (i in 1:1000) {
  V <- B %*% V
  V <- qr.Q(qr(V))
}

lambda <- rep(0, n) 
for (i in 1:n) {
  temp <- V[, i]
  lambda[i] <- t(temp) %*% B %*% temp
}

order_lam <- order(lambda, decreasing = T) 
V <- V[, order_lam]
lambda <- lambda[order_lam]

# eigen(B)

max_ev <- V[, 1]
s <- ifelse(max_ev > 0, 1, -1) 
color <- ifelse(s == 1, "black", "red") 
shape <- ifelse(s == 1, "square", "circle")

plot.igraph(g, vertex.color = color, vertex.shape = shape, vertex.label = NA)

MyKmeans <- function(X, K, centroids, maxIter = 20, epsilon = 1e-5) {
  n <- dim(X)[1]
  for (i in 1:maxIter) {
    dists <- apply(X, 1, function(point)
      sapply(1:nrow(centroids), function(dim)
        dist(rbind(point, centroids[dim, ]))))
    dists <- t(dists)
    clusters <- sapply(1:n, function(x) which.min(dists[x, ]))
    # just for 2-d plot
    plot(X[,1], X[,2], col = clusters + 1)
    points(centroids, col = "black", pch = 19)
    new <- t(sapply(1:K, function(c) colMeans(X[which(clusters == c), ])))
    temp <- sum((new - centroids) ^ 2)
    if (temp > epsilon) {
      centroids <- new
    } else {
      break
    }
  }
  return(list(clusters = clusters, iteration = i))
}

X <- as.matrix(read.csv("synthetic kmeans data.csv", header = T))

plot(X[,1], X[,2])

n <- dim(X)[1]

k <- 2
set.seed(123)
ids <- sample(1:n, k)
centroids <- X[ids, ]

MyKmeans(X, k, centroids)

X <- as.matrix(read.csv("tumor microarray data.csv", header = T))

label <- X[,1]

X <- X[,-1]

for (i in 1:10) {
  result <- kmeans(X, i, nstart = 20)
  cluster <- result$cluster
  print(table(cluster, label))
}