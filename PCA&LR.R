# PCA; Logistic Reg


d <- read.csv("expectancy.csv", header=T)

  
L <- c(1, cumprod(1-d$probability[41:101]))


max_t <- 101 - 40
spline_L <- splinefun(x=0:max_t, y=L)

# produce a bunch of interpolated values to graph
t <- seq(0, max_t, .01)
L_interpolated <- spline_L(t, deriv=0)

plot(40 + 0:max_t, L)
lines(40 + t, L_interpolated, col="red")

final_index <- 12*max_t

# first use the interpolation
total <- rep(0, final_index)
for (i in 1:final_index) {
  if (i==1)
    total[i] <- 200*spline_L(i/12)*exp(-.05*(i)/12)
  else
    total[i] <- total[i-1] + 200*spline_L(i/12)*exp(-.05*(i)/12)
}

# now without the interpolation.  Replace interpolation with
# the probability of living to a given year

raw_total <- rep(0, final_index)
for (i in 1:final_index) {
  # probability of living through the year floor(i/12) gives the number
  # of years after 40, then add 1 since indices are shifted
  p <- L[floor(i/12)+1]
  if (i==1)
    raw_total[i] <- 200*p*exp(-.05*(i)/12)
  else
    raw_total[i] <- raw_total[i-1] + 200*p*exp(-.05*(i)/12)
}

# total with interpolation
total[final_index]
# total without interpolation
raw_total[final_index]

# plot interpolated totals
plot(40+(1:final_index)/12, total/1000, type="l", col="red",
     xlab="year", 
     ylab="total ($1000)")

# PCA

mtrain <- as.matrix(read.csv("mnist_train.csv", header=F))

numbers <- mtrain[,1]

# create response vector
y <- rep(0, nrow(mtrain))
for (i in 1:nrow(mtrain)) {
  if (mtrain[i,1]==3)
    y[i] <- 1
}

X <- mtrain[,2:ncol(mtrain)]
X <- X/250

# subtract off the mean
mu <- rep(0, ncol(X))
for (i in 1:nrow(X)) {
  mu <- mu + X[i,]
}
mu <- mu/nrow(X)

for (i in 1:nrow(X)) {
  X[i,] <- X[i,] - mu
}

# do PCA on only the first 1000 components
Xpca <- X[1:1000,]
ynumberspca <- numbers[1:1000]


Theta <- t(Xpca) %*% Xpca
ev <- eigen(Theta)

# get the first two eigenvectors
v1 <- ev$vectors[,1]
v2 <- ev$vectors[,2]

# compute projections
proj1 <- Xpca %*% v1
proj2 <- Xpca %*% v2

plot(proj1, proj2, col=numbers)

eigenvalues <- ev$values
fraction <- cumsum(eigenvalues)/sum(eigenvalues)
plot(1:length(ev$values), fraction,
     xlab="eigenvalues number", ylab="variance captured")

# looks like roughly 80 captures 90 percent
fraction[80]


show_image <- function(m, oriented=T)
{
  im <- matrix(m, byrow=T, nrow=28)
  
  if (oriented) {
    im_orient <- matrix(0, nrow=28, ncol=28)
    for (i in 1:28)
      im_orient[i,] <- rev(im[,i])
    
    im <- im_orient
  }
  image(im)
}

show_image(X[13,])

# image with 2 principle components
# compute the coefficients
V <- ev$vectors[,1:2]
proj_c <- t(V) %*% X[13,]
# then the projection linearly combines the principle components
Xp <- V %*% proj_c
show_image(Xp)

# image with 10 principle components
# compute the coefficients
V <- ev$vectors[,1:10]
proj_c <- t(V) %*% X[13,]
# then the projection linearly combines the principle components
Xp <- V %*% proj_c
show_image(Xp)

# image with 80 principle components
# compute the coefficients
V <- ev$vectors[,1:80]
proj_c <- t(V) %*% X[13,]
# then the projection linearly combines the principle components
Xp <- V %*% proj_c
show_image(Xp)

# use 80-d PCA
n <- 80
Q <- ev$vectors[,1:n]

# compute the principle component coefficients for all samples in X 

tQ <- t(Q)
pca_coef <- X %*% Q

source("logistic_regression_code.R")
Xpca <- cbind(rep(1, nrow(pca_coef)), pca_coef)
alpha <- rep(0, ncol(Xpca))
alpha_star <- Newton(alpha,Xpca,y)
alpha_star

mtest <- as.matrix(read.csv("mnist_test.csv", header=F))
y_test <- rep(0, nrow(mtest))
for (i in 1:nrow(mtest)) {
  if (mtest[i,1]==3)
    y_test[i] <- 1
}

X_test <- mtest[,2:ncol(mtest)]
X_test <- X_test/250
tQ <- t(Q)
pca_coef_test <- matrix(0, nrow=nrow(X_test), ncol=n)
for (i in 1:nrow(X_test)) {
  pca_coef_test[i,] <- tQ %*% X_test[i,]
}

# generate the pca_coef
X_test_pca <- cbind(rep(1,nrow(pca_coef_test)), pca_coef_test)


# compute sensitivity and specificity.
accuracy <- function(p, alpha_star, X, y)
{
  total_tests <- length(y)
  called_3_right <- 0
  called_not_3_right <- 0
  total_3_tests <- 0
  total_not_3_tests <- 0
  
  for (i in 1:total_tests) {
    prob <- 1/(1 + exp(-sum(alpha_star*X[i,])))
    if (prob >= p)
      pred_y <- 1
    else
      pred_y <- 0
    
    if (y[i] == 1 & pred_y == 1)
      called_3_right <- called_3_right + 1
    if (y[i] == 0 & pred_y == 0)
      called_not_3_right <- called_not_3_right + 1
    
    if (y[i]==1)
      total_3_tests <- total_3_tests + 1
    else
      total_not_3_tests <- total_not_3_tests + 1
  }
  
  sens <- called_3_right/total_3_tests
  spec <- called_not_3_right/total_not_3_tests
  overall <- (called_3_right + called_not_3_right)/total_tests
  
  return (list(sensitivity=sens, specificity=spec, overall=overall))
}

p <- seq(0.05, .96, .05)
specs <- rep(0, length(p))
sens <- rep(0, length(p))
overall <- rep(0, length(p))
for (i in 1:length(p)) {
  acc_out <- accuracy(p[i], alpha_star, X_test_pca, y_test)
  specs[i] <- acc_out$specificity
  sens[i] <- acc_out$sensitivity
  overall[i] <- acc_out$overall
}

df <- data.frame(p=p, sensitivity=sens, specificiy=specs, overall=overall)
df



