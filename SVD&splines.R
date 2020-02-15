library(splines, quietly = T)
bone_mass <- read.table("BoneMassData.txt",
                        header=TRUE,
                        stringsAsFactors = F)
# using dplyr to filter for only female samples
library(dplyr, quietly = T, warn.conflicts = F)
bone_mass <- filter(bone_mass,gender=="female")
x <- bone_mass$age
y <- bone_mass$spnbmd

mygrid <- seq(min(x),max(x),.001)
myknots <- seq(min(x), max(x), length.out = 1000)
Bpp <- splineDesign(knots=myknots,x=mygrid,derivs=2,outer.ok=T)
B <- splineDesign(knots=myknots,x=x,derivs=0,outer.ok=T)

dx <- 1/ncol(B)
omega <- t(Bpp) %*% (dx*Bpp)


fit_and_plot <- function(rho)
{
  alpha <- solve(t(B) %*% B + rho*omega, t(B) %*% y)
  plot(x, y)

  Bplot <- splineDesign(knots=myknots,x=mygrid,derivs=0,outer.ok=T)
  y_pred <- Bplot %*% alpha
  lines(mygrid, y_pred, col="red", lwd=2)
}

# rho=.01
fit_and_plot(.01)
# rho = 1
fit_and_plot(1)
# rho = 100
fit_and_plot(100)

# SVD
norm <- function(x) sqrt(sum(x^2))

A <- matrix(c(1:8,10), nrow=3)
A


V <- eigen(t(A) %*% A)$vectors

S <- matrix(0, nrow=3, ncol=3)
U <- matrix(0, nrow=3, ncol=3)
for (i in 1:3) {
  S[i,i] <- norm(A %*% V[,i])
  U[,i] <- (A %*% V[,i])/S[i,i]
}
V
S
U
# check that A = USV^T
U %*% S %*% t(V)


# compute u3 using gramm-schmidt
u3_start <- c(1,0,0)
u3 <- u3_start - sum(u3_start*U[,1])*U[,1] - sum(u3_start*U[,2])*U[,2]
u3 <- u3/norm(u3)

# add u3 as a column to U and a zero row to S
U <- cbind(U, u3)
S <- rbind(S, c(0,0))


# check that A = USV^T
U %*% S %*% t(V)



