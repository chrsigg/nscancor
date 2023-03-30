x <- matrix(runif(10*5), 10)
y <- matrix(runif(10*5), 10)
z <- matrix(runif(10*5), 10)

xcoef <- matrix(rnorm(2*5), 5)
ycoef <- matrix(rnorm(2*5), 5)
zcoef <- matrix(rnorm(2*5), 5)

# Explained multi-domain correlation
macor(list(x, y, z), list(xcoef, ycoef, zcoef))$cor
