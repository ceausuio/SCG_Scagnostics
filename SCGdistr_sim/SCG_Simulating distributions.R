# Normal Distribution -----------------------------------------------------

n = 1000
x = rnorm(n)
y = rnorm(n)
par(bg = NA)
plot(x, y, pch = ".", col = "blue", main = "Normal")

# Uniform Distribution ----------------------------------------------------

n = 1000
x = runif(n)
y = runif(n)
par(bg = NA)
plot(x, y, pch = ".", col = "blue", main = "Uniform")


# Linear Distribution -----------------------------------------------------

n = 1000
x = runif(n, 0, 10)
b = rnorm(n, 1, 2)
y = x + b
par(bg = NA)
plot(x, y, pch = ".", col = "blue", main = "Linear")

# Heteroscedastic Distributionz -------------------------------------------

n = rep(1:1000, 2)
a = 0
b = 1
sigma2 = n^1.3
eps = rnorm(n, mean = 0, sd = sqrt(sigma2))
y = a + b*n + 3*eps
mod <- lm(y ~ n)
par(bg = NA)
plot(n, y, pch = ".", col = "blue", main = "Heteroscedastic")


# Parable Distribution ----------------------------------------------------

set.seed(1)
n = 1000
x1 = rnorm(n, 1, 25)
x2 = rnorm(n, 4, 85)

a = -0.7
b = 1
c = -100

y1 = a*x1^2 + b*x2 + c
par(bg = NA)
plot(x1, y1, pch = ".", col = "blue", main = "Parabola")


# Clumps Distribution -----------------------------------------------------

set.seed(1)
n = 1000
x1 = rnorm(n, -1, 0.1)
y1 = rnorm(n, -1, 0.1)
par(bg = NA)
plot(x1, y1, pch = ".", xlim = c(-2, 2), ylim = c(-2, 2), col = "blue",
     main = "Clumps")

x2 = rnorm(n, 1, 0.1)
y2 = rnorm(n, -1, 0.1)

points(x2, y2, pch = ".", xlim = c(-2, 2), ylim = c(-2, 2), col = "blue")

x3 = rnorm(n, 0, 0.1)
y3 = rnorm(n, 1, 0.1)

points(x3, y3, pch = ".", xlim = c(-2, 2), ylim = c(-2, 2), col = "blue")


# “Eye” Distribution ------------------------------------------------------

set.seed(1)
n = 1000
Y1 = rnorm(n,0, 0.1)
Y2 = rnorm(n,0, 0.1)
par(bg = NA)
plot(Y1,Y2, pch = ".",ylim=c(-6,6),xlim=c(-6,6), col = "blue", main = "Eye")
rhovec  = c(0)  # direction
dvec    = c(4)           # radius
r   = rhovec[1]
d   = dvec[1]
c   = c(0,0)
c1  = c(1,r)
c2  = c(r,1)
sha = cbind(c1,c2)         # covariance matrix used for the cholesky decomposition within the ellipse function
ucircle = cbind(c(cos(seq(0,360,length=n))),c(sin(seq(0,360,length=n))))
Q       = chol(sha);                 # cholesky decomposition for reshape the unit circle
radius  = d;
D       = radius * (ucircle %*% Q);  # change directions via cholesky decomposition
ellipse = (c + D + cbind(Y1,Y2))
points(ellipse, pch = ".", col = "blue", xlab = "x", ylab = "y")


# AR(1) Process Distribution ----------------------------------------------

n = arima.sim(list(order=c(1,0,0), ar=.5), n= 200) + 10
par(bg = NA)
plot(n, pch = ".", col = "blue", main = "AR(1) process")


# Sinus curve -------------------------------------------------------------

n = 10000
x = rnorm(n)
y = sin(x)
par(bg = NA)
plot(x, y, pch = ".", xlab = "x", ylab = "y", col = "blue", main = "Sinus curve")


# Integers ----------------------------------------------------------------

n = rep(1:1000, 1000)
a = 0
b = 1
sigma2 = n^1.3
eps = rnorm(n, mean = 0, sd = sqrt(sigma2))
y = a + b*n + 2*eps
mod <- lm(y ~ n)
par(bg = NA)
plot(n, y, pch = ".", xlab = "x", ylab = "y", col = "blue",
     xlim = c(0,20), ylim = c(0,20), main = "Integers")

