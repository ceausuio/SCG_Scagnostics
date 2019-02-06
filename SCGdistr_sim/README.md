[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SCGSimulating_Distributions** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml


Name of Quantlet: SCGSimulating_Distributions

Published in: SCG_Scagnostics

Description: Simulates different types of distributions and the respective plots; evaluation of scagnostic characteristics for the different types of distributions for multiple amounts of data points

Keywords: 'normal distribution, uniform distribution, linear distribution, heteroscedastic distribution, parable, clumps, eye distribution, sinus curve, integers, scagnostics'

Author: Ioana Ceausu

Submitted:  Tue, November 6 2018 by Ioana Ceausu


```

### R Code
```r

# Normal Distribution -----------------------------------------------------

n = 1000
x = rnorm(n)
y = rnorm(n)
png(filename = "Normal Distribution", bg = NA)
plot(x, y, pch = ".", col = "blue", main = "Normal")
dev.off()

# Uniform Distribution ----------------------------------------------------

n = 1000
x = runif(n)
y = runif(n)
png(filename = "Uniform Distribution", bg = NA)
plot(x, y, pch = ".", col = "blue", main = "Uniform")
dev.off()

# Linear Distribution -----------------------------------------------------

n = 1000
x = runif(n, 0, 10)
b = rnorm(n, 1, 2)
y = x + b
png(filename = "Linear Distribution", bg = NA)
plot(x, y, pch = ".", col = "blue", main = "Linear")
dev.off()

# Heteroscedastic Distributionz -------------------------------------------

n = rep(1:1000, 2)
a = 0
b = 1
sigma2 = n^1.3
eps = rnorm(n, mean = 0, sd = sqrt(sigma2))
y = a + b*n + 3*eps
mod <- lm(y ~ n)
png(filename = "Heteroscedastic Distribution", bg = NA)
plot(n, y, pch = ".", col = "blue", main = "Heteroscedastic")
dev.off()

# Parable Distribution ----------------------------------------------------

set.seed(1)
n = 1000
x1 = rnorm(n, 1, 25)
x2 = rnorm(n, 4, 85)

a = -0.7
b = 1
c = -100

y1 = a*x1^2 + b*x2 + c
png(filename = "Parable", bg = NA)
plot(x1, y1, pch = ".", col = "blue", main = "Parabola")
dev.off()

# Clumps Distribution -----------------------------------------------------

set.seed(1)
n = 1000
x1 = rnorm(n, -1, 0.1)
y1 = rnorm(n, -1, 0.1)
png(filename = "Clumps", bg = NA)
plot(x1, y1, pch = ".", xlim = c(-2, 2), ylim = c(-2, 2), col = "blue",
     main = "Clumps")

x2 = rnorm(n, 1, 0.1)
y2 = rnorm(n, -1, 0.1)

points(x2, y2, pch = ".", xlim = c(-2, 2), ylim = c(-2, 2), col = "blue")

x3 = rnorm(n, 0, 0.1)
y3 = rnorm(n, 1, 0.1)

points(x3, y3, pch = ".", xlim = c(-2, 2), ylim = c(-2, 2), col = "blue")
dev.off()

# “Eye” Distribution ------------------------------------------------------

set.seed(1)
n = 1000
Y1 = rnorm(n,0, 0.1)
Y2 = rnorm(n,0, 0.1)
png(filename = "Eye Distribution", bg = NA)
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
dev.off()

# AR(1) Process Distribution ----------------------------------------------

n = arima.sim(list(order=c(1,0,0), ar=.5), n= 200) + 10
png(filename = "AR(1) Process", bg = NA)
plot(n, pch = ".", col = "blue", main = "AR(1) process")
dev.off()

# Sinus curve -------------------------------------------------------------

n = 10000
x = rnorm(n)
y = sin(x)
png(filename = "Sinus Curve", bg = NA)
plot(x, y, pch = ".", xlab = "x", ylab = "y", col = "blue", main = "Sinus curve")
dev.off()

# Integers ----------------------------------------------------------------

n = rep(1:1000, 1000)
a = 0
b = 1
sigma2 = n^1.3
eps = rnorm(n, mean = 0, sd = sqrt(sigma2))
y = a + b*n + 2*eps
mod <- lm(y ~ n)
png(filename = "Integers", bg = NA)
plot(n, y, pch = ".", xlab = "x", ylab = "y", col = "blue",
     xlim = c(0,20), ylim = c(0,20), main = "Integers")
dev.off()


```

automatically created on 2019-02-06