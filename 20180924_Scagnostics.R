install.packages("rJava")
install.packages("scagnostics")
library(scagnostics)
scagnostics(boston)
s <- scagnostics(boston)
o <- scagnosticsOutliers(s)
o[o]
g <- scagnosticsGrid(s)
go <- g[o,]
plot(boston[go$x])

e <- scagnosticsExemplars(s)
e[e]
ge = g[e,]
par(mfrow = c(2,2))
for (i in 1:dim(ge)[1])
  plot(boston[[ge$x[i]]], boston[[ge$y[i]]], pch=".", xlab=names(boston)[ge$x[i]], ylab=names(boston)[ge$y[i]])

install.packages("psych")
library(psych)
pairs(boston)
par(bg = NA)
pairs(boston, pch=".", main="Boston Housing Data data by CRIM")

pairs.panels(boston)
pairs.panels(boston,bg=c("red","yellow","blue")[boston$CRIM], pch=21,main="Fisher Iris data by CRIM")
pairs.panels(boston,bg=c("red","yellow","blue")[boston$CRIM], pch=21,main="Boston Housing Data data by CRIM")
pairs.panels(boston,bg=c("red","yellow","blue")[boston$CRIM], pch=21+as.numeric(boston$CRIM),main="Boston Housing Data data by CRIM",hist.col="red")
pairs.panels(boston,bg=c("red","yellow","blue")[boston$CRIM], pch=21+as.numeric(boston$CRIM),main="Boston Housing Data data by CRIM",hist.col="red",stars=TRUE)
pairs.panels(boston,bg=c("red","yellow","blue")[boston$CRIM], pch=21+as.numeric(boston$CRIM),main="Boston Housing Data data by CRIM",hist.col="red",stars=TRUE)

s <- as.matrix(s)
t(s)
s1 = as.matrix(t(s))
pairs(s1)
par(bg = NA)
pairs(s1, pch = ".", col = "blue",main = "Scagnostics of Boston Housing Dataset")
par(bg = NA)
pairs(boston, pch = ".", bg = "grey", col = "blue", main = "Scatterplot matrix (SPLOM) of Boston Housing Dataset")
par(bg = NA)
pairs.panels(t, hist.col="red", stars=TRUE)

#normal distribution
n = 1000
x = rnorm(n)
y = rnorm(n)
par(bg = NA)
plot(x, y, pch = ".", col = "blue", main = "Normal")

#uniform distribution
n = 1000
x = runif(n)
y = runif(n)
par(bg = NA)
plot(x, y, pch = ".", col = "blue", main = "Uniform")

#linear distribution
n = 1000
x = runif(n, 0, 10)
b = rnorm(n, 1, 2)
y = x + b
par(bg = NA)
plot(x, y, pch = ".", col = "blue", main = "Linear")

#heteroscedastic distribution
n = rep(1:1000, 2)
a = 0
b = 1
sigma2 = n^1.3
eps = rnorm(n, mean = 0, sd = sqrt(sigma2))
y = a + b*n + 3*eps
mod <- lm(y ~ n)
par(bg = NA)
plot(n, y, pch = ".", col = "blue", main = "Heteroscedastic")

#Parable Distribution
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

#Clumps distribution
set.seed(1)
n = 1000
x1 = rnorm(n, -1, 0.1)
y1 = rnorm(n, -1, 0.1)
par(bg = NA)
plot(x1, y1, pch = ".", xlim = c(-2, 2), ylim = c(-2, 2), col = "blue", main = "Clumps")

x2 = rnorm(n, 1, 0.1)
y2 = rnorm(n, -1, 0.1)

points(x2, y2, pch = ".", xlim = c(-2, 2), ylim = c(-2, 2), col = "blue")

x3 = rnorm(n, 0, 0.1)
y3 = rnorm(n, 1, 0.1)

points(x3, y3, pch = ".", xlim = c(-2, 2), ylim = c(-2, 2), col = "blue")

#"eye" distribution
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
# plot(cbind(c(d,-d),c(d,-d)),type="n",xlab="", ylab="", main=substitute(list(d==d1,  rho==r1), list(d1=d,r1=sprintf("%0.2f",r))))
points(ellipse, pch = ".", col = "blue", xlab = "x", ylab = "y")

#AR(1) process distribution
n = arima.sim(list(order=c(1,0,0), ar=.5), n= 200) + 10
par(bg = NA)
plot(n, pch = ".", col = "blue", main = "AR(1) process")

#sinus curve
n = 10000
x = rnorm(n)
y = sin(x)
par(bg = NA)
plot(x, y, pch = ".", xlab = "x", ylab = "y", col = "blue", main = "Sinus curve")

#integers
n = rep(1:1000, 1000)
a = 0
b = 1
sigma2 = n^1.3
eps = rnorm(n, mean = 0, sd = sqrt(sigma2))
y = a + b*n + 2*eps
mod <- lm(y ~ n)
par(bg = NA)
plot(n, y, pch = ".", xlab = "x", ylab = "y", col = "blue", xlim = c(0,20), ylim = c(0,20), main = "Integers")

#plot alpha shape 1
x1 = 0
y1 = 2

par(bg = NA)
plot(x1, y1, pch = 19, col = "blue", xlab = "x", ylab = "y", xlim = c(-1, 5), ylim = c(0, 4))

x2 = 1.1
y2 = 1.5

points(x2, y2, pch = 19, col = "blue")

x3 = 1
y3 = 2.7

points(x3, y3, pch = 19, col = "blue")

x4 = 2
y4 = 2.4

points(x4, y4, pch = 19, col = "blue")

x5 = 3
y5 = 1.2

points(x5, y5, pch = 19, col = "blue")

x6 = 3.3
y6 = 2.5

points(x6, y6, pch = 19, col = "blue")


#plot alpha shapes
#alpha shape I - alpha  = 0
## Not run:
# Uniform sample of size n=300 in the annulus B(c, 0.5)\B(c, 0.25)
# with c=(0.5, 0.5).
n = 500
theta = runif(n,0,2*pi)
r = sqrt(runif(n,0.25^2,0.5^2))
x = cbind(0.5+r*cos(theta),0.5+r*sin(theta))
# Value of alpha
alpha = 0
# alpha-shape
ashape.obj = ashape(x, alpha = alpha)
# Plot alpha-shape in blue, sample points in black,
# and Delaunay triangulation in red
par(bg = NA)
plot(ashape.obj, col = c(4, 1, 2), pch = 19)

n = 100
a = rnorm(n, 0.5, 0.05)
b = rnorm(n, 0.5, 0.05)

points (a, b, pch = 19)

#alpha shape I - alpha  = 1
## Not run:
# Uniform sample of size n=300 in the annulus B(c, 0.5)\B(c, 0.25)
# with c=(0.5, 0.5).
n = 500
theta = runif(n,0,2*pi)
r = sqrt(runif(n,0.25^2,0.5^2))
x = cbind(0.5+r*cos(theta),0.5+r*sin(theta))
# Value of alpha
alpha = 1
# alpha-shape
ashape.obj = ashape(x, alpha = alpha)
# Plot alpha-shape in blue, sample points in black,
# and Delaunay triangulation in red
par(bg = NA)
plot(ashape.obj, col = c(4, 1, 2), pch = 19)

n = 100
a = rnorm(n, 0.5, 0.05)
b = rnorm(n, 0.5, 0.05)

points (a, b, pch = 19)

#alpha shape I - alpha  = 1
## Not run:
# Uniform sample of size n=300 in the annulus B(c, 0.5)\B(c, 0.25)
# with c=(0.5, 0.5).
n = 500
theta = runif(n,0,2*pi)
r = sqrt(runif(n,0.25^2,0.5^2))
x = cbind(0.5+r*cos(theta),0.5+r*sin(theta))
# Value of alpha
alpha = 0.2
# alpha-shape
ashape.obj = ashape(x, alpha = alpha)
# Plot alpha-shape in blue, sample points in black,
# and Delaunay triangulation in red
par(bg = NA)
plot(ashape.obj, col = c(4, 1, 2), pch = 19)

n = 100
#a = rnorm(n, 0.5, 0.05)
#b = rnorm(n, 0.5, 0.05)
c = matrix(rnorm(n, 0.5, 0.05), ncol = 2); colnames(c)=c("x","y")
alpha = 1
ashape.obj = ashape(c, alpha = alpha)
#points(ashape.obj$x)
plot(ashape.obj, add=T, pch=19, col = c(4, 1, 2))
#points(c, pch = 19)
