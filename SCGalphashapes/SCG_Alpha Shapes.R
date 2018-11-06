# Plot Alpha Shape 1 ------------------------------------------------------

x1 = 0
y1 = 2

par(bg = NA)
plot(x1, y1, pch = 19, col = "blue", xlab = "x", ylab = "y",
     xlim = c(-1, 5), ylim = c(0, 4))

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


# Plot Alpha Shapes -------------------------------------------------------


#alpha shape I - alpha  = 0

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
c = matrix(rnorm(n, 0.5, 0.05), ncol = 2); colnames(c)=c("x","y")
alpha = 1
ashape.obj = ashape(c, alpha = alpha)

#points(ashape.obj$x)
plot(ashape.obj, add=T, pch=19, col = c(4, 1, 2))


