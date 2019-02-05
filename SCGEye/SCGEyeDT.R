
# Libraries ---------------------------------------------------------------
library(alphahull)
library(RTriangle)


# Calculating the Eye distribution ----------------------------------------

n = 300
theta = runif(n,0,2*pi)
r = sqrt(runif(n,0.25^2,0.5^2))
x = cbind(0.5+r*cos(theta),0.5+r*sin(theta))
# Value of alpha
alpha = 0
# alpha-shape
ashape.obj = ashape(x, alpha = alpha)

# Vizualisation -----------------------------------------------------------

par(bg = NA)
plot(ashape.obj, col = "blue", pch = 19, xlab = "", ylab = "")

plot(df1, pch = 19, col = "blue")

df = as.data.frame(x)

n = 300
a = rnorm(n, 0.5, 0.05)
b = rnorm(n, 0.5, 0.05)
c = cbind(a, b)


df1 = rbind(x, c)

pdf = pslg(df1)

pdf1 = pslg(df1)


plot(pdf1, col = "blue", pch = 19)
tp <- triangulate(pdf)
plot(tp)


