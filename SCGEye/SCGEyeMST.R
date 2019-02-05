## clear history
rm(list = ls(all = TRUE))
graphics.off()

# Libraries ---------------------------------------------------------------
## install and load packages
library(igraph)
library(lattice)
library(alphahull)


# Calculating Eye distribution --------------------------------------------

n = 500
theta = runif(n,0,2*pi)
r = sqrt(runif(n,0.25^2,0.5^2))
x = cbind(0.5+r*cos(theta),0.5+r*sin(theta))
# Value of alpha
alpha = 0
# alpha-shape
ashape.obj = ashape(x, alpha = alpha)


# Calculating MST  --------------------------------------------------------

df = as.data.frame(x)

n = 500
a = rnorm(n, 0.5, 0.05)
b = rnorm(n, 0.5, 0.05)
c = cbind(a, b)

d = rbind(x, c)
df = as.data.frame(d)


colnames(df) = c("theta", "r")

d = dist(cbind(df$theta, df$r))
matrix_of_distances <- as.matrix(d)
G = graph.adjacency(as.matrix(d), weighted = TRUE)
mst = minimum.spanning.tree(G)
edgelist = matrix(data = as.integer(igraph::get.edgelist(mst)), ncol = 2, byrow = F)



# Visualisations ----------------------------------------------------------

plot(df, pch = 16, col = "blue")
for (i in 1:dim(edgelist)[1]){
  from = edgelist[i, 1]
  to = edgelist[i, 2]
  lines(df[c(from, to), 'theta'], df[c(from, to), 'r'])
}

par(mfrow = c(1, 2), bg = NA)
hist(d)
boxplot(d)
