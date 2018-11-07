# Packages and libraries needed -------------------------------------------

library("MASS")
require(stats)
library(igraph)

df = Boston
#plot(df$dis, df$medv)
d = dist(cbind(df$dis, df$medv))
matrix_of_distances <- as.matrix(d)
G = graph.adjacency(as.matrix(d), weighted = TRUE)

## Some graphical parameters
## MST and plot
mst = minimum.spanning.tree(G)
edgelist = matrix(data = as.integer(igraph::get.edgelist(mst)), ncol = 2, byrow = F)

df = Boston[, c('dis', 'medv')]

plot(df, pch = 16, col = "blue")
for (i in 1:dim(edgelist)[1]){
  from = edgelist[i, 1]
  to = edgelist[i, 2]
  lines(df[c(from, to), 'dis'], df[c(from, to), 'medv'])
}
