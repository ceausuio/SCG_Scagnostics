[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SCGMST** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml


Name of Quantlet: SCGMST

Published in: SCG_Scagnostics

Description: Calculates and plots a Minimum Spanning Tree from a Scatter Plot

Keywords: 'minimum spanning tree, mst, scatter plot, scagnostics'

Author: Ioana Ceausu, Elizaveta Zinovyeva

Submitted:  Tue, November 6 2018 by Ioana Ceausu


```

### R Code
```r

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

```

automatically created on 2019-02-06