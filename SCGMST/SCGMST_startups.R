require(stats)
library(igraph)
library(readxl)



#crunchbase database
#EDP_10years = read_excel("/Users/ceausuio/SCG_Startups/20190403_EDP_short_startups.xlsx")

#View(EDP_10years)



edp = as.data.frame(EDP_10years)
stringy = cbind(edp$fins_wg_m1, edp$inv_pln_oe_3y)
d = dist(cbind(edp$fins_wg_m1, edp$inv_pln_oe_3y))

dev.off()
par(bg = NA)
plot(stringy, col = "dark blue")



matrix_of_distances <- as.matrix(d)
G = graph.adjacency(as.matrix(d), weighted = TRUE)
mst = minimum.spanning.tree(G)
edgelist = matrix(data = as.integer(igraph::get.edgelist(mst)), ncol = 2, byrow = F)

df = edp[, c('fins_wg_m1', 'inv_pln_oe_3y')]

par(bg = NA)
plot(df, pch = 16, col = "blue")
for (i in 1:dim(edgelist)[1]){
  from = edgelist[i, 1]
  to = edgelist[i, 2]
  lines(df[c(from, to), 'fins_wg_m1'], df[c(from, to), 'inv_pln_oe_3y'])
}

boxplot(d)

#########


sparse = cbind(edp$strup_age, edp$inv_pln_oe_3y)
d1= dist(cbind(edp$strup_age, edp$inv_pln_oe_3y))

dev.off()
par(bg = NA)
plot(sparse, col = "dark blue")



matrix_of_distances <- as.matrix(d1)
G1 = graph.adjacency(as.matrix(d1), weighted = TRUE)
mst = minimum.spanning.tree(G1)
edgelist = matrix(data = as.integer(igraph::get.edgelist(mst)), ncol = 2, byrow = F)

df1 = edp[, c('strup_age', 'inv_pln_oe_3y')]

par(bg = NA)
plot(df1, pch = 16, col = "blue")
for (i in 1:dim(edgelist)[1]){
  from = edgelist[i, 1]
  to = edgelist[i, 2]
  lines(df1[c(from, to), 'strup_age'], df1[c(from, to), 'inv_pln_oe_3y'])
}

boxplot(d1)
