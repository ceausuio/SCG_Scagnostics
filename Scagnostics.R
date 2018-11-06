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
pairs(boston, pch = ".", bg = "grey", col = "blue", main = "Scatterplot matrix (SPLOM) of Boston Housing Dataset")

s <- as.matrix(s)
t(s)
s1 = as.matrix(t(s))
pairs(s1)
par(bg = NA)
pairs(s1, pch = ".", col = "blue",main = "Scagnostics of Boston Housing Dataset")
par(bg = NA)
pairs.panels(t, hist.col="red", stars=TRUE)
