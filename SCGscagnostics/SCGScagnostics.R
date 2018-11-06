
# Packages needed ---------------------------------------------------------

install.packages("rJava")
install.packages("scagnostics")
install.packages("psych")

# Library -----------------------------------------------------------------

library(scagnostics)
library(psych)


# Calculating the scagnostics  --------------------------------------------

scagnostics(boston)
s = scagnostics(boston)
o = scagnosticsOutliers(s)
o[o]
g = scagnosticsGrid(s)
go = g[o,]
plot(boston[go$x])

e = scagnosticsExemplars(s)
e[e]
ge = g[e,]
par(mfrow = c(2,2))
for (i in 1:dim(ge)[1])
  plot(boston[[ge$x[i]]], boston[[ge$y[i]]], pch=".",
       xlab=names(boston)[ge$x[i]], ylab=names(boston)[ge$y[i]])


# Plotting the SPLOM  -----------------------------------------------------

pairs(boston)
par(bg = NA)

png(filename = "SPLOM BHD.png")
pairs(boston, pch = ".", bg = "grey", col = "blue",
      main = "Scatterplot matrix (SPLOM) of Boston Housing Dataset")
dev.off()


# Plotting the SPLOM for the Scagnostic measures --------------------------

s <- as.matrix(s)
t(s)
s1 = as.matrix(t(s))
pairs(s1)
par(bg = NA)

png(filename = "Scagnostics SPLOM BHD.png")
pairs(s1, pch = ".", col = "blue",main = "Scagnostics of Boston Housing Dataset")
dev.off()

