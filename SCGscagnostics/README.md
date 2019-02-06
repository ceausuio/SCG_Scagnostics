[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SCGScagnostics** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml


Name of Quantlet: SCGScagnostics

Published in: SCG_Scagnostics

Description: "Calculates the scagnostic measures for the dataset and plots the SPLOM, the scagnostics SPLOM and the heat-map of the scagnostic measures"

Keywords: 'scagnostics, scagnostic coefficients, SPLOM, heat-map'

Author: Ioana Ceausu

Submitted:  Tue, November 6 2018 by Ioana Ceausu

Datafile: boston.xls

```

![Picture1](Heatmap%20BHD.png)

![Picture2](Heatmap%20PC%20Factors.png)

![Picture3](SPLOM%20BHD.png)

![Picture4](Scagnostics%20SPLOM%20BHD.png)

![Picture5](Scagnostics%20SPLOM%20PC%20Factors.png)

### R Code
```r


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

png(filename = "SPLOM BHD.png", bg = NA)
pairs(boston, pch = ".", col = "blue",
      main = "Scatterplot matrix (SPLOM) of Boston Housing Dataset")
dev.off()


# Plotting the SPLOM for the Scagnostic measures --------------------------

s = as.matrix(s)
t(s)
s1 = as.matrix(t(s))
pairs(s1)
par(bg = NA)

png(filename = "Scagnostics SPLOM BHD.png")
pairs(s1, pch = ".", col = "blue",main = "Scagnostics of Boston Housing Dataset")
dev.off()


# Example with PC Factors ------------------------------------------------

png(filename = "SPLOM BHD.png", bg = NA)
pairs(factors_scag_gr, pch = ".", col = "green",
      main = "Scatterplot matrix (SPLOM) of PC Factors")
dev.off()

scagnostics(factors_scag_gr)
s = scagnostics(factors_scag_gr)
s1 = as.matrix(t(s))

png(filename = "Scagnostics SPLOM PC Factors.png")
pairs(s1, pch = ".", col = "green", main = "Scagnostics of PC Factors")
dev.off()

```

automatically created on 2019-02-06