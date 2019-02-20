[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **SCGDTscatterplot** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of Quantlet: SCGDTscatterplot

Published in: SCG_Scagnostics

Description: Calculates and plots the Delaunay Triangulation for a Scatter Plot

Keywords: 'scagnostics, scatterplot, Delaunay Triangulation'

Author: Ioana Ceausu

Submitted:  Monday, January 2019 by Ioana Ceausu

Datafile: boston.xls

```

### R Code
```r

install.packages("tripack")
library(tripack)
library(MASS)

par(bg = NA)
plot(Boston$ptratio, Boston$tax, pch = 19, col = "blue", xlab = "zn", ylab = "dis",
     main = "ZN - DIS Scatterplot")

tdf = tri.mesh(Boston$ptratio, Boston$tax, duplicate = "remove")
plot(tdf, main = "")


```

automatically created on 2019-02-20