# Packages needed ---------------------------------------------------------

install.packages("rJava")
install.packages("scagnostics")
install.packages("gplots")
install.packages("MASS")

# Library -----------------------------------------------------------------

library(scagnostics)
library(gplots)
library(MASS)


# BHD -----------------------------------------------------------------

s = scagnostics(Boston)
s = as.matrix(s)


# Heatmap BHD -----------------------------------------------------------------

png(filename = "Heatmap BHD.png", bg = NA)
hm <- heatmap.2(s, Rowv=NA, Colv=NA,
                  dendrogram = "none",  ## to suppress warnings
                  margins=c(5,5), cexRow=0.5, cexCol=1.0, key=TRUE, keysize=1.5,
                  trace="none", scale(s1))

dev.off()
# PC Factors -----------------------------------------------------------------

s = scagnostics(factors_scag_gr)
s = as.matrix(s)


#  Heatmap of PC Factors --------------------------------------------------------

png(filename = "Heatmap PC Factors.png", bg = NA)
hm <- heatmap.2(s, Rowv=NA, Colv=NA,
                dendrogram = "none",  ## to suppress warnings
                margins=c(5,5), cexRow=0.5, cexCol=1.0, key=TRUE, keysize=1.5,
                trace="none", scale(s1))

dev.off()
