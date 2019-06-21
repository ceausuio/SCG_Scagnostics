library(scagnostics)
library(psych)

library(gplots)
library(RColorBrewer)


library(readr)
p2p = read_csv("BUBA/Data P2P data_final_stand.csv")
View(p2p)

p2p = p2p[,2:26]

View(p2p)
p2p = as.data.frame(p2p)

dev.off()
png(
  filename = 'P2P SPLOM1.png',
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
par(
  bg = NA,
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2
)
pairs(p2p, pch = ".", col = "dark red", main = "P2P SPLOM")
dev.off()

s = scagnostics(p2p)
View(s)


s1 = t(s)

dev.off()
png(
  filename = 'P2P Scag SPLOM1.png',
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
par(
  bg = NA,
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2
)
pairs(s1, pch = ".", col = "dark red", main = "P2P Scagnostics SPLOM")
dev.off()


s1 = as.matrix(s1)

write.csv(s1, file = paste0("20190523_Scagnostics of P2P.csv"))

p2p_scag = read_csv("~/BUBA/20190523_Scagnostics of P2P.csv")
View(p2p_scag)
p2p_scag = as.data.frame(p2p_scag)
rownames(p2p_scag) =  p2p_scag$X1
p2p_scag = p2p_scag[,2:10]

p2p_scag_status = p2p_scag[grep("status", rownames(p2p_scag)), ]
p2p_scag_status  = as.matrix(p2p_scag_status)
p2p_scag_status = round(p2p_scag_status, 3)


dev.off()
png(
  filename = 'P2P Status Scag SPLOM2.png',
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
par(
  bg = NA,
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 4,
  cex.lab  = 2
)
pairs(p2p_scag_status, pch = ".", col = "dark red", main = "P2P Status Scagnostics SPLOM")
dev.off()

mypal = colorRampPalette(c("red", "white", "blue"))(n = 299)

dev.off()
graphics.off()
png(
  filename = 'P2P Status Scag Heatmap5.png',
  width     = 4,
  height    = 4,
  units     = "in",
  res       = 1200,
  pointsize = 3
)
par(
  bg = NA,
  mar      = c(5, 5, 5, 5),
  #xaxs     = "i",
  #yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 4
)
heatmap.2(p2p_scag_status,
          col = mypal,
          cellnote = p2p_scag_status,
          scale = "column",
          density.inf = "none",
          dendrogram = "none",
          #sepwidth = c(5, 5),
          notecol = "black",
          cexRow = 1,
          cexCol = 1,
          trace = "none",
          keysize = 0)
dev.off()

dev.off()
png(
  filename = 'ratio 012 vs status.png',
  width     = 4,
  height    = 4,
  units     = "in",
  res       = 1200,
  pointsize = 3
)
par(bg = NA)
plot(p2p$ratio012, p2p$status, pch = ".",
     col = "dark red",
     xlab = "ratio0012",
     ylab = "status",
     main = "ratio 0012 vs. status")
dev.off()

dev.off()
png(
  filename = 'turnover vs status.png',
  width     = 4,
  height    = 4,
  units     = "in",
  res       = 1200,
  pointsize = 3
)
par(bg = NA)
plot(p2p$turnover, p2p$status, pch = ".",
     col = "dark red",
     xlab = "turnover",
     ylab = "status",
     main = "turnover vs. status")
dev.off()

dev.off()
png(
  filename = 'DIO vs status.png',
  width     = 4,
  height    = 4,
  units     = "in",
  res       = 1200,
  pointsize = 3
)
par(bg = NA)
plot(p2p$DIO, p2p$status, pch = ".",
     col = "dark red",
     xlab = "DIO",
     ylab = "status",
     main = "DIO vs. status")
dev.off()
