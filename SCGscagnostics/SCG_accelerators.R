library(scagnostics)
library(psych)

library(gplots)
library(RColorBrewer)
library(MASS)

library(readxl)


EDP_10years = read_excel("Library/Mobile Documents/com~apple~CloudDocs/Humbold/Scagnostics/Data/EDP/20190124_EDP_ApplicationDataRelease-2017_scagnostics.xlsx")
View(EDP_10years)


df_acc = as.data.frame(EDP_10years)
View(df_acc)
colnames(df_acc) = c("age", "rev_sf", "ft_em_m1", "pt_em_m1", "wg_m1", "om_sf",
                     "oe_sf","td_sf", "plan_3y", "fnd1_age", "fnd1_gen")

dev.off()
par(bg = NA)
pairs(df_acc, pch = ".", col = "dark red", main = "SPLOM of accelerated startups")


plot(df_acc$inv_plans_outequity_3years, df_acc$fins_revenues_sincefound,
     pch = 16, col = "dark red", xlab = "Age", ylab = "Investment plan - 3 years" )


s = scagnostics(df_acc)

s1 = t(s)

pairs(s1, pch = ".", col = "dark red", main = "Scagnostics SPLOM of accelerated startups")


heatmap.2(s, scale="column",
          main="Heatmap of Scagnostics of Accelerated Startups",
          density="density",
          dendrogram = "none",
          notecol="black",
          margins=c(10,5), cexRow=1, cexCol=1, trace = "none")


s1 = as.matrix(s1)
View(s1)

write.csv(s1, file = paste0("Scagnostics of startups.csv"))

df_acc_3y = read_excel("Library/Mobile Documents/com~apple~CloudDocs/Humbold/Scagnostics/Data/EDP/20190124_EDP_ApplicationDataRelease-2017_scag plan 3y.xlsx")
View(df_acc_3y)
head(df_acc_3y)

df_acc_3y = df_acc_3y[,2:10]
rownames(df_acc_3y) = c("plan_3y*age",
                        "plan_3y*rev_sf",
                        "plan_3y*ft_em_m1",
                        "plan_3y*pt_em_m1",
                        "plan_3y*wg_m1",
                        "plan_3y*om_sf",
                        "plan_3y*oe_sf",
                        "plan_3y*td_sf",
                        "plan_3y*fnd1_age",
                        "plan_3y*fnd1_gen")


df_acc_3y = as.matrix(df_acc_3y)

View(df_acc_3y)

dev.off()
par(bg = NA)
heatmap.2(df_acc_3y,
          cellnote = df_acc_3y,
          scale="column",
          density="density",
          dendrogram = "none",
          notecol="black",
          margins=c(10,10), cexRow=1, cexCol=1, trace = "none")

