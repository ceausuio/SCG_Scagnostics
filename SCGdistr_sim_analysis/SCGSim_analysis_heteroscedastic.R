rm(list = ls())
options(java.parameters = "-Xmx8000m")

library(scagnostics)
library(alphahull)

# Functions  --------------------------------------------------------------

#function to simulate heteroscedastic distribution
scag_het = function(n){

  N = rep(1:n, 2)
  a = 0
  b = 1
  sigma2 = N^1.3
  eps = rnorm(N, mean = 0, sd = sqrt(sigma2))
  y = a + b*N + 3*eps
  mod <- lm(y ~ N)

  df_het = as.data.frame(cbind(N,y))
  s_het = scagnostics(df_het)

  het_eval = rbind(s_het)
}

#test = scag_het(100)

#function to add to the same list the each simulation iteration
scag_simul = function (sim_settings_n, sim_no_rep){
  #sim_settings_n = 5
  #sim_no_rep = 3
  scag_simul_output = scag_het(sim_settings_n)
  for (i in (1:(sim_no_rep-1))){
    #print(paste0("i",i, "scag_simul_output:",scag_simul_output ))
    added_sim_round = scag_het(sim_settings_n)
    scag_simul_output = rbind(scag_simul_output,added_sim_round)
  } # end of for cycle i
  rownames(scag_simul_output) = paste0("sim",c(1:sim_no_rep))
  return(scag_simul_output)

} #end of scag_simul function

#scag_simul(100, 2) - simno_rep > 2 always

#individual simulation
no_rep = 1000
n = 100
simul_final_output = scag_simul(n,no_rep)

write.csv(simul_final_output, file = paste0("DFScagSim_analysis 100 by 1000.csv"))


#"automation" of simulation for all n
no_rep = 1000
vec_ns = as.list(c(100,500,seq(1000,50000,by=1000)))

simul_final_output_list = lapply(vec_ns,scag_simul, sim_no_rep = no_rep)

write.csv(simul_final_output_list, file = paste0("DFScagSim_analysis uniform full table.csv"))

#Verification of outputs
dim(simul_final_output_list[[52]]) #last "page" in the list
length(simul_final_output_list)

dim(simul_final_output_list[[1]]) #first "page" in the list
length(simul_final_output_list[[1]])

simul_final_output_list[[1]][,1] #first column in the first "page" of the list
list_no1 = simul_final_output_list[[1]] #first "page" in the list
list_no1[,1] #first column in the first "page" of the list


# Creating matrix for each scagnostic measure -----------------------------

#empty matrix for a scagnostic matrix
output_scag_measure = matrix(nrow=no_rep, ncol=length(vec_ns))
dim(output_scag_measure)

column_index = 9
# Column legend: 1 = Outlying; 2 = Skewed; 3 = Clumpy; 4 = Sparse; 5 = Striate;
# 6 = Convex; 7 = Skinny; 8 = Stringy; 9 = Monotonic

for (page in 1:length(vec_ns)){
  working_page = simul_final_output_list[[page]]
  output_scag_measure[,page] =  working_page[,column_index]
}

output_scag_measure #matrix with all 52 columns with the assigned index "column_index"

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_het Skewed.csv"))


# 1 Analysis of scagnostic measures: Outlying----------------------------------------------

#Outlying measure

#preparing data set
osg_het_Outlying = t(output_scag_measure)
dim(osg_het_Outlying)
rownames(osg_het_Outlying) = vec_ns



#Calculating the mean, q25 and q75
#i = 1
mq_het_Outlying = function(x){

  mean_het_Outlying = mean(x)
  q_het_Outlying = quantile(x)
  mq_het_Outlying = as.data.frame(rbind(mean_het_Outlying, q_het_Outlying[2],
                                        q_het_Outlying[4]))
  rownames(mq_het_Outlying) = c("mean", "25%", "75%")
  return (mq_het_Outlying)
}

#test
mq_het_Outlying(osg_het_Outlying[1,])

mean_het_Outlying = mean(osg_het_Outlying[1,])
q_het_Outlying = quantile(osg_het_Outlying[1,])

mq_het_Outlying_output = apply(osg_het_Outlying, 1, mq_het_Outlying)

unlist(mq_het_Outlying_output)
#dim(unlist(mq_het_Outlying_output))

mq_het_Outlying_output_final = matrix(unlist(mq_het_Outlying_output),
                                      ncol = 3, by = TRUE)
#dim(mq_het_Outlying_output_final)

mq_het_Outlying_output_final = cbind(vec_ns,mq_het_Outlying_output_final)
colnames(mq_het_Outlying_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_het_Outlying_output_final) = vec_ns

mq_het_Outlying_output_final


# 1 Visualisation of mean, q25 and q75: Outlying --------------------------------------

#Outlying
dev.off()
par(bg = NA)
plot(mq_het_Outlying_output_final[,1], mq_het_Outlying_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "blue",
     main = "Coefficient: Outlying", sub = "Heteroscedastic distribution")
lines(mq_het_Outlying_output_final[,1], mq_het_Outlying_output_final[,3],
      col = "dark green")
lines(mq_het_Outlying_output_final[,1], mq_het_Outlying_output_final[,4],
      col = "dark red")
dev.off()




# 2 Analysis of scagnostic measures: Skewed ------------------------------------------------------------------

#preparing data set
osg_het_Skewed = t(output_scag_measure)
dim(osg_het_Skewed)
rownames(osg_het_Skewed) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_het_Skewed = function(x){

  mean_het_Skewed = mean(x)
  q_het_Skewed = quantile(x)
  mq_het_Skewed = as.data.frame(rbind(mean_het_Skewed, q_het_Skewed[2],
                                      q_het_Skewed[4]))
  rownames(mq_het_Skewed) = c("mean", "25%", "75%")
  return (mq_het_Skewed)
}

#test
mq_het_Skewed(osg_het_Skewed[1,])

mean_het_Skewed = mean(osg_het_Skewed[1,])
q_het_Skewed = quantile(osg_het_Skewed[1,])

mq_het_Skewed_output = apply(osg_het_Skewed, 1, mq_het_Skewed)

unlist(mq_het_Skewed_output)
#dim(unlist(mq_het_Outlying_output))

mq_het_Skewed_output_final = matrix(unlist(mq_het_Skewed_output),
                                    ncol = 3, by = TRUE)
dim(mq_het_Skewed_output_final)


mq_het_Skewed_output_final = cbind(vec_ns,mq_het_Skewed_output_final)
colnames(mq_het_Skewed_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_het_Skewed_output_final) = vec_ns

mq_het_Skewed_output_final

# 2 Visualisation of mean, q25 and q75: Skewed -----------------------------------

#Skewed

dev.off()
par(bg = NA)
plot(mq_het_Skewed_output_final[,1], mq_het_Skewed_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skewed", sub = "Heteroscedastic Distribution")
lines(mq_het_Skewed_output_final[,1], mq_het_Skewed_output_final[,3],
      col = "dark green")
lines(mq_het_Skewed_output_final[,1], mq_het_Skewed_output_final[,4],
      col = "dark red")
dev.off()




# 3 Analysis of scagnostic measures: Clumpy ------------------------------------------------------------------

#preparing data set
osg_het_Clumpy = t(output_scag_measure)
dim(osg_het_Clumpy)
rownames(osg_het_Clumpy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_het_Clumpy = function(x){

  mean_het_Clumpy = mean(x)
  q_het_Clumpy = quantile(x)
  mq_het_Clumpy = as.data.frame(rbind(mean_het_Clumpy, q_het_Clumpy[2],
                                      q_het_Clumpy[4]))
  rownames(mq_het_Clumpy) = c("mean", "25%", "75%")
  return (mq_het_Clumpy)
}

#test
mq_het_Clumpy(osg_het_Clumpy[1,])

mean_het_Clumpy = mean(osg_het_Clumpy[1,])
q_het_Clumpy = quantile(osg_het_Clumpy[1,])


mq_het_Clumpy_output = apply(osg_het_Clumpy, 1, mq_het_Clumpy)

unlist(mq_het_Clumpy_output)

mq_het_Clumpy_output_final = matrix(unlist(mq_het_Clumpy_output),
                                    ncol = 3, by = TRUE)
dim(mq_het_Clumpy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_het Skewed_mq.csv"))

mq_het_Clumpy_output_final = cbind(vec_ns,mq_het_Clumpy_output_final)
colnames(mq_het_Clumpy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_het_Clumpy_output_final) = vec_ns

mq_het_Clumpy_output_final

# 3 Visualisation of mean, q25 and q75: Clumpy -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_het_Clumpy_output_final[,1], mq_het_Clumpy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Clumpy", sub = "Heteroscedastic distribution")
lines(mq_het_Clumpy_output_final[,1], mq_het_Clumpy_output_final[,3],
      col = "dark green")
lines(mq_het_Clumpy_output_final[,1], mq_het_Clumpy_output_final[,4],
      col = "dark red")
dev.off()




# 4 Analysis of scagnostic measures: Sparse ------------------------------------------------------------------

#preparing data set
osg_het_Sparse = t(output_scag_measure)
dim(osg_het_Sparse)
rownames(osg_het_Sparse) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_het_Sparse = function(x){

  mean_het_Sparse = mean(x)
  q_het_Sparse = quantile(x)
  mq_het_Sparse = as.data.frame(rbind(mean_het_Sparse, q_het_Sparse[2],
                                      q_het_Sparse[4]))
  rownames(mq_het_Sparse) = c("mean", "25%", "75%")
  return (mq_het_Sparse)
}

#test
mq_het_Sparse(osg_het_Sparse[1,])

mean_het_Sparse = mean(osg_het_Sparse[1,])
q_het_Sparse = quantile(osg_het_Sparse[1,])


mq_het_Sparse_output = apply(osg_het_Sparse, 1, mq_het_Sparse)

unlist(mq_het_Sparse_output)

mq_het_Sparse_output_final = matrix(unlist(mq_het_Sparse_output),
                                    ncol = 3, by = TRUE)
dim(mq_het_Sparse_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_het Skewed_mq.csv"))

mq_het_Sparse_output_final = cbind(vec_ns,mq_het_Sparse_output_final)
colnames(mq_het_Sparse_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_het_Sparse_output_final) = vec_ns

mq_het_Sparse_output_final

# 4 Visualisation of mean, q25 and q75: Sparse -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_het_Sparse_output_final[,1], mq_het_Sparse_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "25q, mean, q75", col = "dark blue",
     main = "Coefficient: Sparse", sub = "Heteroscedastic distribution")
lines(mq_het_Sparse_output_final[,1], mq_het_Sparse_output_final[,3],
      col = "dark green")
lines(mq_het_Sparse_output_final[,1], mq_het_Sparse_output_final[,4],
      col = "dark red")
dev.off()


# 5 Analysis of scagnostic measures: Striate ------------------------------------------------------------------

#preparing data set
osg_het_Striate = t(output_scag_measure)
dim(osg_het_Striate)
rownames(osg_het_Striate) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_het_Striate = function(x){

  mean_het_Striate = mean(x)
  q_het_Striate = quantile(x)
  mq_het_Striate = as.data.frame(rbind(mean_het_Striate, q_het_Striate[2],
                                       q_het_Striate[4]))
  rownames(mq_het_Striate) = c("mean", "25%", "75%")
  return (mq_het_Striate)
}

#test
mq_het_Striate(osg_het_Striate[1,])

mean_het_Striate = mean(osg_het_Striate[1,])
q_het_Striate = quantile(osg_het_Striate[1,])


mq_het_Striate_output = apply(osg_het_Striate, 1, mq_het_Striate)

unlist(mq_het_Striate_output)

mq_het_Striate_output_final = matrix(unlist(mq_het_Striate_output),
                                     ncol = 3, by = TRUE)
dim(mq_het_Striate_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_het Skewed_mq.csv"))

mq_het_Striate_output_final = cbind(vec_ns,mq_het_Striate_output_final)
colnames(mq_het_Striate_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_het_Striate_output_final) = vec_ns

mq_het_Striate_output_final

# 5 Visualisation of mean, q25 and q75: Striate -----------------------------------

#Striate
dev.off()
par(bg = NA)
plot(mq_het_Striate_output_final[,1], mq_het_Striate_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Striated", sub = "Heteroscedastic distribution")
lines(mq_het_Striate_output_final[,1], mq_het_Striate_output_final[,3],
      col = "dark green")
lines(mq_het_Striate_output_final[,1], mq_het_Striate_output_final[,4],
      col = "dark red")
dev.off()



# 6 Analysis of scagnostic measures: Convex ------------------------------------------------------------------

#preparing data set
osg_het_Convex = t(output_scag_measure)
dim(osg_het_Convex)
rownames(osg_het_Convex) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_het_Convex = function(x){

  mean_het_Convex = mean(x)
  q_het_Convex = quantile(x)
  mq_het_Convex = as.data.frame(rbind(mean_het_Convex, q_het_Convex[2],
                                      q_het_Convex[4]))
  rownames(mq_het_Convex) = c("mean", "25%", "75%")
  return (mq_het_Convex)
}

#test
mq_het_Convex(osg_het_Convex[1,])

mean_het_Convex = mean(osg_het_Convex[1,])
q_het_Convex = quantile(osg_het_Convex[1,])


mq_het_Convex_output = apply(osg_het_Convex, 1, mq_het_Convex)

unlist(mq_het_Convex_output)

mq_het_Convex_output_final = matrix(unlist(mq_het_Convex_output),
                                    ncol = 3, by = TRUE)
dim(mq_het_Convex_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_het Skewed_mq.csv"))

mq_het_Convex_output_final = cbind(vec_ns, mq_het_Convex_output_final)
colnames(mq_het_Convex_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_het_Convex_output_final) = vec_ns

mq_het_Convex_output_final

# 6 Visualisation of mean, q25 and q75: Convex -----------------------------------

#Convex
dev.off()
par(bg = NA)
plot(mq_het_Convex_output_final[,1], mq_het_Convex_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Convex", sub = "Heteroscedastic distribution")
lines(mq_het_Convex_output_final[,1], mq_het_Convex_output_final[,3],
      col = "dark green")
lines(mq_het_Convex_output_final[,1], mq_het_Convex_output_final[,4],
      col = "dark red")
dev.off()


# 7 Analysis of scagnostic measures: Skinny ------------------------------------------------------------------

#preparing data set
osg_het_Skinny = t(output_scag_measure)
dim(osg_het_Skinny)
rownames(osg_het_Skinny) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_het_Skinny = function(x){

  mean_het_Skinny = mean(x)
  q_het_Skinny = quantile(x)
  mq_het_Skinny = as.data.frame(rbind(mean_het_Skinny, q_het_Skinny[2],
                                      q_het_Skinny[4]))
  rownames(mq_het_Skinny) = c("mean", "25%", "75%")
  return (mq_het_Skinny)
}

#test
mq_het_Skinny(osg_het_Skinny[1,])

mean_het_Skinny = mean(osg_het_Skinny[1,])
q_het_Skinny = quantile(osg_het_Skinny[1,])


mq_het_Skinny_output = apply(osg_het_Skinny, 1, mq_het_Skinny)

unlist(mq_het_Skinny_output)

mq_het_Skinny_output_final = matrix(unlist(mq_het_Skinny_output),
                                    ncol = 3, by = TRUE)
dim(mq_het_Skinny_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_het Skewed_mq.csv"))

mq_het_Skinny_output_final = cbind(vec_ns, mq_het_Skinny_output_final)
colnames(mq_het_Skinny_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_het_Skinny_output_final) = vec_ns

mq_het_Skinny_output_final

# 7 Visualisation of mean, q25 and q75: Skinny -----------------------------------

#Skinny
dev.off()
par(bg = NA)
plot(mq_het_Skinny_output_final[,1], mq_het_Skinny_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skinny", sub = "Heteroscedastic distribution")
lines(mq_het_Skinny_output_final[,1], mq_het_Skinny_output_final[,3],
      col = "dark green")
lines(mq_het_Skinny_output_final[,1], mq_het_Skinny_output_final[,4],
      col = "dark red")
dev.off()


# 8 Analysis of scagnostic measures: Stringy ------------------------------------------------------------------

#preparing data set
osg_het_Stringy = t(output_scag_measure)
dim(osg_het_Stringy)
rownames(osg_het_Stringy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_het_Stringy = function(x){

  mean_het_Stringy = mean(x)
  q_het_Stringy = quantile(x)
  mq_het_Stringy = as.data.frame(rbind(mean_het_Stringy, q_het_Stringy[2],
                                       q_het_Stringy[4]))
  rownames(mq_het_Stringy) = c("mean", "25%", "75%")
  return (mq_het_Stringy)
}

#test
mq_het_Stringy(osg_het_Stringy[1,])

mean_het_Stringy = mean(osg_het_Stringy[1,])
q_het_Stringy = quantile(osg_het_Stringy[1,])


mq_het_Stringy_output = apply(osg_het_Stringy, 1, mq_het_Stringy)

unlist(mq_het_Stringy_output)

mq_het_Stringy_output_final = matrix(unlist(mq_het_Stringy_output),
                                     ncol = 3, by = TRUE)
dim(mq_het_Stringy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_het Skewed_mq.csv"))

mq_het_Stringy_output_final = cbind(vec_ns, mq_het_Stringy_output_final)
colnames(mq_het_Stringy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_het_Stringy_output_final) = vec_ns

mq_het_Stringy_output_final

# 8 Visualisation of mean, q25 and q75: Stringy -----------------------------------

#Stringy
dev.off()
par(bg = NA)
plot(mq_het_Stringy_output_final[,1], mq_het_Stringy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Stringy", sub = "Heteroscedastic distribution")
lines(mq_het_Stringy_output_final[,1], mq_het_Stringy_output_final[,3],
      col = "dark green")
lines(mq_het_Stringy_output_final[,1], mq_het_Stringy_output_final[,4],
      col = "dark red")
dev.off()



# 9 Analysis of scagnostic measures: Monotonic ------------------------------------------------------------------

#preparing data set
osg_het_Monotonic = t(output_scag_measure)
dim(osg_het_Monotonic)
rownames(osg_het_Monotonic) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_het_Monotonic = function(x){

  mean_het_Monotonic = mean(x)
  q_het_Monotonic = quantile(x)
  mq_het_Monotonic = as.data.frame(rbind(mean_het_Monotonic, q_het_Monotonic[2],
                                         q_het_Monotonic[4]))
  rownames(mq_het_Monotonic) = c("mean", "25%", "75%")
  return (mq_het_Monotonic)
}

#test
mq_het_Monotonic(osg_het_Monotonic[1,])

mean_het_Monotonic = mean(osg_het_Monotonic[1,])
q_het_Monotonic = quantile(osg_het_Monotonic[1,])


mq_het_Monotonic_output = apply(osg_het_Monotonic, 1, mq_het_Monotonic)

unlist(mq_het_Monotonic_output)

mq_het_Monotonic_output_final = matrix(unlist(mq_het_Monotonic_output),
                                       ncol = 3, by = TRUE)
dim(mq_het_Monotonic_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_het Skewed_mq.csv"))

mq_het_Monotonic_output_final = cbind(vec_ns, mq_het_Monotonic_output_final)
colnames(mq_het_Monotonic_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_het_Monotonic_output_final) = vec_ns

mq_het_Monotonic_output_final

# 9 Visualisation of mean, q25 and q75: Monotonic --------------------------------

#Monotonic
dev.off()
par(bg = NA)
plot(mq_het_Monotonic_output_final[,1], mq_het_Monotonic_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Monotonic", sub = "Heteroscedastic distribution")
lines(mq_het_Monotonic_output_final[,1], mq_het_Monotonic_output_final[,3],
      col = "dark green")
lines(mq_het_Monotonic_output_final[,1], mq_het_Monotonic_output_final[,4],
      col = "dark red")
dev.off()
