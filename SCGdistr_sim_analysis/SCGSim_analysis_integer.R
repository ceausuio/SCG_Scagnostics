rm(list = ls())
options(java.parameters = "-Xmx10000m")

library(scagnostics)
library(alphahull)

# Functions  --------------------------------------------------------------

#function to simulate Integer distribution
scag_int = function(n){

  N = rep(1:n, n)
  a = 0
  b = 1
  sigma2 = N^1.3
  eps = rnorm(N, mean = 0, sd = sqrt(sigma2))
  y = a + b*N + 2*eps
  mod <- lm(y ~ N)

  df_int = as.data.frame(cbind(N, y))
  s_int = scagnostics(df_int)

  int_eval = rbind(s_int)
}

#test = scag_int(100)

#function to add to the same list the each simulation iteration
scag_simul = function (sim_settings_n, sim_no_rep){
  #sim_settings_n = 5
  #sim_no_rep = 3
  scag_simul_output = scag_int(sim_settings_n)
  for (i in (1:(sim_no_rep-1))){
    #print(paste0("i",i, "scag_simul_output:",scag_simul_output ))
    added_sim_round = scag_int(sim_settings_n)
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

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_int Skewed.csv"))


# 1 Analysis of scagnostic measures: Outlying----------------------------------------------

#Outlying measure

#preparing data set
osg_int_Outlying = t(output_scag_measure)
dim(osg_int_Outlying)
rownames(osg_int_Outlying) = vec_ns



#Calculating the mean, q25 and q75
#i = 1
mq_int_Outlying = function(x){

  mean_int_Outlying = mean(x)
  q_int_Outlying = quantile(x)
  mq_int_Outlying = as.data.frame(rbind(mean_int_Outlying, q_int_Outlying[2],
                                        q_int_Outlying[4]))
  rownames(mq_int_Outlying) = c("mean", "25%", "75%")
  return (mq_int_Outlying)
}

#test
mq_int_Outlying(osg_int_Outlying[1,])

mean_int_Outlying = mean(osg_int_Outlying[1,])
q_int_Outlying = quantile(osg_int_Outlying[1,])

mq_int_Outlying_output = apply(osg_int_Outlying, 1, mq_int_Outlying)

unlist(mq_int_Outlying_output)
#dim(unlist(mq_int_Outlying_output))

mq_int_Outlying_output_final = matrix(unlist(mq_int_Outlying_output),
                                      ncol = 3, by = TRUE)
#dim(mq_int_Outlying_output_final)

mq_int_Outlying_output_final = cbind(vec_ns,mq_int_Outlying_output_final)
colnames(mq_int_Outlying_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_int_Outlying_output_final) = vec_ns

mq_int_Outlying_output_final


# 1 Visualisation of mean, q25 and q75: Outlying --------------------------------------

#Outlying
dev.off()
par(bg = NA)
plot(mq_int_Outlying_output_final[,1], mq_int_Outlying_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "blue",
     main = "Coefficient: Outlying", sub = "Integer distribution")
lines(mq_int_Outlying_output_final[,1], mq_int_Outlying_output_final[,3],
      col = "dark green")
lines(mq_int_Outlying_output_final[,1], mq_int_Outlying_output_final[,4],
      col = "dark red")
dev.off()




# 2 Analysis of scagnostic measures: Skewed ------------------------------------------------------------------

#preparing data set
osg_int_Skewed = t(output_scag_measure)
dim(osg_int_Skewed)
rownames(osg_int_Skewed) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_int_Skewed = function(x){

  mean_int_Skewed = mean(x)
  q_int_Skewed = quantile(x)
  mq_int_Skewed = as.data.frame(rbind(mean_int_Skewed, q_int_Skewed[2],
                                      q_int_Skewed[4]))
  rownames(mq_int_Skewed) = c("mean", "25%", "75%")
  return (mq_int_Skewed)
}

#test
mq_int_Skewed(osg_int_Skewed[1,])

mean_int_Skewed = mean(osg_int_Skewed[1,])
q_int_Skewed = quantile(osg_int_Skewed[1,])

mq_int_Skewed_output = apply(osg_int_Skewed, 1, mq_int_Skewed)

unlist(mq_int_Skewed_output)
#dim(unlist(mq_int_Outlying_output))

mq_int_Skewed_output_final = matrix(unlist(mq_int_Skewed_output),
                                    ncol = 3, by = TRUE)
dim(mq_int_Skewed_output_final)


mq_int_Skewed_output_final = cbind(vec_ns,mq_int_Skewed_output_final)
colnames(mq_int_Skewed_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_int_Skewed_output_final) = vec_ns

mq_int_Skewed_output_final

# 2 Visualisation of mean, q25 and q75: Skewed -----------------------------------

#Skewed

dev.off()
par(bg = NA)
plot(mq_int_Skewed_output_final[,1], mq_int_Skewed_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skewed", sub = "Integer Distribution")
lines(mq_int_Skewed_output_final[,1], mq_int_Skewed_output_final[,3],
      col = "dark green")
lines(mq_int_Skewed_output_final[,1], mq_int_Skewed_output_final[,4],
      col = "dark red")
dev.off()




# 3 Analysis of scagnostic measures: Clumpy ------------------------------------------------------------------

#preparing data set
osg_int_clumpy = t(output_scag_measure)
dim(osg_int_clumpy)
rownames(osg_int_clumpy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_int_clumpy = function(x){

  mean_int_clumpy = mean(x)
  q_int_clumpy = quantile(x)
  mq_int_clumpy = as.data.frame(rbind(mean_int_clumpy, q_int_clumpy[2],
                                      q_int_clumpy[4]))
  rownames(mq_int_clumpy) = c("mean", "25%", "75%")
  return (mq_int_clumpy)
}

#test
mq_int_clumpy(osg_int_clumpy[1,])

mean_int_clumpy = mean(osg_int_clumpy[1,])
q_int_clumpy = quantile(osg_int_clumpy[1,])


mq_int_clumpy_output = apply(osg_int_clumpy, 1, mq_int_clumpy)

unlist(mq_int_clumpy_output)

mq_int_clumpy_output_final = matrix(unlist(mq_int_clumpy_output),
                                    ncol = 3, by = TRUE)
dim(mq_int_clumpy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_int Skewed_mq.csv"))

mq_int_clumpy_output_final = cbind(vec_ns,mq_int_clumpy_output_final)
colnames(mq_int_clumpy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_int_clumpy_output_final) = vec_ns

mq_int_clumpy_output_final

# 3 Visualisation of mean, q25 and q75: clumpy -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_int_clumpy_output_final[,1], mq_int_clumpy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Clumpy", sub = "Integer distribution")
lines(mq_int_clumpy_output_final[,1], mq_int_clumpy_output_final[,3],
      col = "dark green")
lines(mq_int_clumpy_output_final[,1], mq_int_clumpy_output_final[,4],
      col = "dark red")
dev.off()




# 4 Analysis of scagnostic measures: Sparse ------------------------------------------------------------------

#preparing data set
osg_int_Sparse = t(output_scag_measure)
dim(osg_int_Sparse)
rownames(osg_int_Sparse) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_int_Sparse = function(x){

  mean_int_Sparse = mean(x)
  q_int_Sparse = quantile(x)
  mq_int_Sparse = as.data.frame(rbind(mean_int_Sparse, q_int_Sparse[2],
                                      q_int_Sparse[4]))
  rownames(mq_int_Sparse) = c("mean", "25%", "75%")
  return (mq_int_Sparse)
}

#test
mq_int_Sparse(osg_int_Sparse[1,])

mean_int_Sparse = mean(osg_int_Sparse[1,])
q_int_Sparse = quantile(osg_int_Sparse[1,])


mq_int_Sparse_output = apply(osg_int_Sparse, 1, mq_int_Sparse)

unlist(mq_int_Sparse_output)

mq_int_Sparse_output_final = matrix(unlist(mq_int_Sparse_output),
                                    ncol = 3, by = TRUE)
dim(mq_int_Sparse_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_int Skewed_mq.csv"))

mq_int_Sparse_output_final = cbind(vec_ns,mq_int_Sparse_output_final)
colnames(mq_int_Sparse_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_int_Sparse_output_final) = vec_ns

mq_int_Sparse_output_final

# 4 Visualisation of mean, q25 and q75: Sparse -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_int_Sparse_output_final[,1], mq_int_Sparse_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "25q, mean, q75", col = "dark blue",
     main = "Coefficient: Sparse", sub = "Integer distribution")
lines(mq_int_Sparse_output_final[,1], mq_int_Sparse_output_final[,3],
      col = "dark green")
lines(mq_int_Sparse_output_final[,1], mq_int_Sparse_output_final[,4],
      col = "dark red")
dev.off()


# 5 Analysis of scagnostic measures: Striate ------------------------------------------------------------------

#preparing data set
osg_int_Striate = t(output_scag_measure)
dim(osg_int_Striate)
rownames(osg_int_Striate) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_int_Striate = function(x){

  mean_int_Striate = mean(x)
  q_int_Striate = quantile(x)
  mq_int_Striate = as.data.frame(rbind(mean_int_Striate, q_int_Striate[2],
                                       q_int_Striate[4]))
  rownames(mq_int_Striate) = c("mean", "25%", "75%")
  return (mq_int_Striate)
}

#test
mq_int_Striate(osg_int_Striate[1,])

mean_int_Striate = mean(osg_int_Striate[1,])
q_int_Striate = quantile(osg_int_Striate[1,])


mq_int_Striate_output = apply(osg_int_Striate, 1, mq_int_Striate)

unlist(mq_int_Striate_output)

mq_int_Striate_output_final = matrix(unlist(mq_int_Striate_output),
                                     ncol = 3, by = TRUE)
dim(mq_int_Striate_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_int Skewed_mq.csv"))

mq_int_Striate_output_final = cbind(vec_ns,mq_int_Striate_output_final)
colnames(mq_int_Striate_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_int_Striate_output_final) = vec_ns

mq_int_Striate_output_final

# 5 Visualisation of mean, q25 and q75: Striate -----------------------------------

#Striate
dev.off()
par(bg = NA)
plot(mq_int_Striate_output_final[,1], mq_int_Striate_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Striated", sub = "Integer distribution")
lines(mq_int_Striate_output_final[,1], mq_int_Striate_output_final[,3],
      col = "dark green")
lines(mq_int_Striate_output_final[,1], mq_int_Striate_output_final[,4],
      col = "dark red")
dev.off()



# 6 Analysis of scagnostic measures: Convex ------------------------------------------------------------------

#preparing data set
osg_int_Convex = t(output_scag_measure)
dim(osg_int_Convex)
rownames(osg_int_Convex) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_int_Convex = function(x){

  mean_int_Convex = mean(x)
  q_int_Convex = quantile(x)
  mq_int_Convex = as.data.frame(rbind(mean_int_Convex, q_int_Convex[2],
                                      q_int_Convex[4]))
  rownames(mq_int_Convex) = c("mean", "25%", "75%")
  return (mq_int_Convex)
}

#test
mq_int_Convex(osg_int_Convex[1,])

mean_int_Convex = mean(osg_int_Convex[1,])
q_int_Convex = quantile(osg_int_Convex[1,])


mq_int_Convex_output = apply(osg_int_Convex, 1, mq_int_Convex)

unlist(mq_int_Convex_output)

mq_int_Convex_output_final = matrix(unlist(mq_int_Convex_output),
                                    ncol = 3, by = TRUE)
dim(mq_int_Convex_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_int Skewed_mq.csv"))

mq_int_Convex_output_final = cbind(vec_ns, mq_int_Convex_output_final)
colnames(mq_int_Convex_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_int_Convex_output_final) = vec_ns

mq_int_Convex_output_final

# 6 Visualisation of mean, q25 and q75: Convex -----------------------------------

#Convex
dev.off()
par(bg = NA)
plot(mq_int_Convex_output_final[,1], mq_int_Convex_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Convex", sub = "Integer distribution")
lines(mq_int_Convex_output_final[,1], mq_int_Convex_output_final[,3],
      col = "dark green")
lines(mq_int_Convex_output_final[,1], mq_int_Convex_output_final[,4],
      col = "dark red")
dev.off()


# 7 Analysis of scagnostic measures: Skinny ------------------------------------------------------------------

#preparing data set
osg_int_Skinny = t(output_scag_measure)
dim(osg_int_Skinny)
rownames(osg_int_Skinny) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_int_Skinny = function(x){

  mean_int_Skinny = mean(x)
  q_int_Skinny = quantile(x)
  mq_int_Skinny = as.data.frame(rbind(mean_int_Skinny, q_int_Skinny[2],
                                      q_int_Skinny[4]))
  rownames(mq_int_Skinny) = c("mean", "25%", "75%")
  return (mq_int_Skinny)
}

#test
mq_int_Skinny(osg_int_Skinny[1,])

mean_int_Skinny = mean(osg_int_Skinny[1,])
q_int_Skinny = quantile(osg_int_Skinny[1,])


mq_int_Skinny_output = apply(osg_int_Skinny, 1, mq_int_Skinny)

unlist(mq_int_Skinny_output)

mq_int_Skinny_output_final = matrix(unlist(mq_int_Skinny_output),
                                    ncol = 3, by = TRUE)
dim(mq_int_Skinny_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_int Skewed_mq.csv"))

mq_int_Skinny_output_final = cbind(vec_ns, mq_int_Skinny_output_final)
colnames(mq_int_Skinny_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_int_Skinny_output_final) = vec_ns

mq_int_Skinny_output_final

# 7 Visualisation of mean, q25 and q75: Skinny -----------------------------------

#Skinny
dev.off()
par(bg = NA)
plot(mq_int_Skinny_output_final[,1], mq_int_Skinny_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skinny", sub = "Integer distribution")
lines(mq_int_Skinny_output_final[,1], mq_int_Skinny_output_final[,3],
      col = "dark green")
lines(mq_int_Skinny_output_final[,1], mq_int_Skinny_output_final[,4],
      col = "dark red")
dev.off()


# 8 Analysis of scagnostic measures: Stringy ------------------------------------------------------------------

#preparing data set
osg_int_Stringy = t(output_scag_measure)
dim(osg_int_Stringy)
rownames(osg_int_Stringy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_int_Stringy = function(x){

  mean_int_Stringy = mean(x)
  q_int_Stringy = quantile(x)
  mq_int_Stringy = as.data.frame(rbind(mean_int_Stringy, q_int_Stringy[2],
                                       q_int_Stringy[4]))
  rownames(mq_int_Stringy) = c("mean", "25%", "75%")
  return (mq_int_Stringy)
}

#test
mq_int_Stringy(osg_int_Stringy[1,])

mean_int_Stringy = mean(osg_int_Stringy[1,])
q_int_Stringy = quantile(osg_int_Stringy[1,])


mq_int_Stringy_output = apply(osg_int_Stringy, 1, mq_int_Stringy)

unlist(mq_int_Stringy_output)

mq_int_Stringy_output_final = matrix(unlist(mq_int_Stringy_output),
                                     ncol = 3, by = TRUE)
dim(mq_int_Stringy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_int Skewed_mq.csv"))

mq_int_Stringy_output_final = cbind(vec_ns, mq_int_Stringy_output_final)
colnames(mq_int_Stringy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_int_Stringy_output_final) = vec_ns

mq_int_Stringy_output_final

# 8 Visualisation of mean, q25 and q75: Stringy -----------------------------------

#Stringy
dev.off()
par(bg = NA)
plot(mq_int_Stringy_output_final[,1], mq_int_Stringy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Stringy", sub = "Integer distribution")
lines(mq_int_Stringy_output_final[,1], mq_int_Stringy_output_final[,3],
      col = "dark green")
lines(mq_int_Stringy_output_final[,1], mq_int_Stringy_output_final[,4],
      col = "dark red")
dev.off()



# 9 Analysis of scagnostic measures: Monotonic ------------------------------------------------------------------

#preparing data set
osg_int_Monotonic = t(output_scag_measure)
dim(osg_int_Monotonic)
rownames(osg_int_Monotonic) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_int_Monotonic = function(x){

  mean_int_Monotonic = mean(x)
  q_int_Monotonic = quantile(x)
  mq_int_Monotonic = as.data.frame(rbind(mean_int_Monotonic, q_int_Monotonic[2],
                                         q_int_Monotonic[4]))
  rownames(mq_int_Monotonic) = c("mean", "25%", "75%")
  return (mq_int_Monotonic)
}

#test
mq_int_Monotonic(osg_int_Monotonic[1,])

mean_int_Monotonic = mean(osg_int_Monotonic[1,])
q_int_Monotonic = quantile(osg_int_Monotonic[1,])


mq_int_Monotonic_output = apply(osg_int_Monotonic, 1, mq_int_Monotonic)

unlist(mq_int_Monotonic_output)

mq_int_Monotonic_output_final = matrix(unlist(mq_int_Monotonic_output),
                                       ncol = 3, by = TRUE)
dim(mq_int_Monotonic_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_int Skewed_mq.csv"))

mq_int_Monotonic_output_final = cbind(vec_ns, mq_int_Monotonic_output_final)
colnames(mq_int_Monotonic_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_int_Monotonic_output_final) = vec_ns

mq_int_Monotonic_output_final

# 9 Visualisation of mean, q25 and q75: Monotonic --------------------------------

#Monotonic
dev.off()
par(bg = NA)
plot(mq_int_Monotonic_output_final[,1], mq_int_Monotonic_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Monotonic", sub = "Integer distribution")
lines(mq_int_Monotonic_output_final[,1], mq_int_Monotonic_output_final[,3],
      col = "dark green")
lines(mq_int_Monotonic_output_final[,1], mq_int_Monotonic_output_final[,4],
      col = "dark red")
dev.off()
