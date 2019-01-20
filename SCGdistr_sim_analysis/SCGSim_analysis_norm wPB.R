rm(list = ls())
options(java.parameters = "-Xmx8000m")

library(scagnostics)
library(alphahull)

# Functions  --------------------------------------------------------------

#function to simulate normal distribution
scag_norm = function(n){

  #n = 100
  x = rnorm(n)
  y = rnorm(n)

  df_nor = cbind(x,y)
  s_nor = scagnostics(df_nor)

  norm_eval = rbind(s_nor)
}

#test = scag_norm(100)

#function to add to the same list the each simulation iteration
scag_simul = function (sim_settings_n, sim_no_rep){
  #sim_settings_n = 5
  #sim_no_rep = 3
  scag_simul_output = scag_norm(sim_settings_n)
  for (i in (1:(sim_no_rep-1))){
    #print(paste0("i",i, "scag_simul_output:",scag_simul_output ))
    added_sim_round = scag_norm(sim_settings_n)
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

write.csv(simul_final_output_list, file = paste0("DFScagSim_analysis full table.csv"))

#Verification of out put
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

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_norm Skewed.csv"))


# 1 Analysis of scagnostic measures: Outlying----------------------------------------------

#Outlying measure

  #preparing data set
osg_norm_Outlying = t(output_scag_measure)
dim(osg_norm_Outlying)
rownames(osg_norm_Outlying) = vec_ns



  #Calculating the mean, q25 and q75
#i = 1
mq_norm_Outlying = function(x){

  mean_norm_Outlying = mean(x)
  q_norm_Outlying = quantile(x)
  mq_norm_Outlying = as.data.frame(rbind(mean_norm_Outlying, q_norm_Outlying[2],
                                         q_norm_Outlying[4]))
  rownames(mq_norm_Outlying) = c("mean", "25%", "75%")
  return (mq_norm_Outlying)
}

#test
  mq_norm_Outlying(osg_norm_Outlying[1,])

  mean_norm_Outlying = mean(osg_norm_Outlying[1,])
  q_norm_Outlying = quantile(osg_norm_Outlying[1,])

mq_norm_Outlying_output = apply(osg_norm_Outlying, 1, mq_norm_Outlying)

unlist(mq_norm_Outlying_output)
#dim(unlist(mq_norm_Outlying_output))

mq_norm_Outlying_output_final = matrix(unlist(mq_norm_Outlying_output),
                                       ncol = 3, by = TRUE)
#dim(mq_norm_Outlying_output_final)

mq_norm_Outlying_output_final = cbind(vec_ns,mq_norm_Outlying_output_final)
colnames(mq_norm_Outlying_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_norm_Outlying_output_final) = vec_ns

mq_norm_Outlying_output_final


# 1 Visualisation of mean, q25 and q75: Outlying --------------------------------------

#Outlying
dev.off()
par(bg = NA)
plot(mq_norm_Outlying_output_final[,1], mq_norm_Outlying_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "n - no. data points", ylab = "q25, mean, q75", col = "blue",
     main = "Coefficent: Outlying", sub = "Normal distribution")
lines(mq_norm_Outlying_output_final[,1], mq_norm_Outlying_output_final[,3],
      col = "dark green")
lines(mq_norm_Outlying_output_final[,1], mq_norm_Outlying_output_final[,4],
      col = "dark red")
dev.off()




# 2 Analysis of scagnostic measures: Skewed ------------------------------------------------------------------

#preparing data set
osg_norm_Skewed = t(output_scag_measure)
dim(osg_norm_Skewed)
rownames(osg_norm_Skewed) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_norm_Skewed = function(x){

  mean_norm_Skewed = mean(x)
  q_norm_Skewed = quantile(x)
  mq_norm_Skewed = as.data.frame(rbind(mean_norm_Skewed, q_norm_Skewed[2],
                                         q_norm_Skewed[4]))
  rownames(mq_norm_Skewed) = c("mean", "25%", "75%")
  return (mq_norm_Skewed)
}

#test
mq_norm_Skewed(osg_norm_Skewed[1,])

mean_norm_Skewed = mean(osg_norm_Skewed[1,])
q_norm_Skewed = quantile(osg_norm_Skewed[1,])

mq_norm_Skewed_output = apply(osg_norm_Skewed, 1, mq_norm_Skewed)

unlist(mq_norm_Skewed_output)
#dim(unlist(mq_norm_Outlying_output))

mq_norm_Skewed_output_final = matrix(unlist(mq_norm_Skewed_output),
                                       ncol = 3, by = TRUE)
dim(mq_norm_Skewed_output_final)


mq_norm_Skewed_output_final = cbind(vec_ns,mq_norm_Skewed_output_final)
colnames(mq_norm_Skewed_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_norm_Skewed_output_final) = vec_ns

mq_norm_Skewed_output_final

# 2 Visualisation of mean, q25 and q75: Skewed -----------------------------------

#Skewed

dev.off()
par(bg = NA)
plot(mq_norm_Skewed_output_final[,1], mq_norm_Skewed_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "n - no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skewed", sub = "Normal distribution")
lines(mq_norm_Skewed_output_final[,1], mq_norm_Skewed_output_final[,3],
      col = "dark green")
lines(mq_norm_Skewed_output_final[,1], mq_norm_Skewed_output_final[,4],
      col = "dark red")
dev.off()




# 3 Analysis of scagnostic measures: Clumpy ------------------------------------------------------------------

#preparing data set
osg_norm_Clumpy = t(output_scag_measure)
dim(osg_norm_Clumpy)
rownames(osg_norm_Clumpy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_norm_Clumpy = function(x){

  mean_norm_Clumpy = mean(x)
  q_norm_Clumpy = quantile(x)
  mq_norm_Clumpy = as.data.frame(rbind(mean_norm_Clumpy, q_norm_Clumpy[2],
                                       q_norm_Clumpy[4]))
  rownames(mq_norm_Clumpy) = c("mean", "25%", "75%")
  return (mq_norm_Clumpy)
}

#test
mq_norm_Clumpy(osg_norm_Clumpy[1,])

mean_norm_Clumpy = mean(osg_norm_Clumpy[1,])
q_norm_Clumpy = quantile(osg_norm_Clumpy[1,])


mq_norm_Clumpy_output = apply(osg_norm_Clumpy, 1, mq_norm_Clumpy)

unlist(mq_norm_Clumpy_output)

mq_norm_Clumpy_output_final = matrix(unlist(mq_norm_Clumpy_output),
                                     ncol = 3, by = TRUE)
dim(mq_norm_Clumpy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_norm Skewed_mq.csv"))

mq_norm_Clumpy_output_final = cbind(vec_ns,mq_norm_Clumpy_output_final)
colnames(mq_norm_Clumpy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_norm_Clumpy_output_final) = vec_ns

mq_norm_Clumpy_output_final

# 3 Visualisation of mean, q25 and q75: Clumpy -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_norm_Clumpy_output_final[,1], mq_norm_Clumpy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "n - no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Clumpy", sub = "Normal distribution")
lines(mq_norm_Clumpy_output_final[,1], mq_norm_Clumpy_output_final[,3],
      col = "dark green")
lines(mq_norm_Clumpy_output_final[,1], mq_norm_Clumpy_output_final[,4],
      col = "dark red")
dev.off()




# 4 Analysis of scagnostic measures: Sparse ------------------------------------------------------------------

#preparing data set
osg_norm_Sparse = t(output_scag_measure)
dim(osg_norm_Sparse)
rownames(osg_norm_Sparse) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_norm_Sparse = function(x){

  mean_norm_Sparse = mean(x)
  q_norm_Sparse = quantile(x)
  mq_norm_Sparse = as.data.frame(rbind(mean_norm_Sparse, q_norm_Sparse[2],
                                       q_norm_Sparse[4]))
  rownames(mq_norm_Sparse) = c("mean", "25%", "75%")
  return (mq_norm_Sparse)
}

#test
mq_norm_Sparse(osg_norm_Sparse[1,])

mean_norm_Sparse = mean(osg_norm_Sparse[1,])
q_norm_Sparse = quantile(osg_norm_Sparse[1,])


mq_norm_Sparse_output = apply(osg_norm_Sparse, 1, mq_norm_Sparse)

unlist(mq_norm_Sparse_output)

mq_norm_Sparse_output_final = matrix(unlist(mq_norm_Sparse_output),
                                     ncol = 3, by = TRUE)
dim(mq_norm_Sparse_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_norm Skewed_mq.csv"))

mq_norm_Sparse_output_final = cbind(vec_ns,mq_norm_Sparse_output_final)
colnames(mq_norm_Sparse_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_norm_Sparse_output_final) = vec_ns

mq_norm_Sparse_output_final

# 4 Visualisation of mean, q25 and q75: Sparse -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_norm_Sparse_output_final[,1], mq_norm_Sparse_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "n - no. data points", ylab = "25q, mean, q75", col = "dark blue",
     main = "Coefficient: Sparse", sub = "Normal distribution")
lines(mq_norm_Sparse_output_final[,1], mq_norm_Sparse_output_final[,3],
      col = "dark green")
lines(mq_norm_Sparse_output_final[,1], mq_norm_Sparse_output_final[,4],
      col = "dark red")
dev.off()


# 5 Analysis of scagnostic measures: Striate ------------------------------------------------------------------

#preparing data set
osg_norm_Striate = t(output_scag_measure)
dim(osg_norm_Striate)
rownames(osg_norm_Striate) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_norm_Striate = function(x){

  mean_norm_Striate = mean(x)
  q_norm_Striate = quantile(x)
  mq_norm_Striate = as.data.frame(rbind(mean_norm_Striate, q_norm_Striate[2],
                                       q_norm_Striate[4]))
  rownames(mq_norm_Striate) = c("mean", "25%", "75%")
  return (mq_norm_Striate)
}

#test
mq_norm_Striate(osg_norm_Striate[1,])

mean_norm_Striate = mean(osg_norm_Striate[1,])
q_norm_Striate = quantile(osg_norm_Striate[1,])


mq_norm_Striate_output = apply(osg_norm_Striate, 1, mq_norm_Striate)

unlist(mq_norm_Striate_output)

mq_norm_Striate_output_final = matrix(unlist(mq_norm_Striate_output),
                                     ncol = 3, by = TRUE)
dim(mq_norm_Striate_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_norm Skewed_mq.csv"))

mq_norm_Striate_output_final = cbind(vec_ns,mq_norm_Striate_output_final)
colnames(mq_norm_Striate_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_norm_Striate_output_final) = vec_ns

mq_norm_Striate_output_final

# 5 Visualisation of mean, q25 and q75: Striate -----------------------------------

#Striate
dev.off()
par(bg = NA)
plot(mq_norm_Striate_output_final[,1], mq_norm_Striate_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "n - no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Striated", sub = "Normal distribution")
lines(mq_norm_Striate_output_final[,1], mq_norm_Striate_output_final[,3],
      col = "dark green")
lines(mq_norm_Striate_output_final[,1], mq_norm_Striate_output_final[,4],
      col = "dark red")
dev.off()



# 6 Analysis of scagnostic measures: Convex ------------------------------------------------------------------

#preparing data set
osg_norm_Convex = t(output_scag_measure)
dim(osg_norm_Convex)
rownames(osg_norm_Convex) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_norm_Convex = function(x){

  mean_norm_Convex = mean(x)
  q_norm_Convex = quantile(x)
  mq_norm_Convex = as.data.frame(rbind(mean_norm_Convex, q_norm_Convex[2],
                                        q_norm_Convex[4]))
  rownames(mq_norm_Convex) = c("mean", "25%", "75%")
  return (mq_norm_Convex)
}

#test
mq_norm_Convex(osg_norm_Convex[1,])

mean_norm_Convex = mean(osg_norm_Convex[1,])
q_norm_Convex = quantile(osg_norm_Convex[1,])


mq_norm_Convex_output = apply(osg_norm_Convex, 1, mq_norm_Convex)

unlist(mq_norm_Convex_output)

mq_norm_Convex_output_final = matrix(unlist(mq_norm_Convex_output),
                                      ncol = 3, by = TRUE)
dim(mq_norm_Convex_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_norm Skewed_mq.csv"))

mq_norm_Convex_output_final = cbind(vec_ns, mq_norm_Convex_output_final)
colnames(mq_norm_Convex_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_norm_Convex_output_final) = vec_ns

mq_norm_Convex_output_final

# 6 Visualisation of mean, q25 and q75: Convex -----------------------------------

#Convex
dev.off()
par(bg = NA)
plot(mq_norm_Convex_output_final[,1], mq_norm_Convex_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "n - no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Convex", sub = "Normal distribution")
lines(mq_norm_Convex_output_final[,1], mq_norm_Convex_output_final[,3],
      col = "dark green")
lines(mq_norm_Convex_output_final[,1], mq_norm_Convex_output_final[,4],
      col = "dark red")
dev.off()


# 7 Analysis of scagnostic measures: Skinny ------------------------------------------------------------------

#preparing data set
osg_norm_Skinny = t(output_scag_measure)
dim(osg_norm_Skinny)
rownames(osg_norm_Skinny) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_norm_Skinny = function(x){

  mean_norm_Skinny = mean(x)
  q_norm_Skinny = quantile(x)
  mq_norm_Skinny = as.data.frame(rbind(mean_norm_Skinny, q_norm_Skinny[2],
                                       q_norm_Skinny[4]))
  rownames(mq_norm_Skinny) = c("mean", "25%", "75%")
  return (mq_norm_Skinny)
}

#test
mq_norm_Skinny(osg_norm_Skinny[1,])

mean_norm_Skinny = mean(osg_norm_Skinny[1,])
q_norm_Skinny = quantile(osg_norm_Skinny[1,])


mq_norm_Skinny_output = apply(osg_norm_Skinny, 1, mq_norm_Skinny)

unlist(mq_norm_Skinny_output)

mq_norm_Skinny_output_final = matrix(unlist(mq_norm_Skinny_output),
                                     ncol = 3, by = TRUE)
dim(mq_norm_Skinny_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_norm Skewed_mq.csv"))

mq_norm_Skinny_output_final = cbind(vec_ns, mq_norm_Skinny_output_final)
colnames(mq_norm_Skinny_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_norm_Skinny_output_final) = vec_ns

mq_norm_Skinny_output_final

# 7 Visualisation of mean, q25 and q75: Skinny -----------------------------------

#Skinny
dev.off()
par(bg = NA)
plot(mq_norm_Skinny_output_final[,1], mq_norm_Skinny_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "n - no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skinny", sub = "Normal distribution")
lines(mq_norm_Skinny_output_final[,1], mq_norm_Skinny_output_final[,3],
      col = "dark green")
lines(mq_norm_Skinny_output_final[,1], mq_norm_Skinny_output_final[,4],
      col = "dark red")
dev.off()


# 8 Analysis of scagnostic measures: Stringy ------------------------------------------------------------------

#preparing data set
osg_norm_Stringy = t(output_scag_measure)
dim(osg_norm_Stringy)
rownames(osg_norm_Stringy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_norm_Stringy = function(x){

  mean_norm_Stringy = mean(x)
  q_norm_Stringy = quantile(x)
  mq_norm_Stringy = as.data.frame(rbind(mean_norm_Stringy, q_norm_Stringy[2],
                                       q_norm_Stringy[4]))
  rownames(mq_norm_Stringy) = c("mean", "25%", "75%")
  return (mq_norm_Stringy)
}

#test
mq_norm_Stringy(osg_norm_Stringy[1,])

mean_norm_Stringy = mean(osg_norm_Stringy[1,])
q_norm_Stringy = quantile(osg_norm_Stringy[1,])


mq_norm_Stringy_output = apply(osg_norm_Stringy, 1, mq_norm_Stringy)

unlist(mq_norm_Stringy_output)

mq_norm_Stringy_output_final = matrix(unlist(mq_norm_Stringy_output),
                                     ncol = 3, by = TRUE)
dim(mq_norm_Stringy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_norm Skewed_mq.csv"))

mq_norm_Stringy_output_final = cbind(vec_ns, mq_norm_Stringy_output_final)
colnames(mq_norm_Stringy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_norm_Stringy_output_final) = vec_ns

mq_norm_Stringy_output_final

# 8 Visualisation of mean, q25 and q75: Stringy -----------------------------------

#Stringy
dev.off()
par(bg = NA)
plot(mq_norm_Stringy_output_final[,1], mq_norm_Stringy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "n - no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Stringy", sub = "Normal distribution")
lines(mq_norm_Stringy_output_final[,1], mq_norm_Stringy_output_final[,3],
      col = "dark green")
lines(mq_norm_Stringy_output_final[,1], mq_norm_Stringy_output_final[,4],
      col = "dark red")
dev.off()

Analysis of scagnostic measures: Stringy ------------------------------------------------------------------

  #preparing data set
  osg_norm_Stringy = t(output_scag_measure)
dim(osg_norm_Stringy)
rownames(osg_norm_Stringy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_norm_Stringy = function(x){

  mean_norm_Stringy = mean(x)
  q_norm_Stringy = quantile(x)
  mq_norm_Stringy = as.data.frame(rbind(mean_norm_Stringy, q_norm_Stringy[2],
                                        q_norm_Stringy[4]))
  rownames(mq_norm_Stringy) = c("mean", "25%", "75%")
  return (mq_norm_Stringy)
}

#test
mq_norm_Stringy(osg_norm_Stringy[1,])

mean_norm_Stringy = mean(osg_norm_Stringy[1,])
q_norm_Stringy = quantile(osg_norm_Stringy[1,])


mq_norm_Stringy_output = apply(osg_norm_Stringy, 1, mq_norm_Stringy)

unlist(mq_norm_Stringy_output)

mq_norm_Stringy_output_final = matrix(unlist(mq_norm_Stringy_output),
                                      ncol = 3, by = TRUE)
dim(mq_norm_Stringy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_norm Skewed_mq.csv"))

mq_norm_Stringy_output_final = cbind(vec_ns, mq_norm_Stringy_output_final)
colnames(mq_norm_Stringy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_norm_Stringy_output_final) = vec_ns

mq_norm_Stringy_output_final



# 9 Analysis of scagnostic measures: Monotonic ------------------------------------------------------------------

#preparing data set
osg_norm_Monotonic = t(output_scag_measure)
dim(osg_norm_Monotonic)
rownames(osg_norm_Monotonic) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_norm_Monotonic = function(x){

  mean_norm_Monotonic = mean(x)
  q_norm_Monotonic = quantile(x)
  mq_norm_Monotonic = as.data.frame(rbind(mean_norm_Monotonic, q_norm_Monotonic[2],
                                        q_norm_Stringy[4]))
  rownames(mq_norm_Monotonic) = c("mean", "25%", "75%")
  return (mq_norm_Monotonic)
}

#test
mq_norm_Monotonic(osg_norm_Monotonic[1,])

mean_norm_Monotonic = mean(osg_norm_Monotonic[1,])
q_norm_Monotonic = quantile(osg_norm_Monotonic[1,])


mq_norm_Monotonic_output = apply(osg_norm_Monotonic, 1, mq_norm_Monotonic)

unlist(mq_norm_Monotonic_output)

mq_norm_Monotonic_output_final = matrix(unlist(mq_norm_Monotonic_output),
                                      ncol = 3, by = TRUE)
dim(mq_norm_Monotonic_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_norm Skewed_mq.csv"))

mq_norm_Monotonic_output_final = cbind(vec_ns, mq_norm_Monotonic_output_final)
colnames(mq_norm_Monotonic_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_norm_Monotonic_output_final) = vec_ns

mq_norm_Monotonic_output_final

# 9 Visualisation of mean, q25 and q75: Monotonic --------------------------------

#Monotonic
dev.off()
par(bg = NA)
plot(mq_norm_Monotonic_output_final[,1], mq_norm_Monotonic_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "n - no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Monotonic", sub = "Normal distribution")
lines(mq_norm_Monotonic_output_final[,1], mq_norm_Monotonic_output_final[,3],
      col = "dark green")
lines(mq_norm_Monotonic_output_final[,1], mq_norm_Monotonic_output_final[,4],
      col = "dark red")
dev.off()
