rm(list = ls())
options(java.parameters = "-Xmx8000m")

library(scagnostics)
library(alphahull)

# Functions  --------------------------------------------------------------

#function to simulate Clumps distribution
scag_clu = function(n){

  x1 = rnorm(n, -1, 0.3)
  y1 = rnorm(n, -1, 0.3)


  x2 = rnorm(n, 1, 0.3)
  y2 = rnorm(n, -1, 0.3)

  x3 = rnorm(n, 0, 0.3)
  y3 = rnorm(n, 1, 0.3)


  df_clu = cbind(c(x1, x2, x3), c(y1, y2, y3))
  s_clu = scagnostics(df_clu)

  clu_eval = rbind(s_clu)
}

#test = scag_clu(100)

#function to add to the same list the each simulation iteration
scag_simul = function (sim_settings_n, sim_no_rep){
  #sim_settings_n = 5
  #sim_no_rep = 3
  scag_simul_output = scag_clu(sim_settings_n)
  for (i in (1:(sim_no_rep-1))){
    #print(paste0("i",i, "scag_simul_output:",scag_simul_output ))
    added_sim_round = scag_clu(sim_settings_n)
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

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_clu Skewed.csv"))


# 1 Analysis of scagnostic measures: Outlying----------------------------------------------

#Outlying measure

#preparing data set
osg_clu_Outlying = t(output_scag_measure)
dim(osg_clu_Outlying)
rownames(osg_clu_Outlying) = vec_ns



#Calculating the mean, q25 and q75
#i = 1
mq_clu_Outlying = function(x){

  mean_clu_Outlying = mean(x)
  q_clu_Outlying = quantile(x)
  mq_clu_Outlying = as.data.frame(rbind(mean_clu_Outlying, q_clu_Outlying[2],
                                        q_clu_Outlying[4]))
  rownames(mq_clu_Outlying) = c("mean", "25%", "75%")
  return (mq_clu_Outlying)
}

#test
mq_clu_Outlying(osg_clu_Outlying[1,])

mean_clu_Outlying = mean(osg_clu_Outlying[1,])
q_clu_Outlying = quantile(osg_clu_Outlying[1,])

mq_clu_Outlying_output = apply(osg_clu_Outlying, 1, mq_clu_Outlying)

unlist(mq_clu_Outlying_output)
#dim(unlist(mq_clu_Outlying_output))

mq_clu_Outlying_output_final = matrix(unlist(mq_clu_Outlying_output),
                                      ncol = 3, by = TRUE)
#dim(mq_clu_Outlying_output_final)

mq_clu_Outlying_output_final = cbind(vec_ns,mq_clu_Outlying_output_final)
colnames(mq_clu_Outlying_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_clu_Outlying_output_final) = vec_ns

mq_clu_Outlying_output_final


# 1 Visualisation of mean, q25 and q75: Outlying --------------------------------------

#Outlying
dev.off()
par(bg = NA)
plot(mq_clu_Outlying_output_final[,1], mq_clu_Outlying_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "blue",
     main = "Coefficient: Outlying", sub = "Clumps distribution")
lines(mq_clu_Outlying_output_final[,1], mq_clu_Outlying_output_final[,3],
      col = "dark green")
lines(mq_clu_Outlying_output_final[,1], mq_clu_Outlying_output_final[,4],
      col = "dark red")
dev.off()




# 2 Analysis of scagnostic measures: Skewed ------------------------------------------------------------------

#preparing data set
osg_clu_Skewed = t(output_scag_measure)
dim(osg_clu_Skewed)
rownames(osg_clu_Skewed) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_clu_Skewed = function(x){

  mean_clu_Skewed = mean(x)
  q_clu_Skewed = quantile(x)
  mq_clu_Skewed = as.data.frame(rbind(mean_clu_Skewed, q_clu_Skewed[2],
                                      q_clu_Skewed[4]))
  rownames(mq_clu_Skewed) = c("mean", "25%", "75%")
  return (mq_clu_Skewed)
}

#test
mq_clu_Skewed(osg_clu_Skewed[1,])

mean_clu_Skewed = mean(osg_clu_Skewed[1,])
q_clu_Skewed = quantile(osg_clu_Skewed[1,])

mq_clu_Skewed_output = apply(osg_clu_Skewed, 1, mq_clu_Skewed)

unlist(mq_clu_Skewed_output)
#dim(unlist(mq_clu_Outlying_output))

mq_clu_Skewed_output_final = matrix(unlist(mq_clu_Skewed_output),
                                    ncol = 3, by = TRUE)
dim(mq_clu_Skewed_output_final)


mq_clu_Skewed_output_final = cbind(vec_ns,mq_clu_Skewed_output_final)
colnames(mq_clu_Skewed_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_clu_Skewed_output_final) = vec_ns

mq_clu_Skewed_output_final

# 2 Visualisation of mean, q25 and q75: Skewed -----------------------------------

#Skewed

dev.off()
par(bg = NA)
plot(mq_clu_Skewed_output_final[,1], mq_clu_Skewed_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skewed", sub = "Clumps Distribution")
lines(mq_clu_Skewed_output_final[,1], mq_clu_Skewed_output_final[,3],
      col = "dark green")
lines(mq_clu_Skewed_output_final[,1], mq_clu_Skewed_output_final[,4],
      col = "dark red")
dev.off()




# 3 Analysis of scagnostic measures: Clumpy ------------------------------------------------------------------

#preparing data set
osg_clu_Clumpy = t(output_scag_measure)
dim(osg_clu_Clumpy)
rownames(osg_clu_Clumpy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_clu_Clumpy = function(x){

  mean_clu_Clumpy = mean(x)
  q_clu_Clumpy = quantile(x)
  mq_clu_Clumpy = as.data.frame(rbind(mean_clu_Clumpy, q_clu_Clumpy[2],
                                      q_clu_Clumpy[4]))
  rownames(mq_clu_Clumpy) = c("mean", "25%", "75%")
  return (mq_clu_Clumpy)
}

#test
mq_clu_Clumpy(osg_clu_Clumpy[1,])

mean_clu_Clumpy = mean(osg_clu_Clumpy[1,])
q_clu_Clumpy = quantile(osg_clu_Clumpy[1,])


mq_clu_Clumpy_output = apply(osg_clu_Clumpy, 1, mq_clu_Clumpy)

unlist(mq_clu_Clumpy_output)

mq_clu_Clumpy_output_final = matrix(unlist(mq_clu_Clumpy_output),
                                    ncol = 3, by = TRUE)
dim(mq_clu_Clumpy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_clu Skewed_mq.csv"))

mq_clu_Clumpy_output_final = cbind(vec_ns,mq_clu_Clumpy_output_final)
colnames(mq_clu_Clumpy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_clu_Clumpy_output_final) = vec_ns

mq_clu_Clumpy_output_final

# 3 Visualisation of mean, q25 and q75: Clumpy -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_clu_Clumpy_output_final[,1], mq_clu_Clumpy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Clumpy", sub = "Clumps distribution")
lines(mq_clu_Clumpy_output_final[,1], mq_clu_Clumpy_output_final[,3],
      col = "dark green")
lines(mq_clu_Clumpy_output_final[,1], mq_clu_Clumpy_output_final[,4],
      col = "dark red")
dev.off()




# 4 Analysis of scagnostic measures: Sparse ------------------------------------------------------------------

#preparing data set
osg_clu_Sparse = t(output_scag_measure)
dim(osg_clu_Sparse)
rownames(osg_clu_Sparse) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_clu_Sparse = function(x){

  mean_clu_Sparse = mean(x)
  q_clu_Sparse = quantile(x)
  mq_clu_Sparse = as.data.frame(rbind(mean_clu_Sparse, q_clu_Sparse[2],
                                      q_clu_Sparse[4]))
  rownames(mq_clu_Sparse) = c("mean", "25%", "75%")
  return (mq_clu_Sparse)
}

#test
mq_clu_Sparse(osg_clu_Sparse[1,])

mean_clu_Sparse = mean(osg_clu_Sparse[1,])
q_clu_Sparse = quantile(osg_clu_Sparse[1,])


mq_clu_Sparse_output = apply(osg_clu_Sparse, 1, mq_clu_Sparse)

unlist(mq_clu_Sparse_output)

mq_clu_Sparse_output_final = matrix(unlist(mq_clu_Sparse_output),
                                    ncol = 3, by = TRUE)
dim(mq_clu_Sparse_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_clu Skewed_mq.csv"))

mq_clu_Sparse_output_final = cbind(vec_ns,mq_clu_Sparse_output_final)
colnames(mq_clu_Sparse_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_clu_Sparse_output_final) = vec_ns

mq_clu_Sparse_output_final

# 4 Visualisation of mean, q25 and q75: Sparse -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_clu_Sparse_output_final[,1], mq_clu_Sparse_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "25q, mean, q75", col = "dark blue",
     main = "Coefficient: Sparse", sub = "Clumps distribution")
lines(mq_clu_Sparse_output_final[,1], mq_clu_Sparse_output_final[,3],
      col = "dark green")
lines(mq_clu_Sparse_output_final[,1], mq_clu_Sparse_output_final[,4],
      col = "dark red")
dev.off()


# 5 Analysis of scagnostic measures: Striate ------------------------------------------------------------------

#preparing data set
osg_clu_Striate = t(output_scag_measure)
dim(osg_clu_Striate)
rownames(osg_clu_Striate) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_clu_Striate = function(x){

  mean_clu_Striate = mean(x)
  q_clu_Striate = quantile(x)
  mq_clu_Striate = as.data.frame(rbind(mean_clu_Striate, q_clu_Striate[2],
                                       q_clu_Striate[4]))
  rownames(mq_clu_Striate) = c("mean", "25%", "75%")
  return (mq_clu_Striate)
}

#test
mq_clu_Striate(osg_clu_Striate[1,])

mean_clu_Striate = mean(osg_clu_Striate[1,])
q_clu_Striate = quantile(osg_clu_Striate[1,])


mq_clu_Striate_output = apply(osg_clu_Striate, 1, mq_clu_Striate)

unlist(mq_clu_Striate_output)

mq_clu_Striate_output_final = matrix(unlist(mq_clu_Striate_output),
                                     ncol = 3, by = TRUE)
dim(mq_clu_Striate_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_clu Skewed_mq.csv"))

mq_clu_Striate_output_final = cbind(vec_ns,mq_clu_Striate_output_final)
colnames(mq_clu_Striate_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_clu_Striate_output_final) = vec_ns

mq_clu_Striate_output_final

# 5 Visualisation of mean, q25 and q75: Striate -----------------------------------

#Striate
dev.off()
par(bg = NA)
plot(mq_clu_Striate_output_final[,1], mq_clu_Striate_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Striated", sub = "Clumps distribution")
lines(mq_clu_Striate_output_final[,1], mq_clu_Striate_output_final[,3],
      col = "dark green")
lines(mq_clu_Striate_output_final[,1], mq_clu_Striate_output_final[,4],
      col = "dark red")
dev.off()



# 6 Analysis of scagnostic measures: Convex ------------------------------------------------------------------

#preparing data set
osg_clu_Convex = t(output_scag_measure)
dim(osg_clu_Convex)
rownames(osg_clu_Convex) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_clu_Convex = function(x){

  mean_clu_Convex = mean(x)
  q_clu_Convex = quantile(x)
  mq_clu_Convex = as.data.frame(rbind(mean_clu_Convex, q_clu_Convex[2],
                                      q_clu_Convex[4]))
  rownames(mq_clu_Convex) = c("mean", "25%", "75%")
  return (mq_clu_Convex)
}

#test
mq_clu_Convex(osg_clu_Convex[1,])

mean_clu_Convex = mean(osg_clu_Convex[1,])
q_clu_Convex = quantile(osg_clu_Convex[1,])


mq_clu_Convex_output = apply(osg_clu_Convex, 1, mq_clu_Convex)

unlist(mq_clu_Convex_output)

mq_clu_Convex_output_final = matrix(unlist(mq_clu_Convex_output),
                                    ncol = 3, by = TRUE)
dim(mq_clu_Convex_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_clu Skewed_mq.csv"))

mq_clu_Convex_output_final = cbind(vec_ns, mq_clu_Convex_output_final)
colnames(mq_clu_Convex_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_clu_Convex_output_final) = vec_ns

mq_clu_Convex_output_final

# 6 Visualisation of mean, q25 and q75: Convex -----------------------------------

#Convex
dev.off()
par(bg = NA)
plot(mq_clu_Convex_output_final[,1], mq_clu_Convex_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Convex", sub = "Clumps distribution")
lines(mq_clu_Convex_output_final[,1], mq_clu_Convex_output_final[,3],
      col = "dark green")
lines(mq_clu_Convex_output_final[,1], mq_clu_Convex_output_final[,4],
      col = "dark red")
dev.off()


# 7 Analysis of scagnostic measures: Skinny ------------------------------------------------------------------

#preparing data set
osg_clu_Skinny = t(output_scag_measure)
dim(osg_clu_Skinny)
rownames(osg_clu_Skinny) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_clu_Skinny = function(x){

  mean_clu_Skinny = mean(x)
  q_clu_Skinny = quantile(x)
  mq_clu_Skinny = as.data.frame(rbind(mean_clu_Skinny, q_clu_Skinny[2],
                                      q_clu_Skinny[4]))
  rownames(mq_clu_Skinny) = c("mean", "25%", "75%")
  return (mq_clu_Skinny)
}

#test
mq_clu_Skinny(osg_clu_Skinny[1,])

mean_clu_Skinny = mean(osg_clu_Skinny[1,])
q_clu_Skinny = quantile(osg_clu_Skinny[1,])


mq_clu_Skinny_output = apply(osg_clu_Skinny, 1, mq_clu_Skinny)

unlist(mq_clu_Skinny_output)

mq_clu_Skinny_output_final = matrix(unlist(mq_clu_Skinny_output),
                                    ncol = 3, by = TRUE)
dim(mq_clu_Skinny_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_clu Skewed_mq.csv"))

mq_clu_Skinny_output_final = cbind(vec_ns, mq_clu_Skinny_output_final)
colnames(mq_clu_Skinny_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_clu_Skinny_output_final) = vec_ns

mq_clu_Skinny_output_final

# 7 Visualisation of mean, q25 and q75: Skinny -----------------------------------

#Skinny
dev.off()
par(bg = NA)
plot(mq_clu_Skinny_output_final[,1], mq_clu_Skinny_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skinny", sub = "Clumps distribution")
lines(mq_clu_Skinny_output_final[,1], mq_clu_Skinny_output_final[,3],
      col = "dark green")
lines(mq_clu_Skinny_output_final[,1], mq_clu_Skinny_output_final[,4],
      col = "dark red")
dev.off()


# 8 Analysis of scagnostic measures: Stringy ------------------------------------------------------------------

#preparing data set
osg_clu_Stringy = t(output_scag_measure)
dim(osg_clu_Stringy)
rownames(osg_clu_Stringy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_clu_Stringy = function(x){

  mean_clu_Stringy = mean(x)
  q_clu_Stringy = quantile(x)
  mq_clu_Stringy = as.data.frame(rbind(mean_clu_Stringy, q_clu_Stringy[2],
                                       q_clu_Stringy[4]))
  rownames(mq_clu_Stringy) = c("mean", "25%", "75%")
  return (mq_clu_Stringy)
}

#test
mq_clu_Stringy(osg_clu_Stringy[1,])

mean_clu_Stringy = mean(osg_clu_Stringy[1,])
q_clu_Stringy = quantile(osg_clu_Stringy[1,])


mq_clu_Stringy_output = apply(osg_clu_Stringy, 1, mq_clu_Stringy)

unlist(mq_clu_Stringy_output)

mq_clu_Stringy_output_final = matrix(unlist(mq_clu_Stringy_output),
                                     ncol = 3, by = TRUE)
dim(mq_clu_Stringy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_clu Skewed_mq.csv"))

mq_clu_Stringy_output_final = cbind(vec_ns, mq_clu_Stringy_output_final)
colnames(mq_clu_Stringy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_clu_Stringy_output_final) = vec_ns

mq_clu_Stringy_output_final

# 8 Visualisation of mean, q25 and q75: Stringy -----------------------------------

#Stringy
dev.off()
par(bg = NA)
plot(mq_clu_Stringy_output_final[,1], mq_clu_Stringy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Stringy", sub = "Clumps distribution")
lines(mq_clu_Stringy_output_final[,1], mq_clu_Stringy_output_final[,3],
      col = "dark green")
lines(mq_clu_Stringy_output_final[,1], mq_clu_Stringy_output_final[,4],
      col = "dark red")
dev.off()



# 9 Analysis of scagnostic measures: Monotonic ------------------------------------------------------------------

#preparing data set
osg_clu_Monotonic = t(output_scag_measure)
dim(osg_clu_Monotonic)
rownames(osg_clu_Monotonic) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_clu_Monotonic = function(x){

  mean_clu_Monotonic = mean(x)
  q_clu_Monotonic = quantile(x)
  mq_clu_Monotonic = as.data.frame(rbind(mean_clu_Monotonic, q_clu_Monotonic[2],
                                         q_clu_Stringy[4]))
  rownames(mq_clu_Monotonic) = c("mean", "25%", "75%")
  return (mq_clu_Monotonic)
}

#test
mq_clu_Monotonic(osg_clu_Monotonic[1,])

mean_clu_Monotonic = mean(osg_clu_Monotonic[1,])
q_clu_Monotonic = quantile(osg_clu_Monotonic[1,])


mq_clu_Monotonic_output = apply(osg_clu_Monotonic, 1, mq_clu_Monotonic)

unlist(mq_clu_Monotonic_output)

mq_clu_Monotonic_output_final = matrix(unlist(mq_clu_Monotonic_output),
                                       ncol = 3, by = TRUE)
dim(mq_clu_Monotonic_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_clu Skewed_mq.csv"))

mq_clu_Monotonic_output_final = cbind(vec_ns, mq_clu_Monotonic_output_final)
colnames(mq_clu_Monotonic_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_clu_Monotonic_output_final) = vec_ns

mq_clu_Monotonic_output_final

# 9 Visualisation of mean, q25 and q75: Monotonic --------------------------------

#Monotonic
dev.off()
par(bg = NA)
plot(mq_clu_Monotonic_output_final[,1], mq_clu_Monotonic_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Monotonic", sub = "Clumps distribution")
lines(mq_clu_Monotonic_output_final[,1], mq_clu_Monotonic_output_final[,3],
      col = "dark green")
lines(mq_clu_Monotonic_output_final[,1], mq_clu_Monotonic_output_final[,4],
      col = "dark red")
dev.off()
