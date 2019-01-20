rm(list = ls())
options(java.parameters = "-Xmx8000m")

library(scagnostics)
library(alphahull)

# Functions  --------------------------------------------------------------

#function to simulate Parable distribution
scag_par = function(n){

  x1 = rnorm(n, 1, 25)
  x2 = rnorm(n, 4, 85)

  a = -0.7
  b = 1
  c = -100

  y1 = a*x1^2 + b*x2 + c

  df_par = as.data.frame(cbind(x1, y1))
  s_par = scagnostics(df_par)

  par_eval = rbind(s_par)
}

#test = scag_par(100)

#function to add to the same list the each simulation iteration
scag_simul = function (sim_settings_n, sim_no_rep){
  #sim_settings_n = 5
  #sim_no_rep = 3
  scag_simul_output = scag_par(sim_settings_n)
  for (i in (1:(sim_no_rep-1))){
    #print(paste0("i",i, "scag_simul_output:",scag_simul_output ))
    added_sim_round = scag_par(sim_settings_n)
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

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_par Skewed.csv"))


# 1 Analysis of scagnostic measures: Outlying----------------------------------------------

#Outlying measure

#preparing data set
osg_par_Outlying = t(output_scag_measure)
dim(osg_par_Outlying)
rownames(osg_par_Outlying) = vec_ns



#Calculating the mean, q25 and q75
#i = 1
mq_par_Outlying = function(x){

  mean_par_Outlying = mean(x)
  q_par_Outlying = quantile(x)
  mq_par_Outlying = as.data.frame(rbind(mean_par_Outlying, q_par_Outlying[2],
                                        q_par_Outlying[4]))
  rownames(mq_par_Outlying) = c("mean", "25%", "75%")
  return (mq_par_Outlying)
}

#test
mq_par_Outlying(osg_par_Outlying[1,])

mean_par_Outlying = mean(osg_par_Outlying[1,])
q_par_Outlying = quantile(osg_par_Outlying[1,])

mq_par_Outlying_output = apply(osg_par_Outlying, 1, mq_par_Outlying)

unlist(mq_par_Outlying_output)
#dim(unlist(mq_par_Outlying_output))

mq_par_Outlying_output_final = matrix(unlist(mq_par_Outlying_output),
                                      ncol = 3, by = TRUE)
#dim(mq_par_Outlying_output_final)

mq_par_Outlying_output_final = cbind(vec_ns,mq_par_Outlying_output_final)
colnames(mq_par_Outlying_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_par_Outlying_output_final) = vec_ns

mq_par_Outlying_output_final


# 1 Visualisation of mean, q25 and q75: Outlying --------------------------------------

#Outlying
dev.off()
par(bg = NA)
plot(mq_par_Outlying_output_final[,1], mq_par_Outlying_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "blue",
     main = "Coefficient: Outlying", sub = "Parable distribution")
lines(mq_par_Outlying_output_final[,1], mq_par_Outlying_output_final[,3],
      col = "dark green")
lines(mq_par_Outlying_output_final[,1], mq_par_Outlying_output_final[,4],
      col = "dark red")
dev.off()




# 2 Analysis of scagnostic measures: Skewed ------------------------------------------------------------------

#preparing data set
osg_par_Skewed = t(output_scag_measure)
dim(osg_par_Skewed)
rownames(osg_par_Skewed) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_par_Skewed = function(x){

  mean_par_Skewed = mean(x)
  q_par_Skewed = quantile(x)
  mq_par_Skewed = as.data.frame(rbind(mean_par_Skewed, q_par_Skewed[2],
                                      q_par_Skewed[4]))
  rownames(mq_par_Skewed) = c("mean", "25%", "75%")
  return (mq_par_Skewed)
}

#test
mq_par_Skewed(osg_par_Skewed[1,])

mean_par_Skewed = mean(osg_par_Skewed[1,])
q_par_Skewed = quantile(osg_par_Skewed[1,])

mq_par_Skewed_output = apply(osg_par_Skewed, 1, mq_par_Skewed)

unlist(mq_par_Skewed_output)
#dim(unlist(mq_par_Outlying_output))

mq_par_Skewed_output_final = matrix(unlist(mq_par_Skewed_output),
                                    ncol = 3, by = TRUE)
dim(mq_par_Skewed_output_final)


mq_par_Skewed_output_final = cbind(vec_ns,mq_par_Skewed_output_final)
colnames(mq_par_Skewed_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_par_Skewed_output_final) = vec_ns

mq_par_Skewed_output_final

# 2 Visualisation of mean, q25 and q75: Skewed -----------------------------------

#Skewed

dev.off()
par(bg = NA)
plot(mq_par_Skewed_output_final[,1], mq_par_Skewed_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skewed", sub = "Parable Distribution")
lines(mq_par_Skewed_output_final[,1], mq_par_Skewed_output_final[,3],
      col = "dark green")
lines(mq_par_Skewed_output_final[,1], mq_par_Skewed_output_final[,4],
      col = "dark red")
dev.off()




# 3 Analysis of scagnostic measures: Clumpy ------------------------------------------------------------------

#preparing data set
osg_par_Clumpy = t(output_scag_measure)
dim(osg_par_Clumpy)
rownames(osg_par_Clumpy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_par_Clumpy = function(x){

  mean_par_Clumpy = mean(x)
  q_par_Clumpy = quantile(x)
  mq_par_Clumpy = as.data.frame(rbind(mean_par_Clumpy, q_par_Clumpy[2],
                                      q_par_Clumpy[4]))
  rownames(mq_par_Clumpy) = c("mean", "25%", "75%")
  return (mq_par_Clumpy)
}

#test
mq_par_Clumpy(osg_par_Clumpy[1,])

mean_par_Clumpy = mean(osg_par_Clumpy[1,])
q_par_Clumpy = quantile(osg_par_Clumpy[1,])


mq_par_Clumpy_output = apply(osg_par_Clumpy, 1, mq_par_Clumpy)

unlist(mq_par_Clumpy_output)

mq_par_Clumpy_output_final = matrix(unlist(mq_par_Clumpy_output),
                                    ncol = 3, by = TRUE)
dim(mq_par_Clumpy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_par Skewed_mq.csv"))

mq_par_Clumpy_output_final = cbind(vec_ns,mq_par_Clumpy_output_final)
colnames(mq_par_Clumpy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_par_Clumpy_output_final) = vec_ns

mq_par_Clumpy_output_final

# 3 Visualisation of mean, q25 and q75: Clumpy -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_par_Clumpy_output_final[,1], mq_par_Clumpy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Clumpy", sub = "Parable distribution")
lines(mq_par_Clumpy_output_final[,1], mq_par_Clumpy_output_final[,3],
      col = "dark green")
lines(mq_par_Clumpy_output_final[,1], mq_par_Clumpy_output_final[,4],
      col = "dark red")
dev.off()




# 4 Analysis of scagnostic measures: Sparse ------------------------------------------------------------------

#preparing data set
osg_par_Sparse = t(output_scag_measure)
dim(osg_par_Sparse)
rownames(osg_par_Sparse) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_par_Sparse = function(x){

  mean_par_Sparse = mean(x)
  q_par_Sparse = quantile(x)
  mq_par_Sparse = as.data.frame(rbind(mean_par_Sparse, q_par_Sparse[2],
                                      q_par_Sparse[4]))
  rownames(mq_par_Sparse) = c("mean", "25%", "75%")
  return (mq_par_Sparse)
}

#test
mq_par_Sparse(osg_par_Sparse[1,])

mean_par_Sparse = mean(osg_par_Sparse[1,])
q_par_Sparse = quantile(osg_par_Sparse[1,])


mq_par_Sparse_output = apply(osg_par_Sparse, 1, mq_par_Sparse)

unlist(mq_par_Sparse_output)

mq_par_Sparse_output_final = matrix(unlist(mq_par_Sparse_output),
                                    ncol = 3, by = TRUE)
dim(mq_par_Sparse_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_par Skewed_mq.csv"))

mq_par_Sparse_output_final = cbind(vec_ns,mq_par_Sparse_output_final)
colnames(mq_par_Sparse_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_par_Sparse_output_final) = vec_ns

mq_par_Sparse_output_final

# 4 Visualisation of mean, q25 and q75: Sparse -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_par_Sparse_output_final[,1], mq_par_Sparse_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "25q, mean, q75", col = "dark blue",
     main = "Coefficient: Sparse", sub = "Parable distribution")
lines(mq_par_Sparse_output_final[,1], mq_par_Sparse_output_final[,3],
      col = "dark green")
lines(mq_par_Sparse_output_final[,1], mq_par_Sparse_output_final[,4],
      col = "dark red")
dev.off()


# 5 Analysis of scagnostic measures: Striate ------------------------------------------------------------------

#preparing data set
osg_par_Striate = t(output_scag_measure)
dim(osg_par_Striate)
rownames(osg_par_Striate) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_par_Striate = function(x){

  mean_par_Striate = mean(x)
  q_par_Striate = quantile(x)
  mq_par_Striate = as.data.frame(rbind(mean_par_Striate, q_par_Striate[2],
                                       q_par_Striate[4]))
  rownames(mq_par_Striate) = c("mean", "25%", "75%")
  return (mq_par_Striate)
}

#test
mq_par_Striate(osg_par_Striate[1,])

mean_par_Striate = mean(osg_par_Striate[1,])
q_par_Striate = quantile(osg_par_Striate[1,])


mq_par_Striate_output = apply(osg_par_Striate, 1, mq_par_Striate)

unlist(mq_par_Striate_output)

mq_par_Striate_output_final = matrix(unlist(mq_par_Striate_output),
                                     ncol = 3, by = TRUE)
dim(mq_par_Striate_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_par Skewed_mq.csv"))

mq_par_Striate_output_final = cbind(vec_ns,mq_par_Striate_output_final)
colnames(mq_par_Striate_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_par_Striate_output_final) = vec_ns

mq_par_Striate_output_final

# 5 Visualisation of mean, q25 and q75: Striate -----------------------------------

#Striate
dev.off()
par(bg = NA)
plot(mq_par_Striate_output_final[,1], mq_par_Striate_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Striated", sub = "Parable distribution")
lines(mq_par_Striate_output_final[,1], mq_par_Striate_output_final[,3],
      col = "dark green")
lines(mq_par_Striate_output_final[,1], mq_par_Striate_output_final[,4],
      col = "dark red")
dev.off()



# 6 Analysis of scagnostic measures: Convex ------------------------------------------------------------------

#preparing data set
osg_par_Convex = t(output_scag_measure)
dim(osg_par_Convex)
rownames(osg_par_Convex) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_par_Convex = function(x){

  mean_par_Convex = mean(x)
  q_par_Convex = quantile(x)
  mq_par_Convex = as.data.frame(rbind(mean_par_Convex, q_par_Convex[2],
                                      q_par_Convex[4]))
  rownames(mq_par_Convex) = c("mean", "25%", "75%")
  return (mq_par_Convex)
}

#test
mq_par_Convex(osg_par_Convex[1,])

mean_par_Convex = mean(osg_par_Convex[1,])
q_par_Convex = quantile(osg_par_Convex[1,])


mq_par_Convex_output = apply(osg_par_Convex, 1, mq_par_Convex)

unlist(mq_par_Convex_output)

mq_par_Convex_output_final = matrix(unlist(mq_par_Convex_output),
                                    ncol = 3, by = TRUE)
dim(mq_par_Convex_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_par Skewed_mq.csv"))

mq_par_Convex_output_final = cbind(vec_ns, mq_par_Convex_output_final)
colnames(mq_par_Convex_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_par_Convex_output_final) = vec_ns

mq_par_Convex_output_final

# 6 Visualisation of mean, q25 and q75: Convex -----------------------------------

#Convex
dev.off()
par(bg = NA)
plot(mq_par_Convex_output_final[,1], mq_par_Convex_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Convex", sub = "Parable distribution")
lines(mq_par_Convex_output_final[,1], mq_par_Convex_output_final[,3],
      col = "dark green")
lines(mq_par_Convex_output_final[,1], mq_par_Convex_output_final[,4],
      col = "dark red")
dev.off()


# 7 Analysis of scagnostic measures: Skinny ------------------------------------------------------------------

#preparing data set
osg_par_Skinny = t(output_scag_measure)
dim(osg_par_Skinny)
rownames(osg_par_Skinny) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_par_Skinny = function(x){

  mean_par_Skinny = mean(x)
  q_par_Skinny = quantile(x)
  mq_par_Skinny = as.data.frame(rbind(mean_par_Skinny, q_par_Skinny[2],
                                      q_par_Skinny[4]))
  rownames(mq_par_Skinny) = c("mean", "25%", "75%")
  return (mq_par_Skinny)
}

#test
mq_par_Skinny(osg_par_Skinny[1,])

mean_par_Skinny = mean(osg_par_Skinny[1,])
q_par_Skinny = quantile(osg_par_Skinny[1,])


mq_par_Skinny_output = apply(osg_par_Skinny, 1, mq_par_Skinny)

unlist(mq_par_Skinny_output)

mq_par_Skinny_output_final = matrix(unlist(mq_par_Skinny_output),
                                    ncol = 3, by = TRUE)
dim(mq_par_Skinny_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_par Skewed_mq.csv"))

mq_par_Skinny_output_final = cbind(vec_ns, mq_par_Skinny_output_final)
colnames(mq_par_Skinny_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_par_Skinny_output_final) = vec_ns

mq_par_Skinny_output_final

# 7 Visualisation of mean, q25 and q75: Skinny -----------------------------------

#Skinny
dev.off()
par(bg = NA)
plot(mq_par_Skinny_output_final[,1], mq_par_Skinny_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skinny", sub = "Parable distribution")
lines(mq_par_Skinny_output_final[,1], mq_par_Skinny_output_final[,3],
      col = "dark green")
lines(mq_par_Skinny_output_final[,1], mq_par_Skinny_output_final[,4],
      col = "dark red")
dev.off()


# 8 Analysis of scagnostic measures: Stringy ------------------------------------------------------------------

#preparing data set
osg_par_Stringy = t(output_scag_measure)
dim(osg_par_Stringy)
rownames(osg_par_Stringy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_par_Stringy = function(x){

  mean_par_Stringy = mean(x)
  q_par_Stringy = quantile(x)
  mq_par_Stringy = as.data.frame(rbind(mean_par_Stringy, q_par_Stringy[2],
                                       q_par_Stringy[4]))
  rownames(mq_par_Stringy) = c("mean", "25%", "75%")
  return (mq_par_Stringy)
}

#test
mq_par_Stringy(osg_par_Stringy[1,])

mean_par_Stringy = mean(osg_par_Stringy[1,])
q_par_Stringy = quantile(osg_par_Stringy[1,])


mq_par_Stringy_output = apply(osg_par_Stringy, 1, mq_par_Stringy)

unlist(mq_par_Stringy_output)

mq_par_Stringy_output_final = matrix(unlist(mq_par_Stringy_output),
                                     ncol = 3, by = TRUE)
dim(mq_par_Stringy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_par Skewed_mq.csv"))

mq_par_Stringy_output_final = cbind(vec_ns, mq_par_Stringy_output_final)
colnames(mq_par_Stringy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_par_Stringy_output_final) = vec_ns

mq_par_Stringy_output_final

# 8 Visualisation of mean, q25 and q75: Stringy -----------------------------------

#Stringy
dev.off()
par(bg = NA)
plot(mq_par_Stringy_output_final[,1], mq_par_Stringy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Stringy", sub = "Parable distribution")
lines(mq_par_Stringy_output_final[,1], mq_par_Stringy_output_final[,3],
      col = "dark green")
lines(mq_par_Stringy_output_final[,1], mq_par_Stringy_output_final[,4],
      col = "dark red")
dev.off()



# 9 Analysis of scagnostic measures: Monotonic ------------------------------------------------------------------

#preparing data set
osg_par_Monotonic = t(output_scag_measure)
dim(osg_par_Monotonic)
rownames(osg_par_Monotonic) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_par_Monotonic = function(x){

  mean_par_Monotonic = mean(x)
  q_par_Monotonic = quantile(x)
  mq_par_Monotonic = as.data.frame(rbind(mean_par_Monotonic, q_par_Monotonic[2],
                                         q_par_Stringy[4]))
  rownames(mq_par_Monotonic) = c("mean", "25%", "75%")
  return (mq_par_Monotonic)
}

#test
mq_par_Monotonic(osg_par_Monotonic[1,])

mean_par_Monotonic = mean(osg_par_Monotonic[1,])
q_par_Monotonic = quantile(osg_par_Monotonic[1,])


mq_par_Monotonic_output = apply(osg_par_Monotonic, 1, mq_par_Monotonic)

unlist(mq_par_Monotonic_output)

mq_par_Monotonic_output_final = matrix(unlist(mq_par_Monotonic_output),
                                       ncol = 3, by = TRUE)
dim(mq_par_Monotonic_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_par Skewed_mq.csv"))

mq_par_Monotonic_output_final = cbind(vec_ns, mq_par_Monotonic_output_final)
colnames(mq_par_Monotonic_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_par_Monotonic_output_final) = vec_ns

mq_par_Monotonic_output_final

# 9 Visualisation of mean, q25 and q75: Monotonic --------------------------------

#Monotonic
dev.off()
par(bg = NA)
plot(mq_par_Monotonic_output_final[,1], mq_par_Monotonic_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Monotonic", sub = "Parable distribution")
lines(mq_par_Monotonic_output_final[,1], mq_par_Monotonic_output_final[,3],
      col = "dark green")
lines(mq_par_Monotonic_output_final[,1], mq_par_Monotonic_output_final[,4],
      col = "dark red")
dev.off()
