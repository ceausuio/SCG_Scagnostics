rm(list = ls())
options(java.parameters = "-Xmx8000m")

library(scagnostics)
library(alphahull)

# Functions  --------------------------------------------------------------

#function to simulate unial distribution
scag_uni = function(n){

  #n = 100
  x = runif(n)
  y = runif(n)

  df_uni = as.data.frame(cbind(x,y))
  s_uni = scagnostics(df_uni)

  uni_eval = rbind(s_uni)
}

#test = scag_uni(100)

#function to add to the same list the each simulation iteration
scag_simul = function (sim_settings_n, sim_no_rep){
  #sim_settings_n = 5
  #sim_no_rep = 3
  scag_simul_output = scag_uni(sim_settings_n)
  for (i in (1:(sim_no_rep-1))){
    #print(paste0("i",i, "scag_simul_output:",scag_simul_output ))
    added_sim_round = scag_uni(sim_settings_n)
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

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_uni Skewed.csv"))


# 1 Analysis of scagnostic measures: Outlying----------------------------------------------

#Outlying measure

#preparing data set
osg_uni_Outlying = t(output_scag_measure)
dim(osg_uni_Outlying)
rownames(osg_uni_Outlying) = vec_ns



#Calculating the mean, q25 and q75
#i = 1
mq_uni_Outlying = function(x){

  mean_uni_Outlying = mean(x)
  q_uni_Outlying = quantile(x)
  mq_uni_Outlying = as.data.frame(rbind(mean_uni_Outlying, q_uni_Outlying[2],
                                         q_uni_Outlying[4]))
  rownames(mq_uni_Outlying) = c("mean", "25%", "75%")
  return (mq_uni_Outlying)
}

#test
mq_uni_Outlying(osg_uni_Outlying[1,])

mean_uni_Outlying = mean(osg_uni_Outlying[1,])
q_uni_Outlying = quantile(osg_uni_Outlying[1,])

mq_uni_Outlying_output = apply(osg_uni_Outlying, 1, mq_uni_Outlying)

unlist(mq_uni_Outlying_output)
#dim(unlist(mq_uni_Outlying_output))

mq_uni_Outlying_output_final = matrix(unlist(mq_uni_Outlying_output),
                                       ncol = 3, by = TRUE)
#dim(mq_uni_Outlying_output_final)

mq_uni_Outlying_output_final = cbind(vec_ns,mq_uni_Outlying_output_final)
colnames(mq_uni_Outlying_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_uni_Outlying_output_final) = vec_ns

mq_uni_Outlying_output_final


# 1 Visualisation of mean, q25 and q75: Outlying --------------------------------------

#Outlying
dev.off()
par(bg = NA)
plot(mq_uni_Outlying_output_final[,1], mq_uni_Outlying_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "blue",
     main = "Coefficient: Outlying", sub = "Uniform distribution")
lines(mq_uni_Outlying_output_final[,1], mq_uni_Outlying_output_final[,3],
      col = "dark green")
lines(mq_uni_Outlying_output_final[,1], mq_uni_Outlying_output_final[,4],
      col = "dark red")
dev.off()




# 2 Analysis of scagnostic measures: Skewed ------------------------------------------------------------------

#preparing data set
osg_uni_Skewed = t(output_scag_measure)
dim(osg_uni_Skewed)
rownames(osg_uni_Skewed) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_uni_Skewed = function(x){

  mean_uni_Skewed = mean(x)
  q_uni_Skewed = quantile(x)
  mq_uni_Skewed = as.data.frame(rbind(mean_uni_Skewed, q_uni_Skewed[2],
                                       q_uni_Skewed[4]))
  rownames(mq_uni_Skewed) = c("mean", "25%", "75%")
  return (mq_uni_Skewed)
}

#test
mq_uni_Skewed(osg_uni_Skewed[1,])

mean_uni_Skewed = mean(osg_uni_Skewed[1,])
q_uni_Skewed = quantile(osg_uni_Skewed[1,])

mq_uni_Skewed_output = apply(osg_uni_Skewed, 1, mq_uni_Skewed)

unlist(mq_uni_Skewed_output)
#dim(unlist(mq_uni_Outlying_output))

mq_uni_Skewed_output_final = matrix(unlist(mq_uni_Skewed_output),
                                     ncol = 3, by = TRUE)
dim(mq_uni_Skewed_output_final)


mq_uni_Skewed_output_final = cbind(vec_ns,mq_uni_Skewed_output_final)
colnames(mq_uni_Skewed_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_uni_Skewed_output_final) = vec_ns

mq_uni_Skewed_output_final

# 2 Visualisation of mean, q25 and q75: Skewed -----------------------------------

#Skewed

dev.off()
par(bg = NA)
plot(mq_uni_Skewed_output_final[,1], mq_uni_Skewed_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skewed", sub = "Uniform Distribution")
lines(mq_uni_Skewed_output_final[,1], mq_uni_Skewed_output_final[,3],
      col = "dark green")
lines(mq_uni_Skewed_output_final[,1], mq_uni_Skewed_output_final[,4],
      col = "dark red")
dev.off()




# 3 Analysis of scagnostic measures: Clumpy ------------------------------------------------------------------

#preparing data set
osg_uni_Clumpy = t(output_scag_measure)
dim(osg_uni_Clumpy)
rownames(osg_uni_Clumpy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_uni_Clumpy = function(x){

  mean_uni_Clumpy = mean(x)
  q_uni_Clumpy = quantile(x)
  mq_uni_Clumpy = as.data.frame(rbind(mean_uni_Clumpy, q_uni_Clumpy[2],
                                       q_uni_Clumpy[4]))
  rownames(mq_uni_Clumpy) = c("mean", "25%", "75%")
  return (mq_uni_Clumpy)
}

#test
mq_uni_Clumpy(osg_uni_Clumpy[1,])

mean_uni_Clumpy = mean(osg_uni_Clumpy[1,])
q_uni_Clumpy = quantile(osg_uni_Clumpy[1,])


mq_uni_Clumpy_output = apply(osg_uni_Clumpy, 1, mq_uni_Clumpy)

unlist(mq_uni_Clumpy_output)

mq_uni_Clumpy_output_final = matrix(unlist(mq_uni_Clumpy_output),
                                     ncol = 3, by = TRUE)
dim(mq_uni_Clumpy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_uni Skewed_mq.csv"))

mq_uni_Clumpy_output_final = cbind(vec_ns,mq_uni_Clumpy_output_final)
colnames(mq_uni_Clumpy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_uni_Clumpy_output_final) = vec_ns

mq_uni_Clumpy_output_final

# 3 Visualisation of mean, q25 and q75: Clumpy -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_uni_Clumpy_output_final[,1], mq_uni_Clumpy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Clumpy", sub = "Uniform distribution")
lines(mq_uni_Clumpy_output_final[,1], mq_uni_Clumpy_output_final[,3],
      col = "dark green")
lines(mq_uni_Clumpy_output_final[,1], mq_uni_Clumpy_output_final[,4],
      col = "dark red")
dev.off()




# 4 Analysis of scagnostic measures: Sparse ------------------------------------------------------------------

#preparing data set
osg_uni_Sparse = t(output_scag_measure)
dim(osg_uni_Sparse)
rownames(osg_uni_Sparse) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_uni_Sparse = function(x){

  mean_uni_Sparse = mean(x)
  q_uni_Sparse = quantile(x)
  mq_uni_Sparse = as.data.frame(rbind(mean_uni_Sparse, q_uni_Sparse[2],
                                       q_uni_Sparse[4]))
  rownames(mq_uni_Sparse) = c("mean", "25%", "75%")
  return (mq_uni_Sparse)
}

#test
mq_uni_Sparse(osg_uni_Sparse[1,])

mean_uni_Sparse = mean(osg_uni_Sparse[1,])
q_uni_Sparse = quantile(osg_uni_Sparse[1,])


mq_uni_Sparse_output = apply(osg_uni_Sparse, 1, mq_uni_Sparse)

unlist(mq_uni_Sparse_output)

mq_uni_Sparse_output_final = matrix(unlist(mq_uni_Sparse_output),
                                     ncol = 3, by = TRUE)
dim(mq_uni_Sparse_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_uni Skewed_mq.csv"))

mq_uni_Sparse_output_final = cbind(vec_ns,mq_uni_Sparse_output_final)
colnames(mq_uni_Sparse_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_uni_Sparse_output_final) = vec_ns

mq_uni_Sparse_output_final

# 4 Visualisation of mean, q25 and q75: Sparse -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_uni_Sparse_output_final[,1], mq_uni_Sparse_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "25q, mean, q75", col = "dark blue",
     main = "Coefficient: Sparse", sub = "Uniform distribution")
lines(mq_uni_Sparse_output_final[,1], mq_uni_Sparse_output_final[,3],
      col = "dark green")
lines(mq_uni_Sparse_output_final[,1], mq_uni_Sparse_output_final[,4],
      col = "dark red")
dev.off()


# 5 Analysis of scagnostic measures: Striate ------------------------------------------------------------------

#preparing data set
osg_uni_Striate = t(output_scag_measure)
dim(osg_uni_Striate)
rownames(osg_uni_Striate) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_uni_Striate = function(x){

  mean_uni_Striate = mean(x)
  q_uni_Striate = quantile(x)
  mq_uni_Striate = as.data.frame(rbind(mean_uni_Striate, q_uni_Striate[2],
                                        q_uni_Striate[4]))
  rownames(mq_uni_Striate) = c("mean", "25%", "75%")
  return (mq_uni_Striate)
}

#test
mq_uni_Striate(osg_uni_Striate[1,])

mean_uni_Striate = mean(osg_uni_Striate[1,])
q_uni_Striate = quantile(osg_uni_Striate[1,])


mq_uni_Striate_output = apply(osg_uni_Striate, 1, mq_uni_Striate)

unlist(mq_uni_Striate_output)

mq_uni_Striate_output_final = matrix(unlist(mq_uni_Striate_output),
                                      ncol = 3, by = TRUE)
dim(mq_uni_Striate_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_uni Skewed_mq.csv"))

mq_uni_Striate_output_final = cbind(vec_ns,mq_uni_Striate_output_final)
colnames(mq_uni_Striate_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_uni_Striate_output_final) = vec_ns

mq_uni_Striate_output_final

# 5 Visualisation of mean, q25 and q75: Striate -----------------------------------

#Striate
dev.off()
par(bg = NA)
plot(mq_uni_Striate_output_final[,1], mq_uni_Striate_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Striated", sub = "Uniform distribution")
lines(mq_uni_Striate_output_final[,1], mq_uni_Striate_output_final[,3],
      col = "dark green")
lines(mq_uni_Striate_output_final[,1], mq_uni_Striate_output_final[,4],
      col = "dark red")
dev.off()



# 6 Analysis of scagnostic measures: Convex ------------------------------------------------------------------

#preparing data set
osg_uni_Convex = t(output_scag_measure)
dim(osg_uni_Convex)
rownames(osg_uni_Convex) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_uni_Convex = function(x){

  mean_uni_Convex = mean(x)
  q_uni_Convex = quantile(x)
  mq_uni_Convex = as.data.frame(rbind(mean_uni_Convex, q_uni_Convex[2],
                                       q_uni_Convex[4]))
  rownames(mq_uni_Convex) = c("mean", "25%", "75%")
  return (mq_uni_Convex)
}

#test
mq_uni_Convex(osg_uni_Convex[1,])

mean_uni_Convex = mean(osg_uni_Convex[1,])
q_uni_Convex = quantile(osg_uni_Convex[1,])


mq_uni_Convex_output = apply(osg_uni_Convex, 1, mq_uni_Convex)

unlist(mq_uni_Convex_output)

mq_uni_Convex_output_final = matrix(unlist(mq_uni_Convex_output),
                                     ncol = 3, by = TRUE)
dim(mq_uni_Convex_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_uni Skewed_mq.csv"))

mq_uni_Convex_output_final = cbind(vec_ns, mq_uni_Convex_output_final)
colnames(mq_uni_Convex_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_uni_Convex_output_final) = vec_ns

mq_uni_Convex_output_final

# 6 Visualisation of mean, q25 and q75: Convex -----------------------------------

#Convex
dev.off()
par(bg = NA)
plot(mq_uni_Convex_output_final[,1], mq_uni_Convex_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Convex", sub = "Uniform distribution")
lines(mq_uni_Convex_output_final[,1], mq_uni_Convex_output_final[,3],
      col = "dark green")
lines(mq_uni_Convex_output_final[,1], mq_uni_Convex_output_final[,4],
      col = "dark red")
dev.off()


# 7 Analysis of scagnostic measures: Skinny ------------------------------------------------------------------

#preparing data set
osg_uni_Skinny = t(output_scag_measure)
dim(osg_uni_Skinny)
rownames(osg_uni_Skinny) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_uni_Skinny = function(x){

  mean_uni_Skinny = mean(x)
  q_uni_Skinny = quantile(x)
  mq_uni_Skinny = as.data.frame(rbind(mean_uni_Skinny, q_uni_Skinny[2],
                                       q_uni_Skinny[4]))
  rownames(mq_uni_Skinny) = c("mean", "25%", "75%")
  return (mq_uni_Skinny)
}

#test
mq_uni_Skinny(osg_uni_Skinny[1,])

mean_uni_Skinny = mean(osg_uni_Skinny[1,])
q_uni_Skinny = quantile(osg_uni_Skinny[1,])


mq_uni_Skinny_output = apply(osg_uni_Skinny, 1, mq_uni_Skinny)

unlist(mq_uni_Skinny_output)

mq_uni_Skinny_output_final = matrix(unlist(mq_uni_Skinny_output),
                                     ncol = 3, by = TRUE)
dim(mq_uni_Skinny_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_uni Skewed_mq.csv"))

mq_uni_Skinny_output_final = cbind(vec_ns, mq_uni_Skinny_output_final)
colnames(mq_uni_Skinny_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_uni_Skinny_output_final) = vec_ns

mq_uni_Skinny_output_final

# 7 Visualisation of mean, q25 and q75: Skinny -----------------------------------

#Skinny
dev.off()
par(bg = NA)
plot(mq_uni_Skinny_output_final[,1], mq_uni_Skinny_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skinny", sub = "Uniform distribution")
lines(mq_uni_Skinny_output_final[,1], mq_uni_Skinny_output_final[,3],
      col = "dark green")
lines(mq_uni_Skinny_output_final[,1], mq_uni_Skinny_output_final[,4],
      col = "dark red")
dev.off()


# 8 Analysis of scagnostic measures: Stringy ------------------------------------------------------------------

#preparing data set
osg_uni_Stringy = t(output_scag_measure)
dim(osg_uni_Stringy)
rownames(osg_uni_Stringy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_uni_Stringy = function(x){

  mean_uni_Stringy = mean(x)
  q_uni_Stringy = quantile(x)
  mq_uni_Stringy = as.data.frame(rbind(mean_uni_Stringy, q_uni_Stringy[2],
                                        q_uni_Stringy[4]))
  rownames(mq_uni_Stringy) = c("mean", "25%", "75%")
  return (mq_uni_Stringy)
}

#test
mq_uni_Stringy(osg_uni_Stringy[1,])

mean_uni_Stringy = mean(osg_uni_Stringy[1,])
q_uni_Stringy = quantile(osg_uni_Stringy[1,])


mq_uni_Stringy_output = apply(osg_uni_Stringy, 1, mq_uni_Stringy)

unlist(mq_uni_Stringy_output)

mq_uni_Stringy_output_final = matrix(unlist(mq_uni_Stringy_output),
                                      ncol = 3, by = TRUE)
dim(mq_uni_Stringy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_uni Skewed_mq.csv"))

mq_uni_Stringy_output_final = cbind(vec_ns, mq_uni_Stringy_output_final)
colnames(mq_uni_Stringy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_uni_Stringy_output_final) = vec_ns

mq_uni_Stringy_output_final

# 8 Visualisation of mean, q25 and q75: Stringy -----------------------------------

#Stringy
dev.off()
par(bg = NA)
plot(mq_uni_Stringy_output_final[,1], mq_uni_Stringy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Stringy", sub = "Uniform distribution")
lines(mq_uni_Stringy_output_final[,1], mq_uni_Stringy_output_final[,3],
      col = "dark green")
lines(mq_uni_Stringy_output_final[,1], mq_uni_Stringy_output_final[,4],
      col = "dark red")
dev.off()



# 9 Analysis of scagnostic measures: Monotonic ------------------------------------------------------------------

#preparing data set
osg_uni_Monotonic = t(output_scag_measure)
dim(osg_uni_Monotonic)
rownames(osg_uni_Monotonic) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_uni_Monotonic = function(x){

  mean_uni_Monotonic = mean(x)
  q_uni_Monotonic = quantile(x)
  mq_uni_Monotonic = as.data.frame(rbind(mean_uni_Monotonic, q_uni_Monotonic[2],
                                          q_uni_Stringy[4]))
  rownames(mq_uni_Monotonic) = c("mean", "25%", "75%")
  return (mq_uni_Monotonic)
}

#test
mq_uni_Monotonic(osg_uni_Monotonic[1,])

mean_uni_Monotonic = mean(osg_uni_Monotonic[1,])
q_uni_Monotonic = quantile(osg_uni_Monotonic[1,])


mq_uni_Monotonic_output = apply(osg_uni_Monotonic, 1, mq_uni_Monotonic)

unlist(mq_uni_Monotonic_output)

mq_uni_Monotonic_output_final = matrix(unlist(mq_uni_Monotonic_output),
                                        ncol = 3, by = TRUE)
dim(mq_uni_Monotonic_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_uni Skewed_mq.csv"))

mq_uni_Monotonic_output_final = cbind(vec_ns, mq_uni_Monotonic_output_final)
colnames(mq_uni_Monotonic_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_uni_Monotonic_output_final) = vec_ns

mq_uni_Monotonic_output_final

# 9 Visualisation of mean, q25 and q75: Monotonic --------------------------------

#Monotonic
dev.off()
par(bg = NA)
plot(mq_uni_Monotonic_output_final[,1], mq_uni_Monotonic_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Monotonic", sub = "Uniform distribution")
lines(mq_uni_Monotonic_output_final[,1], mq_uni_Monotonic_output_final[,3],
      col = "dark green")
lines(mq_uni_Monotonic_output_final[,1], mq_uni_Monotonic_output_final[,4],
      col = "dark red")
dev.off()
