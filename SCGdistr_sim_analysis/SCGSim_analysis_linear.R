rm(list = ls())
options(java.parameters = "-Xmx8000m")

library(scagnostics)
library(alphahull)

# Functions  --------------------------------------------------------------

#function to simulate liniar distribution
scag_lin = function(n){

  x = runif(n, 0, 10)
  b = rnorm(n, 1, 2)
  y = x + b

  df_lin = as.data.frame(cbind(x,y))
  s_lin = scagnostics(df_lin)

  lin_eval = rbind(s_lin)
}

#test = scag_lin(100)

#function to add to the same list the each simulation iteration
scag_simul = function (sim_settings_n, sim_no_rep){
  #sim_settings_n = 5
  #sim_no_rep = 3
  scag_simul_output = scag_lin(sim_settings_n)
  for (i in (1:(sim_no_rep-1))){
    #print(paste0("i",i, "scag_simul_output:",scag_simul_output ))
    added_sim_round = scag_lin(sim_settings_n)
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

column_index = 6
# Column legend: 1 = Outlying; 2 = Skewed; 3 = Clumpy; 4 = Sparse; 5 = Striate;
# 6 = Convex; 7 = Skinny; 8 = Stringy; 9 = Monotonic

for (page in 1:length(vec_ns)){
  working_page = simul_final_output_list[[page]]
  output_scag_measure[,page] =  working_page[,column_index]
}

output_scag_measure #matrix with all 52 columns with the assigned index "column_index"

write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_lin Skewed.csv"))


# 1 Analysis of scagnostic measures: Outlying----------------------------------------------

#Outlying measure

#preparing data set
osg_lin_Outlying = t(output_scag_measure)
dim(osg_lin_Outlying)
rownames(osg_lin_Outlying) = vec_ns



#Calculating the mean, q25 and q75
#i = 1
mq_lin_Outlying = function(x){

  mean_lin_Outlying = mean(x)
  q_lin_Outlying = quantile(x)
  mq_lin_Outlying = as.data.frame(rbind(mean_lin_Outlying, q_lin_Outlying[2],
                                        q_lin_Outlying[4]))
  rownames(mq_lin_Outlying) = c("mean", "25%", "75%")
  return (mq_lin_Outlying)
}

#test
mq_lin_Outlying(osg_lin_Outlying[1,])

mean_lin_Outlying = mean(osg_lin_Outlying[1,])
q_lin_Outlying = quantile(osg_lin_Outlying[1,])

mq_lin_Outlying_output = apply(osg_lin_Outlying, 1, mq_lin_Outlying)

unlist(mq_lin_Outlying_output)
#dim(unlist(mq_lin_Outlying_output))

mq_lin_Outlying_output_final = matrix(unlist(mq_lin_Outlying_output),
                                      ncol = 3, by = TRUE)
#dim(mq_lin_Outlying_output_final)

mq_lin_Outlying_output_final = cbind(vec_ns,mq_lin_Outlying_output_final)
colnames(mq_lin_Outlying_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_lin_Outlying_output_final) = vec_ns

mq_lin_Outlying_output_final


# 1 Visualisation of mean, q25 and q75: Outlying --------------------------------------

#Outlying
dev.off()
par(bg = NA)
plot(mq_lin_Outlying_output_final[,1], mq_lin_Outlying_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "blue",
     main = "Coefficient: Outlying", sub = "Linear distribution")
lines(mq_lin_Outlying_output_final[,1], mq_lin_Outlying_output_final[,3],
      col = "dark green")
lines(mq_lin_Outlying_output_final[,1], mq_lin_Outlying_output_final[,4],
      col = "dark red")
dev.off()




# 2 Analysis of scagnostic measures: Skewed ------------------------------------------------------------------

#preparing data set
osg_lin_Skewed = t(output_scag_measure)
dim(osg_lin_Skewed)
rownames(osg_lin_Skewed) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_lin_Skewed = function(x){

  mean_lin_Skewed = mean(x)
  q_lin_Skewed = quantile(x)
  mq_lin_Skewed = as.data.frame(rbind(mean_lin_Skewed, q_lin_Skewed[2],
                                      q_lin_Skewed[4]))
  rownames(mq_lin_Skewed) = c("mean", "25%", "75%")
  return (mq_lin_Skewed)
}

#test
mq_lin_Skewed(osg_lin_Skewed[1,])

mean_lin_Skewed = mean(osg_lin_Skewed[1,])
q_lin_Skewed = quantile(osg_lin_Skewed[1,])

mq_lin_Skewed_output = apply(osg_lin_Skewed, 1, mq_lin_Skewed)

unlist(mq_lin_Skewed_output)
#dim(unlist(mq_lin_Outlying_output))

mq_lin_Skewed_output_final = matrix(unlist(mq_lin_Skewed_output),
                                    ncol = 3, by = TRUE)
dim(mq_lin_Skewed_output_final)


mq_lin_Skewed_output_final = cbind(vec_ns,mq_lin_Skewed_output_final)
colnames(mq_lin_Skewed_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_lin_Skewed_output_final) = vec_ns

mq_lin_Skewed_output_final

# 2 Visualisation of mean, q25 and q75: Skewed -----------------------------------

#Skewed

dev.off()
par(bg = NA)
plot(mq_lin_Skewed_output_final[,1], mq_lin_Skewed_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skewed", sub = "Linear Distribution")
lines(mq_lin_Skewed_output_final[,1], mq_lin_Skewed_output_final[,3],
      col = "dark green")
lines(mq_lin_Skewed_output_final[,1], mq_lin_Skewed_output_final[,4],
      col = "dark red")
dev.off()




# 3 Analysis of scagnostic measures: Clumpy ------------------------------------------------------------------

#preparing data set
osg_lin_Clumpy = t(output_scag_measure)
dim(osg_lin_Clumpy)
rownames(osg_lin_Clumpy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_lin_Clumpy = function(x){

  mean_lin_Clumpy = mean(x)
  q_lin_Clumpy = quantile(x)
  mq_lin_Clumpy = as.data.frame(rbind(mean_lin_Clumpy, q_lin_Clumpy[2],
                                      q_lin_Clumpy[4]))
  rownames(mq_lin_Clumpy) = c("mean", "25%", "75%")
  return (mq_lin_Clumpy)
}

#test
mq_lin_Clumpy(osg_lin_Clumpy[1,])

mean_lin_Clumpy = mean(osg_lin_Clumpy[1,])
q_lin_Clumpy = quantile(osg_lin_Clumpy[1,])


mq_lin_Clumpy_output = apply(osg_lin_Clumpy, 1, mq_lin_Clumpy)

unlist(mq_lin_Clumpy_output)

mq_lin_Clumpy_output_final = matrix(unlist(mq_lin_Clumpy_output),
                                    ncol = 3, by = TRUE)
dim(mq_lin_Clumpy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_lin Skewed_mq.csv"))

mq_lin_Clumpy_output_final = cbind(vec_ns,mq_lin_Clumpy_output_final)
colnames(mq_lin_Clumpy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_lin_Clumpy_output_final) = vec_ns

mq_lin_Clumpy_output_final

# 3 Visualisation of mean, q25 and q75: Clumpy -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_lin_Clumpy_output_final[,1], mq_lin_Clumpy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Clumpy", sub = "Linear distribution")
lines(mq_lin_Clumpy_output_final[,1], mq_lin_Clumpy_output_final[,3],
      col = "dark green")
lines(mq_lin_Clumpy_output_final[,1], mq_lin_Clumpy_output_final[,4],
      col = "dark red")
dev.off()




# 4 Analysis of scagnostic measures: Sparse ------------------------------------------------------------------

#preparing data set
osg_lin_Sparse = t(output_scag_measure)
dim(osg_lin_Sparse)
rownames(osg_lin_Sparse) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_lin_Sparse = function(x){

  mean_lin_Sparse = mean(x)
  q_lin_Sparse = quantile(x)
  mq_lin_Sparse = as.data.frame(rbind(mean_lin_Sparse, q_lin_Sparse[2],
                                      q_lin_Sparse[4]))
  rownames(mq_lin_Sparse) = c("mean", "25%", "75%")
  return (mq_lin_Sparse)
}

#test
mq_lin_Sparse(osg_lin_Sparse[1,])

mean_lin_Sparse = mean(osg_lin_Sparse[1,])
q_lin_Sparse = quantile(osg_lin_Sparse[1,])


mq_lin_Sparse_output = apply(osg_lin_Sparse, 1, mq_lin_Sparse)

unlist(mq_lin_Sparse_output)

mq_lin_Sparse_output_final = matrix(unlist(mq_lin_Sparse_output),
                                    ncol = 3, by = TRUE)
dim(mq_lin_Sparse_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_lin Skewed_mq.csv"))

mq_lin_Sparse_output_final = cbind(vec_ns,mq_lin_Sparse_output_final)
colnames(mq_lin_Sparse_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_lin_Sparse_output_final) = vec_ns

mq_lin_Sparse_output_final

# 4 Visualisation of mean, q25 and q75: Sparse -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_lin_Sparse_output_final[,1], mq_lin_Sparse_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "25q, mean, q75", col = "dark blue",
     main = "Coefficient: Sparse", sub = "Linear distribution")
lines(mq_lin_Sparse_output_final[,1], mq_lin_Sparse_output_final[,3],
      col = "dark green")
lines(mq_lin_Sparse_output_final[,1], mq_lin_Sparse_output_final[,4],
      col = "dark red")
dev.off()


# 5 Analysis of scagnostic measures: Striate ------------------------------------------------------------------

#preparing data set
osg_lin_Striate = t(output_scag_measure)
dim(osg_lin_Striate)
rownames(osg_lin_Striate) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_lin_Striate = function(x){

  mean_lin_Striate = mean(x)
  q_lin_Striate = quantile(x)
  mq_lin_Striate = as.data.frame(rbind(mean_lin_Striate, q_lin_Striate[2],
                                       q_lin_Striate[4]))
  rownames(mq_lin_Striate) = c("mean", "25%", "75%")
  return (mq_lin_Striate)
}

#test
mq_lin_Striate(osg_lin_Striate[1,])

mean_lin_Striate = mean(osg_lin_Striate[1,])
q_lin_Striate = quantile(osg_lin_Striate[1,])


mq_lin_Striate_output = apply(osg_lin_Striate, 1, mq_lin_Striate)

unlist(mq_lin_Striate_output)

mq_lin_Striate_output_final = matrix(unlist(mq_lin_Striate_output),
                                     ncol = 3, by = TRUE)
dim(mq_lin_Striate_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_lin Skewed_mq.csv"))

mq_lin_Striate_output_final = cbind(vec_ns,mq_lin_Striate_output_final)
colnames(mq_lin_Striate_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_lin_Striate_output_final) = vec_ns

mq_lin_Striate_output_final

# 5 Visualisation of mean, q25 and q75: Striate -----------------------------------

#Striate
dev.off()
par(bg = NA)
plot(mq_lin_Striate_output_final[,1], mq_lin_Striate_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Striated", sub = "Linear distribution")
lines(mq_lin_Striate_output_final[,1], mq_lin_Striate_output_final[,3],
      col = "dark green")
lines(mq_lin_Striate_output_final[,1], mq_lin_Striate_output_final[,4],
      col = "dark red")
dev.off()



# 6 Analysis of scagnostic measures: Convex ------------------------------------------------------------------

#preparing data set
osg_lin_Convex = t(output_scag_measure)
dim(osg_lin_Convex)
rownames(osg_lin_Convex) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_lin_Convex = function(x){

  mean_lin_Convex = mean(x)
  q_lin_Convex = quantile(x)
  mq_lin_Convex = as.data.frame(rbind(mean_lin_Convex, q_lin_Convex[2],
                                      q_lin_Convex[4]))
  rownames(mq_lin_Convex) = c("mean", "25%", "75%")
  return (mq_lin_Convex)
}

#test
mq_lin_Convex(osg_lin_Convex[1,])

mean_lin_Convex = mean(osg_lin_Convex[1,])
q_lin_Convex = quantile(osg_lin_Convex[1,])


mq_lin_Convex_output = apply(osg_lin_Convex, 1, mq_lin_Convex)

unlist(mq_lin_Convex_output)

mq_lin_Convex_output_final = matrix(unlist(mq_lin_Convex_output),
                                    ncol = 3, by = TRUE)
dim(mq_lin_Convex_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_lin Skewed_mq.csv"))

mq_lin_Convex_output_final = cbind(vec_ns, mq_lin_Convex_output_final)
colnames(mq_lin_Convex_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_lin_Convex_output_final) = vec_ns

mq_lin_Convex_output_final

# 6 Visualisation of mean, q25 and q75: Convex -----------------------------------

#Convex
dev.off()
par(bg = NA)
plot(mq_lin_Convex_output_final[,1], mq_lin_Convex_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Convex", sub = "Linear distribution")
lines(mq_lin_Convex_output_final[,1], mq_lin_Convex_output_final[,3],
      col = "dark green")
lines(mq_lin_Convex_output_final[,1], mq_lin_Convex_output_final[,4],
      col = "dark red")
dev.off()


# 7 Analysis of scagnostic measures: Skinny ------------------------------------------------------------------

#preparing data set
osg_lin_Skinny = t(output_scag_measure)
dim(osg_lin_Skinny)
rownames(osg_lin_Skinny) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_lin_Skinny = function(x){

  mean_lin_Skinny = mean(x)
  q_lin_Skinny = quantile(x)
  mq_lin_Skinny = as.data.frame(rbind(mean_lin_Skinny, q_lin_Skinny[2],
                                      q_lin_Skinny[4]))
  rownames(mq_lin_Skinny) = c("mean", "25%", "75%")
  return (mq_lin_Skinny)
}

#test
mq_lin_Skinny(osg_lin_Skinny[1,])

mean_lin_Skinny = mean(osg_lin_Skinny[1,])
q_lin_Skinny = quantile(osg_lin_Skinny[1,])


mq_lin_Skinny_output = apply(osg_lin_Skinny, 1, mq_lin_Skinny)

unlist(mq_lin_Skinny_output)

mq_lin_Skinny_output_final = matrix(unlist(mq_lin_Skinny_output),
                                    ncol = 3, by = TRUE)
dim(mq_lin_Skinny_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_lin Skewed_mq.csv"))

mq_lin_Skinny_output_final = cbind(vec_ns, mq_lin_Skinny_output_final)
colnames(mq_lin_Skinny_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_lin_Skinny_output_final) = vec_ns

mq_lin_Skinny_output_final

# 7 Visualisation of mean, q25 and q75: Skinny -----------------------------------

#Skinny
dev.off()
par(bg = NA)
plot(mq_lin_Skinny_output_final[,1], mq_lin_Skinny_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skinny", sub = "Linear distribution")
lines(mq_lin_Skinny_output_final[,1], mq_lin_Skinny_output_final[,3],
      col = "dark green")
lines(mq_lin_Skinny_output_final[,1], mq_lin_Skinny_output_final[,4],
      col = "dark red")
dev.off()


# 8 Analysis of scagnostic measures: Stringy ------------------------------------------------------------------

#preparing data set
osg_lin_Stringy = t(output_scag_measure)
dim(osg_lin_Stringy)
rownames(osg_lin_Stringy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_lin_Stringy = function(x){

  mean_lin_Stringy = mean(x)
  q_lin_Stringy = quantile(x)
  mq_lin_Stringy = as.data.frame(rbind(mean_lin_Stringy, q_lin_Stringy[2],
                                       q_lin_Stringy[4]))
  rownames(mq_lin_Stringy) = c("mean", "25%", "75%")
  return (mq_lin_Stringy)
}

#test
mq_lin_Stringy(osg_lin_Stringy[1,])

mean_lin_Stringy = mean(osg_lin_Stringy[1,])
q_lin_Stringy = quantile(osg_lin_Stringy[1,])


mq_lin_Stringy_output = apply(osg_lin_Stringy, 1, mq_lin_Stringy)

unlist(mq_lin_Stringy_output)

mq_lin_Stringy_output_final = matrix(unlist(mq_lin_Stringy_output),
                                     ncol = 3, by = TRUE)
dim(mq_lin_Stringy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_lin Skewed_mq.csv"))

mq_lin_Stringy_output_final = cbind(vec_ns, mq_lin_Stringy_output_final)
colnames(mq_lin_Stringy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_lin_Stringy_output_final) = vec_ns

mq_lin_Stringy_output_final

# 8 Visualisation of mean, q25 and q75: Stringy -----------------------------------

#Stringy
dev.off()
par(bg = NA)
plot(mq_lin_Stringy_output_final[,1], mq_lin_Stringy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Stringy", sub = "Linear distribution")
lines(mq_lin_Stringy_output_final[,1], mq_lin_Stringy_output_final[,3],
      col = "dark green")
lines(mq_lin_Stringy_output_final[,1], mq_lin_Stringy_output_final[,4],
      col = "dark red")
dev.off()



# 9 Analysis of scagnostic measures: Monotonic ------------------------------------------------------------------

#preparing data set
osg_lin_Monotonic = t(output_scag_measure)
dim(osg_lin_Monotonic)
rownames(osg_lin_Monotonic) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_lin_Monotonic = function(x){

  mean_lin_Monotonic = mean(x)
  q_lin_Monotonic = quantile(x)
  mq_lin_Monotonic = as.data.frame(rbind(mean_lin_Monotonic, q_lin_Monotonic[2],
                                         q_lin_Monotonic[4]))
  rownames(mq_lin_Monotonic) = c("mean", "25%", "75%")
  return (mq_lin_Monotonic)
}

#test
mq_lin_Monotonic(osg_lin_Monotonic[1,])

mean_lin_Monotonic = mean(osg_lin_Monotonic[1,])
q_lin_Monotonic = quantile(osg_lin_Monotonic[1,])


mq_lin_Monotonic_output = apply(osg_lin_Monotonic, 1, mq_lin_Monotonic)

unlist(mq_lin_Monotonic_output)

mq_lin_Monotonic_output_final = matrix(unlist(mq_lin_Monotonic_output),
                                       ncol = 3, by = TRUE)
dim(mq_lin_Monotonic_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_lin Skewed_mq.csv"))

mq_lin_Monotonic_output_final = cbind(vec_ns, mq_lin_Monotonic_output_final)
colnames(mq_lin_Monotonic_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_lin_Monotonic_output_final) = vec_ns

mq_lin_Monotonic_output_final

# 9 Visualisation of mean, q25 and q75: Monotonic --------------------------------

#Monotonic
dev.off()
par(bg = NA)
plot(mq_lin_Monotonic_output_final[,1], mq_lin_Monotonic_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Monotonic", sub = "Linear distribution")
lines(mq_lin_Monotonic_output_final[,1], mq_lin_Monotonic_output_final[,3],
      col = "dark green")
lines(mq_lin_Monotonic_output_final[,1], mq_lin_Monotonic_output_final[,4],
      col = "dark red")
dev.off()
