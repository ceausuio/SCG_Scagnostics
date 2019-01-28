rm(list = ls())
options(java.parameters = "-Xmx8000m")

library(scagnostics)
library(alphahull)

# Functions  --------------------------------------------------------------

#function to simulate Arima distribution
scag_ari = function(n){

  N = arima.sim(list(order=c(1, 0, 0), ar=0.93), n)

  df_ari = (cbind(1:n,N))
  s_ari = scagnostics((df_ari))

  ari_eval = rbind(s_ari)
}

#test = scag_ari(100)

#function to add to the same list the each simulation iteration
scag_simul = function (sim_settings_n, sim_no_rep){
  #sim_settings_n = 5
  #sim_no_rep = 3
  scag_simul_output = scag_ari(sim_settings_n)
  for (i in (1:(sim_no_rep-1))){
    #print(paste0("i",i, "scag_simul_output:",scag_simul_output ))
    added_sim_round = scag_ari(sim_settings_n)
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

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_ari Skewed.csv"))


# 1 Analysis of scagnostic measures: Outlying----------------------------------------------

#Outlying measure

#preparing data set
osg_ari_Outlying = t(output_scag_measure)
dim(osg_ari_Outlying)
rownames(osg_ari_Outlying) = vec_ns



#Calculating the mean, q25 and q75
#i = 1
mq_ari_Outlying = function(x){

  mean_ari_Outlying = mean(x)
  q_ari_Outlying = quantile(x)
  mq_ari_Outlying = as.data.frame(rbind(mean_ari_Outlying, q_ari_Outlying[2],
                                        q_ari_Outlying[4]))
  rownames(mq_ari_Outlying) = c("mean", "25%", "75%")
  return (mq_ari_Outlying)
}

#test
mq_ari_Outlying(osg_ari_Outlying[1,])

mean_ari_Outlying = mean(osg_ari_Outlying[1,])
q_ari_Outlying = quantile(osg_ari_Outlying[1,])

mq_ari_Outlying_output = apply(osg_ari_Outlying, 1, mq_ari_Outlying)

unlist(mq_ari_Outlying_output)
#dim(unlist(mq_ari_Outlying_output))

mq_ari_Outlying_output_final = matrix(unlist(mq_ari_Outlying_output),
                                      ncol = 3, by = TRUE)
#dim(mq_ari_Outlying_output_final)

mq_ari_Outlying_output_final = cbind(vec_ns,mq_ari_Outlying_output_final)
colnames(mq_ari_Outlying_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_ari_Outlying_output_final) = vec_ns

mq_ari_Outlying_output_final


# 1 Visualisation of mean, q25 and q75: Outlying --------------------------------------

#Outlying
dev.off()
par(bg = NA)
plot(mq_ari_Outlying_output_final[,1], mq_ari_Outlying_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "blue",
     main = "Coefficient: Outlying", sub = "Arima distribution")
lines(mq_ari_Outlying_output_final[,1], mq_ari_Outlying_output_final[,3],
      col = "dark green")
lines(mq_ari_Outlying_output_final[,1], mq_ari_Outlying_output_final[,4],
      col = "dark red")
dev.off()




# 2 Analysis of scagnostic measures: Skewed ------------------------------------------------------------------

#preparing data set
osg_ari_Skewed = t(output_scag_measure)
dim(osg_ari_Skewed)
rownames(osg_ari_Skewed) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_ari_Skewed = function(x){

  mean_ari_Skewed = mean(x)
  q_ari_Skewed = quantile(x)
  mq_ari_Skewed = as.data.frame(rbind(mean_ari_Skewed, q_ari_Skewed[2],
                                      q_ari_Skewed[4]))
  rownames(mq_ari_Skewed) = c("mean", "25%", "75%")
  return (mq_ari_Skewed)
}

#test
mq_ari_Skewed(osg_ari_Skewed[1,])

mean_ari_Skewed = mean(osg_ari_Skewed[1,])
q_ari_Skewed = quantile(osg_ari_Skewed[1,])

mq_ari_Skewed_output = apply(osg_ari_Skewed, 1, mq_ari_Skewed)

unlist(mq_ari_Skewed_output)
#dim(unlist(mq_ari_Outlying_output))

mq_ari_Skewed_output_final = matrix(unlist(mq_ari_Skewed_output),
                                    ncol = 3, by = TRUE)
dim(mq_ari_Skewed_output_final)


mq_ari_Skewed_output_final = cbind(vec_ns,mq_ari_Skewed_output_final)
colnames(mq_ari_Skewed_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_ari_Skewed_output_final) = vec_ns

mq_ari_Skewed_output_final

# 2 Visualisation of mean, q25 and q75: Skewed -----------------------------------

#Skewed

dev.off()
par(bg = NA)
plot(mq_ari_Skewed_output_final[,1], mq_ari_Skewed_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skewed", sub = "Arima Distribution")
lines(mq_ari_Skewed_output_final[,1], mq_ari_Skewed_output_final[,3],
      col = "dark green")
lines(mq_ari_Skewed_output_final[,1], mq_ari_Skewed_output_final[,4],
      col = "dark red")
dev.off()




# 3 Analysis of scagnostic measures: Clumpy ------------------------------------------------------------------

#preparing data set
osg_ari_clumpy = t(output_scag_measure)
dim(osg_ari_clumpy)
rownames(osg_ari_clumpy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_ari_clumpy = function(x){

  mean_ari_clumpy = mean(x)
  q_ari_clumpy = quantile(x)
  mq_ari_clumpy = as.data.frame(rbind(mean_ari_clumpy, q_ari_clumpy[2],
                                      q_ari_clumpy[4]))
  rownames(mq_ari_clumpy) = c("mean", "25%", "75%")
  return (mq_ari_clumpy)
}

#test
mq_ari_clumpy(osg_ari_clumpy[1,])

mean_ari_clumpy = mean(osg_ari_clumpy[1,])
q_ari_clumpy = quantile(osg_ari_clumpy[1,])


mq_ari_clumpy_output = apply(osg_ari_clumpy, 1, mq_ari_clumpy)

unlist(mq_ari_clumpy_output)

mq_ari_clumpy_output_final = matrix(unlist(mq_ari_clumpy_output),
                                    ncol = 3, by = TRUE)
dim(mq_ari_clumpy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_ari Skewed_mq.csv"))

mq_ari_clumpy_output_final = cbind(vec_ns,mq_ari_clumpy_output_final)
colnames(mq_ari_clumpy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_ari_clumpy_output_final) = vec_ns

mq_ari_clumpy_output_final

# 3 Visualisation of mean, q25 and q75: clumpy -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_ari_clumpy_output_final[,1], mq_ari_clumpy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Clumpy", sub = "Arima distribution")
lines(mq_ari_clumpy_output_final[,1], mq_ari_clumpy_output_final[,3],
      col = "dark green")
lines(mq_ari_clumpy_output_final[,1], mq_ari_clumpy_output_final[,4],
      col = "dark red")
dev.off()




# 4 Analysis of scagnostic measures: Sparse ------------------------------------------------------------------

#preparing data set
osg_ari_Sparse = t(output_scag_measure)
dim(osg_ari_Sparse)
rownames(osg_ari_Sparse) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_ari_Sparse = function(x){

  mean_ari_Sparse = mean(x)
  q_ari_Sparse = quantile(x)
  mq_ari_Sparse = as.data.frame(rbind(mean_ari_Sparse, q_ari_Sparse[2],
                                      q_ari_Sparse[4]))
  rownames(mq_ari_Sparse) = c("mean", "25%", "75%")
  return (mq_ari_Sparse)
}

#test
mq_ari_Sparse(osg_ari_Sparse[1,])

mean_ari_Sparse = mean(osg_ari_Sparse[1,])
q_ari_Sparse = quantile(osg_ari_Sparse[1,])


mq_ari_Sparse_output = apply(osg_ari_Sparse, 1, mq_ari_Sparse)

unlist(mq_ari_Sparse_output)

mq_ari_Sparse_output_final = matrix(unlist(mq_ari_Sparse_output),
                                    ncol = 3, by = TRUE)
dim(mq_ari_Sparse_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_ari Skewed_mq.csv"))

mq_ari_Sparse_output_final = cbind(vec_ns,mq_ari_Sparse_output_final)
colnames(mq_ari_Sparse_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_ari_Sparse_output_final) = vec_ns

mq_ari_Sparse_output_final

# 4 Visualisation of mean, q25 and q75: Sparse -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_ari_Sparse_output_final[,1], mq_ari_Sparse_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "25q, mean, q75", col = "dark blue",
     main = "Coefficient: Sparse", sub = "Arima distribution")
lines(mq_ari_Sparse_output_final[,1], mq_ari_Sparse_output_final[,3],
      col = "dark green")
lines(mq_ari_Sparse_output_final[,1], mq_ari_Sparse_output_final[,4],
      col = "dark red")
dev.off()


# 5 Analysis of scagnostic measures: Striate ------------------------------------------------------------------

#preparing data set
osg_ari_Striate = t(output_scag_measure)
dim(osg_ari_Striate)
rownames(osg_ari_Striate) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_ari_Striate = function(x){

  mean_ari_Striate = mean(x)
  q_ari_Striate = quantile(x)
  mq_ari_Striate = as.data.frame(rbind(mean_ari_Striate, q_ari_Striate[2],
                                       q_ari_Striate[4]))
  rownames(mq_ari_Striate) = c("mean", "25%", "75%")
  return (mq_ari_Striate)
}

#test
mq_ari_Striate(osg_ari_Striate[1,])

mean_ari_Striate = mean(osg_ari_Striate[1,])
q_ari_Striate = quantile(osg_ari_Striate[1,])


mq_ari_Striate_output = apply(osg_ari_Striate, 1, mq_ari_Striate)

unlist(mq_ari_Striate_output)

mq_ari_Striate_output_final = matrix(unlist(mq_ari_Striate_output),
                                     ncol = 3, by = TRUE)
dim(mq_ari_Striate_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_ari Skewed_mq.csv"))

mq_ari_Striate_output_final = cbind(vec_ns,mq_ari_Striate_output_final)
colnames(mq_ari_Striate_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_ari_Striate_output_final) = vec_ns

mq_ari_Striate_output_final

# 5 Visualisation of mean, q25 and q75: Striate -----------------------------------

#Striate
dev.off()
par(bg = NA)
plot(mq_ari_Striate_output_final[,1], mq_ari_Striate_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Striated", sub = "Arima distribution")
lines(mq_ari_Striate_output_final[,1], mq_ari_Striate_output_final[,3],
      col = "dark green")
lines(mq_ari_Striate_output_final[,1], mq_ari_Striate_output_final[,4],
      col = "dark red")
dev.off()



# 6 Analysis of scagnostic measures: Convex ------------------------------------------------------------------

#preparing data set
osg_ari_Convex = t(output_scag_measure)
dim(osg_ari_Convex)
rownames(osg_ari_Convex) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_ari_Convex = function(x){

  mean_ari_Convex = mean(x)
  q_ari_Convex = quantile(x)
  mq_ari_Convex = as.data.frame(rbind(mean_ari_Convex, q_ari_Convex[2],
                                      q_ari_Convex[4]))
  rownames(mq_ari_Convex) = c("mean", "25%", "75%")
  return (mq_ari_Convex)
}

#test
mq_ari_Convex(osg_ari_Convex[1,])

mean_ari_Convex = mean(osg_ari_Convex[1,])
q_ari_Convex = quantile(osg_ari_Convex[1,])


mq_ari_Convex_output = apply(osg_ari_Convex, 1, mq_ari_Convex)

unlist(mq_ari_Convex_output)

mq_ari_Convex_output_final = matrix(unlist(mq_ari_Convex_output),
                                    ncol = 3, by = TRUE)
dim(mq_ari_Convex_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_ari Skewed_mq.csv"))

mq_ari_Convex_output_final = cbind(vec_ns, mq_ari_Convex_output_final)
colnames(mq_ari_Convex_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_ari_Convex_output_final) = vec_ns

mq_ari_Convex_output_final

# 6 Visualisation of mean, q25 and q75: Convex -----------------------------------

#Convex
dev.off()
par(bg = NA)
plot(mq_ari_Convex_output_final[,1], mq_ari_Convex_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Convex", sub = "Arima distribution")
lines(mq_ari_Convex_output_final[,1], mq_ari_Convex_output_final[,3],
      col = "dark green")
lines(mq_ari_Convex_output_final[,1], mq_ari_Convex_output_final[,4],
      col = "dark red")
dev.off()


# 7 Analysis of scagnostic measures: Skinny ------------------------------------------------------------------

#preparing data set
osg_ari_Skinny = t(output_scag_measure)
dim(osg_ari_Skinny)
rownames(osg_ari_Skinny) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_ari_Skinny = function(x){

  mean_ari_Skinny = mean(x)
  q_ari_Skinny = quantile(x)
  mq_ari_Skinny = as.data.frame(rbind(mean_ari_Skinny, q_ari_Skinny[2],
                                      q_ari_Skinny[4]))
  rownames(mq_ari_Skinny) = c("mean", "25%", "75%")
  return (mq_ari_Skinny)
}

#test
mq_ari_Skinny(osg_ari_Skinny[1,])

mean_ari_Skinny = mean(osg_ari_Skinny[1,])
q_ari_Skinny = quantile(osg_ari_Skinny[1,])


mq_ari_Skinny_output = apply(osg_ari_Skinny, 1, mq_ari_Skinny)

unlist(mq_ari_Skinny_output)

mq_ari_Skinny_output_final = matrix(unlist(mq_ari_Skinny_output),
                                    ncol = 3, by = TRUE)
dim(mq_ari_Skinny_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_ari Skewed_mq.csv"))

mq_ari_Skinny_output_final = cbind(vec_ns, mq_ari_Skinny_output_final)
colnames(mq_ari_Skinny_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_ari_Skinny_output_final) = vec_ns

mq_ari_Skinny_output_final

# 7 Visualisation of mean, q25 and q75: Skinny -----------------------------------

#Skinny
dev.off()
par(bg = NA)
plot(mq_ari_Skinny_output_final[,1], mq_ari_Skinny_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skinny", sub = "Arima distribution")
lines(mq_ari_Skinny_output_final[,1], mq_ari_Skinny_output_final[,3],
      col = "dark green")
lines(mq_ari_Skinny_output_final[,1], mq_ari_Skinny_output_final[,4],
      col = "dark red")
dev.off()


# 8 Analysis of scagnostic measures: Stringy ------------------------------------------------------------------

#preparing data set
osg_ari_Stringy = t(output_scag_measure)
dim(osg_ari_Stringy)
rownames(osg_ari_Stringy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_ari_Stringy = function(x){

  mean_ari_Stringy = mean(x)
  q_ari_Stringy = quantile(x)
  mq_ari_Stringy = as.data.frame(rbind(mean_ari_Stringy, q_ari_Stringy[2],
                                       q_ari_Stringy[4]))
  rownames(mq_ari_Stringy) = c("mean", "25%", "75%")
  return (mq_ari_Stringy)
}

#test
mq_ari_Stringy(osg_ari_Stringy[1,])

mean_ari_Stringy = mean(osg_ari_Stringy[1,])
q_ari_Stringy = quantile(osg_ari_Stringy[1,])


mq_ari_Stringy_output = apply(osg_ari_Stringy, 1, mq_ari_Stringy)

unlist(mq_ari_Stringy_output)

mq_ari_Stringy_output_final = matrix(unlist(mq_ari_Stringy_output),
                                     ncol = 3, by = TRUE)
dim(mq_ari_Stringy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_ari Skewed_mq.csv"))

mq_ari_Stringy_output_final = cbind(vec_ns, mq_ari_Stringy_output_final)
colnames(mq_ari_Stringy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_ari_Stringy_output_final) = vec_ns

mq_ari_Stringy_output_final

# 8 Visualisation of mean, q25 and q75: Stringy -----------------------------------

#Stringy
dev.off()
par(bg = NA)
plot(mq_ari_Stringy_output_final[,1], mq_ari_Stringy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Stringy", sub = "Arima distribution")
lines(mq_ari_Stringy_output_final[,1], mq_ari_Stringy_output_final[,3],
      col = "dark green")
lines(mq_ari_Stringy_output_final[,1], mq_ari_Stringy_output_final[,4],
      col = "dark red")
dev.off()



# 9 Analysis of scagnostic measures: Monotonic ------------------------------------------------------------------

#preparing data set
osg_ari_Monotonic = t(output_scag_measure)
dim(osg_ari_Monotonic)
rownames(osg_ari_Monotonic) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_ari_Monotonic = function(x){

  mean_ari_Monotonic = mean(x)
  q_ari_Monotonic = quantile(x)
  mq_ari_Monotonic = as.data.frame(rbind(mean_ari_Monotonic, q_ari_Monotonic[2],
                                         q_ari_Monotonic[4]))
  rownames(mq_ari_Monotonic) = c("mean", "25%", "75%")
  return (mq_ari_Monotonic)
}

#test
mq_ari_Monotonic(osg_ari_Monotonic[1,])

mean_ari_Monotonic = mean(osg_ari_Monotonic[1,])
q_ari_Monotonic = quantile(osg_ari_Monotonic[1,])


mq_ari_Monotonic_output = apply(osg_ari_Monotonic, 1, mq_ari_Monotonic)

unlist(mq_ari_Monotonic_output)

mq_ari_Monotonic_output_final = matrix(unlist(mq_ari_Monotonic_output),
                                       ncol = 3, by = TRUE)
dim(mq_ari_Monotonic_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_ari Skewed_mq.csv"))

mq_ari_Monotonic_output_final = cbind(vec_ns, mq_ari_Monotonic_output_final)
colnames(mq_ari_Monotonic_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_ari_Monotonic_output_final) = vec_ns

mq_ari_Monotonic_output_final

# 9 Visualisation of mean, q25 and q75: Monotonic --------------------------------

#Monotonic
dev.off()
par(bg = NA)
plot(mq_ari_Monotonic_output_final[,1], mq_ari_Monotonic_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Monotonic", sub = "Arima distribution")
lines(mq_ari_Monotonic_output_final[,1], mq_ari_Monotonic_output_final[,3],
      col = "dark green")
lines(mq_ari_Monotonic_output_final[,1], mq_ari_Monotonic_output_final[,4],
      col = "dark red")
dev.off()
