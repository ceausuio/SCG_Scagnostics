rm(list = ls())
options(java.parameters = "-Xmx8000m")

library(scagnostics)
library(alphahull)

# Functions  --------------------------------------------------------------

#function to simulate Eye distribution
scag_eye = function(n){

  theta = runif(n,0,2*pi)
  r = sqrt(runif(n,0.25^2,0.5^2))
  x = cbind(0.5+r*cos(theta),0.5+r*sin(theta))
  # Value of alpha
  alpha = 0
  # alpha-shape
  ashape.obj = ashape(x, alpha = alpha)
  # Plot alpha-shape in blue, sample points in black,
  # and Delaunay triangulation in red


  df = as.data.frame(x)
  set.seed(1)
  # n = 500
  a = rnorm(n, 0.5, 0.05)
  b = rnorm(n, 0.5, 0.05)
  c = cbind(a, b)


  df_eye = rbind(x, c)
  s_eye = scagnostics(df_eye)


  eye_eval = rbind(s_eye)
}

test = scag_eye(100)

#function to add to the same list the each simulation iteration
scag_simul = function (sim_settings_n, sim_no_rep){
  #sim_settings_n = 5
  #sim_no_rep = 3
  scag_simul_output = scag_eye(sim_settings_n)
  for (i in (1:(sim_no_rep-1))){
    #print(paste0("i",i, "scag_simul_output:",scag_simul_output ))
    added_sim_round = scag_eye(sim_settings_n)
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

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_eye Skewed.csv"))


# 1 Analysis of scagnostic measures: Outlying----------------------------------------------

#Outlying measure

#preparing data set
osg_eye_Outlying = t(output_scag_measure)
dim(osg_eye_Outlying)
rownames(osg_eye_Outlying) = vec_ns



#Calculating the mean, q25 and q75
#i = 1
mq_eye_Outlying = function(x){

  mean_eye_Outlying = mean(x)
  q_eye_Outlying = quantile(x)
  mq_eye_Outlying = as.data.frame(rbind(mean_eye_Outlying, q_eye_Outlying[2],
                                        q_eye_Outlying[4]))
  rownames(mq_eye_Outlying) = c("mean", "25%", "75%")
  return (mq_eye_Outlying)
}

#test
mq_eye_Outlying(osg_eye_Outlying[1,])

mean_eye_Outlying = mean(osg_eye_Outlying[1,])
q_eye_Outlying = quantile(osg_eye_Outlying[1,])

mq_eye_Outlying_output = apply(osg_eye_Outlying, 1, mq_eye_Outlying)

unlist(mq_eye_Outlying_output)
#dim(unlist(mq_eye_Outlying_output))

mq_eye_Outlying_output_final = matrix(unlist(mq_eye_Outlying_output),
                                      ncol = 3, by = TRUE)
#dim(mq_eye_Outlying_output_final)

mq_eye_Outlying_output_final = cbind(vec_ns,mq_eye_Outlying_output_final)
colnames(mq_eye_Outlying_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_eye_Outlying_output_final) = vec_ns

mq_eye_Outlying_output_final


# 1 Visualisation of mean, q25 and q75: Outlying --------------------------------------

#Outlying
dev.off()
par(bg = NA)
plot(mq_eye_Outlying_output_final[,1], mq_eye_Outlying_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "blue",
     main = "Coefficient: Outlying", sub = "Eye distribution")
lines(mq_eye_Outlying_output_final[,1], mq_eye_Outlying_output_final[,3],
      col = "dark green")
lines(mq_eye_Outlying_output_final[,1], mq_eye_Outlying_output_final[,4],
      col = "dark red")
dev.off()




# 2 Analysis of scagnostic measures: Skewed ------------------------------------------------------------------

#preparing data set
osg_eye_Skewed = t(output_scag_measure)
dim(osg_eye_Skewed)
rownames(osg_eye_Skewed) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_eye_Skewed = function(x){

  mean_eye_Skewed = mean(x)
  q_eye_Skewed = quantile(x)
  mq_eye_Skewed = as.data.frame(rbind(mean_eye_Skewed, q_eye_Skewed[2],
                                      q_eye_Skewed[4]))
  rownames(mq_eye_Skewed) = c("mean", "25%", "75%")
  return (mq_eye_Skewed)
}

#test
mq_eye_Skewed(osg_eye_Skewed[1,])

mean_eye_Skewed = mean(osg_eye_Skewed[1,])
q_eye_Skewed = quantile(osg_eye_Skewed[1,])

mq_eye_Skewed_output = apply(osg_eye_Skewed, 1, mq_eye_Skewed)

unlist(mq_eye_Skewed_output)
#dim(unlist(mq_eye_Outlying_output))

mq_eye_Skewed_output_final = matrix(unlist(mq_eye_Skewed_output),
                                    ncol = 3, by = TRUE)
dim(mq_eye_Skewed_output_final)


mq_eye_Skewed_output_final = cbind(vec_ns,mq_eye_Skewed_output_final)
colnames(mq_eye_Skewed_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_eye_Skewed_output_final) = vec_ns

mq_eye_Skewed_output_final

# 2 Visualisation of mean, q25 and q75: Skewed -----------------------------------

#Skewed

dev.off()
par(bg = NA)
plot(mq_eye_Skewed_output_final[,1], mq_eye_Skewed_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skewed", sub = "Eye Distribution")
lines(mq_eye_Skewed_output_final[,1], mq_eye_Skewed_output_final[,3],
      col = "dark green")
lines(mq_eye_Skewed_output_final[,1], mq_eye_Skewed_output_final[,4],
      col = "dark red")
dev.off()




# 3 Analysis of scagnostic measures: Clumpy ------------------------------------------------------------------

#preparing data set
osg_eye_clumpy = t(output_scag_measure)
dim(osg_eye_clumpy)
rownames(osg_eye_clumpy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_eye_clumpy = function(x){

  mean_eye_clumpy = mean(x)
  q_eye_clumpy = quantile(x)
  mq_eye_clumpy = as.data.frame(rbind(mean_eye_clumpy, q_eye_clumpy[2],
                                      q_eye_clumpy[4]))
  rownames(mq_eye_clumpy) = c("mean", "25%", "75%")
  return (mq_eye_clumpy)
}

#test
mq_eye_clumpy(osg_eye_clumpy[1,])

mean_eye_clumpy = mean(osg_eye_clumpy[1,])
q_eye_clumpy = quantile(osg_eye_clumpy[1,])


mq_eye_clumpy_output = apply(osg_eye_clumpy, 1, mq_eye_clumpy)

unlist(mq_eye_clumpy_output)

mq_eye_clumpy_output_final = matrix(unlist(mq_eye_clumpy_output),
                                    ncol = 3, by = TRUE)
dim(mq_eye_clumpy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_eye Skewed_mq.csv"))

mq_eye_clumpy_output_final = cbind(vec_ns,mq_eye_clumpy_output_final)
colnames(mq_eye_clumpy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_eye_clumpy_output_final) = vec_ns

mq_eye_clumpy_output_final

# 3 Visualisation of mean, q25 and q75: clumpy -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_eye_clumpy_output_final[,1], mq_eye_clumpy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: clumpy", sub = "Eye distribution")
lines(mq_eye_clumpy_output_final[,1], mq_eye_clumpy_output_final[,3],
      col = "dark green")
lines(mq_eye_clumpy_output_final[,1], mq_eye_clumpy_output_final[,4],
      col = "dark red")
dev.off()




# 4 Analysis of scagnostic measures: Sparse ------------------------------------------------------------------

#preparing data set
osg_eye_Sparse = t(output_scag_measure)
dim(osg_eye_Sparse)
rownames(osg_eye_Sparse) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_eye_Sparse = function(x){

  mean_eye_Sparse = mean(x)
  q_eye_Sparse = quantile(x)
  mq_eye_Sparse = as.data.frame(rbind(mean_eye_Sparse, q_eye_Sparse[2],
                                      q_eye_Sparse[4]))
  rownames(mq_eye_Sparse) = c("mean", "25%", "75%")
  return (mq_eye_Sparse)
}

#test
mq_eye_Sparse(osg_eye_Sparse[1,])

mean_eye_Sparse = mean(osg_eye_Sparse[1,])
q_eye_Sparse = quantile(osg_eye_Sparse[1,])


mq_eye_Sparse_output = apply(osg_eye_Sparse, 1, mq_eye_Sparse)

unlist(mq_eye_Sparse_output)

mq_eye_Sparse_output_final = matrix(unlist(mq_eye_Sparse_output),
                                    ncol = 3, by = TRUE)
dim(mq_eye_Sparse_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_eye Skewed_mq.csv"))

mq_eye_Sparse_output_final = cbind(vec_ns,mq_eye_Sparse_output_final)
colnames(mq_eye_Sparse_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_eye_Sparse_output_final) = vec_ns

mq_eye_Sparse_output_final

# 4 Visualisation of mean, q25 and q75: Sparse -----------------------------------

#Skewed
dev.off()
par(bg = NA)
plot(mq_eye_Sparse_output_final[,1], mq_eye_Sparse_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "25q, mean, q75", col = "dark blue",
     main = "Coefficient: Sparse", sub = "Eye distribution")
lines(mq_eye_Sparse_output_final[,1], mq_eye_Sparse_output_final[,3],
      col = "dark green")
lines(mq_eye_Sparse_output_final[,1], mq_eye_Sparse_output_final[,4],
      col = "dark red")
dev.off()


# 5 Analysis of scagnostic measures: Striate ------------------------------------------------------------------

#preparing data set
osg_eye_Striate = t(output_scag_measure)
dim(osg_eye_Striate)
rownames(osg_eye_Striate) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_eye_Striate = function(x){

  mean_eye_Striate = mean(x)
  q_eye_Striate = quantile(x)
  mq_eye_Striate = as.data.frame(rbind(mean_eye_Striate, q_eye_Striate[2],
                                       q_eye_Striate[4]))
  rownames(mq_eye_Striate) = c("mean", "25%", "75%")
  return (mq_eye_Striate)
}

#test
mq_eye_Striate(osg_eye_Striate[1,])

mean_eye_Striate = mean(osg_eye_Striate[1,])
q_eye_Striate = quantile(osg_eye_Striate[1,])


mq_eye_Striate_output = apply(osg_eye_Striate, 1, mq_eye_Striate)

unlist(mq_eye_Striate_output)

mq_eye_Striate_output_final = matrix(unlist(mq_eye_Striate_output),
                                     ncol = 3, by = TRUE)
dim(mq_eye_Striate_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_eye Skewed_mq.csv"))

mq_eye_Striate_output_final = cbind(vec_ns,mq_eye_Striate_output_final)
colnames(mq_eye_Striate_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_eye_Striate_output_final) = vec_ns

mq_eye_Striate_output_final

# 5 Visualisation of mean, q25 and q75: Striate -----------------------------------

#Striate
dev.off()
par(bg = NA)
plot(mq_eye_Striate_output_final[,1], mq_eye_Striate_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Striated", sub = "Eye distribution")
lines(mq_eye_Striate_output_final[,1], mq_eye_Striate_output_final[,3],
      col = "dark green")
lines(mq_eye_Striate_output_final[,1], mq_eye_Striate_output_final[,4],
      col = "dark red")
dev.off()



# 6 Analysis of scagnostic measures: Convex ------------------------------------------------------------------

#preparing data set
osg_eye_Convex = t(output_scag_measure)
dim(osg_eye_Convex)
rownames(osg_eye_Convex) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_eye_Convex = function(x){

  mean_eye_Convex = mean(x)
  q_eye_Convex = quantile(x)
  mq_eye_Convex = as.data.frame(rbind(mean_eye_Convex, q_eye_Convex[2],
                                      q_eye_Convex[4]))
  rownames(mq_eye_Convex) = c("mean", "25%", "75%")
  return (mq_eye_Convex)
}

#test
mq_eye_Convex(osg_eye_Convex[1,])

mean_eye_Convex = mean(osg_eye_Convex[1,])
q_eye_Convex = quantile(osg_eye_Convex[1,])


mq_eye_Convex_output = apply(osg_eye_Convex, 1, mq_eye_Convex)

unlist(mq_eye_Convex_output)

mq_eye_Convex_output_final = matrix(unlist(mq_eye_Convex_output),
                                    ncol = 3, by = TRUE)
dim(mq_eye_Convex_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_eye Skewed_mq.csv"))

mq_eye_Convex_output_final = cbind(vec_ns, mq_eye_Convex_output_final)
colnames(mq_eye_Convex_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_eye_Convex_output_final) = vec_ns

mq_eye_Convex_output_final

# 6 Visualisation of mean, q25 and q75: Convex -----------------------------------

#Convex
dev.off()
par(bg = NA)
plot(mq_eye_Convex_output_final[,1], mq_eye_Convex_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Convex", sub = "Eye distribution")
lines(mq_eye_Convex_output_final[,1], mq_eye_Convex_output_final[,3],
      col = "dark green")
lines(mq_eye_Convex_output_final[,1], mq_eye_Convex_output_final[,4],
      col = "dark red")
dev.off()


# 7 Analysis of scagnostic measures: Skinny ------------------------------------------------------------------

#preparing data set
osg_eye_Skinny = t(output_scag_measure)
dim(osg_eye_Skinny)
rownames(osg_eye_Skinny) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_eye_Skinny = function(x){

  mean_eye_Skinny = mean(x)
  q_eye_Skinny = quantile(x)
  mq_eye_Skinny = as.data.frame(rbind(mean_eye_Skinny, q_eye_Skinny[2],
                                      q_eye_Skinny[4]))
  rownames(mq_eye_Skinny) = c("mean", "25%", "75%")
  return (mq_eye_Skinny)
}

#test
mq_eye_Skinny(osg_eye_Skinny[1,])

mean_eye_Skinny = mean(osg_eye_Skinny[1,])
q_eye_Skinny = quantile(osg_eye_Skinny[1,])


mq_eye_Skinny_output = apply(osg_eye_Skinny, 1, mq_eye_Skinny)

unlist(mq_eye_Skinny_output)

mq_eye_Skinny_output_final = matrix(unlist(mq_eye_Skinny_output),
                                    ncol = 3, by = TRUE)
dim(mq_eye_Skinny_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_eye Skewed_mq.csv"))

mq_eye_Skinny_output_final = cbind(vec_ns, mq_eye_Skinny_output_final)
colnames(mq_eye_Skinny_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_eye_Skinny_output_final) = vec_ns

mq_eye_Skinny_output_final

# 7 Visualisation of mean, q25 and q75: Skinny -----------------------------------

#Skinny
dev.off()
par(bg = NA)
plot(mq_eye_Skinny_output_final[,1], mq_eye_Skinny_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Skinny", sub = "Eye distribution")
lines(mq_eye_Skinny_output_final[,1], mq_eye_Skinny_output_final[,3],
      col = "dark green")
lines(mq_eye_Skinny_output_final[,1], mq_eye_Skinny_output_final[,4],
      col = "dark red")
dev.off()


# 8 Analysis of scagnostic measures: Stringy ------------------------------------------------------------------

#preparing data set
osg_eye_Stringy = t(output_scag_measure)
dim(osg_eye_Stringy)
rownames(osg_eye_Stringy) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_eye_Stringy = function(x){

  mean_eye_Stringy = mean(x)
  q_eye_Stringy = quantile(x)
  mq_eye_Stringy = as.data.frame(rbind(mean_eye_Stringy, q_eye_Stringy[2],
                                       q_eye_Stringy[4]))
  rownames(mq_eye_Stringy) = c("mean", "25%", "75%")
  return (mq_eye_Stringy)
}

#test
mq_eye_Stringy(osg_eye_Stringy[1,])

mean_eye_Stringy = mean(osg_eye_Stringy[1,])
q_eye_Stringy = quantile(osg_eye_Stringy[1,])


mq_eye_Stringy_output = apply(osg_eye_Stringy, 1, mq_eye_Stringy)

unlist(mq_eye_Stringy_output)

mq_eye_Stringy_output_final = matrix(unlist(mq_eye_Stringy_output),
                                     ncol = 3, by = TRUE)
dim(mq_eye_Stringy_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_eye Skewed_mq.csv"))

mq_eye_Stringy_output_final = cbind(vec_ns, mq_eye_Stringy_output_final)
colnames(mq_eye_Stringy_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_eye_Stringy_output_final) = vec_ns

mq_eye_Stringy_output_final

# 8 Visualisation of mean, q25 and q75: Stringy -----------------------------------

#Stringy
dev.off()
par(bg = NA)
plot(mq_eye_Stringy_output_final[,1], mq_eye_Stringy_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Stringy", sub = "Eye distribution")
lines(mq_eye_Stringy_output_final[,1], mq_eye_Stringy_output_final[,3],
      col = "dark green")
lines(mq_eye_Stringy_output_final[,1], mq_eye_Stringy_output_final[,4],
      col = "dark red")
dev.off()



# 9 Analysis of scagnostic measures: Monotonic ------------------------------------------------------------------

#preparing data set
osg_eye_Monotonic = t(output_scag_measure)
dim(osg_eye_Monotonic)
rownames(osg_eye_Monotonic) = vec_ns


#Calculating the mean, q25 and q75
#i = 1
mq_eye_Monotonic = function(x){

  mean_eye_Monotonic = mean(x)
  q_eye_Monotonic = quantile(x)
  mq_eye_Monotonic = as.data.frame(rbind(mean_eye_Monotonic, q_eye_Monotonic[2],
                                         q_eye_Monotonic[4]))
  rownames(mq_eye_Monotonic) = c("mean", "25%", "75%")
  return (mq_eye_Monotonic)
}

#test
mq_eye_Monotonic(osg_eye_Monotonic[1,])

mean_eye_Monotonic = mean(osg_eye_Monotonic[1,])
q_eye_Monotonic = quantile(osg_eye_Monotonic[1,])


mq_eye_Monotonic_output = apply(osg_eye_Monotonic, 1, mq_eye_Monotonic)

unlist(mq_eye_Monotonic_output)

mq_eye_Monotonic_output_final = matrix(unlist(mq_eye_Monotonic_output),
                                       ncol = 3, by = TRUE)
dim(mq_eye_Monotonic_output_final)

#write.csv(output_scag_measure, file = paste0("DFScagSim_analysis_eye Skewed_mq.csv"))

mq_eye_Monotonic_output_final = cbind(vec_ns, mq_eye_Monotonic_output_final)
colnames(mq_eye_Monotonic_output_final) = c("No. of simulation","mean", "25%", "75%")
rownames(mq_eye_Monotonic_output_final) = vec_ns

mq_eye_Monotonic_output_final

# 9 Visualisation of mean, q25 and q75: Monotonic --------------------------------

#Monotonic
dev.off()
par(bg = NA)
plot(mq_eye_Monotonic_output_final[,1], mq_eye_Monotonic_output_final[,2],
     type = "l", ylim = c(0, 1), xlim = c(0, 50000),
     xlab = "no. data points", ylab = "q25, mean, q75", col = "dark blue",
     main = "Coefficient: Monotonic", sub = "Eye distribution")
lines(mq_eye_Monotonic_output_final[,1], mq_eye_Monotonic_output_final[,3],
      col = "dark green")
lines(mq_eye_Monotonic_output_final[,1], mq_eye_Monotonic_output_final[,4],
      col = "dark red")
dev.off()
