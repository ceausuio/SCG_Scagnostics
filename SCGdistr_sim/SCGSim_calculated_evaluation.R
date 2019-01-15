rm(list = ls())
options(java.parameters = "-Xmx8000m")


library(scagnostics)
library(alphahull)
# n = 500
scag = function(n){
# Normal Distribution -----------------------------------------------------
set.seed(1)

x = rnorm(n)
y = rnorm(n)

df_nor = as.data.frame(cbind(x,y))
s_nor = scagnostics(df_nor)


# Uniform Distribution ----------------------------------------------------
set.seed(1)
#n = 500
x = runif(n)
y = runif(n)

df_uni = as.data.frame(cbind(x,y))
s_uni = scagnostics(df_uni)


# Linear Distribution -----------------------------------------------------
set.seed(1)
# n = 500
x = runif(n, 0, 10)
b = rnorm(n, 1, 2)
y = x + b

df_lin = as.data.frame(cbind(x,y))
s_lin = scagnostics(df_lin)



# Heteroscedastic Distributionz -------------------------------------------
set.seed(1)
N = rep(1:n, 2)
a = 0
b = 1
sigma2 = N^1.3
eps = rnorm(N, mean = 0, sd = sqrt(sigma2))
y = a + b*N + 3*eps
mod <- lm(y ~ N)

df_het = as.data.frame(cbind(N,y))
s_het = scagnostics(df_het)


# Parable Distribution ----------------------------------------------------

set.seed(1)
# n = 500
x1 = rnorm(n, 1, 25)
x2 = rnorm(n, 4, 85)

a = -0.7
b = 1
c = -100

y1 = a*x1^2 + b*x2 + c

df_par = as.data.frame(cbind(x1, y1))
s_par = scagnostics(df_par)


# Clumps Distribution -----------------------------------------------------

set.seed(1)
# n = 500
x1 = rnorm(n, -1, 0.3)
y1 = rnorm(n, -1, 0.3)


x2 = rnorm(n, 1, 0.3)
y2 = rnorm(n, -1, 0.3)

x3 = rnorm(n, 0, 0.3)
y3 = rnorm(n, 1, 0.3)


df_clu = cbind(c(x1, x2, x3), c(y1, y2, y3))
s_clu = scagnostics(df_clu)


# “Eye” Distribution ------------------------------------------------------
set.seed(1)
# n = 500
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


# AR(1) Process Distribution ----------------------------------------------
set.seed(1)
# N = 500
N = arima.sim(list(order=c(1, 0, 0), ar=0.93), n)

df_ari = (cbind(1:n,N))
s_ari = scagnostics((df_ari))


# Sinus curve -------------------------------------------------------------
set.seed(1)
# n = 500
x = rnorm(n)
y = sin(x)

df_sin = as.data.frame(cbind(x, y))
s_sin = scagnostics(df_sin)


# Integers ----------------------------------------------------------------
set.seed(1)
N = rep(1:n, n)
a = 0
b = 1
sigma2 = N^1.3
eps = rnorm(N, mean = 0, sd = sqrt(sigma2))
y = a + b*N + 2*eps
mod <- lm(y ~ N)

df_int = as.data.frame(cbind(N, y))
s_int = scagnostics(df_int)


# Scagnostics of the Simulations ----------------------------------------------------------------

df_scag = rbind(s_nor, s_uni, s_lin, s_het, s_par, s_clu, s_eye, s_ari, s_sin,
                s_int)
#View(df_scag)

write.csv(df_scag, file = paste0("DFScagSim_",n,".csv"))

}
seq.scag = seq(100,10000,length=10) # also with seq.scag = seq(100,100000,length=10)

i=0
for(j in seq.scag){
  i=i+1
  scag(n=j)
  print(i)
}