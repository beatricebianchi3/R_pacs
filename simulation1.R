# simulazione 1 sistemata 
library(MASS)
library(plyr)
library(Matrix)
library(ggplot2)
library(RandomFields)
library(reshape2)
library(viridis)


library(fdaPDE2)

# Mesh generation
unit_square <- MeshUnitSquare(n=10)
Vh <- FunctionalSpace(unit_square, type = "fe", order = 1)
f <- Function(Vh)

# plot of the mesh  
plot(unit_square$nodes, main = "Mesh", asp = 1, pch = 16, col = "blue")
for (i in 1:nrow(unit_square$elements)) {
  polygon(unit_square$nodes[unit_square$elements[i, ], ], border = "black")
}


# Locations generation for beta
n_loc <- 225 # original # 225
locations <- matrix(data=NA,nrow=n_loc,ncol=2)
set.seed(543678)
locations[,1]<-runif(n_loc,0,1)
locations[,2]<-runif(n_loc,0,1)

# Plot locations
points(locations, col='red', pch=16)

epsilon=n_loc*10^-4

# Function generation:
# Function 2 in gamSim form mgcv package [2017]
#library(mgcv)
gamSim_2 <- function(x,y) {
  (0.4*pi^0.3)*(1.2*exp( - ((x-0.2)^2)/(0.3^2) - ((y-0.3)^2)/(0.4^2)) + 0.8*exp( - ((x-0.7)^2)/(0.3^2) - ((y-0.8)^2)/(0.4^2)))
}

exact_data <- as.matrix(gamSim_2(locations[, 1], locations[, 2]), ncol = 1)
# scale true value of f 
exact_data <- scale(exact_data)

# Covariates (random fields) definition
# Gaussian random field with mean=0, scale=0.05
set.seed(1)
model <- RMgauss(scale = 0.05)
cov_1 <- RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data
cov_1 <- unlist(cov_1)
plot(cov_1)

# Matérn random field with nu=1, sigma=2 and scale=0.1
set.seed(2)
model <- RMmatern(nu = 1, var = 2, scale = 0.1)
cov_2 <- RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data
cov_2 <- unlist(cov_2)
plot(cov_2)

# Deterministic field + Gaussian random field with mean=0, scale=0.05
set.seed(3)
model <- RMgauss(scale = 0.05)
fun_3 <- function(x,y){
  cos(5*(x+y))-(2*x-x*y^2)^2
}
add_3 <- RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data
add_3 <- unlist(add_3)
cov_3 <- fun_3(locations[,1],locations[,2]) + add_3
plot(cov_3)

# Deterministic field + Matérn random field with nu=1, sigma=2 and scale=0.1
set.seed(4)
model <- RMmatern(nu = 1, var = 2, scale = 0.1)
fun_4 <- function(x,y){
  cos(5*(x+y))-(2*x-x*y^2)^2
}
add_4 <- RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data
add_4 <- unlist(add_4)
cov_4 <- fun_4(locations[,1],locations[,2]) + add_4
plot(cov_4)

# Setting Parameters for smoothing.R
lambda = 10^seq(-3,3,by=0.25)

# Set beta_0
beta_0 = 0.0

#### Notice: add time estimated when done : 13 s for 20 simulations
# Number of repetitions for each simulation
rep = 10 # 1000 for paper

# List of hypotheses under which data are generated
beta_H1_list <- beta_0+ seq(from = 0, by = 0.01, length.out = 11)

# Define the standard deviation of epsilon_1,...,epsilon_n_loc
sd <- 0.1

# Real simulations: 1
covariates <- cov_1
covariates <- scale(covariates)

res_1_ex <- list()
res_1_non_ex <- list()
rmse_1_fdaPDE <- 0

for (i in 1:length(beta_H1_list)) {
  ps_ex = matrix(data = NA, nrow = 3, ncol = rep)
  ps_non_ex = matrix(data = NA, nrow = 3, ncol = rep)
  
  for (k in 1:rep) {
    set.seed(1)
    rand= rnorm(n_loc, 0, sd = sd)
    model <- RMgauss(scale = 0.05)
    covariates <- scale(unlist(RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data))
    observations <- covariates * beta_H1_list[1] + exact_data + rand # Build observations in H1
    
    data <- data.frame(
      y = observations,
      cov1 = covariates
    )
    
    model <- SRPDE(y ~ f + cov1, data = data)
    model$fit(
      calibration = 1e-2
    )
    
    model$inference("wald", "one_at_the_time", "exact", C = matrix(c(1), nrow = 1, ncol = 1), beta0 = 0)
    ps_ex[1, k] <- model$pvalues
    
    model$inference("speckman","one_at_the_time", "exact", C = matrix(c(1), nrow = 1, ncol = 1), beta0 = 0)
    ps_ex[2, k] <- model$pvalues
    
    model$inference("esf","one_at_the_time", "exact", C = matrix(c(1), nrow = 1, ncol = 1), beta0 = 0)
    ps_ex[3, k] <- model$pvalues
    
    model$inference("wald", "one_at_the_time", "nonexact", C = matrix(c(1), nrow = 1, ncol = 1), beta0 = 0)
    ps_non_ex[1, k] <- model$pvalues
    
    model$inference("speckman","one_at_the_time", "nonexact", C = matrix(c(1), nrow = 1, ncol = 1), beta0 = beta_0)
    ps_non_ex[2, k] <- model$pvalues
    
    model$inference("esf","one_at_the_time", "nonexact", C = matrix(c(1), nrow = 1, ncol = 1), beta0 = beta_0)
    ps_non_ex[3, k] <- model$pvalues
    
    #rmse_1_fdaPDE <- rmse_1_fdaPDE + Local_Solution_ex$solution$rmse[Local_Solution_ex$optimization$lambda_position] + Local_Solution_non_ex$solution$rmse[Local_Solution_non_ex$optimization$lambda_position]
  }
  
  res_1_ex[[i]] <- ps_ex
  res_1_non_ex[[i]] <- ps_non_ex
  
  print(i)
  # Sys.sleep(2)
}

#rmse_1_fdaPDE <- rmse_1_fdaPDE/(2*rep*length(beta_H1_list))





# Real simulations: 2
covariates <- cov_2
covariates <- scale(covariates)

res_2_ex <- list()
res_2_non_ex <- list()
rmse_2_fdaPDE <- 0

for (i in 1:length(beta_H1_list)) {
  ps_ex = matrix(data = NA, nrow = 4, ncol = rep)
  ps_non_ex = matrix(data = NA, nrow = 4, ncol = rep)
  
  for (k in 1:rep) {
    set.seed(k)
    rand= rnorm(n_loc, 0, sd = sd)
    model <- RMmatern(nu = 1, var = 2, scale = 0.1)
    covariates <- scale(unlist(RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data))
    observations <- covariates * beta_H1_list[i] + f + rand # Build observations in H1
    
    Local_Solution_ex <- smooth.FEM(locations = locations, observations=observations,
                                    covariates = covariates,
                                    FEMbasis=FEMbasis, lambda=lambda,
                                    lambda.selection.criterion='grid', lambda.selection.lossfunction = "GCV", DOF.evaluation = "exact",
                                    inference.data.object = Inference_beta_ex)
    
    Local_Solution_non_ex <- smooth.FEM(locations = locations, observations=observations,
                                        covariates = covariates,
                                        FEMbasis=FEMbasis, lambda=lambda,
                                        lambda.selection.criterion='grid', lambda.selection.lossfunction = "GCV", DOF.evaluation = "exact",
                                        inference.data.object = Inference_beta_non_ex)
    
    ps_ex[, k] <- c(Local_Solution_ex$inference$beta$p_values$wald[[1]], Local_Solution_ex$inference$beta$p_values$speckman[[1]], Local_Solution_ex$inference$beta$p_values$eigen_sign_flip[[1]],  Local_Solution_ex$inference$beta$p_values$enh_eigen_sign_flip[[1]])
    ps_non_ex[, k] <- c(Local_Solution_non_ex$inference$beta$p_values$wald[[1]], Local_Solution_non_ex$inference$beta$p_values$speckman[[1]], Local_Solution_non_ex$inferenc$beta$p_values$eigen_sign_flip[[1]],  Local_Solution_non_ex$inference$beta$p_values$enh_eigen_sign_flip[[1]])
    rmse_2_fdaPDE <- rmse_2_fdaPDE + Local_Solution_ex$solution$rmse[Local_Solution_ex$optimization$lambda_position] + Local_Solution_non_ex$solution$rmse[Local_Solution_non_ex$optimization$lambda_position]
  }
  res_2_ex[[i]] <- ps_ex
  res_2_non_ex[[i]] <- ps_non_ex
  
  print(i)
}

rmse_2_fdaPDE <- rmse_2_fdaPDE/(2*rep*length(beta_H1_list))

# Real simulations: 3
covariates <- cov_3
covariates <- scale(covariates)

res_3_ex <- list()
res_3_non_ex <- list()
rmse_3_fdaPDE <- 0

for (i in 1:length(beta_H1_list)) {
  ps_ex = matrix(data = NA, nrow = 4, ncol = rep)
  ps_non_ex = matrix(data = NA, nrow = 4, ncol = rep)
  
  for (k in 1:rep) {
    set.seed(k)
    rand= rnorm(n_loc, 0, sd = sd)
    model <- RMgauss(scale = 0.05)
    add_3 <- unlist(RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data)
    covariates <- scale(fun_3(locations[,1],locations[,2]) + add_3)
    observations <- covariates * beta_H1_list[i] + f + rand # Build observations in H1
    
    
    Local_Solution_ex <- smooth.FEM(locations = locations, observations=observations,
                                    covariates = covariates,
                                    FEMbasis=FEMbasis, lambda=lambda,
                                    lambda.selection.criterion='grid', lambda.selection.lossfunction = "GCV", DOF.evaluation = "exact",
                                    inference.data.object = Inference_beta_ex)
    
    Local_Solution_non_ex <- smooth.FEM(locations = locations, observations=observations,
                                        covariates = covariates,
                                        FEMbasis=FEMbasis, lambda=lambda,
                                        lambda.selection.criterion='grid', lambda.selection.lossfunction = "GCV", DOF.evaluation = "exact",
                                        inference.data.object = Inference_beta_non_ex)
    
    ps_ex[, k] <- c(Local_Solution_ex$inference$beta$p_values$wald[[1]], Local_Solution_ex$inference$beta$p_values$speckman[[1]], Local_Solution_ex$inference$beta$p_values$eigen_sign_flip[[1]],  Local_Solution_ex$inference$beta$p_values$enh_eigen_sign_flip[[1]])
    ps_non_ex[, k] <- c(Local_Solution_non_ex$inference$beta$p_values$wald[[1]], Local_Solution_non_ex$inference$beta$p_values$speckman[[1]], Local_Solution_non_ex$inferenc$beta$p_values$eigen_sign_flip[[1]],  Local_Solution_non_ex$inference$beta$p_values$enh_eigen_sign_flip[[1]])
    rmse_3_fdaPDE <- rmse_3_fdaPDE + Local_Solution_ex$solution$rmse[Local_Solution_ex$optimization$lambda_position] + Local_Solution_non_ex$solution$rmse[Local_Solution_non_ex$optimization$lambda_position]
  }
  res_3_ex[[i]] <- ps_ex
  res_3_non_ex[[i]] <- ps_non_ex
  
  print(i)
  Sys.sleep(2)
}

rmse_3_fdaPDE <- rmse_3_fdaPDE/(2*rep*length(beta_H1_list))

# Real simulations: 4
covariates <- cov_4
covariates <- scale(covariates)

res_4_ex <- list()
res_4_non_ex <- list()
rmse_4_fdaPDE <- 0

for (i in 1:length(beta_H1_list)) {
  ps_ex = matrix(data = NA, nrow = 4, ncol = rep)
  ps_non_ex = matrix(data = NA, nrow = 4, ncol = rep)
  
  for (k in 1:rep) {
    set.seed(k)
    rand= rnorm(n_loc, 0, sd = sd)
    model <- RMmatern(nu = 1, var = 2, scale = 0.1)
    add_4 <- unlist(RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data)
    covariates <- scale(fun_4(locations[,1],locations[,2]) + add_4)
    observations <- covariates * beta_H1_list[i] + f + rand # Build observations in H1
    
    
    # Local_Solution_ex <- smooth.FEM(locations = locations, observations=observations,
    #                                 covariates = covariates,
    #                                 FEMbasis=FEMbasis, lambda=lambda,
    #                                 lambda.selection.criterion='grid', lambda.selection.lossfunction = "GCV", DOF.evaluation = "exact",
    #                                 inference.data.object = Inference_beta_ex)
    
    Local_Solution_non_ex <- smooth.FEM(locations = locations, observations=observations,
                                        covariates = covariates,
                                        FEMbasis=FEMbasis, lambda=lambda,
                                        lambda.selection.criterion='grid', lambda.selection.lossfunction = "GCV", DOF.evaluation = "exact",
                                        inference.data.object = Inference_beta_non_ex)
    
    # ps_ex[, k] <- c(Local_Solution_ex$inference$beta$p_values$wald[[1]], Local_Solution_ex$inference$beta$p_values$speckman[[1]], Local_Solution_ex$inference$beta$p_values$eigen_sign_flip[[1]],  Local_Solution_ex$inference$beta$p_values$enh_eigen_sign_flip[[1]])
    ps_non_ex[, k] <- c(Local_Solution_non_ex$inference$beta$p_values$wald[[1]], Local_Solution_non_ex$inference$beta$p_values$speckman[[1]], Local_Solution_non_ex$inferenc$beta$p_values$eigen_sign_flip[[1]],  Local_Solution_non_ex$inference$beta$p_values$enh_eigen_sign_flip[[1]])
    # rmse_4_fdaPDE <- rmse_4_fdaPDE + Local_Solution_ex$solution$rmse[Local_Solution_ex$optimization$lambda_position] + Local_Solution_non_ex$solution$rmse[Local_Solution_non_ex$optimization$lambda_position]
  }
  # res_4_ex[[i]] <- ps_ex
  res_4_non_ex[[i]] <- ps_non_ex
  
  print(i)
}

rmse_4_fdaPDE <- rmse_4_fdaPDE/(2*rep*length(beta_H1_list))

# Results
# Change apply's input to obtain the different possible results
power_mat <- matrix(NA, nrow = 11, ncol = 3)
for(i in 1:11) {
  power_mat[i, ] <- apply(res_1_ex[[i]], 1, function(x) mean(x < 0.05))[c(1,2,3)]
}
t(power_mat)

#### Plotting ####
# Plotting with ggplot2
# genera n colori
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#genero 3 colori perchè le colonne di power_mat sono 3
cols <- gg_color_hue(ncol(power_mat))
cols <- cols[c(1,2,3)]

colnames(power_mat) <- c("Wald", "Speck", "ESF")
dat <- melt(power_mat) # cambia il formato dell matrice ma è sempre la stessa power_mat
dat$Var1 <- rep(beta_H1_list, length(cols)) # inserisco come prima colonna di dat i valori di beta 1  
colnames(dat) <- c("beta", "Test", "value") # quindi chiama le colonne beta, test e value 

dat$Test <- factor(dat$Test, levels = c("Wald", "Speck", "ESF")) # ome seconda colonna inserisco i nomi 
#dat$Test <- factor(dat$Test, levels = c("SESF","ESF","Wald", "Speck"))

# plotta sulle x i beta1  e sulle y la probabilità di aver commesso errore di 1 primo ovvero se il pvalue è <0.05
# 
plot_3<-ggplot(dat, aes(x = beta, y = value, colour = Test, linetype = Test, shape = Test)) +
  geom_abline(intercept = 0.05, slope = 0, linetype = "dashed") +
  geom_line(size = 1.5) + geom_point(size = 4) +
  theme_bw() +
  scale_shape_manual(values = c(1, 3, 15, 16)) +
  scale_linetype_manual(values=c("dotted", "dotdash","dashed","solid")) +
  scale_color_manual(values = cols) +
  ggtitle("(d)")+
  theme(plot.title = element_text(hjust = 0.5))+
  #theme(legend.key.width = unit(2,"cm")) +
  theme(legend.position = c(0.77, 0.36)) + ylim(c(-0.01, 1.01)) +
  #scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1))
  theme(legend.title = element_text(size = 20, face = "bold"), legend.text=element_text(size=20)) +
  theme(plot.title = element_text(size = 20, face = "bold")) + #theme(legend.position = "top") +
  theme(axis.text=element_text(size = 12), axis.title=element_text(size=20)) +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) +
  xlab(expr(beta)) + ylab("Power") #+ annotate("text", x = 0.19, y = 0.03, label = "0.05")

plot_3

require(gridExtra)
grid.arrange( plot_1, plot_2, plot_3, plot_4,ncol=2, nrow=2)

### SMALL  ####
# Results
# Change apply's input to obtain the different possible results

power_mat <- matrix(NA, nrow = 11, ncol = 4)
for(i in 1:11) {
  power_mat[i, ] <- apply(res_4_ex[[i]], 1, function(x) mean(x < 0.05))[1:4]
}
t(power_mat)

# Plotting with ggplot2
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <- gg_color_hue(ncol(power_mat))
cols <- cols[c(1,4,2,3)]

colnames(power_mat) <- c("Wald", "Speck", "ESF", "PESF")
dat <- melt(power_mat)
dat$Var1 <- rep(beta_H1_list, length(cols))
colnames(dat) <- c("beta", "Test", "value")

dat$Test <- factor(dat$Test, levels = c("Wald", "Speck", "ESF", "PESF"))

plot_3<-ggplot(dat, aes(x = beta, y = value, colour = Test, linetype = Test, shape = Test)) +
  geom_abline(intercept = 0.05, slope = 0, linetype = "dashed") +
  geom_line(size = 0.8) + geom_point(size = 2) +
  theme_bw() + 
  scale_shape_manual(values = c(1, 3, 15, 16)) +
  scale_linetype_manual(values=c("dotted", "dotdash","dashed","solid")) +
  scale_color_manual(values = cols) +
  scale_color_manual(values = cols) +
  ggtitle("(d)")+
  theme(plot.title = element_text(hjust = 0.5))+
  #theme(legend.key.width = unit(2,"cm")) +
  theme(legend.position = c(0.88, 0.30)) + ylim(c(-0.01, 1.01)) +
  #scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1))
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text=element_text(size=9)) +
  theme(plot.title = element_text(size = 16, face = "bold")) + #theme(legend.position = "top") +
  theme(axis.text=element_text(size = 16), axis.title=element_text(size=16)) +
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6)) +
  xlab(expr(beta)) + ylab("Power") #+ annotate("text", x = 0.19, y = 0.03, label = "0.05")

plot_3

require(gridExtra)
grid.arrange( plot_1, plot_2, plot_3, plot_4,ncol=2, nrow=2)

#### Medium #### 
#In paper
# Results
# Change apply's input to obtain the different possible results

power_mat <- matrix(NA, nrow = 11, ncol = 4)
for(i in 1:11) {
  power_mat[i, ] <- apply(res_3_non_ex[[i]], 1, function(x) mean(x < 0.05))[1:4]
}
t(power_mat)

# Plotting with ggplot2
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <- gg_color_hue(ncol(power_mat))
cols <- cols[c(1,4,2,3)]

colnames(power_mat) <- c("Wald", "Speck", "ESF", "PESF")
dat <- melt(power_mat)
dat$Var1 <- rep(beta_H1_list, length(cols))
colnames(dat) <- c("beta", "Test", "value")

dat$Test <- factor(dat$Test, levels = c("Wald", "Speck", "ESF", "PESF"))

plot_4<-ggplot(dat, aes(x = beta, y = value, colour = Test, linetype = Test, shape = Test)) +
  geom_abline(intercept = 0.05, slope = 0, linetype = "dashed") +
  geom_line(size = 1.2) + geom_point(size = 3) +
  theme_bw() + 
  scale_shape_manual(values = c(1, 3, 15, 16)) +
  scale_linetype_manual(values=c("dotted", "dotdash","dashed","solid")) +
  scale_color_manual(values = cols) +
  ggtitle("(c) - FSPAI")+
  theme(plot.title = element_text(hjust = 0.5))+
  #theme(legend.key.width = unit(2,"cm")) +
  theme(legend.position = c(0.87, 0.27)) + ylim(c(-0.01, 1.01)) +
  #scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1))
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text=element_text(size=14)) +
  theme(plot.title = element_text(size = 20, face = "bold")) + #theme(legend.position = "top") +
  theme(axis.text=element_text(size = 16), axis.title=element_text(size=16)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  xlab(expr(beta)) + ylab("Power") #+ annotate("text", x = 0.19, y = 0.03, label = "0.05")

plot_4

require(gridExtra)
grid.arrange( plot_1, plot_2, plot_3, plot_4,ncol=2, nrow=2)

#### Simulation 1 with GAM ####
library(mgcv)

### Supposing that the previous part has already been done and we can use the same data
rep=200 # 1000 for paper

loc1 <- as.numeric(locations[,1])
loc2 <- as.numeric(locations[,2])
# Sim mgcv
# Real simulations: 1
start <- Sys.time()
covariates <- cov_1
covariates <- scale(unlist(covariates))

res_1_mgcv <- list()
rmse_1_mgcv <- 0

for (i in 1:length(beta_H1_list)) {
  ps_mgcv = matrix(data = NA, nrow = 1, ncol = rep)
  
  for (k in 1:rep) {
    set.seed(k)
    rand= rnorm(n_loc, 0, sd = sd)
    model <- RMgauss(scale = 0.05)
    covariates <- as.numeric(scale(unlist(RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data)))
    observations <- covariates * beta_H1_list[i] + f + rand # Build observations in H1
    
    Local_Solution_mgcv <- gam(observations ~ te(loc1,loc2,k=30,d=2,bs="tp")+covariates)
    
    ps_mgcv[, k] <- summary(Local_Solution_mgcv)$p.table[2,4]
    pred_data <- as.data.frame(cbind(loc1,loc2,covariates))
    colnames(pred_data) <- c("loc1","loc2","covariates")
    rmse_1_mgcv <- rmse_1_mgcv + sqrt(mean((as.numeric(observations) - as.numeric(predict.gam(Local_Solution_mgcv, newdata = pred_data)))^2))
  }
  res_1_mgcv[[i]] <- ps_mgcv
  
  print(i)
}

rmse_1_mgcv <- rmse_1_mgcv/(rep*length(beta_H1_list))

# Real simulations: 2
covariates <- cov_2
covariates <- scale(covariates)

res_2_mgcv <- list()
rmse_2_mgcv <- 0

for (i in 1:length(beta_H1_list)) {
  ps_mgcv = matrix(data = NA, nrow = 1, ncol = rep)
  
  for (k in 1:rep) {
    set.seed(k)
    rand= rnorm(n_loc, 0, sd = sd)
    model <- RMmatern(nu = 1, var = 2, scale = 0.1)
    covariates <- as.numeric(scale(unlist(RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data)))
    observations <- covariates * beta_H1_list[i] + f + rand # Build observations in H1
    
    Local_Solution_mgcv <- gam(observations ~ te(loc1,loc2,k=30,d=2,bs="tp")+covariates)
    
    ps_mgcv[, k] <- summary(Local_Solution_mgcv)$p.table[2,4]
    pred_data <- as.data.frame(cbind(loc1,loc2,covariates))
    colnames(pred_data) <- c("loc1","loc2","covariates")
    rmse_2_mgcv <- rmse_2_mgcv + sqrt(mean((as.numeric(observations) - as.numeric(predict.gam(Local_Solution_mgcv, newdata = pred_data)))^2))
  }
  res_2_mgcv[[i]] <- ps_mgcv
  
  print(i)
}

rmse_2_mgcv <- rmse_2_mgcv/(rep*length(beta_H1_list))

# Real simulations: 3
covariates <- cov_3
covariates <- scale(covariates)

res_3_mgcv <- list()
rmse_3_mgcv <- 0

for (i in 1:length(beta_H1_list)) {
  ps_mgcv = matrix(data = NA, nrow = 1, ncol = rep)
  
  for (k in 1:rep) {
    set.seed(k)
    rand= rnorm(n_loc, 0, sd = sd)
    model <- RMgauss(scale = 0.05)
    add_3 <- unlist(RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data)
    covariates <- as.numeric(scale(fun_3(locations[,1],locations[,2]) + add_3))
    observations <- covariates * beta_H1_list[i] + f + rand # Build observations in H1
    
    Local_Solution_mgcv <- gam(observations ~ te(loc1,loc2,k=30,d=2,bs="tp")+covariates)
    
    ps_mgcv[, k] <- summary(Local_Solution_mgcv)$p.table[2,4]
    pred_data <- as.data.frame(cbind(loc1,loc2,covariates))
    colnames(pred_data) <- c("loc1","loc2","covariates")
    rmse_3_mgcv <- rmse_3_mgcv + sqrt(mean((as.numeric(observations) - as.numeric(predict.gam(Local_Solution_mgcv, newdata = pred_data)))^2))
  }
  res_3_mgcv[[i]] <- ps_mgcv
  
  print(i)
}

rmse_3_mgcv <- rmse_3_mgcv/(rep*length(beta_H1_list))

# Real simulations: 4
covariates <- cov_4
covariates <- scale(covariates)

res_4_mgcv <- list()
rmse_4_mgcv <- 0

for (i in 1:length(beta_H1_list)) {
  ps_mgcv = matrix(data = NA, nrow = 1, ncol = rep)
  
  for (k in 1:rep) {
    set.seed(k)
    rand= rnorm(n_loc, 0, sd = sd)
    model <- RMmatern(nu = 1, var = 2, scale = 0.1)
    add_4 <- unlist(RandomFields::RFsimulate(model, x = SpatialPoints(locations), n = 1)@data)
    covariates <- as.numeric(scale(fun_4(locations[,1],locations[,2]) + add_4))
    observations <- covariates * beta_H1_list[i] + f + rand # Build observations in H1
    
    Local_Solution_mgcv <- gam(observations ~ te(loc1,loc2,k=20,d=2,bs="tp")+covariates)
    
    ps_mgcv[, k] <- summary(Local_Solution_mgcv)$p.table[2,4]
    pred_data <- as.data.frame(cbind(loc1,loc2,covariates))
    colnames(pred_data) <- c("loc1","loc2","covariates")
    rmse_4_mgcv <- rmse_4_mgcv + sqrt(mean((as.numeric(observations) - as.numeric(predict.gam(Local_Solution_mgcv, newdata = pred_data)))^2))
  }
  res_4_mgcv[[i]] <- ps_mgcv
  
  print(i)
}

rmse_4_mgcv <- rmse_4_mgcv/(rep*length(beta_H1_list))

### Power matrix
power_mat <- matrix(NA, nrow = 11, ncol = 1)
for(i in 1:length(beta_H1_list)) {
  power_mat[i, ] <- apply(res_1_mgcv[[i]], 1, function(x) mean(x < 0.05))[1]
}
t(power_mat)
### Supposing that the previous part has already been done and we can use the same data

comparison_mat <- matrix(NA, nrow = 11,ncol = 2)
for(i in 1:length(beta_H1_list)) {
  comparison_mat[i, 1] <- apply(res_4_mgcv[[i]], 1, function(x) mean(x < 0.05))[[1]]
  comparison_mat[i, 2] <- apply(res_4_ex[[i]], 1, function(x) mean(x < 0.05))[[4]]
}
t(comparison_mat)

# Plotting comparison
# Plotting with ggplot2
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <- gg_color_hue(4)
cols <- cols[c(1,4)]

colnames(comparison_mat) <- c("GAM-Wald", "Enh-ESF")
dat <- melt(comparison_mat)
dat$Var1 <- rep(beta_H1_list, length(cols))
colnames(dat) <- c("beta", "Test", "value")

dat$Test <- factor(dat$Test, levels = c("GAM-Wald", "Enh-ESF"))

plot_4<-ggplot(dat, aes(x = beta, y = value, colour = Test, linetype = Test, shape = Test)) +
  geom_abline(intercept = 0.05, slope = 0, linetype = "dashed") +
  geom_line(size = 1.2) + geom_point(size = 3) +
  theme_bw() + #scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash")) +
  scale_color_manual(values = cols) +
  ggtitle("(d) - FSPAI")+
  theme(plot.title = element_text(hjust = 0.5))+
  #theme(legend.key.width = unit(2,"cm")) +
  theme(legend.position = c(0.84, 0.27)) + ylim(c(-0.01, 1.01)) +
  #scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1))
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text=element_text(size=14)) +
  theme(plot.title = element_text(size = 20, face = "bold")) + #theme(legend.position = "top") +
  theme(axis.text=element_text(size = 16), axis.title=element_text(size=16)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  xlab(expr(beta)) + ylab("Power") #+ annotate("text", x = 0.19, y = 0.03, label = "0.05")

plot_4

require(gridExtra)
grid.arrange( plot_1, plot_2, plot_3, plot_4,ncol=2, nrow=2)