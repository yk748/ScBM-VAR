rm(list=ls())

library(superheat)
library(igraph)
library(MASS)
library(plotly)
library(combinat)
library(networkD3)
library(dplyr)
library(ggsankey)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(reshape2)
library(pheatmap)

source('ScBM_library.R')

#####################################################
# Employees on nonfarm payrolls by industry sector and selected industry detail
# https://www.bls.gov/news.release/empsit.t17.htm

dir_employ <- "./Employees"
date <- mapply(t=1:121,function(t)paste0(rep(c(1990:2020),each=4)[t],
                                         rep(c("-Q1","-Q2","-Q3","-Q4"),31)[t]))
# Total private
## Goods-producing
### 1. Mining and logging
dat_mine <- read.csv(paste0(dir_employ,"/Mining and Logging.csv"),header=TRUE)
colnames(dat_mine) <- c("date","value")
dat_mine <- dat_mine[which(dat_mine$date=="1990-01-01"):
                       which(dat_mine$date=="2020-03-01"),]
dat_mineQ <- data.frame(date=date,
                        value=rep(NA,121))
for(t in 1:121){dat_mineQ$value[t] <- sum(log(dat_mine$value[(3*(t-1)+1):(3*t)]))}
dat_minQ <- data.frame(date=date[-1],value=diff(dat_mineQ$value))

### 2. Construction
dat_cons <- read.csv(paste0(dir_employ,"/Construction.csv"),header=TRUE)
colnames(dat_cons) <- c("date","value")
dat_cons <- dat_cons[which(dat_cons$date=="1990-01-01"):
                       which(dat_cons$date=="2020-03-01"),]
dat_consQ <- data.frame(date=date,
                        value=rep(NA,121))
for(t in 1:121){dat_consQ$value[t] <- sum(log(dat_cons$value[(3*(t-1)+1):(3*t)]))}
dat_consQ <- data.frame(date=date[-1],value=diff(dat_consQ$value))

### Manufacturing
#### 3. Durable goods
dat_dur <- read.csv(paste0(dir_employ,"/Durable goods.csv"),header=TRUE)
colnames(dat_dur) <- c("date","value")
dat_dur <- dat_dur[c(which(dat_dur$date=="1990-01-01"):
                       which(dat_dur$date=="2020-03-01")),]
dat_durQ <- data.frame(date=date,
                       value=rep(NA,121))
for(t in 1:121){dat_durQ$value[t] <- sum(log(dat_dur$value[(3*(t-1)+1):(3*t)]))}
dat_durQ <- data.frame(date=date[-1],value=diff(dat_durQ$value))

#### 4. Nondurable goods
dat_nondur <- read.csv(paste0(dir_employ,"/Nondurable goods.csv"),header=TRUE)
colnames(dat_nondur) <- c("date","value")
dat_nondur <- dat_nondur[c(which(dat_nondur$date=="1990-01-01"):
                             which(dat_nondur$date=="2020-03-01")),]
dat_nondurQ <- data.frame(date=date,
                          value=rep(NA,121))
for(t in 1:121){dat_nondurQ$value[t] <- sum(log(dat_nondur$value[(3*(t-1)+1):(3*t)]))}
dat_nondurQ <- data.frame(date=date[-1],value=diff(dat_nondurQ$value))

## Private service-providing
### Trade, transportation, and utilities
#### 5. Wholesale trade
dat_whole <- read.csv(paste0(dir_employ,"/Wholesale Trade.csv"),header=TRUE)
colnames(dat_whole) <- c("date","value")
dat_whole <- dat_whole[c(which(dat_whole$date=="1990-01-01"):
                           which(dat_whole$date=="2020-03-01")),]
dat_wholeQ <- data.frame(date=date,
                         value=rep(NA,121))
for(t in 1:121){dat_wholeQ$value[t] <- sum(log(dat_whole$value[(3*(t-1)+1):(3*t)]))}
dat_wholeQ <- data.frame(date=date[-1],value=diff(dat_wholeQ$value))

#### 6. Retail trade
dat_retail <- read.csv(paste0(dir_employ,"/Retail Trade.csv"),header=TRUE)
colnames(dat_retail) <- c("date","value")
dat_retail <- dat_retail[c(which(dat_retail$date=="1990-01-01"):
                             which(dat_retail$date=="2020-03-01")),]
dat_retailQ <- data.frame(date=date,
                          value=rep(NA,121))
for(t in 1:121){dat_retailQ$value[t] <- sum(log(dat_retail$value[(3*(t-1)+1):(3*t)]))}
dat_retailQ <- data.frame(date=date[-1],value=diff(dat_retailQ$value))

#### 7. Transportation and warehousing
dat_trans <- read.csv(paste0(dir_employ,"/Transportation and Warehousing.csv"),header=TRUE)
colnames(dat_trans) <- c("date","value")
dat_trans <- dat_trans[c(which(dat_trans$date=="1990-01-01"):
                           which(dat_trans$date=="2020-03-01")),]
dat_transQ <- data.frame(date=date,
                         value=rep(NA,121))
for(t in 1:121){dat_transQ$value[t] <- sum(log(dat_trans$value[(3*(t-1)+1):(3*t)]))}
dat_transQ <- data.frame(date=date[-1],value=diff(dat_transQ$value))

#### 8. Utilities
dat_util <- read.csv(paste0(dir_employ,"/Utilities.csv"),header=TRUE)
colnames(dat_util) <- c("date","value")
dat_util <- dat_util[c(which(dat_util$date=="1990-01-01"):
                         which(dat_util$date=="2020-03-01")),]
dat_utilQ <- data.frame(date=date,
                        value=rep(NA,121))
for(t in 1:121){dat_utilQ$value[t] <- sum(log(dat_util$value[(3*(t-1)+1):(3*t)]))}
dat_utilQ <- data.frame(date=date[-1],value=diff(dat_utilQ$value))

### 9. Information
dat_info <- read.csv(paste0(dir_employ,"/Information.csv"),header=TRUE)
colnames(dat_info) <- c("date","value")
dat_info <- dat_info[c(which(dat_info$date=="1990-01-01"):
                         which(dat_info$date=="2020-03-01")),]
dat_infoQ <- data.frame(date=date,
                        value=rep(NA,121))
for(t in 1:121){dat_infoQ$value[t] <- sum(log(dat_info$value[(3*(t-1)+1):(3*t)]))}
dat_infoQ <- data.frame(date=date[-1],value=diff(dat_infoQ$value))

### Financial activities
#### 10. Finance and insurance
dat_fin <- read.csv(paste0(dir_employ,"/Finance and insurance.csv"),header=TRUE)
colnames(dat_fin) <- c("date","value")
dat_fin <- dat_fin[c(which(dat_fin$date=="1990-01-01"):
                       which(dat_fin$date=="2020-03-01")),]
dat_finQ <- data.frame(date=date,
                       value=rep(NA,121))
for(t in 1:121){dat_finQ$value[t] <- sum(log(dat_fin$value[(3*(t-1)+1):(3*t)]))}
dat_finQ <- data.frame(date=date[-1],value=diff(dat_finQ$value))

#### 11. Real estate and rental and leasing
dat_real <- read.csv(paste0(dir_employ,"/Real Estate and Rental and Leasing.csv"),header=TRUE)
colnames(dat_real) <- c("date","value")
dat_real <- dat_real[c(which(dat_real$date=="1990-01-01"):
                         which(dat_real$date=="2020-03-01")),]
dat_realQ <- data.frame(date=date,
                        value=rep(NA,121))
for(t in 1:121){dat_realQ$value[t] <- sum(log(dat_real$value[(3*(t-1)+1):(3*t)]))}
dat_realQ <- data.frame(date=date[-1],value=diff(dat_realQ$value))

### Professional and business services
#### 12. Professional, scientific, and technical services
dat_prof <- read.csv(paste0(dir_employ,"/Professional, Scientific, and Technical Services.csv"),header=TRUE)
colnames(dat_prof) <- c("date","value")
dat_prof <- dat_prof[c(which(dat_prof$date=="1990-01-01"):
                         which(dat_prof$date=="2020-03-01")),]
dat_profQ <- data.frame(date=date,
                        value=rep(NA,121))
for(t in 1:121){dat_profQ$value[t] <- sum(log(dat_prof$value[(3*(t-1)+1):(3*t)]))}
dat_profQ <- data.frame(date=date[-1],value=diff(dat_profQ$value))

#### 13. Management of companies and enterprises
dat_manage <- read.csv(paste0(dir_employ,"/Management of Companies and Enterprises.csv"),header=TRUE)
colnames(dat_manage) <- c("date","value")
dat_manage <- dat_manage[c(which(dat_manage$date=="1990-01-01"):
                             which(dat_manage$date=="2020-03-01")),]
dat_manageQ <- data.frame(date=date,
                          value=rep(NA,121))
for(t in 1:121){dat_manageQ$value[t] <- sum(log(dat_manage$value[(3*(t-1)+1):(3*t)]))}
dat_manageQ <- data.frame(date=date[-1],value=diff(dat_manageQ$value))

#### 14. Administrative and support and waste management and remediation services
dat_admin <- read.csv(paste0(dir_employ,"/Administrative and Support Services.csv"),header=TRUE)
colnames(dat_admin) <- c("date","value")
dat_admin <- dat_admin[c(which(dat_admin$date=="1990-01-01"):
                           which(dat_admin$date=="2020-03-01")),]
dat_adminQ <- data.frame(date=date,
                         value=rep(NA,121))
for(t in 1:121){dat_adminQ$value[t] <- sum(log(dat_admin$value[(3*(t-1)+1):(3*t)]))}
dat_adminQ <- data.frame(date=date[-1],value=diff(dat_adminQ$value))

### Private education and health services
#### 15. Private educational services
dat_edu <- read.csv(paste0(dir_employ,"/Private Educational Services.csv"),header=TRUE)
colnames(dat_edu) <- c("date","value")
dat_edu <- dat_edu[c(which(dat_edu$date=="1990-01-01"):
                       which(dat_edu$date=="2020-03-01")),]
dat_eduQ <- data.frame(date=date,
                       value=rep(NA,121))
for(t in 1:121){dat_eduQ$value[t] <- sum(log(dat_edu$value[(3*(t-1)+1):(3*t)]))}
dat_eduQ <- data.frame(date=date[-1],value=diff(dat_eduQ$value))

#### 16. Health care and social assistance
dat_health <- read.csv(paste0(dir_employ,"/Health Care and Social Assistance.csv"),header=TRUE)
colnames(dat_health) <- c("date","value")
dat_health <- dat_health[c(which(dat_health$date=="1990-01-01"):
                             which(dat_health$date=="2020-03-01")),]
dat_healthQ <- data.frame(date=date,
                          value=rep(NA,121))
for(t in 1:121){dat_healthQ$value[t] <- sum(log(dat_health$value[(3*(t-1)+1):(3*t)]))}
dat_healthQ <- data.frame(date=date[-1],value=diff(dat_healthQ$value))

### Leisure and hospitality
#### 17. Arts, entertainment, and recreation
dat_art <- read.csv(paste0(dir_employ,"/Arts, Entertainment, and Recreation.csv"),header=TRUE)
colnames(dat_art) <- c("date","value")
dat_art <- dat_art[c(which(dat_art$date=="1990-01-01"):
                       which(dat_art$date=="2020-03-01")),]
dat_artQ <- data.frame(date=date,
                       value=rep(NA,121))
for(t in 1:121){dat_artQ$value[t] <- sum(log(dat_art$value[(3*(t-1)+1):(3*t)]))}
dat_artQ <- data.frame(date=date[-1],value=diff(dat_artQ$value))

#### 18. Accommodation and food services
dat_acc <- read.csv(paste0(dir_employ,"/Accommodation and Food Services.csv"),header=TRUE)
colnames(dat_acc) <- c("date","value")
dat_acc <- dat_acc[c(which(dat_acc$date=="1990-01-01"):
                       which(dat_acc$date=="2020-03-01")),]
dat_accQ <- data.frame(date=date,
                       value=rep(NA,121))
for(t in 1:121){dat_accQ$value[t] <- sum(log(dat_acc$value[(3*(t-1)+1):(3*t)]))}
dat_accQ <- data.frame(date=date[-1],value=diff(dat_accQ$value))

### 19. Other services
dat_other <- read.csv(paste0(dir_employ,"/Other Services.csv"),header=TRUE)
colnames(dat_other) <- c("date","value")
dat_other <- dat_other[c(which(dat_other$date=="1990-01-01"):
                           which(dat_other$date=="2020-03-01")),]
dat_otherQ <- data.frame(date=date,
                         value=rep(NA,121))
for(t in 1:121){dat_otherQ$value[t] <- sum(log(dat_other$value[(3*(t-1)+1):(3*t)]))}
dat_otherQ <- data.frame(date=date[-1],value=diff(dat_otherQ$value))

# Government
## 20. Federal
dat_fed <- read.csv(paste0(dir_employ,"/Federal.csv"),header=TRUE)
colnames(dat_fed) <- c("date","value")
dat_fed <- dat_fed[c(which(dat_fed$date=="1990-01-01"):
                       which(dat_fed$date=="2020-03-01")),]
dat_fedQ <- data.frame(date=date,
                       value=rep(NA,121))
for(t in 1:121){dat_fedQ$value[t] <- sum(log(dat_fed$value[(3*(t-1)+1):(3*t)]))}
dat_fedQ <- data.frame(date=date[-1],value=diff(dat_fedQ$value))

## 21. State government
dat_state <- read.csv(paste0(dir_employ,"/State Government.csv"),header=TRUE)
colnames(dat_state) <- c("date","value")
dat_state <- dat_state[c(which(dat_state$date=="1990-01-01"):
                           which(dat_state$date=="2020-03-01")),]
dat_stateQ <- data.frame(date=date,
                         value=rep(NA,121))
for(t in 1:121){dat_stateQ$value[t] <- sum(log(dat_state$value[(3*(t-1)+1):(3*t)]))}
dat_stateQ <- data.frame(date=date[-1],value=diff(dat_stateQ$value))

## 22. Local government
dat_local <- read.csv(paste0(dir_employ,"/Local Government.csv"),header=TRUE)
colnames(dat_local) <- c("date","value")
dat_local <- dat_local[c(which(dat_local$date=="1990-01-01"):
                           which(dat_local$date=="2020-03-01")),]
dat_localQ <- data.frame(date=date,
                         value=rep(NA,121))
for(t in 1:121){dat_localQ$value[t] <- sum(log(dat_local$value[(3*(t-1)+1):(3*t)]))}
dat_localQ <- data.frame(date=date[-1],value=diff(dat_localQ$value))


data <- cbind(dat_minQ$value,
              dat_consQ$value,
              dat_durQ$value,
              dat_nondurQ$value,
              dat_wholeQ$value,
              dat_retailQ$value,
              dat_transQ$value,
              dat_utilQ$value,
              dat_infoQ$value,
              dat_finQ$value,
              dat_realQ$value,
              dat_profQ$value,
              dat_manageQ$value,
              dat_adminQ$value,
              dat_eduQ$value,
              dat_healthQ$value,
              dat_artQ$value,
              dat_accQ$value,
              dat_otherQ$value,
              dat_fedQ$value,
              dat_stateQ$value,
              dat_localQ$value)
colnames(data) <- c(1:22)

#####################################################
# Figure 1:
Y_rep <- data[,c(1:8)]
colnames(Y_rep) <- c("Mining","Construction","Duration","Nonduration",
                     "Wholesale","Retail","Transportation","Utility")

pdf("timeplot_PVAR.pdf", width = 12, height = 9)
par(mar = c(5, 4, 4, 2) + 1)
plot.ts(Y_rep,"single",main="")
dev.off()
pdf("acf_PVAR.pdf", width = 9, height = 9)
par(mar = c(5, 4, 4, 2) + 1)
plot(acf(Y_rep[,c(1,4,5,6)],lag.max=20,plot=FALSE),ylab="")
dev.off()
pdf("pacf_PVAR.pdf", width = 9, height = 9)
par(mar = c(5, 4, 4, 2) + 1)
plot(pacf(Y_rep[,c(1,4,5,6)],lag.max=20,plot=FALSE),ylab="")
dev.off()


#####################################################
# Estimation:
Yt <- t(data)
TT <- dim(Yt)[2]
q <- dim(Yt)[1]
Est_PVAR_ols <- PVAR_ols(Yt,s=4,p=1)
superheat(Est_PVAR_ols$Phi_hat)

Phi_series <- Est_PVAR_ols$Phi_hat

#####################################################
# Figure C1:
# Check the number of communities:
svd_1 <- svd(Phi_series[,1:q])
svd_2 <- svd(Phi_series[,(q+1):(2*q)])
svd_3 <- svd(Phi_series[,(2*q+1):(3*q)])
svd_4 <- svd(Phi_series[,(3*q+1):(4*q)])


pdf("rank_PVAR.pdf", width = 13, height = 9)
par(mfrow=c(2,2))
# season 1
plot(1:q,
     mapply(k1=1:q,function(k1)norm(svd_1$u[,c(1:k1)] %*% diag(svd_1$d[c(1:k1)],k1) %*% t(svd_1$v[,c(1:k1)]),"F")^2/norm(Phi_series[,1:q],"F")^2),
     type="p",ylab="Proportion of variation",xlab="Mode index",main="season 1",
     ylim=c(0,1), pch=3,col="blue",cex=1.5)
lines(mapply(k1=1:q,function(k1)norm(svd_1$u[,c(1:k1)] %*% diag(svd_1$d[c(1:k1)],k1) %*% t(svd_1$v[,c(1:k1)]),"F")^2/norm(Phi_series[,1:q],"F")^2),
      col="blue",lty="dashed",cex=1.5)
points(1:q,mapply(k1=1:q,function(k1)norm(svd_1$u[,k1] %*% diag(svd_1$d[k1],1) %*% t(svd_1$v[,k1]),"F")^2/norm(Phi_series[,1:q],"F")^2),
       col="red",pch=1,cex=1.5)
lines(mapply(k1=1:q,function(k1)norm(svd_1$u[,k1] %*% diag(svd_1$d[k1],1) %*% t(svd_1$v[,k1]),"F")^2/norm(Phi_series[,1:q],"F")^2),
      col="red",cex=1.5)
abline(h=c(0,0.7),lty="dotted")

# season 2
plot(1:q,
     mapply(k2=1:q,function(k2)norm(svd_2$u[,c(1:k2)] %*% diag(svd_2$d[c(1:k2)],k2) %*% t(svd_2$v[,c(1:k2)]),"F")^2/norm(Phi_series[,(q+1):(2*q)],"F")^2),
     type="p",ylab="Proportion of variation",xlab="Mode index",main="season 2",
     ylim=c(0,1), pch=3,col="blue",cex=1.5)
lines(mapply(k2=1:q,function(k2)norm(svd_2$u[,c(1:k2)] %*% diag(svd_2$d[c(1:k2)],k2) %*% t(svd_2$v[,c(1:k2)]),"F")^2/norm(Phi_series[,(q+1):(2*q)],"F")^2),
      col="blue",lty="dashed",cex=1.5)
points(1:q,mapply(k2=1:q,function(k2)norm(svd_2$u[,k2] %*% diag(svd_2$d[k2],1) %*% t(svd_2$v[,k2]),"F")^2/norm(Phi_series[,(q+1):(2*q)],"F")^2),
       col="red",pch=1,cex=1.5)
lines(mapply(k2=1:q,function(k2)norm(svd_2$u[,k2] %*% diag(svd_2$d[k2],1) %*% t(svd_2$v[,k2]),"F")^2/norm(Phi_series[,(q+1):(2*q)],"F")^2),
      col="red",cex=1.5)
abline(h=c(0,0.7),lty="dotted")

# season 3
plot(1:q,
     mapply(k3=1:q,function(k3)norm(svd_3$u[,c(1:k3)] %*% diag(svd_3$d[c(1:k3)],k3) %*% t(svd_3$v[,c(1:k3)]),"F")^2/norm(Phi_series[,(2*q+1):(3*q)],"F")^2),
     type="p",ylab="Proportion of variation",xlab="Mode index",main="season 3",
     ylim=c(0,1), pch=3,col="blue",cex=1.5)
lines(mapply(k3=1:q,function(k3)norm(svd_3$u[,c(1:k3)] %*% diag(svd_3$d[c(1:k3)],k3) %*% t(svd_3$v[,c(1:k3)]),"F")^2/norm(Phi_series[,(2*q+1):(3*q)],"F")^2),
      col="blue",lty="dashed",cex=1.5)
points(1:q,mapply(k3=1:q,function(k3)norm(svd_3$u[,k3] %*% diag(svd_3$d[k3],1) %*% t(svd_3$v[,k3]),"F")^2/norm(Phi_series[,(2*q+1):(3*q)],"F")^2),
       col="red",pch=1,cex=1.5)
lines(mapply(k3=1:q,function(k3)norm(svd_3$u[,k3] %*% diag(svd_3$d[k3],1) %*% t(svd_3$v[,k3]),"F")^2/norm(Phi_series[,(2*q+1):(3*q)],"F")^2),
      col="red",cex=1.5)
abline(h=c(0,0.7),lty="dotted")

# season 4
plot(1:q,
     mapply(k4=1:q,function(k4)norm(svd_4$u[,c(1:k4)] %*% diag(svd_4$d[c(1:k4)],k4) %*% t(svd_4$v[,c(1:k4)]),"F")^2/norm(Phi_series[,(3*q+1):(4*q)],"F")^2),
     type="p",ylab="Proportion of variation",xlab="Mode index",main="season 4",
     ylim=c(0,1), pch=3,col="blue",cex=1.5)
lines(mapply(k4=1:q,function(k4)norm(svd_4$u[,c(1:k4)] %*% diag(svd_4$d[c(1:k4)],k4) %*% t(svd_4$v[,c(1:k4)]),"F")^2/norm(Phi_series[,(3*q+1):(4*q)],"F")^2),
      col="blue",lty="dashed",cex=1.5)
points(1:q,mapply(k4=1:q,function(k4)norm(svd_4$u[,k4] %*% diag(svd_4$d[k4],1) %*% t(svd_4$v[,k4]),"F")^2/norm(Phi_series[,(3*q+1):(4*q)],"F")^2),
       col="red",pch=1,cex=1.5)
lines(mapply(k4=1:q,function(k4)norm(svd_4$u[,k4] %*% diag(svd_4$d[k4],1) %*% t(svd_4$v[,k4]),"F")^2/norm(Phi_series[,(3*q+1):(4*q)],"F")^2),
      col="red",cex=1.5)
abline(h=c(0,0.7),lty="dotted")
dev.off()


# seasons determined:
k1 <- 2; k2 <- 3; k3 <- 4; k4 <- 2

#####################################################
# Run the algorithm:
set.seed(123)
n_comm <- c(2,3,3,4,4,4,4,2)
alpha_hat <- alpha_cv(Phi_series,q,p=1,n_comm=n_comm,type="PVAR",clst="kmeans",fold=5)

PisCES_output <- PisCES(Phi_series,q,p=1,n_comm=n_comm,s=4,alpha=alpha_hat$alpha_hat)
PVAR_4R1L <- Blockbuster(cbind(PisCES_output$VR_bar[[4]],
                               PisCES_output$VL_bar[[1]]),q,2,clst="kmeans")
PVAR_1R2L <- Blockbuster(cbind(PisCES_output$VR_bar[[1]],
                               PisCES_output$VL_bar[[2]]),q,3,clst="kmeans")
PVAR_2R3L <- Blockbuster(cbind(PisCES_output$VR_bar[[2]],
                               PisCES_output$VL_bar[[3]]),q,4,clst="kmeans")
PVAR_3R4L <- Blockbuster(cbind(PisCES_output$VR_bar[[3]],
                               PisCES_output$VL_bar[[4]]),q,4,clst="kmeans")

# Refine the labels:
result_list <- list(PVAR_4R1L,PVAR_1R2L,PVAR_2R3L,PVAR_3R4L)
PVAR_4R1L <- result_list[[1]]
PVAR_1R2L <- result_list[[2]]
PVAR_2R3L <- result_list[[3]]
PVAR_3R4L <- result_list[[4]]

flow <- mapply(x=c(1,2,2,3,3,4,4,1),function(x)result_list[[x]]$group_1q)
rownames(flow) <-  c("Mining","Construction","Durable","Nondurable","Wholesale",
                     "Retail","Transportation","Utility","Information",
                     "Finance","RealEstate","Professional","Management","Administration",
                     "Education","HealthCare","Arts","Accommodation","Other",
                     "Federal","State","Local")
colnames(flow) <- c("Q1.send","Q1.receive",
                    "Q2.send","Q2.receive",
                    "Q3.send","Q3.receive",
                    "Q4.send","Q4.receive")

#####################################################
# Figure 2:
# Construct data flows:
data_flow <- data.frame(x = c(rep("Q1.sending",q),rep("Q1.receiving",q),
                              rep("Q2.sending",q),rep("Q2.receiving",q),
                              rep("Q3.sending",q),rep("Q3.receiving",q),
                              rep("Q4.sending",q),rep("Q4.receiving",q)), 
                        node = c(flow[,1],flow[,2],flow[,3],flow[,4],flow[,5],
                                 flow[,6],flow[,7],flow[,8]),
                        next_x = c(rep("Q1.receiving",q),rep("Q2.sending",q),
                                   rep("Q2.receiving",q),rep("Q3.sending",q),
                                   rep("Q3.receiving",q),rep("Q4.sending",q),
                                   rep("Q4.receiving",q),rep(NA,q)),
                        next_node = c(flow[,2],flow[,3],flow[,4],flow[,5],flow[,6],
                                      flow[,7],flow[,8],flow[,1]))

data_flow$x <- factor(data_flow$x,levels=c("Q1.sending","Q1.receiving",
                                           "Q2.sending","Q2.receiving",
                                           "Q3.sending","Q3.receiving",
                                           "Q4.sending","Q4.receiving"))
data_flow$next_x <- factor(data_flow$next_x,levels=c("Q1.sending","Q1.receiving",
                                                     "Q2.sending","Q2.receiving",
                                                     "Q3.sending","Q3.receiving",
                                                     "Q4.sending","Q4.receiving"))

data_flow$node <- factor(data_flow$node,levels=c(1,2,3,4))
data_flow$next_node <- factor(data_flow$next_node,levels=c(1,2,3,4))

colorful <- c(brewer.pal(12,"Paired"),brewer.pal(10,"Set3"))

# Plotting:
pdf("sanky_PVAR.pdf", width = 13, height = 9)
ggplot(data_flow, aes(x = x, next_x = next_x, 
                      node = node, next_node = next_node, 
                      fill = node, label = node)) +
  geom_sankey(flow.alpha = 0.3) +
  theme_sankey(base_size = 15) +
  labs(x = NULL) + labs(y = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  #Q1.in, group 1
  annotate("text", x = 1, y = c(-3:(-2-length(names(which(flow[,1]==1))))), 
           color=colorful[which(flow[,1]==1)],
           size=4,label=names(which(flow[,1]==1)),fontface=2)  + 
  #Q1.in, group 2
  annotate("text", x = 1, y = c(1:(length(names(which(flow[,1]==2))))), 
           color=colorful[which(flow[,1]==2)],
           size=4,label=names(which(flow[,1]==2)),fontface=2)  +
  #Q2.in, group 1
  annotate("text", x = 3, y = c(-7:(-6-length(names(which(flow[,3]==1))))), 
           color=colorful[which(flow[,3]==1)],
           size=4,label=names(which(flow[,3]==1)),fontface=2)  + 
  #Q2.in, group 2
  annotate("text", x = 3, y = c(-3.5:(-4.5+length(names(which(flow[,3]==2))))), 
           color=colorful[which(flow[,3]==2)], 
           size=4,label=names(which(flow[,3]==2)),fontface=2) + 
  #Q2.in, group 3
  annotate("text", x = 3, y = c(7.5:(6.5+length(names(which(flow[,3]==3))))), 
           color=colorful[which(flow[,3]==3)], 
           size=4,label=names(which(flow[,3]==3)),fontface=2)+
  #Q3.in, group 1
  annotate("text", x = 5, y = c(-11.5:(-10.5-length(names(which(flow[,5]==1))))), 
           color=colorful[which(flow[,5]==1)],
           size=4,label=names(which(flow[,5]==1)),fontface=2)  + 
  #Q3.in, group 2
  annotate("text", x = 5, y = c(-8:(-9+length(names(which(flow[,5]==2))))), 
           color=colorful[which(flow[,5]==2)],
           size=4,label=names(which(flow[,5]==2)),fontface=2) + 
  #Q3.in, group 3
  annotate("text", x = 5, y = c(1:(length(names(which(flow[,5]==3))))), 
           color=colorful[which(flow[,5]==3)],
           size=4,label=names(which(flow[,5]==3)),fontface=2) + 
  #Q3.in, group 4
  annotate("text", x = 5, y = c(9.5:(8.5+length(names(which(flow[,5]==4))))), 
           color=colorful[which(flow[,5]==4)],
           size=4,label=names(which(flow[,5]==4)),fontface=2)+
  #Q4.in, group 1
  annotate("text", x = 7, y = c(-8.5:(-7.5-length(names(which(flow[,7]==1))))), 
           color=colorful[which(flow[,7]==1)],
           size=4,label=names(which(flow[,7]==1)),fontface=2) + 
  #Q4.in, group 2
  annotate("text", x = 7, y = c(-5:(-6+length(names(which(flow[,7]==2))))), 
           color=colorful[which(flow[,7]==2)],
           size=4,label=names(which(flow[,7]==2)),fontface=2) + 
  #Q4.in, group 3
  annotate("text", x = 7, y = c(3:(2+length(names(which(flow[,7]==3))))), 
           color=colorful[which(flow[,7]==3)],
           size=4,label=names(which(flow[,7]==3)),fontface=2)  + 
  #Q4.in, group 4
  annotate("text", x = 7, y = c(10.5:(9.5+length(names(which(flow[,7]==4))))), 
           color=colorful[which(flow[,7]==4)],
           size=4,label=names(which(flow[,7]==4)),fontface=2)  
dev.off()


#####################################################
# Figure 3:
# Construct discrepancy matrix
disc_table <- matrix(NA,22,22)
for (i in 1:22){
  disc_table[i,] <- rowSums(t(mapply(j=1:22,function(j)abs(flow[j,] != flow[i,]))))
}
disc_table <- disc_table/2

rownames(disc_table) <- colnames(disc_table) <- c("Mining","Construction","Durable",
                                                  "Nondurable","Wholesale",
                                                  "Retail","Transportation",
                                                  "Utility","Information",
                                                  "Finance","RealEstate",
                                                  "Professional","Management",
                                                  "Administration",
                                                  "Education","HealthCare",
                                                  "Arts","Accommodation","Other",
                                                  "Federal","State","Local")

pl_heat <- pheatmap::pheatmap(as.matrix(disc_table),
         clustering_distance_rows = as.dist(disc_table),
         clustering_distance_cols = as.dist(disc_table),
         clustering_method = "average",
         fontsize_row = 14,      
         fontsize_col = 14,     
         angle_col = "45", legend = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(5))

pdf("heatmap_PVAR.pdf", width = 13, height = 9)
print(pl_heat)
grid.text("Discrepancy (All seasons)", x = 0.93, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()



