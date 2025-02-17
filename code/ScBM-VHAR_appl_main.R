rm(list=ls())

library(superheat)
library(igraph)
library(fastRG)
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

source('ScBM-LRVAR_library.R')

#####################################################
# Data loading:
load("./logrv-2010-2017.Rdata")

colnames(logrv) <- c("SPX","FTSE","N225","GDAXI","RUT",
                     "AORD","DJI","NDX","FCHI","HSI",
                     "KS11","AEX","SSMI","IBEX","NSEI",
                     "MXX","BVSP","GSPTSE","STOXX50E","FTSEMIB.MI")

#####################################################
# Figure in supplementary
Y_rep <- logrv[,c(2,3,11,16)]

pdf("timeplot_VHAR.pdf", width = 12, height = 9)
par(mar = c(5, 4, 4, 2) + 1)
plot.ts(Y_rep,"single",main="")
dev.off()
pdf("acf_VHAR.pdf", width = 9, height = 9)
par(mar = c(5, 4, 4, 2) + 1)
plot(acf(Y_rep,lag.max=20,plot=FALSE),ylab="")
dev.off()
pdf("pacf_VHAR.pdf", width = 9, height = 9)
par(mar = c(5, 4, 4, 2) + 1)
plot(pacf(Y_rep,lag.max=20,plot=FALSE),ylab="")
dev.off()


#####################################################
# Estimation:
set.seed(123)
Yt <- t(logrv)
TT <- dim(Yt)[2]
q <- dim(Yt)[1]

Est_VHAR_ols <- VHAR_ols(Yt)
Phi_series <- Est_VHAR_ols$Phi_hat
superheat(Phi_series)

#####################################################
# Figure in supplementary
# Check the number of communities:
svd_1 <- svd(Phi_series[,1:q])
svd_2 <- svd(Phi_series[,(q+1):(2*q)])
svd_3 <- svd(Phi_series[,(2*q+1):(3*q)])


pdf("rank_VHAR.pdf", width = 13, height = 4)
par(mfrow=c(1,3))
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
dev.off()

# seasons determined:
k1 <- 3; k2 <- 4; k3 <- 3;

#####################################################
# Run the algorithm:
set.seed(123)
n_comm <- c(3,4,4,4,4,3)
alpha_hat <- alpha_cv(Est_VHAR_ols$Phi_hat,q,n_comm=n_comm,fold=5)

PisCES_output <- PisCES(as.matrix(Phi_series),q,n_comm=n_comm,s=3,p=1,alpha=alpha_hat$alpha_hat)
PVAR_1L <- blockbuster(PisCES_output$VL_bar[[1]],q,3)
PVAR_1R2L <- blockbuster(cbind(PisCES_output$VR_bar[[1]],
                               PisCES_output$VL_bar[[2]]),q,4)
PVAR_2R3L <- blockbuster(cbind(PisCES_output$VR_bar[[2]],
                               PisCES_output$VL_bar[[3]]),q,4)
PVAR_3R <- blockbuster(PisCES_output$VR_bar[[3]],q,3)

# Refine the labels:
result_list <- list(PVAR_1L,PVAR_1R2L,PVAR_2R3L,PVAR_3R)
refined_list <- Refiner(result_list)
PVAR_1L <- refined_list[[1]]
PVAR_1R2L <- refined_list[[2]]
PVAR_2R3L <- refined_list[[3]]
PVAR_3R <- refined_list[[4]]

#####################################################
# Main Figure
# Construct data flows:
flow <- mapply(x=c(1,2,2,3,3,4),function(x)refined_list[[x]]$group_1q)
rownames(flow) <- c("SPX","FTSE","N225","GDAXI","RUT",
                    "AORD","DJI","NDX","FCHI","HSI",
                    "KS11","AEX","SSMI","IBEX","NSEI",
                    "MXX","BVSP","GSPTSE","STOXX50E","FTSEMIB.MI")
colnames(flow) <- c("daily.sending","daily.receiving","weekly.sending","weekly.receiving","monthly.sending","monthly.receiving")

data_flow <- data.frame(x = c(rep("daily.sending",q),rep("daily.receiving",q),rep("weekly.sending",q),
                              rep("weekly.receiving",q),rep("monthly.sending",q),rep("monthly.receiving",q)), 
                        node = c(flow[,1],flow[,2],flow[,3],flow[,4],flow[,5],flow[,6]),
                        next_x = c(rep("daily.receiving",q),rep("weekly.sending",q),rep("weekly.receiving",q),
                                   rep("monthly.sending",q),rep("monthly.receiving",q),rep(NA,q)),
                        next_node = c(flow[,2],flow[,3],flow[,4],flow[,5],flow[,6],rep(NA,q)))

data_flow$x <- factor(data_flow$x,levels=c("daily.sending","daily.receiving",
                                           "weekly.sending","weekly.receiving",
                                           "monthly.sending","monthly.receiving"))
data_flow$next_x <- factor(data_flow$next_x,levels=c("daily.sending","daily.receiving",
                                                     "weekly.sending","weekly.receiving",
                                                     "monthly.sending","monthly.receiving"))

data_flow$node <- factor(data_flow$node,levels=c(1,2,3,4))
data_flow$next_node <- factor(data_flow$next_node,levels=c(1,2,3,4))

color20 <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Set3"))

# Plotting:
pdf("sanky_VHAR.pdf", width = 13, height = 9)
ggplot(data_flow, aes(x = x, next_x = next_x, 
                      node = node, next_node = next_node, 
                      fill = node, label = node)) +
  geom_sankey(flow.alpha = 0.3) +
  theme_sankey(base_size = 15) +
  labs(x = NULL) + labs(y = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  #daily.in, group 1
  annotate("text", x = 1, y = c(-3:(-2-length(names(which(flow[,1]==1))))), 
           color=color20[which(flow[,1]==1)],
           size=4,label=names(which(flow[,1]==1)),fontface=2) +
  #daily.in, group 2
  annotate("text", x = 1, y = c(-0.5:(-1.5+length(names(which(flow[,1]==2))))), 
           color=color20[which(flow[,1]==2)],
           size=4,label=names(which(flow[,1]==2)),fontface=2) +
  #daily.in, group 3
  annotate("text", x = 1, y = c(7:(6+length(names(which(flow[,1]==3))))), 
           color=color20[which(flow[,1]==3)],
           size=4,label=names(which(flow[,1]==3)),fontface=2) +
  #daily.out, group 1
  annotate("text", x = 2, y = c(-6:(-5-length(names(which(flow[,2]==1))))), 
           color=color20[which(flow[,2]==1)],
           size=4,label=names(which(flow[,2]==1)),fontface=2) +
  #daily.out, group 2
  annotate("text", x = 2, y = c(-0.5:(0.5-length(names(which(flow[,2]==2))))), 
           color=color20[which(flow[,2]==2)],
           size=4,label=names(which(flow[,2]==2)),fontface=2) +
  #daily.out, group 3
  annotate("text", x = 2, y = c(2.5:(1.5+length(names(which(flow[,2]==3))))), 
           color=color20[which(flow[,2]==3)],
           size=4,label=names(which(flow[,2]==3)),fontface=2) +
  #daily.out, group 4
  annotate("text", x = 2, y = c(9:(8+length(names(which(flow[,2]==4))))), 
           color=color20[which(flow[,2]==4)],
           size=4,label=names(which(flow[,2]==4)),fontface=2) +
  #weekly.out, group 1
  annotate("text", x = 4, y = c(-9:(-8-length(names(which(flow[,4]==1))))), 
           color=color20[which(flow[,4]==1)],
           size=4,label=names(which(flow[,4]==1)),fontface=2) +
  #weekly.out, group 2
  annotate("text", x = 4, y = c(-3.5:(-2.5-length(names(which(flow[,4]==2))))), 
           color=color20[which(flow[,4]==2)],
           size=4,label=names(which(flow[,4]==2)),fontface=2) +
  #weekly.out, group 3
  annotate("text", x = 4, y = c(-0.5:(-1.5+length(names(which(flow[,4]==3))))), 
           color=color20[which(flow[,4]==3)],
           size=4,label=names(which(flow[,4]==3)),fontface=2) +
  #weekly.out, group 4
  annotate("text", x = 4, y = c(6:(5+length(names(which(flow[,4]==4))))), 
           color=color20[which(flow[,4]==4)],
           size=4,label=names(which(flow[,4]==4)),fontface=2) +
  #monthly.out, group 1
  annotate("text", x = 6, y = c(-5:(-4-length(names(which(flow[,6]==1))))), 
           color=color20[which(flow[,6]==1)],
           size=4,label=names(which(flow[,6]==1)),fontface=2) +
  #monthly.out, group 2
  annotate("text", x = 6, y = c(-2.5:(-3.5+length(names(which(flow[,6]==2))))), 
           color=color20[which(flow[,6]==2)],
           size=4,label=names(which(flow[,6]==2)),fontface=2) +
  #monthly.out, group 3
  annotate("text", x = 6, y = c(5:(4+length(names(which(flow[,6]==3))))), 
           color=color20[which(flow[,6]==3)],
           size=4,label=names(which(flow[,6]==3)),fontface=2)
dev.off()



#####################################################
# Main Figure
# Construct discrepancy matrix
for (i in 1:20){
  disc_table_daily <- mapply(j=1:20,function(j)1*(flow[j,1] != flow[,1]) 
                                 + 1*(flow[j,2] != flow[,2]))
  disc_table_weekly <- mapply(j=1:20,function(j)1*(flow[j,3] != flow[,3]) 
                                 + 1*(flow[j,4] != flow[,4]))
  disc_table_monthly <- mapply(j=1:20,function(j)1*(flow[j,5] != flow[,5]) 
                                 + 1*(flow[j,6] != flow[,6]))
}

rownames(disc_table_daily) <- rownames(disc_table_weekly) <- rownames(disc_table_monthly) <- c("SPX","FTSE","N225","GDAXI","RUT",
                          "AORD","DJI","NDX","FCHI","HSI",
                          "KS11","AEX","SSMI","IBEX","NSEI",
                          "MXX","BVSP","GSPTSE","STOXX50E","FTSEMIB.MI")
colnames(disc_table_daily) <- colnames(disc_table_weekly) <- colnames(disc_table_monthly) <- c("SPX","FTSE","N225","GDAXI","RUT",
                          "AORD","DJI","NDX","FCHI","HSI",
                          "KS11","AEX","SSMI","IBEX","NSEI",
                          "MXX","BVSP","GSPTSE","STOXX50E","FTSEMIB.MI")
count_df_daily <- melt(disc_table_daily)
count_df_weekly <- melt(disc_table_weekly)
count_df_monthly <- melt(disc_table_monthly)


pdf("heatmap_VHAR.pdf", width = 13, height = 7)
heat_daily <- ggplot(count_df_daily , aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "white", breaks=c(0,1,2)) +
  scale_y_discrete(limits = rev(levels(factor(count_df_daily$Var1)))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "", fill = "Discrepancy") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  theme(legend.position = "none")
heat_weekly <- ggplot(count_df_weekly , aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "white", breaks=c(0,1,2)) +
  scale_y_discrete(limits = rev(levels(factor(count_df_weekly$Var1)))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "", fill = "Discrepancy") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  theme(legend.position = "none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
heat_monthly <- ggplot(count_df_monthly , aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "white", breaks=c(0,1,2)) +
  scale_y_discrete(limits = rev(levels(factor(count_df_monthly$Var1)))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "", fill = "Discrepancy") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  theme(legend.position = "none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend_comm <- get_legend(heat_daily + theme(legend.position = "bottom"))
grid.draw(arrangeGrob(arrangeGrob(heat_daily,heat_weekly,heat_monthly,ncol = 3),
           legend_comm,
           nrow = 2,
           heights = c(10, 1)))
dev.off()


