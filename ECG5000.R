# Librerie ---------------------------------------------------------------------
library(readr)
library(fda)
library(fda.usc)
library(fclust)
library(mclust)
library(clue)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(viridis)  
library(stats)
library(aricode)
library(dendextend)
library(cluster)
library(grc)
library(fgrc)
library(funHDDC)
library(clue)
library(factoextra)
library(fossil)
source("FKM.mahal.R")

# Dati di training  ------------------------------------------------------------
ds <- read_csv("ECG5000_TRAIN.csv", col_names = FALSE)

# Analisi descrittiva ----------------------------------------------------------

# Figura 4.1
lbls <- c("N","R-on-T PVC","PVC","SPB" ,"UB")
pct <- prop.table(table(ds$X1)) * 100
lbls <- paste(lbls, pct)
lbls <- paste(lbls,"%",sep="") 
pie(x = pct,labels = lbls, edges = 2000, col = brewer.pal(5, "Accent"), 
    border = FALSE)

# Figura 4.2
lbls <- c("N","R-on-T PVC","PVC","SPB" ,"UB")
me1 <- apply(X = ds[ds$X1 == 1,-1], MARGIN = 2, FUN = mean)
se1 <- apply(X = ds[ds$X1 == 1,-1], MARGIN = 2, FUN = sd)
me2 <- apply(X = ds[ds$X1 == 2,-1], MARGIN = 2, FUN = mean)
se2 <- apply(X = ds[ds$X1 == 2,-1], MARGIN = 2, FUN = sd)
me3 <- apply(X = ds[ds$X1 == 3,-1], MARGIN = 2, FUN = mean)
se3 <- apply(X = ds[ds$X1 == 3,-1], MARGIN = 2, FUN = sd)
me4 <- apply(X = ds[ds$X1 == 4,-1], MARGIN = 2, FUN = mean)
se4 <- apply(X = ds[ds$X1 == 4,-1], MARGIN = 2, FUN = sd)
me5 <- apply(X = ds[ds$X1 == 5,-1], MARGIN = 2, FUN = mean)
se5 <- apply(X = ds[ds$X1 == 5,-1], MARGIN = 2, FUN = sd)
# 1
DF <- data.frame(X = 1:140, Y = me1) 
plt1 <- ggplot(DF, aes(X, Y)) +                                      
  geom_line(color = "royalblue4", size = 1) + 
  geom_ribbon(aes(ymin = Y + se1, ymax = Y - se1), alpha=0.2, 
              fill = "royalblue1", color = "black", linetype = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(title = lbls[1], x ="Time (samples)", y = "Amplitude")
# 2
DF <- data.frame(X = 1:140, Y = me2) 
plt2 <- ggplot(DF, aes(X, Y)) +                                      
  geom_line(color = "royalblue4", size = 1) + 
  geom_ribbon(aes(ymin = Y + se2, ymax = Y - se2), alpha=0.2, 
              fill = "royalblue1", color = "black", linetype = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(title = lbls[2], x ="Time (samples)", y = "Amplitude")
# 3
DF <- data.frame(X = 1:140, Y = me3) 
plt3 <- ggplot(DF, aes(X, Y)) +                                      
  geom_line(color = "royalblue4", size = 1) + 
  geom_ribbon(aes(ymin = Y + se3, ymax = Y - se3), alpha=0.2, 
              fill = "royalblue1", color = "black", linetype = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(title = lbls[3], x ="Time (samples)", y = "Amplitude")
# 4
DF <- data.frame(X = 1:140, Y = me4) 
plt4 <- ggplot(DF, aes(X, Y)) +                                      
  geom_line(color = "royalblue4", size = 1) + 
  geom_ribbon(aes(ymin = Y + se4, ymax = Y - se4), alpha=0.2, 
              fill = "royalblue1", color = "black", linetype = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(title = lbls[4], x ="Time (samples)", y = "Amplitude")
# 5
DF <- data.frame(X = 1:140, Y = me5) 
plt5 <- ggplot(DF, aes(X, Y)) +                                      
  geom_line(color = "royalblue4", size = 1) + 
  geom_ribbon(aes(ymin = Y + se5, ymax = Y - se5), alpha=0.2, 
              fill = "royalblue1", color = "black", linetype = 0) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(title = lbls[5], x ="Time (samples)", y = "Amplitude")
#
ggarrange(plt1, plt2, plt3, plt4, plt5, ncol = 2, nrow = 3)


# Aggregazione Classi ----------------------------------------------------------
og <- ds$X1 # classi originarie
ds$X1[ds$X1 != 1] <- 2

# Figura 4.3
clr <- ifelse(ds$X1 == 1, "forestgreen", "firebrick1") 
matplot(t(ds[,-1]), type = "l", col = clr, xlab = "Time (samples)", 
        ylab = "Amplitude")
legend("top", legend = c("N","A"), col = c("forestgreen", "firebrick1"), 
       lty = 1, lwd = 2, ncol = 2, title = "Class", inset = 0.01, bty = "n")


# Base ottimale ----------------------------------------------------------------
l <- 10^seq(-5, 5, 0.05)
set.seed(100)
optbasis <- optim.basis(fdataobj = as.matrix(ds[,-1]), lambda = l, 
                        numbasis = seq(10,100,5))
optbasis$numbasis.opt # L
optbasis$lambda.opt # lambda
optbasis$gcv.opt # funzione di perdita

# Figura 4.4
windows()
par(oma=c(0, 0, 0, 5))
matplot(x = l, y = t(optbasis$gcv), type = 'l', xlab = expression(lambda), 
        ylab = "GCV", xlim = c(0, 100), ylim = c(0,0.25), col = viridis(19), 
        lty = 1:3)
legend(105, 0.23, bty='n', xpd=NA, legend = as.character(seq(10,100,5)), 
       ncol = 1, lty = 1:3, col = viridis(19), title = expression(L))

# Approssimazione nelle B-spline -----------------------------------------------
Bsplines <- create.bspline.basis(rangeval = c(1,140), nbasis = 90)
fdParobj <- fdPar(fdobj = Bsplines, Lfdobj = 2, lambda = 10^(-2.75))
fdata <- smooth.basis(argvals = 1:140, y = t(ds[,-1]), fdParobj = fdParobj)$fd

# FPCA -------------------------------------------------------------------------
FPCA <- pca.fd(fdata, nharm = 10)

# Figura 4.5
s <- 1:10
perc <- round(FPCA$values[s]/sum(FPCA$values)*100,2)
cumperc <- round((cumsum(FPCA$values)[s])/sum(FPCA$values)*100,2)
lab <- paste0(cumperc, "%")
lab_perc <- paste0(perc, "%")
df <- data.frame(
  x = s, 
  cum = cumperc,
  lab = lab,
  perc = perc,
  lab_perc = lab_perc
)
plt2 <- ggplot(df, aes(x, cum)) +
  geom_point(color = "royalblue1") +
  geom_text(aes(label = lab, y = cum + 1), 
            vjust = 0) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  scale_y_continuous(breaks = seq(60, 100, by = 10), limits = c(57,100)) +
  geom_line(color = "royalblue1") +
  labs(x ="Number of functional principal components", 
       y = "Cumulative percentage of explained variability") + 
  theme_bw() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt1 <- ggplot(data=df, aes(x=x, y=perc)) + 
  geom_bar(position = 'dodge', stat='identity', fill = "royalblue1") +
  geom_text(aes(label=lab_perc), vjust = -0.5) + 
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0, 60, by = 10), limits = c(0,60)) +
  labs(x ="Functional principal component", 
       y = "Percentage of explained variability") +
  theme_bw() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())
windows()
ggarrange(plt1, plt2, ncol = 2, nrow = 1)

# Figura 4.6
eigen_names <- paste("fPC", 1:10, sep = "")
windows()
par(oma=c(0, 0, 0, 5))
plot(FPCA$harmonics, col = viridis(ncol(FPCA$scores)), lty = 1:3, 
     xlab = "Time (samples)", ylab = "Amplitude")
legend(par("usr")[2]+1, par("usr")[4]-0.25, bty = "n", xpd = NA, 
       legend = eigen_names, ncol = 1, lty = 1:3, title = "Eigenfunctions",
       col = viridis(ncol(FPCA$scores)))

# Figura 4.7
cov <- var.fd(fdata)
cov_values <- eval.bifd(1:140, 1:140, cov)
windows()
persp(1:140, 1:140, cov_values, xlab = "Time", ylab = "Time", 
      zlab = "Covariance", theta = 35, phi = 20, expand = 0.7, 
      col = "royalblue1",shade = 0.09, ticktype = "detailed", 
      nticks = 6, border = "royalblue4")

# Figura 4.8
df <- data.frame(x = factor(ds$X1), y = FPCA$scores[,1])
windows()
ggplot(df, aes(x, y, fill = x)) + 
  geom_violin(trim = FALSE, color = NA) + theme(legend.position='right', t) + 
  labs(x = expression(X[1]), y = 'fPC1', fill = "Class") +  
  scale_fill_manual(values = c("forestgreen","firebrick1"),labels = c("N", "A")) +  
  geom_boxplot(width = 0.2,color="black") + coord_flip() +
  theme_bw()

# Battiti N (A) con fPC1 > -2 (fPC1 < -2)
s <- -2
round(sum(ds$X1 == 1 & FPCA$scores[,1] > s)/sum(ds$X1 == 1)*100,2)
round(sum(ds$X1 == 2 & FPCA$scores[,1] < s)/sum(ds$X1 == 2)*100,2)


# Hard clustering --------------------------------------------------------------
change_labels <- function(true, pred){
  count <- table(pred, true)
  match <- solve_LSAP(count, maximum = TRUE)
  pred <- match[pred]
  return(pred)
}

FMI <- function(true, pred){
  TP <- sum(true == 2 & pred == 2)
  TN <- sum(true == 1 & pred == 1)
  FP <- sum(true == 1 & pred == 2)
  FN <- sum(true == 2 & pred == 1)
  return(sqrt(TP/(TP+FP) * TP/(TP+FN)))
}

# KMFD 
res <- kmeans(x = t(fdata$coefs), centers = 2, nstart = 50)
cl <- change_labels(ds$X1, res$cluster)
t <- table(cl, ds$X1)
round(sum(diag(t))/sum(t) * 100, 2) # CCR
round(adjustedRandIndex(cl, ds$X1), 3) # ARI
round(NMI(cl, ds$X1), 3) # NMI
round(FMI(ds$X1, cl),3) # FMI

# KM-FPCA 
res <- kmeans(x = scale(FPCA$scores[,1:4]), centers = 2, nstart = 50)
cl <- change_labels(ds$X1, res$cluster)
t <- table(cl, ds$X1)
round(sum(diag(t))/sum(t) * 100, 2) # CCR
round(adjustedRandIndex(cl, ds$X1), 3) # ARI
round(NMI(cl, ds$X1), 3) # NMI
round(FM_index(cl, ds$X1)[1], 3) # FMI

# HAC-d0
D <- semimetric.basis(fdata, nderiv = 0, type.basis1 = Bsplines)
res <- cutree(agnes(D, method = "ward", diss = TRUE) , k = 2)
cl <- change_labels(ds$X1, res)
t <- table(cl, ds$X1)
round(sum(diag(t))/sum(t) * 100, 2) # CCR
round(adjustedRandIndex(cl, ds$X1), 3) # ARI
round(NMI(cl, ds$X1), 3) # NMI
round(FMI(ds$X1, cl),3) # FMI

# Figura 4.9
ag <- as.hclust(agnes(D, method = "ward", diss = TRUE))
windows()
fviz_dend(ag, cex = 0.5, k = 2, color_labels_by_k = TRUE, 
          k_colors = c("forestgreen", "firebrick1"), show_labels = FALSE,
          main = "", ggtheme = theme_minimal())

# Tabella 4.3
print(t)

# FigurA 4.10
pt <- round(prop.table(table(cl, og), margin = 2) * 100, 2)
windows()
bp <- barplot(pt, 
              names.arg = lbls, col = c("forestgreen", "firebrick1"), 
              border = FALSE, xlab = "Class", ylab = "Percentage", 
              ylim = c(0,115), yaxt = "n", beside = TRUE, space = c(0.3,1.2))
axis(2, at = seq(0,100,20))
legend("top", legend = c("Cluster 1", "Cluster 2"), pch = 15, 
       col = c("forestgreen", "firebrick1"), bty = "n", ncol = 2)
text(bp, pt+3, paste(pt,"%",sep=""))
rm(pt)

# Figura 4.11
windows()
par(mfrow = c(1,2))
plot(fdata[cl == 1, ], col = 'lightgreen', lty = 1, href = FALSE, 
     xlab = "Time (samples)", ylab = "Amplitude", main = "Cluster 1",
     ylim = c(-6, 4))
lines(mean.fd(fdata), lty = 2)
lines(mean.fd(fdata[cl == 1,]), col = "forestgreen", lty = 1, lwd = 3)
legend("top", legend = c("Mean amplitude","Cluster mean amplitude"), 
       col = c(1, "forestgreen"), lty = c(2,1), lwd = c(1,3), ncol = 1, 
       inset = 0.01, bty = "n", cex = 0.8)
title(main = "Cluster 1")
plot(fdata[cl == 2, ], col = 'lightpink', lty = 1, href = FALSE,
     xlab = "Time (samples)", ylab = "Amplitude", main = "Cluster 2",
     ylim = c(-6,4), main = "Cluster 2")
lines(mean.fd(fdata), lty = 2)
lines(mean.fd(fdata[cl == 2,]), col = "firebrick1", lty = 1, lwd = 3)
legend("top", legend = c("Mean amplitude","Cluster mean amplitude"), 
       col = c(1, "firebrick1"), lty = c(2,1), lwd = c(1,3), ncol = 1, 
       inset = 0.01, bty = "n", cex = 0.8)
title(main = "Cluster 2")
par(mfrow = c(1,1))

# HAC-d1
D <- semimetric.basis(fdata, nderiv = 1, type.basis1 = Bsplines)
res <- cutree(agnes(D, method = "ward", diss = TRUE) , k = 2)
cl <- change_labels(ds$X1, res)
t <- table(cl, ds$X1)
round(sum(diag(t))/sum(t) * 100, 2) # CCR
round(adjustedRandIndex(cl, ds$X1), 3) # ARI
round(NMI(cl, ds$X1), 3) # NMI
round(FMI(ds$X1, cl),3) # FMI

# FPCK
res <- FPCK(X = as.matrix(ds[,-1]), N.comp = 2, N.clust = 2, N.random = 50, 
            Basis = 'Bsp', par.Bsp = c(98, 4), lambda = 10^{-2.75})
cl <- change_labels(ds$X1, res$cluster)
t <- table(cl, ds$X1)
round(sum(diag(t))/sum(t) * 100, 2) # CCR
round(adjustedRandIndex(cl, ds$X1), 3) # ARI
round(NMI(cl, ds$X1), 3) # NMI
round(FMI(ds$X1, cl),3) # FMI

# Fuzzy Clustering -------------------------------------------------------------

set.seed(100)
# FKMFD
res1 <- fclust::FKM(X = t(fdata$coefs), k = 2:5, RS = 10)
# FKMedFD
res2 <- FKM.med(X = t(fdata$coefs), k = 2:5, RS = 50)
# FKM-FPCA
res3 <- fclust::FKM(X = FPCA$scores[,1:4], k = 2:5, RS = 10)
# FKMed-FPCA
res4 <- FKM.med(X = FPCA$scores[,1:4], k = 2:5, RS = 10)
# FKM-FPCA + Mahalanobis
set.seed(100)
res5 <- FKM.mahal.opt(X = FPCA$scores[,1:4], k = 2, prior = ds$X1, 
                      ntau = 10000, alpha = 0.1)

# Figura 4.12
tmp <- viridis(5)
windows()
par(oma=c(0, 0, 1, 0))
plot(2:5, res1$criterion, type = 'o', pch = 16, col = tmp[1], xaxt = "n",
     xlab = expression(k), ylab = expression(FS(k)), ylim = c(0.38,0.8))
lines(2:5, res2$criterion, type = 'o', pch = 16, col = tmp[2])
lines(2:5, res3$criterion, type = 'o', pch = 16, col = tmp[3])
lines(2:5, res4$criterion, type = 'o', pch = 16, col = tmp[4])
lines(2, SIL.F(FPCA$scores[,1:4], res5$FKM$U), type = "p", pch = 16, col = tmp[5])
axis(1, at = 2:5)
legend(par("usr")[1]-0.18, par("usr")[3]+0.51, bty = "n", xpd = NA, 
       legend = c("FKMFD", "FKMedFD", "FKM-FPCA", "FKMed-FPCA", 
                  "FKM-FPCA + Mahalanobis"), 
       ncol = 3, pch = 16, col = tmp, cex = 0.9)
rm(tmp)


# FS(K)
round(rbind(max(res1$criterion), max(res2$criterion), max(res3$criterion), 
            max(res4$criterion), SIL.F(FPCA$scores[,1:4], res5$FKM$U)), 3)
# ARI.F
round(rbind(ARI.F(ds$X1, res1$U), ARI.F(ds$X1, res2$U), 
            ARI.F(ds$X1, res3$U), ARI.F(ds$X1, res4$U),
            ARI.F(ds$X1, res5$FKM$U)), 3)
# AI (%)
u.max <- apply(res5$FKM$U, MARGIN = 1, FUN = max)
round(rbind(mean(res1$clus[,2] < .55), mean(res2$clus[,2] < .55),
            mean(res3$clus[,2] < .55), mean(res4$clus[,2] < .55),
            mean(u.max < .55)) * 100, 2)
# d(C^*, C)
round(rbind(dist_cluster(ds$X1, res1$U)$d, dist_cluster(ds$X1, res3$U)$d, 
            dist_cluster(ds$X1, res4$U)$d, res5$d), 3)

# Figura 4.13
pchs <- ifelse(ds$X1 == 1, 21, 24)
windows()
par(oma=c(2, 0, 0, 0))
par(mfrow = c(1,2))
plot(FPCA$scores[,1], FPCA$scores[,2], 
     col = ifelse(change_labels(ds$X1, res3$clus[,1]) == 1, 
                  "forestgreen", "firebrick1"), pch = pchs,
     xlab = "fPC1", ylab = "fPC2", main = "FKM-FPCA")
for(i in unique(res3$clus[,1])) {
  cluster_points <- FPCA$scores[change_labels(ds$X1, res3$clus[,1]) == i, 1:2]
  hull_indices <- chull(cluster_points)
  hull_indices <- c(hull_indices, hull_indices[1])
  lines(cluster_points[hull_indices, ], col = unique(clr)[i])
}
lines(mean(FPCA$scores[change_labels(ds$X1, res3$clus[,1]) == 1,1]),
      mean(FPCA$scores[change_labels(ds$X1, res3$clus[,1]) == 1,2]), type = "p",
      pch = 8, cex = 2)
lines(mean(FPCA$scores[change_labels(ds$X1, res3$clus[,1]) == 2,1]),
      mean(FPCA$scores[change_labels(ds$X1, res3$clus[,1]) == 2,2]), type = "p",
      pch = 8, cex = 2)

plot(FPCA$scores[,1], FPCA$scores[,2], 
     col = ifelse(res5$FKM$clus == 1, "forestgreen", "firebrick1"), pch = pchs,
     xlab = "fPC1", ylab = "fPC2", main = "FKM-FPCA + Mahalanobis")
for(i in unique(res5$FKM$clus)) {
  cluster_points <- FPCA$scores[res5$FKM$clus == i, 1:2]
  hull_indices <- chull(cluster_points)
  hull_indices <- c(hull_indices, hull_indices[1])
  lines(cluster_points[hull_indices, ], col = unique(clr)[i])
}
lines(mean(FPCA$scores[res5$FKM$clus == 1,1]),
      mean(FPCA$scores[res5$FKM$clus == 1,2]), type = "p",
      pch = 8, cex = 2)
lines(mean(FPCA$scores[res5$FKM$clus == 2,1]),
      mean(FPCA$scores[res5$FKM$clus == 2,2]), type = "p",
      pch = 8, cex = 2)
par(mfrow = c(1,1))
legend(par("usr")[1]+5, par("usr")[3]-4, bty = "n", xpd = NA, 
       legend = c("N", "A"), pch = c(21,24), ncol = 2, title = "Class")
legend(par("usr")[1]+8.5, par("usr")[3]-4, bty = "n", xpd = NA, 
       legend = c("1", "2"), lty = 1, col = unique(clr), ncol = 2, 
       title = "Cluster", seg.len = 0.9)


# Figura 4.14
windows()
hist(res5$FKM$U, border = FALSE, col = "royalblue1", xlab = expression(u[ik]), 
     main = "", ylim = c(0,500))

# Figura 4.15
idx <- which(u.max < 0.55)
tmp <- c(expression(x[32]), expression(x[416]), expression(x[494]), 
         expression(x[500]))
windows()
par(mfrow = c(2,2))
for(i in 1:length(idx)){
  plot(fdata[idx[i],], href = FALSE, xlab = "Time (samples)", 
       ylab = "Amplitude", ylim = c(-4, 3))
  lines(mean.fd(fdata[res5$FKM$clus==1,]), col = "forestgreen", lwd = 2)
  lines(mean.fd(fdata[res5$FKM$clus==2,]), col = "firebrick1", lwd = 2)
  title(main = tmp[i])
}

# misure di bontà come partizione hard
cl <- res5$FKM$clus
t <- table(cl, ds$X1)
round(sum(diag(t))/sum(t) * 100, 2) # CCR
round(adjustedRandIndex(cl, ds$X1), 3) # ARI
round(NMI(cl, ds$X1), 3) # NMI
round(FMI(ds$X1, cl),3) # FMI


# Model-based clustering -------------------------------------------------------
set.seed(100)
res <- funHDDC(data = fdata, K = 2,
               model = c("AkjBkQkDk", "AkjBQkDk",
                         "AkBkQkDk", "ABkQkDk", "AkBQkDk", "ABQkDk"))

# misure di bontà come partizione hard
cl <- change_labels(ds$X1, res$class)
t <- table(cl, ds$X1)
round(sum(diag(t))/sum(t) * 100, 2) # CCR
round(adjustedRandIndex(cl, ds$X1), 3) # ARI
round(NMI(cl, ds$X1), 3) # NMI
round(FMI(ds$X1, cl),3) # FMI


# Confronto con la letteratura -------------------------------------------------

# unisco training e test
tmp <- read_csv("ECG5000_TEST.csv", col_names = FALSE)
ds <- rbind(ds, tmp)
rm(tmp)

# Ricalcolo base ottimale ------------------------------------------------------
l <- 10^seq(-5, 5, 0.05)
set.seed(100)
optbasis <- optim.basis(fdataobj = as.matrix(ds[,-1]), lambda = l, 
                        numbasis = seq(10,100,5))
optbasis$numbasis.opt # L
optbasis$lambda.opt # lambda
optbasis$gcv.opt # funzione di perdita

# Approssimazione nelle B-spline
Bsplines <- create.bspline.basis(rangeval = c(1,140), nbasis = 90)
fdParobj <- fdPar(fdobj = Bsplines, Lfdobj = 2, lambda = 10^(-2.70))
fdata <- smooth.basis(argvals = 1:140, y = t(ds[,-1]), 
                      fdParobj = fdParobj)$fd

# FPCA
FPCA <- pca.fd(fdata, nharm = 10)

# HAC-d0
D <- semimetric.basis(fdata, nderiv = 0, type.basis1 = Bsplines)
cl <- cutree(agnes(D, method = "ward", diss = TRUE) , k = 5)
round(mean(rand.index(cl, ds$X1)), 4) # RI
round(mean(NMI(cl, ds$X1)), 4) # NMI

# FKM-FPCA + Mahalanobis
set.seed(100)
res <- FKM.mahal.opt(X = FPCA$scores[,1:4], k = 5, prior = ds$X1, ntau = 1000, 
                     alpha = 0.1)
round(rand.index(ds$X1, res$FKM$clus), 4)
round(NMI(ds$X1, res$FKM$clus), 4)
