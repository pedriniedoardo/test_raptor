# AIM ---------------------------------------------------------------------
# build the reference to estimate the ageing variable

# -------------------------------------------------------------------------
# source("scr/R/00_tutorial_load_libraries.R")

# read in the data --------------------------------------------------------
load("../../data/test/dschen2016.RData")

fnpx <- "aging_hom_brain_n500_"
ftype <- "pdf"
# read in the test data

# wrangling ---------------------------------------------------------------
# quantile normalization of the matrices
g47 <- log1p(limma::normalizeBetweenArrays(dschen2016$g47, method = "quantile"))
g11 <- log1p(limma::normalizeBetweenArrays(dschen2016$g11, method = "quantile"))

# check the distributions
boxplot(dschen2016$g47[,1:10])
boxplot(g47[,1:10])

# compute monotony
mon_g47 <- abs(apply(g47, 1, cor, y=dschen2016$p47$age, method="spearman"))
mon_g11 <- abs(apply(g11, 1, cor, y=dschen2016$p11$age, method="spearman"))

sel47 <- mon_g47>sqrt(.25)
sel11 <- mon_g11>sqrt(.25)

pca_47 <- summary(prcomp(t(g47[sel47,]), center = T, scale. = F, rank. = 20))
pca_11 <- summary(prcomp(t(g11[sel11,]), center = T, scale. = F, rank. = 20))

# stats function
sts <- function(x, y, sel=T){
  x <- x[sel]
  y <- y[sel]
  
  lmi <- lm(y~x)
  res <- residuals(lmi)
  return(list(lm=lmi, 
              r2 = summary(lmi)$adj.r.squared, 
              rho = cor(x,y, method="spearman"),
              rmse = sqrt(sum(res**2)/length(res))))
}

# build the reference -----------------------------------------------------
# build references with all the samples
selsamp <- c(T,F)

# build the model for a subset of g47 samples
m47 <- ge_im(X = g47[sel47,selsamp], p = dschen2016$p47[selsamp,],
             formula = "X ~ s(age, bs = 'cr')",
             dim_red = "pca", nc = 1)

# build an even rank of 500 steps in ages from the min to the max
n.inter <- 500
ndat <- data.frame(age = seq(min(dschen2016$p47$age),
                             max(dschen2016$p47$age), l=n.inter))

# list containint both predicted expression values and relative ages
r47 <- list(interpGE=predict(m47, ndat),
            time.series=ndat$age)
boxplot(r47$interpGE[,1:10])

# do the same as above for g11 dataset
m11 <- ge_im(X = g11[sel11,selsamp], p = dschen2016$p11[selsamp,],
             formula = "X ~ s(age, bs = 'cr')",
             dim_red = "pca", nc = 1)

ndat2 <- data.frame(age=seq(min(dschen2016$p11$age),
                            max(dschen2016$p11$age),l=n.inter))

r11 <- list(interpGE=predict(m11, ndat2),
            time.series=ndat2$age)



# stage samples on BA47 & BA11 references
# use r47 as reference
ae_g47_r47 <- ae(g47, r47$interpGE, r47$time.series)
ae_g11_r47 <- ae(g11, r47$interpGE, r47$time.series)

# use r11 as reference
ae_g47_r11 <- ae(g47, r11$interpGE, r11$time.series)
ae_g11_r11 <- ae(g11, r11$interpGE, r11$time.series)

# add the estimated ages to the original dataset as metadata
dschen2016$p47$ae47 <- ae_g47_r47$age.estimates[, 1]
dschen2016$p11$ae47 <- ae_g11_r47$age.estimates[, 1]
dschen2016$p47$ae11 <- ae_g47_r11$age.estimates[, 1]
dschen2016$p11$ae11 <- ae_g11_r11$age.estimates[, 1]

# plotting ----------------------------------------------------------------
# plots
fig_custom(figname = paste0(fnpx, "BA47"),
           path = "../../out/image/",
           output = ftype, fig.width = 5*3, fig.height = 5)
## BA47
par(mfrow = c(1,3), pty ='s', bty='l')
# ae47 vs chron
x <- "age"
y <- "ae47"
plot(dschen2016$p47[, c(x,y)], main = paste("BA47 -",y, "v",x), 
     col = transp(col.palette[selsamp +1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p47[,x], dschen2016$p47[,y], sel=selsamp)
sv <- sts(dschen2016$p47[,x], dschen2016$p47[,y], sel=!selsamp)
# test to see if ref. samples and val. samples have bias in age-ae
# t.test(residuals(sv$lm), residuals(sr$lm)) # non-signif
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[2])
abline(sv$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref,\nR2=",round(sr$r2, 3), "\np=", round(sr$rho, 3)),
                 paste0("val,\nR2=",round(sv$r2, 3), "\np=", round(sv$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[2:1])

# ae11 vs chron
x <- "age"
y <- "ae11"
plot(dschen2016$p47[, c(x,y)], main = paste("BA47 -",y, "v",x), 
     col = transp(col.palette[1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p47[,x], dschen2016$p47[,y], sel=T)
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref, R2=",round(sr$r2, 3), "\np=", round(sr$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[1])

# ae47 vs ae11
x <- "ae47"
y <- "ae11"
plot(dschen2016$p47[, c(x,y)], main = paste("BA47 -",y, "v",x), 
     col = transp(col.palette[selsamp +1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p47[,x], dschen2016$p47[,y], sel=T)
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref, R2=",round(sr$r2, 3), "\np=", round(sr$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[1])

dev.off()


## BA11
fig_custom(paste0(fnpx, "BA11"), 
           path = "../../out/image/",
           output = ftype, fig.width = 5*3, fig.height = 5)
par(mfrow = c(1,3), pty ='s', bty='l')
# ae11 vs chron
x <- "age"
y <- "ae11"
plot(dschen2016$p11[, c(x,y)], main = paste("BA11 -",y, "v",x), 
     col = transp(col.palette[selsamp +1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p11[,x], dschen2016$p11[,y], sel=selsamp)
sv <- sts(dschen2016$p11[,x], dschen2016$p11[,y], sel=!selsamp)
# test to see if ref. samples and val. samples have bias in age-ae
# t.test(residuals(sv$lm), residuals(sr$lm)) # non-signif

abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[2])
abline(sv$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref,\nR2=",round(sr$r2, 3), "\np=", round(sr$rho, 3)),
                 paste0("val,\nR2=",round(sv$r2, 3), "\np=", round(sv$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[2:1])

# ae47 vs chron
x <- "age"
y <- "ae47"
plot(dschen2016$p11[, c(x,y)], main = paste("BA11 -",y, "v",x), 
     col = transp(col.palette[1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p11[,x], dschen2016$p11[,y], sel=T)
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref, R2=",round(sr$r2, 3), "\np=", round(sr$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[1])

# ae47 vs ae11
x <- "ae11"
y <- "ae47"
plot(dschen2016$p11[, c(x,y)], main = paste("BA11 -",y, "v",x), 
     col = transp(col.palette[selsamp +1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p11[,x], dschen2016$p11[,y], sel=T)
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref, R2=",round(sr$r2, 3), "\np=", round(sr$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[1])

dev.off()

# -------------------------------------------------------------------------
# compare age estimates of 47 and 11 tissues in a patient on the 47 ref.
its <- intersect(dschen2016$p47$snum, dschen2016$p11$snum)
i47_11 <- match(its, dschen2016$p47$snum)
i11_47 <- match(its, dschen2016$p11$snum)
plot(dschen2016$p47[i47_11,"ae47"], dschen2016$p11[i11_47,"ae47"])
sr <- sts(dschen2016$p47[i47_11,"ae47"], dschen2016$p11[i11_47,"ae47"], sel=T)
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref, R2=",round(sr$r2, 3), "\np=", round(sr$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[1])

# -------------------------------------------------------------------------
# build references with other 50% samples only 200 samples
selsamp2 <- c(F,T)

m47_alt <- ge_im(X = g47[sel47, selsamp2], p = dschen2016$p47[selsamp2,],
             formula = "X ~ s(age, bs = 'cr')",
             dim_red = "pca", nc = 1)
ndat_alt <- data.frame(age=seq(min(dschen2016$p47$age),
                               max(dschen2016$p47$age), l=200))
r47_alt <- list(interpGE=predict(m47_alt, ndat_alt), time.series=ndat_alt$age)

m11_alt <- ge_im(X = g11[sel11,selsamp2], p = dschen2016$p11[selsamp2,],
             formula = "X ~ s(age, bs = 'cr')",
             dim_red = "pca", nc = 1)
ndat2_alt <- data.frame(age=seq(min(dschen2016$p11$age),
                                max(dschen2016$p11$age), l=200))
r11_alt <- list(interpGE=predict(m11_alt, ndat2_alt), time.series=ndat2_alt$age)


# stage samples on BA47 & BA11 references
ae_g47_r47_alt <- ae(g47, r47_alt$interpGE, r47_alt$time.series)
ae_g11_r47_alt <- ae(g11, r47_alt$interpGE, r47_alt$time.series)

ae_g47_r11_alt <- ae(g47, r11_alt$interpGE, r11_alt$time.series)
ae_g11_r11_alt <- ae(g11, r11_alt$interpGE, r11_alt$time.series)

dschen2016$p47$ae47_2 <- ae_g47_r47_alt$age.estimates[, 1]
dschen2016$p11$ae47_2 <- ae_g11_r47_alt$age.estimates[, 1]
dschen2016$p47$ae11_2 <- ae_g47_r11_alt$age.estimates[, 1]
dschen2016$p11$ae11_2 <- ae_g11_r11_alt$age.estimates[, 1]

# -------------------------------------------------------------------------
# plots
fig_custom(paste0(fnpx, "BA47_2"),
           path = "../../out/image/",
           output = ftype, fig.width = 2.5*3, fig.height = 3)
## BA47
par(mfrow = c(1,3), pty ='s', bty='l')
# ae47 vs chron
x <- "age"
y <- "ae47_2"
plot(dschen2016$p47[, c(x,y)], main = paste("BA47 -",y, "v",x), 
     col = transp(col.palette[selsamp2 +1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p47[,x], dschen2016$p47[,y], sel=selsamp2)
sv <- sts(dschen2016$p47[,x], dschen2016$p47[,y], sel=!selsamp2)
# test to see if ref. samples and val. samples have bias in age-ae
# t.test(residuals(sv$lm), residuals(sr$lm)) # non-signif
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[2])
abline(sv$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref,\nR2=",round(sr$r2, 3), "\np=", round(sr$rho, 3)),
                 paste0("val,\nR2=",round(sv$r2, 3), "\np=", round(sv$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[2:1])

# ae11 vs chron
x <- "age"
y <- "ae11_2"
plot(dschen2016$p47[, c(x,y)], main = paste("BA47 -",y, "v",x), 
     col = transp(col.palette[1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p47[,x], dschen2016$p47[,y], sel=T)
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref, R2=",round(sr$r2, 3), "\np=", round(sr$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[1])

# ae47 vs ae11
x <- "ae47_2"
y <- "ae11_2"
plot(dschen2016$p47[, c(x,y)], main = paste("BA47 -",y, "v",x), 
     col = transp(col.palette[selsamp2 +1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p47[,x], dschen2016$p47[,y], sel=T)
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref, R2=",round(sr$r2, 3), "\np=", round(sr$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[1])

dev.off()


## BA11
fig_custom(paste0(fnpx, "BA11_2"),
           path = "../../out/image/",
           output = ftype, fig.width = 2.5*3, fig.height = 3)
par(mfrow = c(1,3), pty ='s', bty='l')
# ae11 vs chron
x <- "age"
y <- "ae11_2"
plot(dschen2016$p11[, c(x,y)], main = paste("BA11 -",y, "v",x), 
     col = transp(col.palette[selsamp2 +1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p11[,x], dschen2016$p11[,y], sel=selsamp2)
sv <- sts(dschen2016$p11[,x], dschen2016$p11[,y], sel=!selsamp2)
# test to see if ref. samples and val. samples have bias in age-ae
# t.test(residuals(sv$lm), residuals(sr$lm)) # non-signif

abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[2])
abline(sv$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref,\nR2=",round(sr$r2, 3), "\np=", round(sr$rho, 3)),
                 paste0("val,\nR2=",round(sv$r2, 3), "\np=", round(sv$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[2:1])

# ae47 vs chron
x <- "age"
y <- "ae47_2"
plot(dschen2016$p11[, c(x,y)], main = paste("BA11 -",y, "v",x), 
     col = transp(col.palette[1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p11[,x], dschen2016$p11[,y], sel=T)
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref, R2=",round(sr$r2, 3), "\np=", round(sr$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[1])

# ae47 vs ae11
x <- "ae11_2"
y <- "ae47_2"
plot(dschen2016$p11[, c(x,y)], main = paste("BA11 -",y, "v",x), 
     col = transp(col.palette[selsamp2 +1]), lwd=2, xaxt="n", yaxt="n"); twoTicks()
sr <- sts(dschen2016$p11[,x], dschen2016$p11[,y], sel=T)
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref, R2=",round(sr$r2, 3), "\np=", round(sr$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[1])

dev.off()

# compare age estimates of 47 and 11 tissues in a patient on the 47 ref.
its2 <- intersect(dschen2016$p47$snum, dschen2016$p11$snum)
i47_11_2 <- match(its2, dschen2016$p47$snum)
i11_47_2 <- match(its2, dschen2016$p11$snum)
plot(dschen2016$p47[i47_11_2,"ae47_2"], dschen2016$p11[i11_47_2,"ae47_2"])
sr <- sts(dschen2016$p47[i47_11_2,"ae47_2"], dschen2016$p11[i11_47_2,"ae47_2"], sel=T)
abline(a=0, b=1, lty=2)
abline(sr$lm, lwd=2, col = col.palette[1])
legend("bottomright", 
       legend=c( paste0("ref, R2=",round(sr$r2, 3), "\np=", round(sr$rho, 3))),
       bty='n', lwd=2, lty=1, col=col.palette[1])
