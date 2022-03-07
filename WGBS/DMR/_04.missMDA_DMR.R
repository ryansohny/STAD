library(missMDA)

######### ALL DNA methylation #########
mat <- as.matrix(read.table("/mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/DMR/DMR_min55/DMR_abs15_ALL.txt", sep="\t", header=TRUE, row.names=1))

t_mat <- t(mat)
n <- rep(c("Normal"), each=84)
t <- rep(c("Tumor"), each=84)
Type <- c(n,t)
t_mat_type <- data.frame(t_mat, Type)

nb <- estim_ncpPCA(t_mat_type[1:2197], ncp.max=50) # 총 2197개의 DMR (abs(T-N) > 15%)
nb$ncp
# nb$ncp를 통해서 몇 개의 ncp를 사용해야하는 지 확인
res.comp <- imputePCA(t_mat_type[1:2197], ncp=28)
t_mat_imputed_type <- data.frame(res.comp$completeObs, Type)
write.csv(t_mat_imputed_type, "DMR_mat_imputed.csv")

######### Only Normal DNA methylation #########
library(missMDA)
mat <- as.matrix(read.table("/mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/DMR/DMR_min55/DMR_abs15_Normal.txt", sep="\t", header=TRUE, row.names=1))
n_mat <- data.frame(t(mat))
nb <- estim_ncpPCA(n_mat[1:2197], ncp.max=50)
res.comp <- imputePCA(n_mat[1:2197], ncp=nb$ncp)
n_mat_imputed <- data.frame(res.comp$completeObs)
nb$ncp # 6
write.csv(n_mat_imputed, "DMR_mat_normal_imputed.csv")

# Only Tumor DNA methylation
library(missMDA)
mat <- as.matrix(read.table("/mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/DMR/DMR_min55/DMR_abs15_Tumor.txt", sep="\t", header=TRUE, row.names=1))

t_mat <- data.frame(t(mat))

nb <- estim_ncpPCA(t_mat[1:2197], ncp.max=50) # 총 2197개의 DMR (abs(T-N) > 15%)
nb$ncp #14
# nb$ncp를 통해서 몇 개의 ncp를 사용해야하는 지 확인
res.comp <- imputePCA(t_mat[1:2197], ncp=nb$ncp)
t_mat_imputed <- data.frame(res.comp$completeObs)
write.csv(t_mat_imputed, "DMR_mat_tumor_imputed.csv")

######### Tumor-Normal DNA methylation #########
library(missMDA)
mat <- as.matrix(read.table("/mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/DMR/DMR_min55/DMR_abs15_Tumor-Normal.txt", sep="\t", header=TRUE, row.names=1))
tn_mat <- data.frame(t(mat))
nb <- estim_ncpPCA(tn_mat[1:2197], ncp.max=50)
res.comp <- imputePCA(tn_mat[1:2197], ncp=nb$ncp)
tn_mat_imputed <- data.frame(res.comp$completeObs)
nb$ncp # 8
write.csv(tn_mat_imputed, "DMR_mat_tumor-normal_imputed.csv")
