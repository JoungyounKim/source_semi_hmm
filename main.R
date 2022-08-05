

rm(list=ls())
origin<-"C:\\Users\\megam\\Desktop\\KJY\\"
direc_in<-origin
direc_out<-origin

source(paste(origin, "HMM-FDR\\","rdata1.hmm.R.txt",sep=""))
source(paste(origin, "HMM-FDR\\","rdata.hmm.R.txt",sep=""))
source(paste(origin, "HMM-FDR\\","bwfw1.hmm.R.txt",sep=""))
source(paste(origin, "HMM-FDR\\","bwfw.hmm.R.txt",sep=""))
source(paste(origin, "HMM-FDR\\","em1.hmm.R.txt",sep=""))
source(paste(origin, "HMM-FDR\\","em.hmm.R.txt",sep=""))
source(paste(origin, "HMM-FDR\\","mt.hmm.R.txt",sep=""))
source(paste(origin, "HMM-FDR\\","mt.gw.R.txt",sep=""))
source(paste(origin, "HMM-FDR\\","mt.bh.R.txt",sep=""))

library("Iso")
library("fdrtool")
library("spatstat")
setwd(direc_in)
source("ftn.R")

L<-1
data_original<-NULL
data<-NULL
## parameters of f_0 for sun&cai
mu0<-0 
sigma0<-0
outfile<-NULL
### Read Data ###

###########################################################################
####  ILI Epidemic data                                                 ####
####  Downloaded from https://websenti.u707.jussieu.fr/sentiweb/?page=table
###########################################################################

ili<-read.csv(paste(direc_in, "ILI-France.csv",sep=""), skip=1, as.is=T)
ili$inc<-as.numeric(ili$inc)
ili$inc100<-as.numeric(ili$inc100)
ili$week_new<-nrow(ili):1
tmpWhich<-which(is.na(ili$inc))

	

data<-log(1+ili$inc100)
data<-data[length(data):1]
	
data_original<-data
	
if(L==1){
	mu0<-2.50
	sigma0<-0.81
	outfile<-"ILI_L1"

}
if(L==2){
	mu0<-2.37
	sigma0<-0.76
	outfile<-"ILI_L2"
	}

## Normalizing ##
data<-(data-mu0)/sigma0
data<-data[!is.na(data)]


####################################
#####         Estimation       #####
####################################
#data<-data[1:1000]

n<-length(data)
q=0.1# FDR level

##################
#### Peeling  ####
##################
est<-estimate(data=data, 1, 1, 0.5, 0.5)

res_peeling<-mt.hmm(1-est$p1, q)


######################
#### Sun & Cai    ####
######################

em.res1<-NULL
em.res2<-NULL
em.res3<-NULL
res_EM1<-NULL
res_EM2<-NULL
res_EM3<-NULL


####################
# Sun&Cai EM algorithm
####################
L=1
em.res1<-em1.hmm(data, maxiter=500)      
res_EM1<-mt.hmm(em.res1$lf, q)
L=2
em.res2<-em.hmm(data,L=L, maxiter=500)
res_EM2<-mt.hmm(em.res2$lf, q)

L=3
em.res3<-em.hmm(data,L=L, maxiter=500)
res_EM3<-mt.hmm(em.res3$lf, q)

em_res_mat<-matrix(0, nrow=7, ncol=2)

em_res_mat[1,]<-em.res1$f1 #mu1, sigma1

em_res_mat[2:3,]<-as.matrix(em.res2$f1) #mu1, sigma1,mu2, sigma2
em_res_mat[4,]<-as.matrix(em.res2$pc)

em_res_mat[5:7,]<-as.matrix(em.res3$f1) #mu1, sigma1,mu2, sigma2,mu3, sigma3

write.table(em_res_mat, paste(direc_out, outfile, "_em_res.txt",sep=""),sep=",", col.names=c("mu", "sd"), row.names=F, quote=F, append=F )



######################
##  Method by LIU
##########################

print("### LIU_02  ###")
est_liu_02<-liu(data=data, 1, 1, 0.5, 0.5, lambda=0.2)
res_liu02<-mt.hmm(1-est_liu_02$p1, q)

print("### LIU_05  ###")
est_liu_05<-liu(data=data, 1, 1, 0.5, 0.5, lambda=0.5)
res_liu05<-mt.hmm(1-est_liu_05$p1, q)

print("### LIU_08  ###")
est_liu_08<-liu(data=data, 1, 1, 0.5, 0.5, lambda=0.8)
res_liu08<-mt.hmm(1-est_liu_08$p1, q)


#############################
#### Save Resultes in TXT ###
#############################

tmp<-data.frame(round(cbind(data, 1-est$p1, res_peeling$de, em.res1$lf, res_EM1$de,  em.res2$lf, res_EM2$de, em.res3$lf, res_EM3$de, 
1-est_liu_02$p1, res_liu02$de,1-est_liu_05$p1, res_liu05$de,1-est_liu_08$p1, res_liu08$de ), digits=10))
names(tmp)<-c("data","p0_peeling","de_peeling","p0_sc_1","res_sc_1","p0_sc_2","res_sc_2","p0_sc_3","res_sc_3", "p0_liu_02", "res_liu_02","p0_liu_05", "res_liu_05","p0_liu_08", "res_liu_08" )
write.table(tmp, paste(direc_out, outfile, "_res.txt",sep=""), col.names=T, row.names=F, quote=F, append=F)
##
tmp_pdf<-est$pdf_kernel
write.table(cbind(tmp_pdf$x, tmp_pdf$y), paste(direc_out, outfile, "_pdf.txt",sep=""), col.names=c("x", "y"), row.names=F, quote=F, append=F)

#############################
#### Save Results in PLOTS ###
#############################

pdf(paste(direc_out, outfile, "_P1.pdf",sep=""), width=10, height=10)

tmpLim<-range(c(est$p1,1-em.res1$lf, 1-em.res2$lf, 1-em.res3$lf ))
par(mfrow=c(5,1))
par( mar=c(1,2,2,1))
plot(data, type="l", main="Plot of Pr(Z=1|Data)", lty=1 )
par( mar=c(1,2,2,1))

plot(est$p1, type="l", lty=1,ylim=tmpLim)
par( mar=c(1,2,2,1))

plot(1-em.res1$lf, type="l", lty=1,ylim=tmpLim)
par( mar=c(1,2,2,1))

plot(1-em.res2$lf, type="l", lty=1,ylim=tmpLim)
par( mar=c(1,2,2,1))
plot(1-em.res3$lf,type="l",  lty=1, xlab="Index",ylim=tmpLim)

dev.off()
pdf(paste(direc_out, outfile, "_Z.pdf",sep=""), width=10, height=10)
tmpLim<-range(c(res_peeling$de,res_EM1$de, res_EM2$de, res_EM3$de ))
par(mfrow=c(5,1))
par( mar=c(1,2,2,1))
plot(data, type="l", main=expression(hat(Z)), lty=1)
par( mar=c(1,2,2,1))
plot(res_peeling$de, type="l", lty=1, ylim=tmpLim)
par( mar=c(1,2,2,1))
plot(res_EM1$de, type="l",lty=1, ylim=tmpLim)
par( mar=c(1,2,2,1))
plot(res_EM2$de, type="l",lty=1, ylim=tmpLim)
par( mar=c(1,2,2,1))
plot(res_EM3$de, type="l",lty=1, ylim=tmpLim)

 
dev.off()

###################################################
############## PDF output  ########################
###################################################
outfile<-NULL
my_y_lim<-NULL
my_x_at<-NULL
my_y_at<-NULL
	outfile<-"ILI_hist_ver_4.pdf"
	my_y_lim<-c(0, 0.5)
	my_x_at<--3.6
	my_y_at<-0.47
	

outfile4<-paste(direc_out, outfile, sep="")
pdf(outfile4, width=9, height=4)
tmpMatrix<-diag(1, nrow=2, ncol=2)
for(i in 1:300){
	tmpMatrix=tmpMatrix%*%matrix(c(est$a0, 1-est$a0, 1-est$a1, est$a1), byrow=T, nrow=2)
}

my_cex_axis=1.5
my_cex_legend=1.7
layout(mat = matrix(c(1,2,3), nrow=1))
par(mar = c(2,2.2,1,0.3))
hist(data, freq=F, ylim=my_y_lim, main="SEMI-HMM", xlab=NA,breaks=20, cex.main=1.5, cex=1.3, cex.axis=my_cex_axis)

lines(est$x, tmpMatrix[1,1]*dnorm(est$x, 0, 1), lwd=1.3, lty=4)

lines(est$x, tmpMatrix[1,2]*est$y, lwd=1)

## 2 
hist(data, freq=F, ylim=my_y_lim, main="SUN&CAI", xlab=NA,breaks=20, cex.main=1.5, cex=1.3, cex.axis=my_cex_axis)
lines(est$x, tmpMatrix[1,1]*dnorm(est$x, 0, 1), lwd=1.3, lty=4)
tmpMatrix<-diag(1, nrow=2, ncol=2)
for(i in 1:300){
	tmpMatrix=tmpMatrix%*%em.res1$A
}
lines(est$x, tmpMatrix[1,2]*dnorm(est$x,em.res1$f1[1], em.res1$f1[2]),col="red", lwd=1.3, lty=1)
tmpMatrix<-diag(1, nrow=2, ncol=2)
for(i in 1:300){
	tmpMatrix=tmpMatrix%*%em.res2$A
}
lines(est$x, (tmpMatrix[1,2])*(em.res2$pc[1]*dnorm(est$x,em.res2$f1[1,1], em.res2$f1[1,2])+em.res2$pc[2]*dnorm(est$x,em.res2$f1[2,1], em.res2$f1[2,2])), col="red", lwd=1,lty=2)
legend(my_x_at, my_y_at, c("L=1", "L=2"), lty=c(1,2), col="red", lwd=1.3, cex=my_cex_legend)

## 3
hist(data, freq=F, ylim=my_y_lim, main="LIU", xlab=NA,breaks=20, cex.main=1.5, cex=1.3,, cex.axis=my_cex_axis)
lines(est$x, (mean(c(est_liu_02$p0, est_liu_05$p0, est_liu_08$p0)))*dnorm(est$x, 0, 1), lwd=1, lty=4)

lines(est_liu_02$x, (1-est_liu_02$p0)*est_liu_02$y, lty=1, lwd=1.3, col="blue")
lines(est_liu_05$x, (1-est_liu_05$p0)*est_liu_05$y, lty=2, lwd=1.3, col="blue")
lines(est_liu_08$x, (1-est_liu_08$p0)*est_liu_08$y, lty=3, lwd=1.3, col="blue")

legend(my_x_at, my_y_at, c(expression(paste(lambda,"=0.2")), expression(paste(lambda,"=0.5")), expression(paste(lambda,"=0.8"))), lty=c(1,2,3), col="blue", lwd=1.3, cex=my_cex_legend)

dev.off()

