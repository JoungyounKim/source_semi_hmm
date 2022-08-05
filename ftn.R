
#####################################################
## Function to return sum  of elements in vector a ##
#####################################################
carefulSum<-function(a){

	a<-a[a>0]
	minA<-min(a)
	tmp<-minA*sum(exp(log(a/minA)))
	return(tmp)
}

carefulSum_in_log<-function(log_a){

	#a<-a[a>0]
	log_minA<-min(log_a)

	tmp<-log(sum(exp(log_a-log_minA)))+log_minA
	
	#tmp<-minA*sum(exp(log(a/minA)))
	return(tmp)
}




carefulSum<-function(a){

	a<-a[a>0]
	minA<-min(a)
	tmp<-minA*sum(exp(log(a/minA)))
	return(tmp)
}


##########################################################
## Function to return transition probability Pr(z1->z2) ##
##########################################################
transit.prob<-function(z1, z2, a0, a1){


	
	if(z1==0 & z2==0){
		if(a0==0)
			return (1*10^{-20})
		return (a0)
	}
	if(z1==0 & z2==1){
		if(a0==1)
			return (1*10^{-20})
		return (1-a0)
	}
	if(z1==1 & z2==0){
		if(a1==1)
			return (1*10^{-20})
		return (1-a1)
	}
	if(z1==1 & z2==1){
		if(a1==0)
			return (1*10^{-20})
		return (a1)
	}
}

####################################################################
# Functio to return pdf from Alternative model (Assumed under H1)  #
####################################################################
#
#  x0: point to evaluate pdf
#  pdf: Kernel estimates for PDF
#
#####################################################################
pdf_H1_new<-function(x0,x,pdf){

	
	if(x0<=pdf$x[1]){
		if(pdf$y[1]<10^{-20})
			return(log(10^{-20}))
		return(log(pdf$y[1]))
		
	}
	if(x0>=pdf$x[length(pdf$x)]){
		if(pdf$y[length(pdf$x)]<10^{-20})
			return(log(10^{-20}))
		return(pdf$y[length(pdf$x)])
	}


	tmpWhich<-which(pdf$x==x0)
	if(length(tmpWhich)==1)
		return(pdf$y[tmpWhich])

	#print("1. In pdf_H1_new")
      at0<-which.min(abs(pdf$x-x0))

	if(pdf$x[at0]-x0<0){
		at1=at0
		at2=at0+1		
	}
	if(pdf$x[at0]-x0>=0){	
		at1=at0-1
		at2=at0
	}
	
	a=x0-pdf$x[at1]
	b=pdf$x[at2]-x0
	a<-a/(a+b)
	tmp<-(1-a)*pdf$y[at1]+a*pdf$y[at2]

	
	if(tmp==0)
		return (log(1*10^{-20}))

	
	return(log(tmp))
}



###############################################
##    Function to compute Prob(Z_k=1 | Data) ##      
###############################################


peeling_for_one<-function(index, data, a0, a1, cdf_tmp,p0,p1){

       
	z_i<-c(0,1)
	z_i_r<-c(0,1)
	z_i_l<-c(0,1)

	f_i<-c(0,0)
	f_i_r<-c(0,0)
	f_i_l<-c(0,0)
	### Forward

	if(index==1){
		f_i_l[1]<-log(p0)
		f_i_l[2]<-log(1-p0)
		}
	if(index>=2){
		f_i_l[1]<-dnorm(data[1], mean=mu0, sd=1, log=T)+log(p0)
		f_i_l[2]<-pdf_H1_new(data[1],data,cdf_tmp)+log(1-p0)

		}

	if(index>2){
		for(i in 2:(index-1)){
			f_tmp<-rep(0, 4)
			f_tmp[1]<-f_i_l[1]+log(transit.prob(0,0,a0, a1))+dnorm(data[i], mean=mu0, sd=1, log=T)
			f_tmp[2]<-f_i_l[1]+log(transit.prob(0,1,a0, a1))+pdf_H1_new(data[i],data,cdf_tmp)


			f_tmp[3]<-f_i_l[2]+log(transit.prob(1,0,a0, a1))+dnorm(data[i], mean=mu0, sd=1, log=T)
			
			f_tmp[4]<-f_i_l[2]+log(transit.prob(1,1,a0, a1))+pdf_H1_new(data[i],data,cdf_tmp)


			f_i_l[1]<-carefulSum_in_log(c(f_tmp[1],f_tmp[3]))
			f_i_l[2]<-carefulSum_in_log(c(f_tmp[2], f_tmp[4]))
			
		}#end_for	
	}#end_if

	
### At index

	f_i[1]<-dnorm(data[index], mean=mu0, sd=1, log=T)
	f_i[2]<-pdf_H1_new(data[index],data,cdf_tmp)


	

###  Backward

	if(index==length(data)){
		f_i_r[1]<-log(1/2)
		f_i_r[2]<-log(1/2)		
	}
	if(index<length(data)){
		f_i_r[1]<-dnorm(data[length(data)], mean=mu0, sd=1, log=T)
		f_i_r[2]<-pdf_H1_new(data[length(data)],data,cdf_tmp)


	}
	if(index<(length(data)-1)){

		for(i in (length(data)-1):(index+1)){
			f_tmp<-rep(0, 4)
			f_tmp[1]<-f_i_r[1]+log(transit.prob(0,0,a0, a1))+dnorm(data[i], mean=mu0, sd=1, log=T)		
			f_tmp[2]<-f_i_r[1]+log(transit.prob(1,0,a0, a1))+pdf_H1_new(data[i],data,cdf_tmp)
			f_tmp[3]<-f_i_r[2]+log(transit.prob(0,1,a0, a1))+dnorm(data[i], mean=mu0, sd=1, log=T)
			f_tmp[4]<-f_i_r[2]+log(transit.prob(1,1,a0, a1))+pdf_H1_new(data[i],data,cdf_tmp)
			f_i_r[1]<-carefulSum_in_log(c(f_tmp[1], f_tmp[3]))
			f_i_r[2]<-carefulSum_in_log(c(f_tmp[2], f_tmp[4]))	
		}#end_for
	}#end_if

	

	tmp_matrix<-matrix(0, nrow=8, ncol=4)
	count=0
	for(i in 1:2){
		for(j in 1:2){
			for(k in 1:2){
				count=count+1		
				tmp=f_i_l[i]+f_i[j]+f_i_r[k]
				if(index>1)
					tmp=tmp+log(transit.prob(i-1, j-1, a0, a1))
				if(index<length(data))
					tmp=tmp+log(transit.prob(j-1, k-1, a0, a1))
				tmp_matrix[count, ]=c(i-1, j-1, k-1, tmp)
			}#endFOR
		}#endFOR
	}#endFOR
	
	den<-carefulSum_in_log(tmp_matrix[,4])
	num<-carefulSum_in_log(tmp_matrix[c(3,4,7,8),4])

	p00<-carefulSum_in_log(tmp_matrix[1:2,4])
	p01<-carefulSum_in_log(tmp_matrix[3:4,4])
	p10<-carefulSum_in_log(tmp_matrix[5:6,4])
	p11<-carefulSum_in_log(tmp_matrix[7:8,4])

       
	if(is.na(den) ||is.nan(den)){
		print(c("index = ", index))
		print(tmp_matrix)
		print(c(num, den))
		stop("ERROR")
	}
	p1<-exp(num-den)
	p00<-exp(p00-den)
	p01<-exp(p01-den)
	p10<-exp(p10-den)
	p11<-exp(p11-den)

	if(index==1){
		p00=p01=p10=p11=0
	}
	if(is.na(p1) || is.nan(p1)){
        
		stop("ERROR")
        }
	
	return(list(p1=p1, p00=p00, p01=p01, p10=p10, p11=p11))
}#end_of_peeling_for_one



lis<-function(my_prob, alpha){
	print("in LIS")	
	my_prob_sort<-sort(my_prob)
	my_cum<-cumsum(my_prob_sort)
	my_cum<-my_cum/(1:length(my_prob))
	print(my_prob)
	print(my_cum)
	

	if(min(my_cum)>alpha){
		return(rep(0,length(my_prob)))
	}
	at.which<-max(which(my_cum<=alpha))
	lis_at_which<-my_prob_sort[at.which]
	tmp_pred<-ifelse(my_prob<lis_at_which,1,0)
	
	return(tmp_pred)
}


#######################################################
##  function to estimate p(z=1),a0, and a1
#######################################################
estimate<-function(data=data, case, mu1, a0_true, a1_true){

	Max_count<-20
	
	a0<-a0_true
	a1<-a1_true
	
	count_conv<-0
	tol<-0.0001
	x_final<-seq(-15, 10, by=0.01)

			
	w<-rep(1, length(data))
	w<-w/sum(w)
	cdf_tmp<-NULL
	pdf_tmp<-NULL
	pdf_kernel<-NULL
	cdf_alpha<-NULL
		
	
	cut_off<-1.0
	

	###################################
	##  Add random noise (to Data) to avoid ties
	###################################
	tmpCount<-0
	while(TRUE){
		data_orig<-sort(unique(data))
		tmpCount<-tmpCount+1
		print("$$$$$")
		print(tmpCount)
		print("$$$$$")
		tmpTable<-table(data)
		multiple_value<-as.numeric(names(tmpTable)[which(tmpTable[]>1)])
		if(length(multiple_value)==0)
			break
		for(j in 1:length(multiple_value)){
			#print(j)
			tmpWhich<-which(abs(data-multiple_value[j])<10^(-10))
			data[tmpWhich]<-duplicateData(data_orig, multiple_value[j], length(tmpWhich))
		}
		
	}


	##########################
	###  Added  on 20170413
	##########################
	pdf_tmp<-rep(0, length(data))
	pdf_tmp[which(abs(data)>cut_off)]<-1
	pdf_kernel<-density(data, weights=pdf_tmp/sum(pdf_tmp))
	
	
	w<-rep(1/n,n)

	
	pdf_true<-NULL
	cdf_true<-NULL
	cdf_tmp_old<-NULL
	prob1_old<-rep(0, length(data))

	a0_trace<-NULL
	a1_trace<-NULL
	p0<-NULL
	p1<-NULL
	cdf_alpha_old<-0


 



	while(TRUE){				
		print("@@@@@@@@@@@@@@@@@@@@@@@@@@")
		print(count_conv)
		print("@@@@@@@@@@@@@@@@@@@@@@@@@@")

		cdf_tmp_old<-cdf_tmp
		#----------------------------#
		## Update density estimation##
		#----------------------------#
	
		
		w[which(w>0)]<-exp(log(w[which(w>0)])-log(carefulSum(w)))
		w<-ifelse(w<1/n, 1/n, w)		
		w[which(w>0)]<-exp(log(w[which(w>0)])-log(carefulSum(w)))


		##########################
		###  Blocked on 20170411
		##########################

		tmp_cdf_finall<-cdf_finall(sort(data), w[order(data)], length(data)/10)
		cdf_tmp<-tmp_cdf_finall$tmp
		tmp_diff<-diff(cdf_tmp)

		
		print(paste("alpha0 = ", tmp_cdf_finall$alpha))

		
		cdf_alpha<-max(min(tmp_cdf_finall$alpha, 0.95), 0.05)## proportion of Non-Null

		
		pdf_tmp<-diff(c(0, cdf_tmp))
		pdf_tmp[pdf_tmp<0]<-0

		

		pdf_tmp_weight<-pdf_tmp
		pdf_tmp_weight[pdf_tmp_weight>0]<-exp(log(pdf_tmp_weight[pdf_tmp_weight>0])-log(carefulSum(pdf_tmp_weight[pdf_tmp_weight>0])))
		pdf_kernel<-density(sort(data), weights=pdf_tmp_weight)####
		
		
		tmp_x<-unique(sort(data))

		

		pdf_true<-rep(0, length(pdf_kernel$x))
		cdf_true<-rep(0, length(pdf_kernel$x))
		if(case==1){
			pdf_true<-dnorm(pdf_kernel$x, mean=mu1, sd=1)
			cdf_true<-pnorm(pdf_kernel$x, mean=mu1, sd=1)
		}
		if(case==2){
			pdf_true<-0.5*dnorm(pdf_kernel$x, mean=mu1, sd=1)+0.5*dnorm(pdf_kernel$x, mean=2, sd=1)
			cdf_true<-0.5*pnorm(pdf_kernel$x, mean=mu1, sd=1)+0.5*pnorm(pdf_kernel$x, mean=2, sd=1)	
		}

		if(case==3){
			pdf_true<-0.4*dnorm(pdf_kernel$x, mean=mu1, sd=1)+0.3*dnorm(pdf_kernel$x, mean=1, sd=1)+0.3*dnorm(pdf_kernel$x, mean=3, sd=1)
			cdf_true<-0.4*pnorm(pdf_kernel$x, mean=mu1, sd=1)+0.3*pnorm(pdf_kernel$x, mean=1, sd=1)+0.3*pnorm(pdf_kernel$x, mean=3, sd=1)
	
		}
		if(case==4){
			
			mixture_ratio<-c(0.1, 0.4, 0.5)
			pdf_true<-mixture_ratio[1]*dnorm(pdf_kernel$x, mean=2, sd=1)+mixture_ratio[2]*dnorm(pdf_kernel$x, mean=(2+mu1), sd=1)+mixture_ratio[3]*dnorm(pdf_kernel$x, mean=(2+4*mu1), sd=1)
			cdf_true<-mixture_ratio[1]*pnorm(pdf_kernel$x, mean=2, sd=1)+mixture_ratio[2]*pnorm(pdf_kernel$x, mean=(2+mu1), sd=1)+mixture_ratio[3]*pnorm(pdf_kernel$x, mean=(2+4*mu1), sd=1)
		}
		if(case==5){
			pdf_true<-0.5*dnorm(pdf_kernel$x, mean=mu1, sd=1)+0.5*dgamma(pdf_kernel$x, shape=-mu1, scale=2)
			cdf_true<-0.5*pnorm(pdf_kernel$x, mean=mu1, sd=1)+0.5*pgamma(pdf_kernel$x, shape=-mu1, scale=2)	

		}

		if(case==6){
			pdf_true<-dunif(pdf_kernel$x, min=-mu1, max=mu1)
			cdf_true<-punif(pdf_kernel$x, min=-mu1, max=mu1)

		}

		count_conv=count_conv+1								
		a0_old<-a0
		a1_old<-a1

		
		if(count_conv==1){
			p1=0.5	
			p0=1-p1
		}

		#----------------------------------------------------#
		### Estimate P(Z_k=1 |Data), Pr(Z_{k-1}, Z_k |Data)###
		#----------------------------------------------------#
		prob1<-rep(0, n)
		
		prob2<-matrix(0, nrow=n, ncol=4)
		prob2_old<-matrix(0, nrow=n, ncol=4)

		

		tmp_check<-NULL
		time_old<-NULL
		time_new<-NULL
		time_old_total<-0
		time_new_total<-0
		count=0
               
		for(index in 1:length(data)){
			
			windowSize<-5
			lowerIndex<-max(index-windowSize, 1)
			upperIndex<-min(index+windowSize, length(data))
			tmpData<-data[lowerIndex:upperIndex]
			tmpIndex<-windowSize+1
			if(index<=windowSize+1)
				tmpIndex=index
	
			aa<-data.frame(cbind(sort(data),pdf_tmp))
			names(aa)<-c("x", "y")

			tmp<-peeling_for_one(tmpIndex, tmpData, a0, a1, pdf_kernel,p0, p1) # modified on 2019/03/14
			

			
			prob1[index]<-tmp$p1 ## Prob to be Z=1
			
			prob2[index,]<-c(tmp$p00, tmp$p01, tmp$p10, tmp$p11)

			

		}#end_for_index

             

		if(length(which(prob1<0)>0))
			stop("Error: prob cannot be less than 0")	

		prob2<-prob2[-1,]

	
		
		a0<-sum(prob2[,1])/sum(1-prob1[-n])
		a1<-sum(prob2[,4])/sum(prob1[-n])


		tmpDiag<-diag(2)
		for(i in 1:100){

			tmpDiag<-tmpDiag%*%matrix(c(a0, 1-a0, 1-a1, a1), byrow=T, nrow=2)
		}
		p0<-tmpDiag[1,1]
		p1<-tmpDiag[1,2]		

		
		p0<-ifelse(p0>=1,1-1e-10,p0)
		p0<-ifelse(p0<0,1e-10,p0)
		p1<-ifelse(p1>=1, 1-1e-10, p1)
		p1<-ifelse(p1<0, 1e-10, p1)
		print(c(p0, p1))
		a0_trace<-c(a0_trace, a0)
		a1_trace<-c(a1_trace, a1)

		a0<-ifelse(a0>=1, 1-1e-10, a0)
		a0<-ifelse(a0<0, 1e-10, a0)
		a1<-ifelse(a1>=1, 1-1e-10, a1)
		a1<-ifelse(a1<0, 1e-10, a1)
		

		w=prob1
		
                if(is.nan(mean(cdf_tmp))){
                        print(cdf_tmp)
                        break
                }
		
		
		if(count_conv>1 && mean(abs(cdf_tmp_old-cdf_tmp), na.rm=T)<tol || abs(cdf_alpha_old-cdf_alpha)<tol)
			break
	
		
		if(count_conv==Max_count)
			break

		
		cdf_tmp_old<-cdf_tmp
		cdf_alpha_old<-cdf_alpha
		
	}#end_while
	print(count_conv)


	return(list(p1=prob1, a0=a0, a1=a1, pdf_kernel=pdf_kernel, x=pdf_kernel$x, y=pdf_kernel$y,alpha=cdf_alpha))
}## end of function::estimate


liu<-function(data=data, case, mu1, a0_true, a1_true, lambda){



	### compute p-val
	p_val<-(1-pnorm(abs(data), mean=0, sd=1))*2
	w_lambda<-sum(p_val>lambda)


	###  Determine p0
	p0<-max(0.01, min(w_lambda/(1-lambda)/length(data),0.99))




	###  get Kernel CDF of data  #####
	pdf_kernel<-density(sort(data),kernel="epanechnikov" )


	### get F1
	   f1<-(pdf_kernel$y-p0*dnorm(pdf_kernel$x))/(1-p0)
	   f1<-ifelse(f1<0, 10^{-20}, f1)
         #f1<-ifelse(f1>1, 1, f1)
		
         f1<-list(x=pdf_kernel$x, y=f1)


	### compute LIS
				
	
		a0<-0.5
		a1<-0.5

		Max_count<-50
		tol<-0.0001

		count_conv<-0
		while(TRUE){	

			count_conv=count_conv+1

			a0_old<-a0
			a1_old<-a1

			p1=1-p0

			#----------------------------------------------------#
			### Estimate P(Z_k=1 |Data), Pr(Z_{k-1}, Z_k |Data)###
			#----------------------------------------------------#
			prob1<-rep(0, n)
			prob2<-matrix(0, nrow=n, ncol=4)
			prob2_old<-matrix(0, nrow=n, ncol=4)


			tmp_check<-NULL
			time_old<-NULL
			time_new<-NULL
			time_old_total<-0
			time_new_total<-0
			count=0
			for(index in 1:length(data)){
			
				windowSize<-5
				lowerIndex<-max(index-windowSize, 1)
				upperIndex<-min(index+windowSize, length(data))
				tmpData<-data[lowerIndex:upperIndex]
				tmpIndex<-windowSize+1
				if(index<=windowSize+1)
					tmpIndex=index
	
				tmp<-peeling_for_one(tmpIndex, tmpData, a0, a1, f1,p0, p1) # modified on 2019/03/14

				
				prob1[index]<-tmp$p1
			
				prob2[index,]<-c(tmp$p00, tmp$p01, tmp$p10, tmp$p11)


			}#end_for_index
		

			if(length(which(prob1<0)>0))
				stop("Error: prob cannot be less than 0")	

			prob2<-prob2[-1,]


			a0<-max(0,min(1,sum(prob2[,1])/sum(1-prob1[-n])))
			a1<-max(0,min(1,sum(prob2[,4])/sum(prob1[-n])))
				
			

			if(a0<1e-20)
				a0=1e-20
			if(a1<1e-20)
				a1=1e-20
		
										
			if(count_conv>1 && max(abs(c(a0-a0_old, a1-a1_old)))<tol)
				break

			if(count_conv==Max_count)
				break
		}# end_while

		
	
		print(count_conv)


		return(list(p1=prob1, p0=p0, x=f1$x, y=f1$y))
	

}## end of function::semi

duplicateData<-function(data, x, n){
	tmpData<-sort(unique(data))
	lower<-NULL
	upper<-NULL
	p<-1/4
	diff<-NULL
	crit<-10^(-10)
	if(abs(x-tmpData[1])>crit && abs(x-tmpData[length(tmpData)])>crit){
		tmpWhich<-which(abs(x-tmpData)<crit)
		lower<-p*tmpData[tmpWhich-1]+(1-p)*x
		upper<-p*tmpData[tmpWhich+1]+(1-p)*x
	}
	if(abs(x-tmpData[1])<crit){
		lower<-x
		upper<-p*tmpData[2]+(1-p)*x
	}
	if(abs(x-tmpData[length(tmpData)])<crit){
		lower<-p*tmpData[length(tmpData)-1]+(1-p)*x
		upper<-x
	}
	
	if(is.na(lower) || is.na(upper)){
		stop("Error")	
	}
	tmpX<-runif(n, min=lower, max=upper)
	return(tmpX)
}


#############################################
#### Original code is in "code_from_lee"  ###
#############################################
# For Weighted ECDF, use ewcdf() under library(spatstat)

library(spatstat)
library(Iso)
library(fdrtool)

# This function calculates the distance $\gamma  \
#d_n(\hat{F}_{s,n}^{\gamma},\check{F}_{s,n}^\gamma)$
# for grid of gamma values in [0,1].
# Input data is a numeric vector containing observations from the mixture model.
# Input w is the vector of weights.
# Bigger gridsize  gives more  accurate estimates of alpha.

EstMixMdl <- function(data, w, gridsize)
{
n <- length(data)
w<-w[order(data)]
data <- sort(data)
#plot(data, w)
data.1 <- unique(data)
Fn <- ewcdf(data, w)
Fn.1 <- Fn(data.1)
## Calculate the known F_b at the data points
## Use Standard Normal CDF
Fb <- pnorm(data.1)
## Compute the weights (= frequency/n) of the unique data values, i.e., dF_n
# Note that Fn.1 is already weighted
Freq <- diff(c(0,Fn.1))
distance <- rep(0,gridsize)
##distance[0]<- sqrt(t((Fn.1-Fb)^2)%*%Freq)
distance[0]<- sqrt(sum((Fn.1-Fb)^2*Freq, na.rm=T))
for(i in 1:gridsize)
{
  a <- i/gridsize               ## Assumes a value of the mixing proportion
  F.hat <- (Fn.1-(1-a)*Fb)/a     ## Computes the naive estimator of F_s
  F.is=pava(F.hat,Freq,decreasing=FALSE) ## Computes the Isotonic Estimator of F_s
  F.is[which(F.is<=0)]=0
  F.is[which(F.is>=1)]=1
  F.is[is.na(F.is)] <- 0
  ##distance[i] <- a*sqrt(t((F.hat-F.is)^2)%*%Freq);
  distance[i] <- a*sqrt(sum((F.hat-F.is)^2*Freq, na.rm=T))
}
  return(list("dist"=distance, "Fn.1"=Fn.1, "Fb"=Fb))
}

# The following function evaluates the numerical second derivative of any function

Comp_2ndDer <- function(dist.alpha, gridsize)
  {
  dder <- diff(dist.alpha)    ## Computes the 1st order differences
  dder <- diff(dder)      ## Computes the 2nd order differences
  dder <- c(0,0,dder)       ## The numerical double derivative vector

  return(dder)
}

# Compute CDF

CDFEst <- function(Fn.1, Fb, Est)
{
## Compute the weights (= frequency/n) of the unique data values, i.e., dF_n
Freq <- diff(c(0,Fn.1))
## Computes the naive estimator of F_s
Est.CDF.naive <- (Fn.1-(1-Est)*Fb)/Est
## Computes the Isotonic Estimator of F_s
Est.CDF=pava(Est.CDF.naive,Freq,decreasing=FALSE)
Est.CDF[which(Est.CDF<=0)]=0
Est.CDF[which(Est.CDF>=1)]=1
Est.CDF[is.na(Est.CDF)]=0
#plot(Est.CDF.naive, main="naive")
#plot(Est.CDF, main="final")
return(cbind(Est.CDF.naive,Est.CDF))
}

# Compute Density
# This was also written by Patra and Sen
# Not sure if we need this.

DensEst <- function(data, Fn.1, Fb, Est)
{
F.hat <- (Fn.1-(1-Est)*Fb)/Est
Freq <- diff(c(0,Fn.1))
F.is <- pava(F.hat,Freq,decreasing=FALSE)
F.is[which(F.is<=0)] <- 0
F.is[which(F.is>=1)] <- 1
F.check <- F.is
data.1 <- unique(sort(data))
x <- data.1
y <- F.check
ll <- gcmlcm(x,y, type="lcm")
xtemp=rep(ll$x.knots,each=2)
ytemp=c(0,rep(ll$slope.knots,each=2),0)
ans<-rbind(t(xtemp),t(ytemp))
return(ans)
}

cdf_finall<-function(data, w, gridsize){

	estMix<-EstMixMdl(data, w, gridsize)
	
	dist.alpha <- estMix$dist
	Fn.1 <- estMix$Fn.1
	Fb <- estMix$Fb

	dder <- Comp_2ndDer(dist.alpha, gridsize)
	Est <- which.max(dder)/gridsize

	tmp<-CDFEst (Fn.1, Fb, Est)

	##############
	### Added on 20170407
	############

	my_unique_data<-unique(sort(data))
	
	tmp<-tmp[,2]
	#plot(tmp)
	return(list("alpha"=Est, "tmp"=tmp))
}


pdf_finall<-function(data, w, gridsize){

	estMix<-EstMixMdl(data, w, gridsize)
	
	dist.alpha <- estMix$dist
	Fn.1 <- estMix$Fn.1
	Fb <- estMix$Fb

	dder <- Comp_2ndDer(dist.alpha, gridsize)
	Est <- which.max(dder)/gridsize

	#tmp<-CDFEst (Fn.1, Fb, Est)
	tmp<-DensEst(data, Fn.1, Fb, Est)

	##############
	### Added on 20170407
	############

	#my_unique_data<-unique(sort(data))
	

	return(tmp)
}








