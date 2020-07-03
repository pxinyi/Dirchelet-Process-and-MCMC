# PACKAGE NEEDED ############################################################################################
library(LaplacesDemon)
library(MASS)
library(mvtnorm)
library(CholWishart)
library(ggplot2)
library(reshape2)



# generating a sample from G0 ########################################################################
generate_location_process=function(){      
  Sigma_sample1=matrix(rInvWishart(n=1, df=nu, Sigma=Chi),ncol = 2,byrow = TRUE)
  mu_sample1=mvrnorm(1,mu0,Sigma_sample1/lambda)
  return(list("mu"=mu_sample1,"Sigma"=Sigma_sample1))
}


# generating weights and locations ###################################################################
stick_breaking_process = function(num_weights, alpha) { 
  betas = rbeta(num_weights, 1, alpha)
  remaining_stick_lengths = c(1, cumprod(1 - betas))[1:num_weights]
  weights = remaining_stick_lengths * betas
  return(weights)
}
max_num_weights=200
G_weights=stick_breaking_process(max_num_weights,alpha)





# FUNCITONS TO SIMPLIFY MCMC PROCESS ########################################################################

dpost_NIM_single=function(xi,nu,lambda,mu0,Chi){
    new_nu=nu+1
    new_lambda=lambda+1
    new_mu0=as.numeric((lambda*mu0+xi)/(lambda+1))
    new_Chi=Chi+(lambda/(lambda+1))*as.numeric((t(t(xi-mu0))%*%t(xi-mu0)))
    
    sample_sigma=matrix(rInvWishart(n=1,df=new_nu,Sigma=new_Chi),ncol=2,byrow = TRUE)
    sample_mu=mvrnorm(1,new_mu0,(sample_sigma/new_lambda))
    return(list("mu"=sample_mu,"Sigma"=sample_sigma))
}


dpost_NIM_multiple=function(datai,nu,lambda,mu0,Chi){
  
  n=dim(datai)[1]
  xbar=colMeans(datai)
  new_nu=nu+n
  new_lambda=lambda+n
  new_mu0=(lambda*mu0+n*xbar)/new_lambda
  
  temp1=sweep(datai, 2, colMeans(datai))
  S=t(temp1)%*%temp1
  temp2=(lambda*n/new_lambda)*matrix(xbar-mu0,ncol = 1)%*%matrix(xbar-mu0,nrow = 1)
  new_Chi=S+Chi+temp2
  sample_sigma=matrix(rInvWishart(n=1, df=new_nu, Sigma=new_Chi),ncol = 2,byrow = TRUE)
  sample_mu=mvrnorm(1,new_mu0,(sample_sigma/new_lambda))
  return(list("mu"=sample_mu,"Sigma"=sample_sigma))
}







# FUNCTIONS FOR MARGINAL DENSITY OF AN OBSERVATION BELONGING TO A NEW CLUSTER ##################################

log_marginal_data=function(data,nu,lambda,mu0,Chi){			
  new_nu=nu+1			
  new_lambda=lambda+1			
  
  temp1=data-mu0
  log_new_Chi_det=rep(0,dim(data)[1])
  for(i in 1:dim(data)[1]){
    temp2=temp1[i,]
    temp3=(lambda/new_lambda)*(matrix(temp2,ncol = 1)%*%matrix(temp2,nrow = 1))
    new_Chi=Chi+temp3
    log_new_Chi_det[i]=as.numeric(determinant(new_Chi)$modulus)
  }
  
  ltemp1=-log(pi)	
  ltemp2=lmvgamma(x=(new_nu/2), p=2)-lmvgamma(x=(nu/2), p=2)
  ltemp3=(nu/2)*as.numeric(determinant(Chi)$modulus)-(new_nu/2)*log_new_Chi_det
  ltemp4=log(lambda/new_lambda)
  
  ltemp=ltemp1+ltemp2+ltemp4
  rlt=ltemp+ltemp3
  
  return(rlt)			
}



log_marginal_x=function(xi,nu,lambda,mu0,Chi){			
  new_nu=nu+1			
  new_lambda=lambda+1			
  new_Chi=Chi+(lambda/(lambda+1))*as.numeric((t(t(xi-mu0))%*%t(xi-mu0)))	
  
  ltemp1=-log(pi)	
  ltemp2=lmvgamma(x=(new_nu/2), p=2)-lmvgamma(x=(nu/2), p=2)
  ltemp3=(nu/2)*as.numeric(determinant(Chi)$modulus)-(new_nu/2)*as.numeric(determinant(new_Chi)$modulus)
  ltemp4=log(lambda/new_lambda)			
  
  rlt=ltemp1+ltemp2+ltemp3+ltemp4
  
  return(rlt)			
}





marginal_x=function(xi,nu,lambda,mu0,Chi){
    new_nu=nu+1
    new_lambda=lambda+1
    new_mu0=as.numeric((lambda*mu0+xi)/(lambda+1))
    new_Chi=Chi+(lambda/(lambda+1))*as.numeric((t(t(xi-mu0))%*%t(xi-mu0)))
    
    temp1=1/pi
    temp2=mvgamma(x=(new_nu/2), p=2)/mvgamma(x=(nu/2), p=2)
    temp3=((det(Chi)**(nu/2))/(det(new_Chi)**(new_nu/2)))
    temp4=lambda/new_lambda
  
    return(temp1*temp2*temp3*temp4)
}




# FUNCTIONS FOR MCMC ITERATION ###################################################################################
mciter=function(iC,iNC,inum,imu,isigma,iM){
  C=iC
  NC=iNC
  Num_clusters=inum
  theta_mu=imu
  theta_sig=isigma
  logxi_M=iM
  
  
  for(indexi in 1:dim(X)[1]){
    #step2&3
    labeli=C[indexi]
    if(NC[labeli]==1){
      
      theta_mu[[labeli]]=theta_mu[[Num_clusters]]
      theta_mu[[Num_clusters]]=NA
      theta_sig[[labeli]]=theta_sig[[Num_clusters]]
      theta_sig[[Num_clusters]]=NA
      C[C==Num_clusters]=labeli
      logxi_M[,labeli]=logxi_M[,Num_clusters]
      logxi_M[,Num_clusters]=NA
      NC[labeli]=NC[Num_clusters]
      NC[Num_clusters]=NA
      Num_clusters=Num_clusters-1
    }else{
      NC[labeli]=NC[labeli]-1
    }
    C[indexi]=NA
    
    #step4
    xi=X[indexi,]
    logp_list=log(NC[1:Num_clusters])+logxi_M[indexi,1:Num_clusters]
    log_new=log(alpha)+log_marginal_dataM[indexi]
    mm=max(logp_list,log_new)
    tr_logp_list=exp(logp_list-mm)
    tr_log_new=exp(log_new-mm)
    prob_list=c(tr_log_new,tr_logp_list)/sum(tr_log_new,tr_logp_list)
    ci=sample(c(-1,1:Num_clusters),size = 1,prob = prob_list)
    
    #step5
    if(ci==-1){
      #print(paste("current index is",as.character(indexi),"and log_new is"))
      #print(log_new)
      #print("##############")
      
      new_label=Num_clusters+1
      C[indexi]=new_label
      new_theta=dpost_NIM_single(xi,nu,lambda,mu0,Chi)
      theta_mu[[new_label]]=new_theta$mu
      theta_sig[[new_label]]=new_theta$Sigma
      NC[new_label]=1
      t=dmvnorm(x=X, mean =new_theta$mu, sigma = new_theta$Sigma, log = TRUE)
      if(new_label > dim(logxi_M)[2]){
        logxi_M=cbind(logxi_M,t)
      }else{
        logxi_M[,new_label]=t
      }
      
      Num_clusters=Num_clusters+1
    } else {
      C[indexi]=ci
      NC[ci]=NC[ci]+1
    }
  }
  
  logxi_M=matrix(NA,nrow = dim(X)[1],ncol = Num_clusters)
  for(label in 1:Num_clusters){
    datai=X[(C==label),,drop=FALSE]
    new_theta=dpost_NIM_multiple(datai,nu,lambda,mu0,Chi)
    theta_mu[[label]]=new_theta$mu
    theta_sig[[label]]=new_theta$Sigma
    logxi_M[,label]=dmvnorm(x=X, mean =new_theta$mu, sigma = new_theta$Sigma, log = TRUE)
  }
  
  return(list("C"=C,"NC"=NC,"Num_clusters"=Num_clusters,"theta_mu"=theta_mu,"theta_sig"=theta_sig,"logxi_M"=logxi_M))
}





mciter_8=function(iC,iNC,inum,imu,isigma,iM,km=3){
  C=iC
  NC=iNC
  Num_clusters=inum
  theta_mu=imu
  theta_sig=isigma
  logxi_M=iM
  
  
  log_marginal_possi_theta=matrix(nrow = dim(X)[1],ncol = km*10)
  possi_theta=list()
  for(i in 1:(km*10)){
    temp1=generate_location_process()
    possi_theta[[i]]=temp1
    log_marginal_possi_theta[,i]=dmvnorm(X, mean =temp1$mu, sigma = temp1$Sigma, log = TRUE)
  }
  
  
  for(indexi in 1:dim(X)[1]){
    #step2&3
    labeli=C[indexi]
    if(NC[labeli]==1){
      theta_mu[[labeli]]=theta_mu[[Num_clusters]]
      theta_mu[[Num_clusters]]=NA
      theta_sig[[labeli]]=theta_sig[[Num_clusters]]
      theta_sig[[Num_clusters]]=NA
      C[C==Num_clusters]=labeli
      logxi_M[,labeli]=logxi_M[,Num_clusters]
      logxi_M[,Num_clusters]=NA
      NC[labeli]=NC[Num_clusters]
      NC[Num_clusters]=NA
      Num_clusters=Num_clusters-1
    }else{
      NC[labeli]=NC[labeli]-1
    }
    C[indexi]=NA
    
    
    #step4
    xi=X[indexi,]
    logp_list=log(NC[1:Num_clusters])+logxi_M[indexi,1:Num_clusters]
    indexci=sample(1:(10*km),km,replace = FALSE)
    new_logp=log(alpha/km)+log_marginal_possi_theta[indexi,indexci]
    mm=max(logp_list,new_logp)
    tr_logp_list=exp(logp_list-mm)
    tr_log_new=exp(new_logp-mm)
    prob_list=c(tr_log_new,tr_logp_list)/sum(tr_log_new,tr_logp_list)
    ci=sample(c(-1*indexci,1:Num_clusters),size = 1,prob = prob_list)
    
    #step5
    if(ci< 0){
      
      new_label=Num_clusters+1
      C[indexi]=new_label
      new_theta=possi_theta[[(-ci)]]
      theta_mu[[new_label]]=new_theta$mu
      theta_sig[[new_label]]=new_theta$Sigma
      NC[new_label]=1
      t=dmvnorm(x=X, mean =new_theta$mu, sigma = new_theta$Sigma, log = TRUE)
      if(new_label > dim(logxi_M)[2]){
        logxi_M=cbind(logxi_M,t)
      }else{
        logxi_M[,new_label]=t
      }
      Num_clusters=Num_clusters+1
      
      new_para=generate_location_process()
      possi_theta[[(-ci)]]=new_para
      log_marginal_possi_theta[,(-ci)]=dmvnorm(X, mean =new_para$mu, sigma = new_para$Sigma, log = TRUE)
      
      
      print(paste("the current indexi is",as.character(indexi),"with parameter as"))
      print(possi_theta[[(-ci)]])
      print(paste("log_density",as.character(log_marginal_possi_theta[indexi,(-ci)])))
      print("###########")
      
    } else {
      C[indexi]=ci
      NC[ci]=NC[ci]+1
    }
  }
  
  
  logxi_M=matrix(NA,nrow = dim(X)[1],ncol = Num_clusters)
  for(label in 1:Num_clusters){
    datai=X[(C==label),,drop=FALSE]
    new_theta=dpost_NIM_multiple(datai,nu,lambda,mu0,Chi)
    theta_mu[[label]]=new_theta$mu
    theta_sig[[label]]=new_theta$Sigma
    logxi_M[,label]=dmvnorm(x=X, mean =new_theta$mu, sigma = new_theta$Sigma, log = TRUE)
  }
  
  return(list("C"=C,"NC"=NC,"Num_clusters"=Num_clusters,"theta_mu"=theta_mu,"theta_sig"=theta_sig,"logxi_M"=logxi_M))
}


    
   







# FUNCTION TO COMPUTE POSTERIOR PREDICTED DENSITY FOR ANY PREDICTED DATA #####################################

marginal_density=function(idata,iC,inum,imu,isigma){
    
    freq=as.matrix(table(iC))/(length(iC)+alpha)
    
    md1=rep(0,dim(idata)[1])
    for(i in 1:inum){
        mu=imu[(2*i-1):(2*i)]
        sigma=matrix(isigma[(4*i-3):(4*i)],ncol = 2,byrow = TRUE)
        pre_density=dmvnorm(x=idata, mean =mu, sigma =sigma, log = FALSE)*freq[i]
        md1=md1+pre_density
    }
    md2=apply(idata,1,FUN = function(i){marginal_x(i,nu,lambda,mu0,Chi)})
    md2=(alpha/(length(iC)+alpha))*md2
    
    return(md1+md2)
}








# FUNCTION TO SIMULTE DATASET BASED ON CLUSTER CENTER AND COVARIANCE MATRIX ###################################

simulation=function(mudf,sigmadf,Num_simu){
    datalist=list()
    N=dim(mudf)[1]
    
    for (i in 1:N){
        mui=as.numeric(mudf[i,])
        sigmai=matrix(as.numeric(sigmadf[i,]/9),ncol = 2,byrow = TRUE)
        samplei=mvrnorm(n=Num_simu,mu=mui,Sigma=sigmai)
        dfi=data.frame(samplei)
        colnames(dfi)=c("x1","x2")
        dfi$C=i
        datalist[[i]]=dfi
    }
    
    output=do.call(rbind, datalist)
    return(output)
}









