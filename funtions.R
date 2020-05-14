#####useful package
#library(LaplacesDemon)
library(MASS)
library(mvtnorm)
library(CholWishart)
library(ggplot2)
library(reshape2)


####useful function


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
    datai=as.matrix(datai)
    n=dim(datai)[1]
    xbar=as.numeric(apply(datai,2,mean))
    new_nu=nu+n
    new_lambda=lambda+n
    new_mu0=(lambda*mu0+n*xbar)/(lambda+n)
    temp1=apply(datai,1,FUN = function(xi){xi-xbar})
    S=temp1%*%t(temp1)
    new_Chi=Chi+(lambda*n/(lambda+n))*t(t(xbar-mu0))%*%t(xbar-mu0)+S
    
    sample_sigma=matrix(rInvWishart(n=1, df=new_nu, Sigma=new_Chi),ncol = 2,byrow = TRUE)
    sample_mu=mvrnorm(1,new_mu0,(sample_sigma/new_lambda))
    return(list("mu"=sample_mu,"Sigma"=sample_sigma))
}


log_marginal_x=function(xi,nu,lambda,mu0,Chi){
    new_nu=nu+1
    new_lambda=lambda+1
    new_mu0=as.numeric((lambda*mu0+xi)/(lambda+1))
    new_Chi=Chi+(lambda/(lambda+1))*as.numeric((t(t(xi-mu0))%*%t(xi-mu0)))
    
    temp1=1/pi
    temp2=mvgamma(x=(new_nu/2), p=2)/mvgamma(x=(nu/2), p=2)
    temp3=((det(Chi)**(nu/2))/(det(new_Chi)**(new_nu/2)))
    temp4=lambda/new_lambda
    
    return(log(temp1*temp2*temp3*temp4))
}




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
        xi=as.numeric(X[indexi,])
        logp_list=log(NC[1:Num_clusters])+as.numeric(logxi_M[indexi,1:Num_clusters])
        log_new=log(alpha)+log_marginal_x(xi,nu,lambda,mu0,Chi)
        mm=max(logp_list,log_new)
        tr_logp_list=exp(logp_list-rep(mm,Num_clusters))
        tr_log_new=exp(log_new-mm)
        prob_list=c(tr_log_new,tr_logp_list)/sum(tr_log_new,tr_logp_list)
        ci=sample(c(-1,1:Num_clusters),size = 1,prob = prob_list)
        
        #step5
        if(ci==-1){
            new_label=Num_clusters+1
            C[indexi]=new_label
            new_theta=dpost_NIM_single(xi,nu,lambda,mu0,Chi)
            theta_mu[[new_label]]=new_theta$mu
            theta_sig[[new_label]]=new_theta$Sigma
            NC[new_label]=1
            logxi_M[,new_label]=dmvnorm(x=X, mean =new_theta$mu, sigma = new_theta$Sigma, log = TRUE)
            Num_clusters=Num_clusters+1
        } else {
            C[indexi]=ci
            NC[ci]=NC[ci]+1
        }
    }
    
    logxi_M=data.frame(matrix(NA,nrow = dim(X)[1],ncol = Num_clusters))
    for(label in 1:Num_clusters){
        datai=X[(C==label),]
        new_theta=dpost_NIM_multiple(datai,nu,lambda,mu0,Chi)
        theta_mu[[label]]=new_theta$mu
        theta_sig[[label]]=new_theta$Sigma
        logxi_M[,label]=dmvnorm(x=X, mean =new_theta$mu, sigma = new_theta$Sigma, log = TRUE)
    }
    
    return(list("C"=C,"NC"=NC,"Num_clusters"=Num_clusters,"theta_mu"=theta_mu,"theta_sig"=theta_sig,"logxi_M"=logxi_M))
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






marginal_density=function(idata,iC,inum,imu,isigma){
    pre_data=idata
    C=iC
    Num_clusters=inum
    theta_mu=imu
    theta_sig=isigma
    
    freq=as.matrix(table(C))/(length(C)+alpha)
    
    md1=rep(0,dim(idata)[1])
    for(i in 1:Num_clusters){
        pre_density=dmvnorm(x=pre_data, mean =theta_mu[[i]], sigma =theta_sig[[i]], log = FALSE)*freq[i]
        md1=md1+pre_density
    }

    md2=apply(pre_data,1,FUN = function(i){marginal_x(i,nu,lambda,mu0,Chi)})
    md2=(alpha/(length(C)+alpha))*md2
    
    return(md1+md2)
}








mciter_8=function(iC,iNC,inum,imu,isigma,iM,m=3){
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
        xi=as.numeric(X[indexi,])
        logp_list=log(NC[1:Num_clusters])+as.numeric(logxi_M[indexi,1:Num_clusters])
        possi_theta=list()
        for(i in 1:m){
            possi_theta[[i]]=generate_location_process()
        }
        new_logp=sapply(1:m, FUN = function(i){dmvnorm(xi, mean =possi_theta[[i]]$mu, sigma = possi_theta[[i]]$Sigma, log = TRUE)})
        mm=max(logp_list,new_logp)
        tr_logp_list=exp(logp_list-rep(mm,Num_clusters))
        tr_log_new=exp(new_logp-mm)
        prob_list=c(tr_log_new,tr_logp_list)/sum(tr_log_new,tr_logp_list)
        ci=sample(c(-1*seq(1,m),1:Num_clusters),size = 1,prob = prob_list)
        
        #step5
        if(ci<0){
            new_label=Num_clusters+1
            C[indexi]=new_label
            new_theta=possi_theta[[-ci]]
            theta_mu[[new_label]]=new_theta$mu
            theta_sig[[new_label]]=new_theta$Sigma
            NC[new_label]=1
            logxi_M[,new_label]=dmvnorm(x=X, mean =new_theta$mu, sigma = new_theta$Sigma, log = TRUE)
            Num_clusters=Num_clusters+1
        } else {
            C[indexi]=ci
            NC[ci]=NC[ci]+1
        }
    }
    
    logxi_M=data.frame(matrix(NA,nrow = dim(X)[1],ncol = Num_clusters))
    for(label in 1:Num_clusters){
        datai=X[(C==label),]
        new_theta=dpost_NIM_multiple(datai,nu=4,lambda=0.01,mu0=c(0,0),Chi=matrix(c(2,0,0,4),ncol = 2))
        theta_mu[[label]]=new_theta$mu
        theta_sig[[label]]=new_theta$Sigma
        logxi_M[,label]=dmvnorm(x=X, mean =new_theta$mu, sigma = new_theta$Sigma, log = TRUE)
    }
    
    return(list("C"=C,"NC"=NC,"Num_clusters"=Num_clusters,"theta_mu"=theta_mu,"theta_sig"=theta_sig,"logxi_M"=logxi_M))
}



simulation=function(Num_clusters,Mu,Sigma,Num_simu){
    datalist=list()
    mu=matrix(Mu[1:(2*Num_clusters)],ncol = 2,byrow = TRUE)
    sigma=matrix(Sigma[1:(4*Num_clusters)],ncol = 4,byrow = TRUE)
    
    for (i in 1:Num_clusters){
        mui=as.numeric(mu[i,])
        sigmai=matrix(sigma[i,]/9,ncol = 2,byrow = TRUE)
        samplei=mvrnorm(Num_simu,mui,sigmai)
        dfi=data.frame(samplei)
        dfi$C=i
        datalist[[i]]=dfi
    }
    output=do.call(rbind, datalist)
    return(output)
}









