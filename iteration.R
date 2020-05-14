##########initialization
#data=fakedata
X=data[,1:2]
theta_mu=list()
theta_sig=list()
Num_clusters=8
C=sample.int(Num_clusters, size = dim(X)[1], replace =TRUE)
NC=as.numeric(table(C))
for(i in 1:max(C)){
    Xi=X[C==i,]
    theta_mu[[i]]=apply(Xi,2,mean)
    theta_sig[[i]]=cov(Xi)
}
m=sapply(1:Num_clusters, FUN=function(i) {dmvnorm(x=X, mean =theta_mu[[i]], sigma =theta_sig[[i]], log = TRUE)})
logxi_M=data.frame(matrix(m,byrow=FALSE,nrow = dim(X)[1]))



########## the prediction points
xrange=seq(min(data$x1)*1.1,max(data$x1)*1.1,0.5)
yrange=seq(min(data$x2)*1.1,max(data$x2)*1.1,0.5)
pre_data1=rep(xrange,times=length(yrange))
pre_data2=rep(yrange,each=length(xrange))
pre_data=cbind(pre_data1,pre_data2)




############run iteration and get trace of C
Num_iteration=500
Num_cut=Num_iteration*0.1
Num_slice=10
Num_trace=floor((Num_iteration-Num_cut)/Num_slice)

C_trace=matrix(NA,nrow = Num_iteration,ncol = dim(X)[1])
NumClusters_trace=rep(0,Num_iteration)
mu_trace=matrix(NA,nrow=Num_iteration,ncol=100)
sigma_trace=matrix(NA,nrow=Num_iteration,ncol=200)
pre_marginal_trace=matrix(NA,nrow=Num_trace,ncol=dim(pre_data)[1])


temp=mciter(C,NC,Num_clusters,theta_mu,theta_sig,logxi_M)
j=1
for(i in 1:Num_iteration){
    temp=mciter(temp$C,temp$NC,temp$Num_clusters,temp$theta_mu,temp$theta_sig,temp$logxi_M)
    
    NumClusters_trace[i]=temp$Num_clusters
    C_trace[i,]=temp$C
    mul=unlist(temp$theta_mu)
    sigl=unlist(temp$theta_sig)
    mu_trace[i,1:length(mul)]=mul
    sigma_trace[i,1:length(sigl)]=sigl
    print(paste0(i,"done"))
    
    if(i>Num_cut && i%%Num_slice==0){
        m=marginal_density(idata=pre_data,iC=temp$C,inum=temp$Num_clusters,imu=temp$theta_mu,isigma=temp$theta_sig)
        #j=i/Num_slice-1
        pre_marginal_trace[j,]=m
        print(paste0(j,"get into trace"))
        j=j+1
    }
}





#############analyze the results

# original data
ggplot(data, aes(x = x1, y = x2))+
  geom_point(aes(shape=as.factor(C),colour = as.factor(C)),size=1.5)+
  scale_shape_manual(values=(1:20))+
  ggtitle("Original cluster assignment")




#initial cluster assignment
ggplot(X, aes(x = x1, y = x2, shape=as.factor(C),colour = as.factor(C)))+
  geom_point()+
  ggtitle("initial clusters assignment of MCMC")



# MCMC results
iter=500
ggplot(data, aes(x = data$x1, y = data$x2))+
    geom_point(aes(shape= as.factor(C_trace[iter,]),colour = as.factor(C_trace[iter,])))+
    scale_shape_manual(values=(1:20))+
    ggtitle("MCMC clusters assignment")


# number of clusters
plot(NumClusters_trace,main = "Num_clusters in MCMC")


# contour plot
pre_marginal=apply(pre_marginal_trace, 2, sum)
iter=500
df <- data.frame("x"=pre_data1,"y"=pre_data2,"z"=pre_marginal)
mudf=data.frame(matrix(mu_trace[iter,1:(2*NumClusters_trace[iter])],ncol = 2,byrow = TRUE))
sigmadf=data.frame(matrix(sigma_trace[iter,1:(4*NumClusters_trace[iter])],ncol = 4,byrow = TRUE))
ecidf=simulation(Num_clusters=NumClusters_trace[iter],Mu=mu_trace[iter,],Sigma=sigma_trace[iter,],Num_simu=100)

ggplot(data, aes(x = data$x1, y = data$x2))+
  geom_point(aes(shape= as.factor(C_trace[iter,]),colour = as.factor(C_trace[iter,])),size=1,alpha=0.6)+
  scale_shape_manual(values=(1:20))+
  geom_point(data=mudf,aes(x=X1,y=X2,colour=as.factor(seq(1,NumClusters_trace[iter])),shape=as.factor(seq(1,NumClusters_trace[iter]))),size=3,alpha=1)+
  stat_contour(data=df,geom="polygon",aes(x=x,y=y,z=z,fill=stat(level)),alpha=0.5)+
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  stat_ellipse(data=ecidf,aes(x=X1, y=X2,color=as.factor(C)),type = "norm")+
  ggtitle("final cluster assignments and post-predicted density contour plot")







# contour plot and MCMC for a certain iteration
iter=350
con_iter=(iter-50)/10
df <- data.frame("x"=pre_data1,"y"=pre_data2,"z"=pre_marginal_trace[con_iter,])

mudf=data.frame(matrix(mu_trace[iter,1:(2*NumClusters_trace[iter])],ncol = 2,byrow = TRUE))
sigmadf=data.frame(matrix(sigma_trace[iter,1:(4*NumClusters_trace[iter])],ncol = 4,byrow = TRUE))

ecidf=simulation(Num_clusters=NumClusters_trace[iter],Mu=mu_trace[iter,],Sigma=sigma_trace[iter,],Num_simu=100)

ggplot(data, aes(x = data$x1, y = data$x2))+
  geom_point(aes(shape= as.factor(C_trace[iter,]),colour = as.factor(C_trace[iter,])),size=1.3,alpha=0.6)+
  scale_shape_manual(values=(1:20))+
  geom_point(data=mudf,aes(x=X1,y=X2,colour=as.factor(seq(1,NumClusters_trace[iter])),shape=as.factor(seq(1,NumClusters_trace[iter]))),size=3,alpha=1)+
  stat_contour(data=df,geom="polygon",aes(x=x,y=y,z=z,fill=stat(level)),alpha=0.3)+
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  stat_ellipse(data=ecidf,aes(x=X1, y=X2,color=as.factor(C)),type = "norm")+
  ggtitle(paste("Clusters/Cluster Assignments/Contourplot of iteration",as.character(iter)))






pdf(file = "My Plot.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) 

dev.off()











