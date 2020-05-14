# fake data ----- the data is designed to verify the MCMC process in iteration.R

# the distance object here controls the how spread out the clusters and the large distance, the more spread out  
distance=5

# set the center of clusters
mu1=c(-5,2)*distance
mu2=c(3,5)*distance
mu3=c(-2,1)*distance
mu4=c(0,4)*distance
mu5=c(-2.5,-1)*distance
mu6=c(1.3,3)*distance

# set the corvariance matrix
sigma1=matrix(c(4,0.4,0.4,1),ncol = 2)
sigma2=matrix(c(1.5,0.3,0.3,4),ncol = 2)
sigma3=matrix(c(3,-0.2,-0.2,2.5),ncol = 2)
sigma4=matrix(c(9,2.5,2.5,7),ncol = 2)
sigma5=matrix(c(8,-3.8,-3.8,3),ncol = 2)
sigma6=matrix(c(5.6,3.9,3.9,4.4),ncol = 2)


# generate data for each cluster based on the random sample size
sample_size=sample(20:100,1)
s1=mvrnorm(sample_size,mu1,sigma1)
df1=data.frame(s1)
df1$C=1

sample_size=sample(20:100,1)
s2=mvrnorm(sample_size,mu2,sigma2)
df2=data.frame(s2)
df2$C=2

sample_size=sample(20:100,1)
s3=mvrnorm(sample_size,mu3,sigma3)
df3=data.frame(s3)
df3$C=3

sample_size=sample(20:100,1)
s4=mvrnorm(sample_size,mu4,sigma4)
df4=data.frame(s4)
df4$C=4

sample_size=sample(20:100,1)
s5=mvrnorm(sample_size,mu5,sigma5)
df5=data.frame(s5)
df5$C=5

sample_size=sample(20:100,1)
s6=mvrnorm(sample_size,mu6,sigma6)
df6=data.frame(s6)
df6$C=6



# create the fake data as dataframe
fakedatal=list(df1,df2,df3,df4,df5,df6)
fakedata=do.call(rbind, fakedatal)
colnames(fakedata)=c("x1","x2","C")



pdf(file = "My Plot_5.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) 


ggplot(fakedata, aes(x = x1, y = x2))+
    geom_point(aes(shape=as.factor(C),colour = as.factor(C)),size=1)+
    scale_shape_manual(values=(1:20))+
    stat_ellipse(aes(x=x1, y=x2,color=as.factor(C)),type = "norm",alpha=0.7)
    ggtitle("Fake_data----Original cluster assignment and covariance matrix")



pre_marginal=apply(pre_marginal_trace, 2, sum)
iter=500
df <- data.frame("x"=pre_data1,"y"=pre_data2,"z"=pre_marginal)
mudf=data.frame(matrix(mu_trace[iter,1:(2*NumClusters_trace[iter])],ncol = 2,byrow = TRUE))
sigmadf=data.frame(matrix(sigma_trace[iter,1:(4*NumClusters_trace[iter])],ncol = 4,byrow = TRUE))
ecidf=simulation(Num_clusters=NumClusters_trace[iter],Mu=mu_trace[iter,],Sigma=sigma_trace[iter,],Num_simu=100)

ggplot(data, aes(x = x1, y = x2))+
    geom_point(aes(shape= as.factor(C_trace[iter,]),colour = as.factor(C_trace[iter,])),size=1,alpha=0.7)+
    scale_shape_manual(values=(1:20))+
    geom_point(data=mudf,aes(x=X1,y=X2,colour=as.factor(seq(1,NumClusters_trace[iter])),shape=as.factor(seq(1,NumClusters_trace[iter]))),size=3,alpha=1)+
    stat_contour(data=df,geom="polygon",aes(x=x,y=y,z=z,fill=stat(level)),alpha=0.6)+
    scale_fill_distiller(palette = "Spectral", direction = -1)+
    stat_ellipse(data=ecidf,aes(x=X1, y=X2,color=as.factor(C)),type = "norm")+
    ggtitle("Fake_data----final cluster assignments and post-predicted density contour plot")

dev.off()

