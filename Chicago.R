# Chicago Crime Map
Chicago=read.csv(file = "Crimes_Map.csv")
Chicago=na.omit(Chicago)

# to reduce the size of the dataset, choose only the crime map for ASSAULT.
Chicagodata=na.omit(Chicago[Chicago$PRIMARY.DESCRIPTION=="ASSAULT",])


ggplot(data=Chicago)+
    geom_point(aes(x=LATITUDE,y=LONGITUDE),size=0.5,alpha=0.1)


ggplot(data=Chicagodata)+
    geom_point(aes(x=LATITUDE,y=LONGITUDE),size=0.5,alpha=0.2)


index_set=sample(dim(Chicagodata)[1],500)
ggplot(data=Chicagodata[index_set,])+
    geom_point(aes(x=LATITUDE,y=LONGITUDE),size=0.5,alpha=0.2)



# set parameters 
alpha=3
p=2
nu=4
Chi=0.1*matrix(c(3,0,0,5),2,2)
mu0=c(0,0)
lambda=0.1


# get the dataset of size 500 and normalize it
data=Chicagodata[index_set,c("LATITUDE","LONGITUDE")]
colnames(data)=c("x1","x2")
data$x1=30*(data$x1-mean(data$x1))
data$x2=30*(data$x2-mean(data$x2))
ggplot(data=data)+
    geom_point(aes(x=x1,y=x2),size=0.5,alpha=0.2)+
    ggtitle("normalized sample data for assualts in Chicago 2020")



# initialization
X=data
theta_mu=list()
theta_sig=list()
Num_clusters=sample(5:10,1)
C=sample.int(Num_clusters, size = dim(X)[1], replace =TRUE)
NC=as.numeric(table(C))
for(i in 1:max(C)){
    Xi=X[C==i,]
    theta_mu[[i]]=apply(Xi,2,mean)
    theta_sig[[i]]=cov(Xi)
}
m=sapply(1:Num_clusters, FUN=function(i) {dmvnorm(x=X, mean =theta_mu[[i]], sigma =theta_sig[[i]], log = TRUE)})
logxi_M=data.frame(matrix(m,byrow=FALSE,nrow = dim(X)[1]))


xrange=seq(min(data$x1),max(data$x1),length.out = 50)
yrange=seq(min(data$x2),max(data$x2),length.out = 50)
pre_data1=rep(xrange,times=length(yrange))
pre_data2=rep(yrange,each=length(xrange))
pre_data=cbind(pre_data1,pre_data2)



Num_iteration=200
Num_cut=Num_iteration*0.1
Num_slice=10
Num_trace=floor((Num_iteration-Num_cut)/Num_slice)

C_trace=matrix(NA,nrow = Num_iteration,ncol = dim(X)[1])
NumClusters_trace=rep(0,Num_iteration)
mu_trace=matrix(NA,nrow=Num_iteration,ncol=200)
sigma_trace=matrix(NA,nrow=Num_iteration,ncol=400)
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





#initial cluster assignment
ggplot(X, aes(x = x1, y = x2, shape=as.factor(C),colour = as.factor(C)))+
    geom_point()+
    scale_shape_manual(values=(1:20))+
    ggtitle("initial clusters assignment of MCMC")


iter=100
ggplot(data, aes(x = data$x1, y = data$x2))+
    geom_point(aes(shape= as.factor(C_trace[iter,]),colour = as.factor(C_trace[iter,])))+
    scale_shape_manual(values=(1:20))+
    ggtitle("MCMC clusters assignment of iteration 100")



pre_marginal=apply(pre_marginal_trace, 2, sum)
iter=200
df <- data.frame("x"=pre_data1,"y"=pre_data2,"z"=pre_marginal)
mudf=data.frame(matrix(mu_trace[iter,1:(2*NumClusters_trace[iter])],ncol = 2,byrow = TRUE))
sigmadf=data.frame(matrix(sigma_trace[iter,1:(4*NumClusters_trace[iter])],ncol = 4,byrow = TRUE))
ecidf=simulation(Num_clusters=NumClusters_trace[iter],Mu=mu_trace[iter,],Sigma=sigma_trace[iter,],Num_simu=100)

ggplot(data, aes(x = data$x1, y = data$x2))+
    geom_point(aes(shape= as.factor(C_trace[iter,]),colour = as.factor(C_trace[iter,])),size=1,alpha=0.6)+
    scale_shape_manual(values=(1:25))+
    geom_point(data=mudf,aes(x=X1,y=X2,colour=as.factor(seq(1,NumClusters_trace[iter])),shape=as.factor(seq(1,NumClusters_trace[iter]))),size=3,alpha=1)+
    stat_contour(data=df,geom="polygon",aes(x=x,y=y,z=z,fill=stat(level)),alpha=0.5)+
    scale_fill_distiller(palette = "Spectral", direction = -1)+
    stat_ellipse(data=ecidf,aes(x=X1, y=X2,color=as.factor(C)),type = "norm")+
    ggtitle("final cluster assignments and post-predicted density contour plot")



pdf(file = "ChicagoCrime__contour plot.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5)

dev.off()



library(plotly)
contour_range= list(end = max(df$z), size =10 , start = min(df$z))
fig <- plot_ly(
    type = 'contour',
    x=df$x,
    y=df$y,
    z=df$z,
    contours = list(showlabels = TRUE,coloring = 'heatmap')
)
fig
