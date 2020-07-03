# load the HOMICIDE crime between year 2012 and 2017
HOMICIDE=read.csv(file = "code/Crimes_HOMICIDE.csv")
Chicago=data.frame(HOMICIDE[,c("Latitude","Longitude")])
colnames(Chicago)=c("x1","x2")            


# Rescale the data ---- data3 is the full dataset and data4 is randomly selected subset
library(scales)
data3=data.frame(Chicago)
ratio=(max(Chicago$x2)-min(Chicago$x2))/(max(Chicago$x1)-min(Chicago$x1))
data3$x1=rescale(Chicago$x2,to=c(-1,1))
data3$x2=rescale(Chicago$x1,to=c(-ratio,ratio))
ggplot(data=data3)+
    geom_point(aes(x=x1,y=x2),size=1,alpha=0.3)+
    coord_fixed(ratio=1)+
    ggtitle("rescaled HOMICIDE data")


data4=data3[sample(1:3000,100),]
ggplot(data=data4)+
    geom_point(aes(x=x1,y=x2),size=1,alpha=0.3)+
    coord_fixed(ratio=1)+
    ggtitle("rescaled HOMICIDE data (subset)")


# large alpha in order to allow more clusters
alpha=10
p=2
nu=4
Chi=matrix(c(0.01,0,0,0.01),2,2)
mu0=c(0,0)
lambda=0.01


data=data4



# INITIALIZATION ###################################################################################
X=as.matrix(data[,1:2],ncol=2)
theta_mu=list()
theta_sig=list()
Num_clusters=sample(seq(2,5),1)
C=sample.int(Num_clusters, size = dim(X)[1], replace =TRUE)
NC=as.numeric(table(C))
for(i in 1:max(C)){
    Xi=X[C==i,]
    theta_mu[[i]]=apply(Xi,2,mean)
    theta_sig[[i]]=cov(Xi)
}
logxi_M=matrix(ncol = Num_clusters,nrow = dim(X)[1])
for(i in 1:Num_clusters){
    logxi_M[,i]=dmvnorm(x=X, mean =theta_mu[[i]], sigma =theta_sig[[i]], log = TRUE)
}






# MCMC ###########################################################################################
Num_iteration=500
Max_NumClusters=100
C_trace=matrix(NA,nrow = Num_iteration,ncol = dim(X)[1])
NumClusters_trace=rep(0,Num_iteration)
NC_trace=matrix(NA,nrow=Num_iteration,ncol=Max_NumClusters)
mu_trace=matrix(NA,nrow=Num_iteration,ncol=Max_NumClusters*2)
sigma_trace=matrix(NA,nrow=Num_iteration,ncol=Max_NumClusters*4)


ptm=proc.time()
temp=mciter_8(C,NC,Num_clusters,theta_mu,theta_sig,logxi_M,10)
for(i in 1:Num_iteration){
    temp=mciter_8(temp$C,temp$NC,temp$Num_clusters,temp$theta_mu,temp$theta_sig,temp$logxi_M,10)
    t=temp$Num_clusters
    NumClusters_trace[i]=t
    NC_trace[i,1:t]=temp$NC[1:t]
    C_trace[i,]=temp$C
    mul=unlist(temp$theta_mu)
    sigl=unlist(temp$theta_sig)
    mu_trace[i,1:length(mul)]=mul
    sigma_trace[i,1:length(sigl)]=sigl
    print(paste0("iteration",i,"done"))
    print("***************************************************")
}
proc.time() - ptm


# COMPUTE POSTERIOR PREDICTED DENSITY ##########################################################
Num_grid=100
xrange=seq(min(data$x1),max(data$x1),length.out = Num_grid)
yrange=seq(min(data$x2),max(data$x2),length.out = Num_grid)
x1=rep(xrange,times=length(yrange))
x2=rep(yrange,each=length(xrange))
pre_data=cbind(x1,x2)


Num_cut=Num_iteration*0.1
Num_slice=10
Seq_index=seq(Num_cut,Num_iteration,Num_slice)
pre_marginal_trace=matrix(NA,nrow=length(Seq_index),ncol=dim(pre_data)[1])


j=1
for(i in Seq_index){
    Num_clusters=NumClusters_trace[i]
    C=C_trace[i,]
    theta_mu=mu_trace[i,1:(2*Num_clusters)]
    theta_sigma=sigma_trace[i,1:(4*Num_clusters)]
    pre_marginal=marginal_density(pre_data,C,Num_clusters,theta_mu,theta_sigma)
    pre_marginal_trace[j,]=pre_marginal
    j=j+1
}








# PLOTS #######################################################################################
pdf(file = "Chicago(subset)---A8--k=3.pdf") 


ggplot(data, aes(x = x1, y = x2))+
    coord_fixed(ratio=1)+
    geom_point(size=1,alpha=0.5)

# plot1---trace of the number of clusters
plot(NumClusters_trace,main = "Num_clusters in MCMC")





# plot2---cluster assignment of each iteration
iter=400
iter_num=NumClusters_trace[iter]
mudf=data.frame(matrix(mu_trace[iter,1:(2*iter_num)],ncol = 2,byrow = TRUE))
sigmadf=data.frame(matrix(sigma_trace[iter,1:(4*iter_num)],ncol = 4,byrow = TRUE))
ecidf=simulation(mudf,sigmadf,200)

ggplot(data, aes(x = x1, y = x2))+
    geom_point(aes(shape= as.factor(C_trace[iter,]),colour = as.factor(C_trace[iter,])))+
    coord_fixed(ratio=1)+
    scale_shape_manual(values=(1:25))+
    geom_point(data=mudf,aes(x=X1,y=X2,colour=as.factor(seq(1,iter_num)),shape=as.factor(seq(1,iter_num))),size=3,alpha=1)+
    stat_ellipse(data=ecidf,aes(x=x1, y=x2,color=as.factor(C)),type = "norm")+
    ggtitle(paste("Clusters assignments of iteration",as.character(iter)))





# plot3---the contour plot
pre_marginal=apply(pre_marginal_trace, 2, sum)
df <- data.frame(pre_data)
df$z=pre_marginal
ggplot(data, aes(x = x1, y = x2))+
    geom_point(size=1,alpha=1)+
    coord_fixed(ratio=1)+
    stat_contour(data=df,geom="polygon",aes(x=x1,y=x2,z=z,fill=stat(level)),alpha=0.5)+
    scale_fill_distiller(palette = "Spectral", direction = -1)+
    ggtitle("contour plot of predicted posterior density")


library(plotly)
size0=(max(df$z)-min(df$z))/20
contour_range= list(end = max(df$z), size =size0, start = min(df$z),showlabels = TRUE)
fig <- plot_ly(
    type = 'contour',
    x=df$x1,
    y=df$x2,
    z=df$z,
    contours = contour_range
)
fig<-fig %>% layout(title = '(1)contour plot of of predicted posterior density')
fig

fig1<- plot_ly(x = xrange, y = yrange, z = matrix(df$z,ncol=Num_grid,byrow=TRUE)) %>% 
    add_surface(contours = list(z = list(show=TRUE,usecolormap=TRUE,highlightcolor="#ff0000",project=list(z=TRUE))))

fig1 <- fig1 %>% layout(title = '(2)contour plot of of predicted posterior density',
                        scene = list(
                            camera=list(
                                eye = list(x=-1.5, y=-1.5, z=0.5)
                            )
                        )
)
fig1


# final plot
iter=500
iter_num=NumClusters_trace[iter]
mudf=data.frame(matrix(mu_trace[iter,1:(2*iter_num)],ncol = 2,byrow = TRUE))
sigmadf=data.frame(matrix(sigma_trace[iter,1:(4*iter_num)],ncol = 4,byrow = TRUE))
ecidf=simulation(mudf,sigmadf,200)

ggplot(data, aes(x = x1, y = x2))+
    geom_point(aes(shape= as.factor(C_trace[iter,]),colour = as.factor(C_trace[iter,])),size=1.3,alpha=0.6)+
    scale_shape_manual(values=(1:20))+
    coord_fixed(ratio=1)+
    geom_point(data=mudf,aes(x=X1,y=X2,colour=as.factor(seq(1,NumClusters_trace[iter])),shape=as.factor(seq(1,NumClusters_trace[iter]))),size=3,alpha=1)+
    stat_ellipse(data=ecidf,aes(x=x1, y=x2,color=as.factor(C)),type = "norm")+
    stat_contour(data=df,geom="polygon",aes(x=x1,y=x2,z=z,fill=stat(level)),alpha=0.3)+
    scale_fill_distiller(palette = "Spectral", direction = -1)+
    ggtitle("final plot")

dev.off()

