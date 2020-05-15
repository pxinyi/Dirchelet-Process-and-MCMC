# reorganize the
library(pdist)


maxc=max(NumClusters_trace)
re=data.frame(matrix(ncol=2*maxc,nrow=Num_iteration))
names=NULL
for(i in 1:maxc){
    name=c(paste0("x",as.character(i)),paste0("y",as.character(i)))
    names=c(names,name)
}
colnames(re)=names


n1=NumClusters_trace[1]
re[1,1:(2*n1)]=mu_trace[1,1:(2*n1)]
for (i in 2:500){
    temp=match_cluster((i-1),i)
    n=dim(temp)[1]
    re[i,1:(2*n)]=as.numeric(t(temp))
}





# example of show the trace of the 1st cluster centers

# distribution
ggplot(data=re)+
    geom_point(aes(x=x1,y=y1),size=0.6)

# traceplot in 2D
fig <- plot_ly(data = re, x = ~x1, y = ~y1, type = 'scatter', mode = 'lines') 
fig


pdf(file = "diagnostic.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5)

# histogram for each dimension
par(mfrow=c(3,4))
hist(re$x1)
hist(re$y1)
hist(re$x2)
hist(re$y3)
hist(re$x3)
hist(re$y3)
hist(re$x4)
hist(re$y4)
hist(re$x5)
hist(re$y5)
hist(re$x6)
hist(re$y6)

# autocorrelation for each dimension
par(mfrow=c(3,4))
acf(x=re$x1)
acf(x=re$y1)
acf(x=re$x2)
acf(x=re$y2)
acf(x=re$x3)
acf(x=re$y3)
acf(x=re$x4)
acf(x=re$y4)
acf(x=re$x5)
acf(x=re$y5)
acf(x=re$x6)
acf(x=re$y6)

dev.off()





match_cluster=function(iter1,iter2){
    N1=NumClusters_trace[iter1]
    N2=NumClusters_trace[iter2]
    
    centerM1=matrix(mu_trace[iter1,1:(2*N1)],ncol=2,byrow = TRUE)
    centerM2=matrix(mu_trace[iter2,1:(2*N2)],ncol=2,byrow = TRUE)
    newcenterM2=matrix(rep(0,(2*N2)),ncol = 2)
    dist=as.matrix(pdist(centerM1, centerM2))
    
    Aval=seq(1,N2)
    for (index1 in 1:min(N1,N2)){
        order=sort(dist[index1,],index.return=TRUE)$ix
        for(j in 1:N2){
            if (is.element(order[j],Aval)){
                index2=order[j]
                break
            }else{
                j=j+1
            }
        }
        Aval[index2]=NA
        newcenterM2[index1,]=centerM2[index2,]
    }
    if(N1<N2){
        Aval_rest=Aval[!is.na(Aval)]
        for (i in 1:length(Aval_rest)){
            newcenterM2[(N1+i),]=centerM2[Aval_rest[i],]
        }
    }
    return(newcenterM2) 
}




















iter1=1
iter2=2

N1=NumClusters_trace[iter1]
N2=NumClusters_trace[iter2]

centerM1=matrix(mu_trace[iter1,1:(2*N1)],ncol=2,byrow = TRUE)
centerM2=matrix(mu_trace[iter2,1:(2*N2)],ncol=2,byrow = TRUE)
newcenterM2=matrix(rep(0,(2*N2)),ncol = 2)
dist=as.matrix(pdist(centerM1, centerM2))


Aval=seq(1,N2)
for (index1 in 1:min(N1,N2)){
    order=sort(dist[index1,],index.return=TRUE)$ix
    print(Aval)
    for(j in 1:N2){
        if (is.element(order[j],Aval)){
            index2=order[j]
            break
        }else{
            j=j+1
        }
    }
    Aval[index2]=NA
    print(index2)
    print(Aval)
    print("********")
    newcenterM2[index1,]=centerM2[index2,]
}

if(N1<N2){
    Aval_rest=Aval[!is.na(Aval)]
    for (i in 1:length(Aval_rest)){
        newcenterM2[(N1+i),]=centerM2[Aval_rest[i],]
    }
}






