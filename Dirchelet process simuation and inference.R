library(LaplacesDemon)
library(MASS)
library(mvtnorm)
library(CholWishart)
library(ggplot2)
library(reshape2)
library(ggforce)


# setting hyperparamaters

# 1st time paramater
"
alpha=3
p=2
nu=4
Chi=matrix(c(3,0,0,5),2,2)
mu0=c(0,0)
lambda=0.05
"



# 2nd time parameter------ increase Chi and lambda so that the cluster are closer 
"
alpha=3
p=2
nu=4
Chi=3*matrix(c(3,0,0,5),2,2)
mu0=c(0,0)
lambda=0.1
"

# 3rd time parameter------ increase Chi and lambda so that the cluster are closer 
"
alpha=3
p=2
nu=4
Chi=3*matrix(c(3,0,0,5),2,2)
mu0=c(0,0)
lambda=1
"


# 4th time parameter------ increase Chi and lambda so that the cluster are closer 
alpha=3
p=2
nu=4
Chi=5*matrix(c(3,0,0,5),2,2)
mu0=c(0,0)
lambda=0.1






# generating a sample from G0
generate_location_process=function(){      
    Sigma_sample1=matrix(rInvWishart(n=1, df=nu, Sigma=Chi),ncol = 2,byrow = TRUE)
    mu_sample1=mvrnorm(1,mu0,Sigma_sample1/lambda)
    return(list("mu"=mu_sample1,"Sigma"=Sigma_sample1))
}


# generating weights and locations
stick_breaking_process = function(num_weights, alpha) { 
    betas = rbeta(num_weights, 1, alpha)
    remaining_stick_lengths = c(1, cumprod(1 - betas))[1:num_weights]
    weights = remaining_stick_lengths * betas
    return(weights)
}
max_num_weights=200
G_weights=stick_breaking_process(max_num_weights,alpha)



# generating the index vector C
num_observation=200
C1=sample.int(length(G_weights), size = num_observation, replace =TRUE, prob = G_weights)
temp1=unique(C1)
num_clusters=length(temp1)


# generate data accordingly
datalist=list() # list of dataframes
true_mu=list()
true_sig=list()
for(i in temp1){
    location=generate_location_process()
    m=location$mu
    s=location$Sigma
    true_mu[[i]]=m
    true_sig[[i]]=s
    Num_sample=sum(C1==i)
    samplei=matrix(mvrnorm(Num_sample,as.numeric(m),s),ncol=2)
    dfi=data.frame(samplei)
    colnames(dfi)=c("x1","x2")
    dfi$C=i
    datalist[[i]]=dfi
}
data=do.call(rbind, datalist)


# plot the simulated data
ggplot(data, aes(x = data$x1, y = data$x2, colour = as.factor(data$C)))+
    geom_point()+
    scale_color_discrete(name="clusters")
                                     
                                     





