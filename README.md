# Dirchelet-Process-and-MCMC

This folders target at the simulation and inference of direchlet mixed guassian model.

The Normal-inverse-Wishart distribution is used as the prior and the data is located in Rˆ2 with determined parameters. 
The standard Gibbs sampling method in the article is used for inference.
The deisgned data is used to verify the MCMC process.
The Chicago crime map data is used as an application of the model.


The simulation.R file generates the data by first sampling each cluster center and covariance matrix from the direchlet process 
and then sampling the data under the guassian model based on a random sample size.  The dataset contains location 
and truly cluster assignment of each data point.

The function.R file includes all the function needed for Gibbs sampling iteration.

The iteration.R file includes the predicted data designed, iteration process and the plots of MCMC results.

The designed_data.R file includes the designed data with 6 clusters and their distance can be changed.

The Chicago****.R file includes the subsets of Chicago crime data and their predicted posterior contour plot.

The match_cluster.R file switches the labels of clusters in each iteration based on the pairwise L2 distance in cluster centers between 
current iteration and previous iteration. Denote
$\alpha$

