# SVRCD-Algorithm
A Stochastic Variance Reduced Coordinate Descent Algorithm for Learning Sparse Bayesian Network from Discrete High-Dimensional Data

SVRCD stands for Stochastic Variance-Reduced Coordinate Descent which is a Bayesian Network (BN) structure learning algorithm. SVRCD learns a sparse BN using High-Dimensional Discrete Data. We have implemented SVRCD in R. We have used discretecdAlgorithm package in R and modify the optimization and the score function in this package to implement SVRCD. You can look at our proposed algorithm in discretecdAlgorithm\src\dCD.cpp. CDOnePoint function have been replaced with SVRCD code.

# How to run SVRCD

The algorithm have been tested using sinthetic data sampled from bipartite graphs, scale-free graphs, and random graphs. There are three R file in this repository which you can run them to sample data and test the algorithm using data.

# Data

bipartite.R samples data from bipartite graph and then calls the algorithm using sampled data. scale-free.R bipartite.R samples data from scale-free network and then calls the algorithm using sampled data. dc.R bipartite.R samples data from random graph and then calls the algorithm using sampled data.
