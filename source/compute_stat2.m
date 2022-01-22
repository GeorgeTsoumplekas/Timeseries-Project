function [s,a,x_pred] = compute_stat2(X,K,n,m,T,training_size,coeff)
%compute_stat2.m Function that predicts the next T values of the timeseries
%using a non-linear LAP with K neighbors model applied on a specific learning set. It also
%calculates the mean of absolute values of the prediction errors and the value a
%with which it is going to be compared with to determine whether we have a
%change point or not
%Inputs:
%   X: the timeseries examined
%   K: number of neighbors
%   n: last value of the timeseries that we consider we know
%   m: embedding dimension
%   T: number of steps for which we are making predictions
%   training_size: size of training set
%   coeff: coefficient used for the calculation of threshold a
%Outputs:
%   s: mean of the abolute values of the prediction errors
%   a: threshold of s to determine change points
%   x_pred: vector containing the predicted values returned by the LAP
%   model

    tau = 1;
    q = 0;
    x_learning = X(n-(training_size-1)-m:n);
    length(x_learning);
    a = coeff*std(x_learning);
    
    x_pred = localpredictmultistep(x_learning,training_size,tau,m,T,K,q);
    s = sum(abs(X(n+1:n+T)-x_pred(:)))/T;
end