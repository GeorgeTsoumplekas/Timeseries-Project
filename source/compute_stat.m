function [s,a,x_pred]= compute_stat(X, p, q, n, T, training_size,coeff)
%compute_stat.m Function that predicts the next T values of the timeseries
%using an ARMA(p,q) model applied on a specific learning set. It also
%calculates the mean of absolute values of the prediction errors and the value a
%with which it is going to be compared with to determine whether we have a
%change point or not
%Inputs:
%   X: the timeseries examined
%   p: order of the AR part of the ARMA model
%   q: order of the MA part of the ARMA model
%   n: last value of the timeseries that we consider we know
%   T: number of steps for which we are making predictions
%   training_size: size of training set
%   coeff: coefficient used for the calculation of threshold a
%Outputs:
%   s: mean of the abolute values of the prediction errors
%   a: threshold of s to determine change points
%   x_pred: vector containing the predicted values returned by the ARMA
%   model

    x_learning=X(n-(training_size-1):n);
    a=coeff*std(x_learning);
    
    x_pred=predictARMAmultistep(x_learning, training_size, p, q, T, '');
    s = sum(abs(X(n+1:n+T)-x_pred(:)))/T;
end