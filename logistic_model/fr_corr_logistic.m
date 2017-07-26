function [mu,rho] = fr_corr_logistic(xth,k,yth);
% Return expected firing rate and correlation from pair of homogeneous
% neurons with logistic input-output function parameterized by threshold
% x-value xth and slope parameter k. Both of these are defined relative to
% the zero mean and unit standard deviation input distribution, assumed
% gaussian. yth is the output probability at xth e.g. 1% or 5%, etc.

alpha = get_alpha_logistic(k,xth,yth); % Get centre x value of logistic

% logistic
% logistic_fun = @(x,alpha,k) 1./(1+(exp(-k*x + alpha)));

% Function to integrate to calculate firing rate
fr_fun = @(x,alpha,k) (exp(-0.5*(x.^2))/sqrt(2*pi) ) * 1./(1+(exp(-k*x - alpha)));

% Function to integrate to calculate correlation
x2_fun = @(x,alpha,k) (exp(-0.5*(x.^2))/sqrt(2*pi) ) .* ((1./(1+(exp(-k*x - alpha)))).^2);

% Expected firing rate
mu = integral(@(x)fr_fun(x,alpha,k),-Inf,Inf);

% Expected correlation
Ex2 = integral(@(x)x2_fun(x,alpha,k),-Inf,Inf);
cov_x = Ex2 - mu.^2;
rho = cov_x/(mu*(1-mu));

