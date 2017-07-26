function [chisq] = errfunc(y,yhat,sigma);
% Return sum of squares error.
% y = target values
% yhat = model estimate values
% sigma = estimates of variability in each variable

chisq = sum(((y-yhat).^2)./(2*sigma.^2));

