function alpha = get_alpha_logistic(beta,xth,yth);
% alpha = get_alpha_logistic(beta,xth,yth)
% Returns alpha parameter of logistic function, given the beta parameter
% (slope) and location of a point on the curve (xth,yth).

alpha = log(yth/(1-yth)) - beta*xth;