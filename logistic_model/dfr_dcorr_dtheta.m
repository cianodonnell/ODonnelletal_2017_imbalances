function [dy_dthresh,dy_dslope] = dfr_dcorr_dtheta(theta,delta_theta)
% [dy_dthresh,dy_dslope] = DFR_CORR_DTHETA(theta,delta_theta);
% Numerically compute the gradient in firing rate and correlation with
% respect to the threshold and slope of the logistic model.
% Theta is the parameter vec at which to compute the gradient
% delta_theta is the fractional change in theta parameters to compute
% gradient with.

yth = 0.01;

% Record baseline
theta0 = theta;

% d_dthresh
theta = theta0; % reset
theta(1) = theta(1) + delta_theta*abs(theta(1));
y_dthresh_pos = stats_from_logistic_params(theta,yth);
theta = theta0; % reset
theta(1) = theta(1) - delta_theta*abs(theta(1));
y_dthresh_neg = stats_from_logistic_params(theta,yth);
dy_dthresh = ( y_dthresh_pos - y_dthresh_neg ) / (2*delta_theta*abs(theta(1)));

% d_dslope
theta = theta0; % reset
theta(3) = theta(3) + delta_theta*abs(theta(3));
y_dslope_pos = stats_from_logistic_params(theta,yth);
theta = theta0; % reset
theta(3) = theta(3) - delta_theta*abs(theta(3));
y_dslope_neg = stats_from_logistic_params(theta,yth);
dy_dslope = ( y_dslope_pos - y_dslope_neg ) / (2*delta_theta*abs(theta(3)));