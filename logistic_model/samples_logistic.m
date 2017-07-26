function [S] = samples_logistic(theta,yth,nneuron,nsample);

% Generate binary sample data from 5-parameter logistic model.
% y = 1/( 1+exp(-slope*x -thresh) )

% Extract parameters from input argument
thresh_mean = theta(1); % Mean threshold
thresh_sd = theta(2); % S.D. of thresholds
slope_mean = theta(3); % Mean slope
slope_sd = theta(4); % S.D. of slopes
paramcorr = theta(5); % Correlation of parameters in population

% Draw parameters
parammeanvec = [thresh_mean slope_mean];
paramcov_mat = zeros(2);
paramcov_mat(1,1) = thresh_sd^2;
paramcov_mat(2,2) = slope_sd^2;
paramcov_mat(1,2) = paramcorr*thresh_sd*slope_sd;
paramcov_mat(2,1) = paramcov_mat(1,2);
R = mvnrnd(parammeanvec,paramcov_mat,nneuron);

% Find centre point corresponding to that threshold
for i = 1:nneuron
    R(i,1) = get_alpha_logistic(R(i,1),R(i,2),yth);
end

% Draw samples
S = sparse(logical(zeros(nsample,nneuron)));
for i = 1:nsample
    pvec = (1./(1+exp(-R(:,2).*randn - R(:,1) )))';
    S(i,:) = logical(floor(rand(1,nneuron)+pvec));
end