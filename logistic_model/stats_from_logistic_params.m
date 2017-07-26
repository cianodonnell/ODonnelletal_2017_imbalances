function [y] = stats_from_logistic_params(theta,yth);

% Calculate the 3 target stats (firing rate mean and s.d., pairwise
% covariance mean) from logistic model parameters.
% y = 1/(1+exp(-k*x - alpha))
% thresh = log(yth/(1-yth)) - k*xth
% Do numerical integration over probability distributions

nd = 40; % Number of discretization points
nlim = 4; % Number of sigma to go away from mean

nd_gauss = 1e3;
nlim_gauss = 4;

nsampcorr = 1e4;

% Extract parameters from input argument
thresh_mean = theta(1); % Mean threshold
thresh_sd = theta(2); % Variance of thresholds
slope_mean = theta(3); % Mean slope
slope_sd = theta(4); % Variance of slopes
paramcorr = theta(5); % Correlation of parameters in population

% Build parameters
parammeanvec = [thresh_mean slope_mean];
paramcov_mat = zeros(2);
paramcov_mat(1,1) = thresh_sd^2;
paramcov_mat(2,2) = slope_sd^2;
paramcov_mat(1,2) = paramcorr*thresh_sd*slope_sd;
paramcov_mat(2,1) = paramcov_mat(1,2);

% Discretizaton vectors
threshvec = linspace(thresh_mean-nlim*thresh_sd,thresh_mean+nlim*thresh_sd,nd);
slopevec = linspace(max(slope_mean-nlim*slope_sd,0),slope_mean+nlim*slope_sd,nd);

% 2D pdf
[X1,X2] = meshgrid(threshvec',slopevec');
X = [X1(:) X2(:)];
p = mvnpdf(X, parammeanvec, paramcov_mat);
slopethresh_pdfmat = reshape(p,nd,nd);
slopethresh_pdfmat = slopethresh_pdfmat./sum(slopethresh_pdfmat(:));

% Input pdf (normal distribution)
invec = linspace(-nlim_gauss,nlim_gauss,nd_gauss);
outpdf = normpdf(invec);
outpdf = outpdf./sum(outpdf);

% Initialize result vector
y = zeros(1,3);

% Calculate firing rate mean, s.d.
fr_mean_given_ts = zeros(nd);
for i = 1:nd
    for j = 1:nd;
        thresh = X1(i,j);
        slope = X2(i,j);
        alpha = get_alpha_logistic(slope,thresh,yth);
        
        fr_mean_given_ts(i,j) = sum(outpdf .* 1./(1+(exp(-slope*invec - alpha))) );
    end
end
fr_mean_mat = slopethresh_pdfmat.*fr_mean_given_ts;
fr_mean = sum(fr_mean_mat(:));
y(1) = fr_mean;
fr_var_mat = slopethresh_pdfmat.*((fr_mean_given_ts-fr_mean).^2);
y(2) = sqrt(sum(fr_var_mat(:)));

% Calculate mean pairwise correlation
xvec_mat = zeros(nd,nd,nd_gauss);
for i = 1:nd
    for j = 1:nd;
        thresh = X1(i,j);
        slope = X2(i,j);
        alpha = get_alpha_logistic(slope,thresh,yth);
        
        xvec_mat(i,j,:) = 1./(1+(exp(-slope*invec - alpha)));
    end
end
var_x = fr_mean_given_ts.*(1-fr_mean_given_ts);

corrval = zeros(nsampcorr,1);
for n = 1:nsampcorr;
    a = mvnrnd(parammeanvec,paramcov_mat,2);
    [~,i] = min(abs(threshvec-a(1)));
    [~,j] = min(abs(slopevec-a(3)));
    [~,ii] = min(abs(threshvec-a(2)));
    [~,jj] = min(abs(slopevec-a(4)));

    xvec = squeeze(xvec_mat(i,j,:));
    yvec = squeeze(xvec_mat(ii,jj,:));
    cov_xy = sum(outpdf.*(xvec.*yvec)') - fr_mean_given_ts(i,j)*fr_mean_given_ts(ii,jj);
    corrval(n) = cov_xy./sqrt(var_x(i,j)*var_x(ii,jj));
end
y(3) = mean(corrval);