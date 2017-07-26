function [theta_hat,y_hat,y] = mcmc_logistic_3param(movie)
% Fit the parameters of a logistic model to three statistics of the
% data.
% 
% 3 parameters to match:
% 1) mean firing rate
% 2) s.d. of firing rates
% 3) mean pairwise correlation
%
% 5 parameters to fit. The threshold-linear function for each neuron in the population
% is drawn from a 2-D multivariate gaussian. Parameters include
% 1) threshold mean
% 2) threshold s.d.
% 3) slope mean
% 4) slope s.d.
% 5) threshold-slope correlation

% PARAMETERS

dt = 4; % Timebin size
nsteps = 2e3; % Number of steps for MCMC
burntime = 0; % Number steps for burntime
%avec = 0.01*ones(1,5); % S.d. of random MCMC parameter change size
avec = [0.025 0.005 0.02 0.01 0.01]; % Variance of candidate samples for each parameter
yth = 0.01; % Threshold output at which to calculate input threshold

plotstep = 1; % plot every xth step

% Get basic stats from movie file
%[~,~,y] = get_basicstats_bin(movie,dt);
[~,~,y] = get_basicstats_prob_3param(movie,dt);

% Initial guess for model parameters
theta1 = [-2 0.1 1 0.1 0];

% Get initial yhat
yhat1 = stats_from_logistic_params(theta1,yth);

% Error
chi1 = errfunc(y,yhat1,0.1*ones(1,3));

% Do MCMC
theta_save = zeros(nsteps-burntime,5);
yhat_save = zeros(nsteps-burntime,3);
chi_save = zeros(nsteps-burntime,1);
for i = 1:nsteps
    
%     if mod(i-burntime,100)==0
%         fprintf('%1.0f/%1.0f steps\n',i-burntime,nsteps-burntime )
%         
%         figure(1)
%         clf
%         plot([1:nsteps-burntime],chi_save)
%         xlabel('Timestep')
%         ylabel('Error')
%         ylim([0 1.1*max(chi_save)]);
%         hold off
%         drawnow
%     end
    
    theta2 = theta1 + avec.*randn(1,5); % New parameter guess
    theta2(2:4) = max(theta2(2:4),1e-3);
    yhat2 = stats_from_logistic_params(theta2,yth); % New yhats
    chi2 = errfunc(y,yhat2,0.1*ones(1,3)); % Error on new
    ratio = exp(-chi2+chi1); % Error ratio
    
%     if rand < ratio
    if ratio > 1
        theta1 = theta2;
        chi1 = chi2;
    end
    
    if i > burntime
        theta_save(i-burntime,:) = theta1;
        yhat_save(i-burntime,:) = yhat2;
        chi_save(i-burntime) = chi1;
    end
    
end

% Get best fit parameters
[~,minind] = min(chi_save);
theta_hat = theta_save(minind,:);

% Get target stats
y_hat = stats_from_logistic_params(theta_hat,yth);

%************
% PLOT
%************
% figure(1)
% semilogy([1:nsteps-burntime],chi_save)
% xlabel('Timestep')
% ylabel('Error')
% ylim([0 1.1*max(chi_save)]);
% 
% 
% % Plot target trajectories
% figure(2);
% 
% subplot(221), plot(yhat_save(1:plotstep:end,1),yhat_save(1:plotstep:end,2));
% hold on
% plot(y(1),y(2),'ro')
% plot(yhat_save(1,1),yhat_save(1,2),'go')
% plot(y_hat(1),y_hat(2),'co')
% hold off
% xlabel('Rate mean')
% ylabel('Rate s.d.')
% 
% subplot(222), plot(yhat_save(1:plotstep:end,1),yhat_save(1:plotstep:end,3));
% hold on
% plot(y(1),y(3),'ro')
% plot(yhat_save(1,1),yhat_save(1,3),'go')
% plot(y_hat(1),y_hat(3),'co')
% hold off
% xlabel('Rate mean')
% ylabel('Corr mean')
% 
% subplot(224), plot(yhat_save(1:plotstep:end,2),yhat_save(1:plotstep:end,3));
% hold on
% plot(y(2),y(3),'ro')
% plot(yhat_save(1,2),yhat_save(1,3),'go')
% plot(y_hat(2),y_hat(3),'co')
% hold off
% xlabel('Rate s.d.')
% ylabel('Corr mean')
% 
% % figure(3)
% % plot3(yhat_save(:,1),yhat_save(:,2),yhat_save(:,3));
% % xlabel('Rate mean')
% % ylabel('Rate s.d.')
% % zlabel('Corr mean')
% % hold on
% % plot3(y(1),y(2),y(3),'ro')
% % plot3(yhat_save(1,1),yhat_save(1,2),yhat_save(1,3),'go')
% % plot3(y_hat(1),y_hat(2),y_hat(3),'co')
% % hold off
% 
% figure(4)
% clf
% subplot(441), plot(theta_save(:,1),theta_save(:,2))
% hold on
% plot(theta_save(1,1),theta_save(1,2),'go')
% plot(theta_save(end,1),theta_save(end,2),'co')
% subplot(442), plot(theta_save(:,1),theta_save(:,3))
% hold on
% plot(theta_save(1,1),theta_save(1,3),'go')
% plot(theta_save(end,1),theta_save(end,3),'co')
% subplot(443), plot(theta_save(:,1),theta_save(:,4))
% hold on
% plot(theta_save(1,1),theta_save(1,4),'go')
% plot(theta_save(end,1),theta_save(end,4),'co')
% subplot(444), plot(theta_save(:,1),theta_save(:,5))
% hold on
% plot(theta_save(1,1),theta_save(1,5),'go')
% plot(theta_save(end,1),theta_save(end,5),'co')
% subplot(446), plot(theta_save(:,2),theta_save(:,3))
% hold on
% plot(theta_save(1,2),theta_save(1,3),'go')
% plot(theta_save(end,2),theta_save(end,3),'co')
% subplot(447), plot(theta_save(:,2),theta_save(:,4))
% hold on
% plot(theta_save(1,2),theta_save(1,4),'go')
% plot(theta_save(end,2),theta_save(end,4),'co')
% subplot(448), plot(theta_save(:,2),theta_save(:,5))
% hold on
% plot(theta_save(1,2),theta_save(1,5),'go')
% plot(theta_save(end,2),theta_save(end,5),'co')
% subplot(4,4,11), plot(theta_save(:,3),theta_save(:,4))
% hold on
% plot(theta_save(1,3),theta_save(1,4),'go')
% plot(theta_save(end,3),theta_save(end,4),'co')
% subplot(4,4,12), plot(theta_save(:,3),theta_save(:,5))
% hold on
% plot(theta_save(1,3),theta_save(1,5),'go')
% plot(theta_save(end,3),theta_save(end,5),'co')
% subplot(4,4,16), plot(theta_save(:,4),theta_save(:,5))
% hold on
% plot(theta_save(1,4),theta_save(1,5),'go')
% plot(theta_save(end,4),theta_save(end,5),'co')