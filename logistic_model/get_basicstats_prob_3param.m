function [pspike,corrmat,y] = get_basicstats_prob_3param(movie,dt)

% Calculates the mean firing rate, s.d. in firing rates and mean pairwise correlation in data. 

%%%%%%%%
% PARAMETERS
%%%%%%%%

nsamp = 1e3; % Number of samples per timestep to test n_on

%%%%%%%%
% INITIAL PROCESSING
%%%%%%%%

% Get number of cells and timesteps
[nt,nunits] = size(movie); 

% Rebin time
tindsnew = [1:dt:nt];
ntnew = length(tindsnew);
newmovie = zeros(ntnew,nunits);
for t = 1:ntnew-1;
    nspikesmovieseg = 0.256*sum(movie(tindsnew(t):tindsnew(t+1)-1,:)); % Expected number of spikes in those time bins
    newmovie(t,:) = 1-poisspdf(zeros(1,nunits),nspikesmovieseg)';
end
nspikesmovieseg = 0.256*sum(movie(tindsnew(end):nt,:));
newmovie(end,:) = 1-poisspdf(zeros(1,nunits),nspikesmovieseg)';

%%%%%%%%%%%%%
% MAKE PSEUDO MOVIE
%%%%%%%%%%%%%
ntpseudo = nsamp*ntnew;
pseudomovie = zeros(ntpseudo,nunits);
for t = 1:ntnew
    pseudomovie(1+(t-1)*nsamp:nsamp+(t-1)*nsamp,:) = floor(rand(nsamp,nunits)+repmat(newmovie(t,:),nsamp,1));
end

%%%%%%%%%
% FIRING RATES
%%%%%%%%%
pspike = mean(newmovie,1);

%%%%%%%%
% Correlations
%%%%%%%%
corrmat = corr(pseudomovie);
corrlist = corrmat;
corrlist(1:nunits+1:end)=[];

%********
% WRAP SUMMARY STATS
%********
y(1) = mean(pspike);
y(2) = std(pspike);
y(3) = mean(corrlist);