function [concentration]=generate_logNdistribution(UPb_mean,UPb_std,nbin,nage)

% generate a lognormal distribution of U-Pb ages, assuming that nages ages
% were measured

mu = log((UPb_mean^2)/sqrt(UPb_std+UPb_mean^2));
sigma = sqrt(log(UPb_std/(UPb_mean^2)+1));
age_distribution = lognrnd(mu,sigma,[1,nage]);

% estimate the relative frequency and plot it.

binsize=(200-0)/(nbin-1);
edges=0:binsize:200;

counts=histcounts(age_distribution,edges);

%convert each bin into a concentration within the whole rock
%counts=h.Values;
concentration=counts./nage*2e-4;

% figure
% bar(edges(1:end-1)+binsize/2,concentration,1)
% axis square
% xlabel('age (Ma)')
% ylabel('concentration')
% title(UPb_mean)
% %axis([0 max(edges) 0 max(concentration)*1.25])
% drawnow
% pause(2)