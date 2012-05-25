function [logMu, logSigma] = meanStd2logMeanLogStd(mu,sigma)
% [logMu, logSigma] = meanStd2logMeanLogStd(mu,sigma)
%
% convert mean and variance of data to parameters for lognormal dist
% works if data consistent with lognormal
%
% Input
%   mu, data mean (i.e., mean(data))
%   sigma, data std (i.e. std(data))
%
% Output
%   logMu,logSigma, parameters for lognormal distribution (mean and std)
logSigma = sqrt( ( log(1 + (sigma./mu).^2) ));
logMu = -0.5*(logSigma.^2) + log(mu);

return
% test

logData = 2*randn(1000,1)+5;
logSigma_true = std(logData)
logMu_true = mean(logData)

data = exp(logData);

mu_est = mean(data)
sigma_est = std(data)


[logMu_est, logSigma_est] = meanVar2logMeanLogVar(mu_est,sigma_est)

[logMu_true, logMu_est
    logSigma_true logSigma_est]

