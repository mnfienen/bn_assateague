function [Su_sfc,Su_sfc_logstd, suLogMean,suLogStd, pdfRange,pdf] = findSuModel(z,su, Su_slope, max_depth_into_bed);
% [Su_sfc,Su_sfc_logstd, suLogMean,suLogStd, pdfRange,pdf] = findSuModel(z,su, Su_slope, max_depth_into_bed);
%
% compute shear stress stats consistent with MBES model
%
% Input
%   z, depth into bed
%   Su, shear strength estimates
%   Su_slope, the ASSUMED rate of increase in shear strength with depth: Su(z) = Su_sfc + Su_slope*z
%   max_depth_into_bed, optional truncation of the model at this depth, default is 1m
%                           
% Output
%    Su_sfc, best estimate of parameter in the ASSUMED model
%    Su_sfc_logstd, estimated log-transformed standard deviation
%      use it this way: pdf = lognpdf(pdfRange,log(Su_sfc),Su_sfc_logstd);
%    suLogMean,suLogStd,  statistics of the log transformed data
%    pdfRanges,pdf, log-normal distribution for Su_sfc
%

% first, get some stats
logSu = log(su);
suLogMean = mean(su);
suLogStd = std(su);

% use only these
MAXDEPTHINTOBED = 1;
if(~exist('max_depth_into_bed'))
    max_depth_into_bed=MAXDEPTHINTOBED;
end
id = find(z<max_depth_into_bed);

%  Su_model.m spits out log(su)
% set some upper and lower bounds
LB =  0;
UB = 30; % kPa
[Su_sfc,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN]=lsqcurvefit('Su_model',su(1),z(id),logSu(id),[LB],[UB],[],Su_slope);
if(isnan(Su_sfc))
    Su_sfc = su(1);
end
% check
Su_est = exp(Su_model(Su_sfc, z, Su_slope));

% assume residuals uniformly distributed in z, (thus apply to z0)
Su_sfc_logstd = std(log(Su_est)-log(su));
pdfRange = 0.1+[0:100]*30/100;
pdf = lognpdf(pdfRange,log(Su_sfc),Su_sfc_logstd);

return
% here is a test
data = [
0.05		2.82
0.15		4.45
0.25		6.09
0.35		7.72
0.5		10.17
0.7		13.44
0.9		16.71
1.25		19.00
1.75		19.00
3.5		19.00];
z = data(:,1);
su = data(:,2);

% some slopes
Su_slope =1.7
Su_slope = 12
Su_slope = 0

% use this function
[Su_sfc,Su_sfc_logstd, suLogMean,suLogStd, pdfRange,pdf] = findSuModel(z,su, Su_slope, 1);
Su_est = exp(Su_model(Su_sfc, z, Su_slope));

% assume residuals uniformly distributed in z, (thus apply to z0)
plot(z,su,'.', z,Su_est, pdf,pdfRange)
plot(z,log(su),'.', z,log(Su_est), pdf,log(pdfRange))
