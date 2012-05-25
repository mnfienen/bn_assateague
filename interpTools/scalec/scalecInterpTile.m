function [ZI, MSEI, NMSEI, MSRI] = scalecInterpTile(x, z, s, xi, lx, filtername, nmseitol,WB);% [ZI, MSEI, NMSEI, MSRI] = scalecInterpTile(x, z, s, xi, lx, filtername, nmseitol,WB);%% remove a plane-trend and pass to scalecInterpTilePerturbations (which does not remove any trend)% % Input%   x, the nxm location of the data- repeated indices not handled well (ignored)%   z, the observations%   s, the the observation errors (i.e., standard deviations, rms, etc.)%      s is used to weight the observations as 1/s%      choose s=0 if no information is available about the weights %   xi, the interpolation locations%   lx, the interpolation weighting length scales%   filtername, the name of a filter to analyze:%      'quadloess'%      'linloess'%      'hanning'%      'boxcar'%   nmseitol, a maximum error level, if exceeded causes doubling of smoothing scales%       NOTE: if nmseitol=1 then this means we accept result with input scales%             if nmseitol<1 then this means interpolation will successive doubling of scales to reach desired noise reduction%   WB, a flag to use the waitbar to show progress. WB=1 will show waitbar, missing,empty, or other value won't% % Output%   zi, the estimate%   msei, the mean square interpolation error estimate (units of z)%   nmsei, the normalized mean square error%   msri, the mean square residuals% modifications by updating the trend only where valid perturbations were interpolated% catch  input[Ni,mi] = size(xi);[N,m] = size(x);if(~exist('WB'))    WB = 0;elseif(isempty(WB))    WB = 0;end% need some consistent weightswtol = 0.01;[wt, var_z] = consistentWeight(z, s.^2, wtol);% do regression to remove a norm fieldbtrend = zeros(mi+1,1);bi = btrend;% only compute against variables with variancestd_x = std(x);id =find(std_x==0);if(length(id)>0) % catch variables with zero variance (e.g., a profile)    std_x(id) = 1+0*id;endvarid = find(std_x>0);if(length(varid>0))    % this is just for getting the data ready, so it is meant to be bullet proof, not statistically pure!    [btrend([1,varid+1]),bi([1,varid+1])] = regr_xzw([ones(N,1),x(:,varid)], z, wt);    fprintf('removed order %d polynomial\n', 1);end% regression failed if nanif(any(isnan(btrend)))    % pad with zero    btrend = zeros(length(btrend),1);end% remove trendztrend = [ones(N,1),x]*btrend;z = z - ztrend; % send to interp tile perturbations[ZI, MSEI, NMSEI, MSRI] = scalecInterpTilePerturbations(x, z, s, xi, lx, filtername,nmseitol,WB);% replace the trend - except where uselessztrend = [ones(Ni,1),xi]*btrend;ztrend = reshape(ztrend, size(ZI));idgood = NMSEI~=1 & NMSEI~=0;ZI(idgood) = ZI(idgood) + ztrend(idgood);