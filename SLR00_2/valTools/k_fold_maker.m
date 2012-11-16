% k_fold index maker
% a m!ke@usgs joint
% mnfienen@usgs.gov
% function allfolds = k_fold_maker(n,k)
% input:
%   n is the length of the sequence of indices
%   k is the number of folds to split it into
%
% returns a struct called allfolds with k elements
% allfolds(i).leftout is the left out indices of fold i
% allfolds(i).retained is the kept indices of fold i

function allfolds = k_fold_maker(n,k)

allinds = [1:n];
inds_per_fold = floor(n/k);
dingleberry = rem(n,k);

for i = 1:k-1
    currn = length(allinds);
    currinds = randsample(currn,inds_per_fold);
    allfolds(i).leftout = allinds(currinds);
    allfolds(i).retained = setdiff([1:n], allinds(currinds));
    % reduce the size of allinds by the currently leftout 
    allinds = setdiff(allinds,allfolds(i).leftout);
end % for

allfolds(i+1).leftout = allinds;
allfolds(i+1).retained = setdiff([1:n],allinds);