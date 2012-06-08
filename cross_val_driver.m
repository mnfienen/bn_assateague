% CROSS VALIDATION DRIVER

function skillz = cross_val_driver(baseNetname,voodooPar,baseCASEfilename,nFolds)
% skillz = K_fold_CV_driver(baseNETname,voodooPar,baseCASEfilename,nFolds) 
% a m!ke@usgs joint
% 
% a function to perform K-fold cross validation of a Netica Bayes Net
% dependencies:
%   valTools\k_fold_maker.m
%   valTools\compile_net.m
%   valTools\ mikeppwd.txt  (contains Netica/Norsys API password as text)
%   interpTools\
%   bayesTools\cas2dat.m
%
% INPUT:
%       baseNETname --> base NETA file name (includes '.neta' extension)
%       voodooPar -- > voodoo parameter for building CPTs (float)
%       baseCASEfilename --> casefilename containing all cases used to
%             build the original net. includes '.cas' extension.
%             Subsets of this file will be created with trailing numbers
%             for each fold.
%       nFolds --> integer number of folds to evenly divide the initial
%             dataset into
%
% OUTPUT:
%       skillz --> a MATLAB structure with skill values for both training
%       sets and for validation (left out) sets as well as some metadata
%% set some paths for dependencies 
addpath('C:\BN_ASSATEAGUE\bayesTools');
addpath('C:\BN_ASSATEAGUE\interpTools');
addpath('C:\BN_ASSATEAGUE\valTools');

%% make the netica stuff global
global NETICA_ENV NETICA_NET

%% generate the set of folds from indices corresponding to the length 
%  of the number of cases
cas_dat_root = baseCASEfilename(1:end-4);
% make a .dat file out of the cas file to be able to load in the case
% information
if ~exist([cas_dat_root,'.dat'],'file')
    cas2dat(cas_dat_root);
end %if

% open the .dat-ified .cas file for data
data = load([cas_dat_root,'.dat']); % THESE NAMES MUST BE IN COLUMN ORDER
% set up the data names for all columns in the dat file
dataNames = { 'model_row' 'islandwidth' 'islandelev' 'maxWT' 'mean_rch' ...
    'max_WT_loc' 'west_index' 'east_index' ...
    'max_WT_col' 'mean_DTW' 'max_DTW'};

% set up the test cases
% MNF --> need to document the setup of this cell array
allTestCases = {
            'updateWT_checkWT'      ...
            {'maxWT','mean_DTW'}  {'maxWT','mean_DTW'}
            'updateAll_predWT'       ...
            {'islandwidth'   'islandelev' 'mean_rch'} {'maxWT','mean_DTW'}
            };
N = length(data);

% make the k_folds
allfolds = k_fold_maker(N,nFolds);

%% loop over the folds to calculate skillz
for cfold = 1:nFolds
    disp(['Fitting to data for fold ', num2str(cfold)])
    % first subset the calibration data
    data_cal = data(allfolds(cfold).retained);
    cfold_cas_calib_file = [cas_dat_root, '_', num2str(cfold), '.cas'];
    ofp = fopen(cfold_cas_calib_file,'w');
    outfilename = [baseNetname(1:end-5), '_', num2str(cfold), '.neta'];
    % Write out the cas file
    for i = 1:length(dataNames)
        fprintf(ofp,'%s ',dataNames{i});
    end
    fprintf(ofp,'\n');
    for i = 1:length(allfolds(cfold).retained)
        fprintf(ofp,'%f ',data(allfolds(cfold).retained(i),:));
        fprintf(ofp,'\n');
    end
    fclose(ofp)
    
    % now compile and relearn the net with the current fold's casefile
    compile_net(baseNetname,cfold_cas_calib_file,voodooPar,outfilename)
       
    
    
    
    
end %"for cfold" loop

