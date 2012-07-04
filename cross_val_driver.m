% CROSS VALIDATION DRIVER

function cross_val_driver(baseNetname,voodooPar,baseCASEfilename,nFolds)
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
%       savedStats --> a MATLAB structure with skill values for both training
%       sets and for validation (left out) sets as well as some metadata.
%       Saved to .mat files
%% make the netica stuff global --> kinda kludgy but required...
global NETICA_ENV NETICA_NET NETICA_NETNAME NETICA_NODESTRUCT NETICA_STREAMER NETICA_ALL_NODES

%% set some paths for dependencies 
addpath('C:\BN_ASSATEAGUE\bayesTools');
addpath('C:\BN_ASSATEAGUE\interpTools');
addpath('C:\BN_ASSATEAGUE\valTools');


%% generate the set of folds from indices corresponding to the length 
%  of the number of cases
cas_dat_root = baseCASEfilename(1:end-4);
% make a directory for figures if there isn't one yet
figdir = [cas_dat_root,'_figures'];
if ~exist(figdir,'dir')
    mkdir(figdir);
end % if

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
% set this to select data
did = 1:1:N; % decimate data to test things --> if 1:1:N; all is analyzed at once
% make the k_folds
allfolds = k_fold_maker(N,nFolds);

%% initalize the output structure
savedStats = struct; % store results

%% calculate skill, etc. for the entire data set priot to cross-validation
    for ctcase = 1:length(allTestCases(:,1))
    % get run parameters
        testcasename = allTestCases{ctcase,1}; % string
        NodeNamesIn = allTestCases{ctcase,2}; % cell
        NodeNamesOut = allTestCases{ctcase,3};
        [ALLpred,ALLDataIn, ALLDataOut, ALLbh] = runPredictBayesMNF(baseNetname, ...
            NodeNamesIn, NodeNamesOut, data, dataNames, ['ALL', testcasename]);
        close (ALLbh);
        
        savedStats.casname = cas_dat_root;
        savedStats.netaname = baseNetname;
        
        savedStats.ALLsavedStats =  PredictBayesPostProc(baseNetname, NodeNamesOut, ALLpred, dataNames, ...
                data,  figdir,['ALL', testcasename]);
        currfilename = [cas_dat_root, '_ct-', num2str(ctcase), '_NO_CV.mat'];
        save(currfilename,'savedStats');
        savedStats = struct; % store results
    end %ctcase
    
%% loop over the folds to calculate skillz
for cfold = 1:nFolds
       close('all');
%% First, fit to the left out data
    disp(['Fitting to data for fold ', num2str(cfold)])
    % first subset the calibration data
    data_cal = data(allfolds(cfold).retained,:);
    data_val = data(allfolds(cfold).leftout,:);
    cfold_cas_calib_file = [cas_dat_root, '_', num2str(cfold), '.cas'];
    ofp = fopen(cfold_cas_calib_file,'w');
    CurrFoldNet = [baseNetname(1:end-5), '_', num2str(cfold), '.neta'];

    % Write out the cas file
    for i = 1:length(dataNames)
        fprintf(ofp,'%s ',dataNames{i});
    end
    fprintf(ofp,'\n');
    for i = 1:length(data_cal)
        fprintf(ofp,'%f ',data_cal(i,:));
        fprintf(ofp,'\n');
    end
    fclose(ofp);
    
    % now compile and relearn the net with the current fold's casefile
    compile_net(baseNetname,cfold_cas_calib_file,voodooPar,CurrFoldNet)

    
%% Next, calculate skills for both training and validation data
        savedStats.casname = cfold_cas_calib_file;
        savedStats.netaname = CurrFoldNet;

for ctcase = 1:length(allTestCases(:,1))
    % get run parameters
        testcasename = allTestCases{ctcase,1}; % string
        NodeNamesIn = allTestCases{ctcase,2}; % cell
        NodeNamesOut = allTestCases{ctcase,3};
    % run the datasets through the Bayes prediction code to get residuals
    % and calculate/plot skills for validation and calibration
        % validation set
        [valpred,valDataIn, valDataOut, valbh] = runPredictBayesMNF(CurrFoldNet, ...
            NodeNamesIn, NodeNamesOut, data_val, dataNames, ['val', testcasename]);
        close (valbh);

        savedStats.valsavedStats =  PredictBayesPostProc(CurrFoldNet, NodeNamesOut, valpred, dataNames, ...
                data_val,  figdir,['val', testcasename]);
        
        
        % calibration set
        [calpred,calDataIn, calDataOut, calbh] = runPredictBayesMNF(CurrFoldNet, ...
            NodeNamesIn, NodeNamesOut, data_cal, dataNames, ['cal', testcasename]);
        close (calbh);

        savedStats.calsavedStats =  PredictBayesPostProc(CurrFoldNet, NodeNamesOut, calpred, dataNames, ...
                data_cal,  figdir,['val', testcasename]);
            
        % dump out the saved stats to disk
        currfilename = [cfold_cas_calib_file(1:end-4), '_fold-', num2str(cfold) '_ct-', num2str(ctcase), '.mat'];
        save(currfilename,'savedStats');
        savedStats = struct;
    end %for ctcase
    
    
end %"for cfold" loop
i=1
