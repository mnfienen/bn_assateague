% groundwater nets
% v1 - a mock up that we used as a whiteboard, includes PRECIP and SLR
% v2-3b, trained on numbers, but the RECHARGE was bogus
% v3c - fixed recharge
% v3d - added depth to watertable, what we ended up with on 16 August 2011
%     - modify ZWT bins to add resolution

%% paths
addpath('C:\BN_ASSATEAGUE\bayesTools')
addpath('C:\BN_ASSATEAGUE\interpTools')
addpath('C:\BN_ASSATEAGUE\valTools')

%% Bayes and filenames
netRoot = 'SLR00_EM10';
netName = [netRoot,'.neta'];
cas_dat_root = 'SEAWAT2NETICA_SLR_00'; % root for either cas file or dat file
[netstuff] = netInfo(netName);
nodeNames = cellstr(char(netstuff.name));
% if there's not already a .dat file coresponding to the .cas file, make
% one
if ~exist([cas_dat_root,'.dat'],'file')
    cas2dat(cas_dat_root);
end %if

%% make a directory for figures if there isn't one yet
figdir = [cas_dat_root,'_figures']
if ~exist(figdir,'dir')
    mkdir(figdir);
end % if

%% get data
data = load([cas_dat_root,'.dat']); % THESE NAMES MUST BE IN THE ORDER OF THE COLUMNS
%dataNames = { 'islandwidth' 'islandelev' 'maxWT' 'mean_rch' 'max_WT_loc' 'mean_DTW' 'max_DTW'}
dataNames = { 'model_row' 'islandwidth' 'islandelev' 'maxWT' 'mean_rch' 'max_WT_loc' 'west_index' 'east_index' 'max_WT_col' 'mean_DTW' 'max_DTW'}
N = length(data)

% set this to select data
did = 1:1:N; % decimate data to test things


%%  define cases
[netNameD,netNameF]=fileparts(netName)
switch netNameF
    case {netRoot}
        allTestCases = {
            'updateWT_checkWT'      {'maxWT','mean_DTW'}  {'maxWT','mean_DTW'}
            'updateAll_predWT'       {'islandwidth'   'islandelev' 'mean_rch'} {'maxWT','mean_DTW'}
            }
    %add cases if we need to treat different nets (differnt variables) differently    
 end

%% include data errors?
USEERROR=0;
if(USEERROR)
%     % suggest typical rms error (or about a bin-width?)
%     widthError = 5; % m
%     velocityError = 0.0001; % m/s
%     dischargeError = 1; %m^2
%     depthError = 0.001; % m
%     slopeError = 1e-3;
%     dataError = 0*data;
%     % scalor input for each variable will be interpreted as gaussian error
%     % structure input can be used to give specific error distributoins (e.g., log-normal, or ouput from previous simulations,...)
%     dataError(:,1) = widthError;
%     dataError(:,2) = dischargeError;
%     dataError(:,3) = velocityError;
%     dataError(:,4) = depthError;
%     dataError(:,5) = slopeError;
else
    clear dataError
end

%% additional sensitivity studies
% if(0)
%     % list cases, REMOVE one variable
%     inputvars =  {'slope','sinuosity','rcurve','isBraided',   'hydraulic_width', 'hydraulic_velocity' 'hydraulic_discharge'}
%     cnt =length(allTestCases(:,1))
%     for i=1:length(inputvars)
%         cnt = cnt+1;
%         allTestCases{cnt,1} = sprintf('remove_%s', inputvars{i});
%         j=[1:(i-1),(i+1):length(inputvars)];
%         allTestCases{cnt,2} = {inputvars{j}};
%         allTestCases{cnt,3} = {'hydraulic_depth'};
%     end
%     
%     % list cases, USE only one variable
%     for i=1:length(inputvars)
%         cnt = cnt+1;
%         allTestCases{cnt,1} = sprintf('only_%s', inputvars{i});
%         allTestCases{cnt,2} = {inputvars{i}};
%         allTestCases{cnt,3} = {'hydraulic_depth'};
%     end
%     
%     % list cases, pick any 2 variables
%     for i=1:length(inputvars)
%         for j=(i+1):length(inputvars)
%             cnt = cnt+1;
%             allTestCases{cnt,1} = sprintf('2vars_%s+%s', inputvars{i},inputvars{j});
%             allTestCases{cnt,2} = {inputvars{i},inputvars{j}};
%             allTestCases{cnt,3} = {'hydraulic_depth'};
%         end
%     end
% end
% 
%% loop throght the cases
cnt = 0;
savedStats = struct; % store results
for i=1:1:length(allTestCases(:,1))
    % get run parameters
    testcasename = allTestCases{i,1} % string
    nodeNamesIn = allTestCases{i,2} % cell
    nodeNamesOut = allTestCases{i,3}
    
    %% run
    if(exist('dataError','var'))
        %  if data have errors:
        fprintf('using data errors\n');
        testcasename=[testcasename,'_useErr'];
        [pred, dataIn,dataOut, bh] = runPredictBayes(netName, nodeNamesIn, nodeNamesOut, data(did,:), dataNames, testcasename, dataError(did,:));
    else
        [pred, dataIn,dataOut, bh] = runPredictBayes(netName, nodeNamesIn, nodeNamesOut, data(did,:), dataNames, testcasename);
    end
    close(bh)
    
    %% evaluate
    VARNAMES = nodeNamesOut
    for j=1:length(VARNAMES)
        cnt = cnt +1;
        VARNAME = VARNAMES{j}
        varid = strmatch(VARNAME, dataNames)
        parid = strmatch(VARNAME,nodeNamesOut)
        pdf = getfield(pred.pdf,VARNAME)+eps;
        pdfRanges = getfield(pred.ranges,VARNAME);
        N = length(pdf(:,1));
        [stats] = PDF2Stats(pdfRanges, pdf, 0.5);
        py975 = getPy(pdfRanges, pdf, .975);
        py95 = getPy(pdfRanges, pdf, .95);
        py75 = getPy(pdfRanges, pdf, .75);
        py25 = getPy(pdfRanges, pdf, .25);
        py05 = getPy(pdfRanges, pdf, .05);
        py025 = getPy(pdfRanges, pdf, .025);
        LR = sum(pred.logLikelihoodRatio(:,parid))
        
        % regression
        [bML,brmse1,skML] = regr_xzw([ones(N,1) stats.bMostProb(1:1:N) ],data(did,varid)-nanmean(data(did,varid)), 1./(stats.bStd(1:1:N)+eps))
        [bM,brmse1,skM] = regr_xzw([ones(N,1) stats.bMean(1:1:N) ],data(did,varid)-nanmean(data(did,varid)), 1./(stats.bStd(1:1:N)+eps))
        rmseM = sqrt(nanmean((stats.bMean(1:1:N)-data(did,varid)).^2))
        meaneM = nanmean((stats.bMean(1:1:N)-data(did,varid)))
        rmseML = sqrt(nanmean((stats.bMostProb(1:1:N)-data(did,varid)).^2))
        meaneML = nanmean((stats.bMostProb(1:1:N)-data(did,varid)))
        
        % save stats
        savedStats(cnt).testcasename = testcasename;
        savedStats(cnt).VARNAME = VARNAME;
        savedStats(cnt).LR = LR;
        savedStats(cnt).skM = skM;
        savedStats(cnt).meaneM = meaneM;
        savedStats(cnt).rmseM = rmseM;
        savedStats(cnt).bM = bM;
        savedStats(cnt).skML = skML;
        savedStats(cnt).meaneML = meaneML;
        savedStats(cnt).rmseML = rmseML;
        savedStats(cnt).bML = bML;
        
        % plot
        figure
        tp = did(:); % these went into analysis
        tid = 1:N;
        subplot 211
        hp=patch([tp;flipud(tp)],[py025;flipud(py975)],[0,0,0]+0.75); set(hp, 'linestyle','none')
        hold on
        hp=patch([tp;flipud(tp)],[py05;flipud(py95)],[0,0,0]+0.5); set(hp, 'linestyle','none')
        hp=patch([tp;flipud(tp)],[py25;flipud(py75)],[0,0,0]+.25); set(hp, 'linestyle','none')
        % add the predicted data
        plot(tp,data(did,varid),'+')
        grid on
        [fp,netNamePrint] = fileparts(netName)
        title([VARNAME,' - ',testcasename,'-', netNamePrint],'interpreter','none')
        % put up other stuff
        plot(tp, stats.bMean(1:N), 'r')
        
        
        subplot 212
        plot(stats.bMean(1:1:N),data(did,varid),'.');
        axis equal
        axis square
        hold on
        plot(xlim,xlim,'k')
        hold off
        grid on
        xlabel([VARNAME, ' predicted'], 'interpreter','none')
        ylabel([VARNAME, ' observed'], 'interpreter','none')
        title(sprintf('skill=%.2f, meane=%.2f, rmse=%.2f, LR=%.0f', skM,meaneM,rmseM,LR))
        print('-dpng', [figdir,'\prediction_',VARNAME,' - ',testcasename,'-',netNamePrint,'.png'])
        hgsave([figdir,'\prediction_',VARNAME,' - ',testcasename,'-',netNamePrint,'.fig'])
    end
    % close all
end