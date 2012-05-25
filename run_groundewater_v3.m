% groundwater nets
% v1 - a mock up that we used as a whiteboard, includes PRECIP and SLR
% v2-3b, trained on numbers, but the RECHARGE was bogus
% v3c - fixed recharge
% v3d - added depth to watertable, what we ended up with on 16 August 2011
%     - modify ZWT bins to add resolution

%% paths
addpath('C:\BN_ASSATEAGUE\bayesTools')
addpath('C:\BN_ASSATEAGUE\interpTools')
%% Bayes
%netName = 'groundwater_v3d.neta'
netName = 'groundwater_v3dEM.neta' % EM train
[netstuff] = netInfo(netName)
nodeNames = cellstr(char(netstuff.name))


%% get data
data = load('SEAWAT2NETICA+dwt.dat'); % RECHARGE fixed, added depth to water table
dataNames = {'modelrow'   'islandwidth'   'islandelev'   'zwatertable'      'recharge'       'xwatertable'       'westindex'   'eastindex'  'maxWTCol' 'dwatertable'  'maxDTW'
}
N = length(data)

% set this to select data
did = 1:1:N; % decimate data to test things


%%  define cases
[netNameD,netNameF]=fileparts(netName)
switch netNameF
    case {'groundwater_v3d','groundwater_v3dEM'}
        allTestCases = {
            'updateZWT_checkZWT'       {'zwatertable','xwatertable','dwatertable'}  {'zwatertable','xwatertable','dwatertable'}
            'updateAll_predZWT'       {'islandwidth'   'islandelev' 'recharge'} {'zwatertable','xwatertable','dwatertable'}
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
    if(exist('dataError'))
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
        print('-dpng', ['prediction_',VARNAME,' - ',testcasename,'-',netNamePrint,'.png'])
        hgsave(['prediction_',VARNAME,' - ',testcasename,'-',netNamePrint,'.fig'])
    end
    % close all
end

return


%% format the stat output
dostats = {'LR', 'skM', 'bM', 'meaneM', 'rmseM'}
dotests =  unique(cellstr(char(savedStats.testcasename)))
dovars = unique(cellstr(char(savedStats.VARNAME)))
for i=1:length(dostats)
    fprintf('%s\n', dostats{i});
    for j=1:length(dotests)
        fprintf('%s,', dotests{j});
        id1 = strmatch(dotests{j},char(savedStats.testcasename));
        for k=1:length(dovars)
            id2 = strmatch(dovars{k},char(savedStats(id1).VARNAME));
            fprintf('%.3f, ', savedStats(id1(id2(end))).(dostats{i}));
        end
        fprintf('\n')
    end
end


%% show priors
sid = [2:9] % order for plots
dataNames(sid)
sum(sum(~isnan(data(:,sid))))

figure
for i=1:length(sid)
    subplot(length(sid),1,i)
    pdf = getfield(pred.pdf,['prior_',dataNames{sid(i)}])+eps;
    pdfRanges = getfield(pred.ranges,dataNames{sid(i)});
    plotPDF2(pdfRanges,pdf,1)
    xlabel(dataNames{sid(i)}, 'interpreter','none')
end
print('-dpng', ['priors','- ',netName,'.png'])
hgsave(['priors','- ',netName,'.fig'])


%% TEST CRITICAL DEPTH PREDICTIONS
% load a data set
load riverES_USGSstream_v6_newindx28jan2011_LR=3937_only_hydraulic_width.734669.mat
load riverES_USGSstream_v6_newindx28jan2011_LR=4924_2vars_slope+hydraulic_velocity.734669.mat
% get prediction
pdf = getfield(pred.pdf,'hydraulic_depth')+eps;
pdfRanges = getfield(pred.ranges,'hydraulic_depth');
dataOut = dataOut(:,1);
likelihoodRatio = sum(pred.logLikelihoodRatio(:,1));
ti = [1:N]';

% get stats
[stats] = PDF2Stats(pdfRanges, pdf, 0.5);
py975 = getPy(pdfRanges, pdf, .975);
py95 = getPy(pdfRanges, pdf, .95);
py75 = getPy(pdfRanges, pdf, .75);
py25 = getPy(pdfRanges, pdf, .25);
py05 = getPy(pdfRanges, pdf, .05);
py025 = getPy(pdfRanges, pdf, .025);
[bML,brmse1,skML] = regr_xzw([ones(N,1) stats.bMostProb(1:1:N) ],dataOut-nanmean(dataOut(1:N)), 1./stats.bStd(1:1:N));
meaneML = nanmean((stats.bMostProb(1:1:N)-dataOut(1:N)));
rmseML = sqrt(nanmean((stats.bMostProb(1:1:N)-dataOut).^2));
[bM,brmse1,skM] = regr_xzw([ones(N,1) (stats.bMean(1:1:N)- nanmean(dataOut))],dataOut-nanmean(dataOut(1:N)), 1./stats.bStd(1:1:N));
meaneM = nanmean((stats.bMean(1:1:N)-dataOut(1:N)));
rmseM = sqrt(nanmean((stats.bMean(1:1:N)-dataOut).^2));
Bskill=1-nansum(nansum((pred.pdf.data_hydraulic_depth-pdf).^2))/N;
Vskill=1-nanvar(stats.bMean(1:1:N)-dataOut)/nanvar(dataOut) ;
fprintf('meane=%.2f, rmse=%.2f, skM=%.2f, skB=%.2f, skV=%.2f, LR=%.0f -- %s\n', meaneM,rmseM,skM,Bskill,Vskill,sum(likelihoodRatio),testcasename)

% plot regression
figure
plot(stats.bMean(1:1:N),dataOut,'.b', stats.bMostProb(1:1:N),dataOut,'.r', [min(pdfRanges) max(pdfRanges)],[min(pdfRanges) max(pdfRanges)])
xlabel('predicted')
ylabel('observed')
legend('mean', 'most prob.')
title('water depth (m)')

% plot confidence
figure
hp=patch([ti;flipud(ti)],[py025;flipud(py975)],[0,0,0]+0.75); set(hp, 'linestyle','none')
hold on
hp=patch([ti;flipud(ti)],[py05;flipud(py95)],[0,0,0]+0.5); set(hp, 'linestyle','none')
hp=patch([ti;flipud(ti)],[py25;flipud(py75)],[0,0,0]+.25); set(hp, 'linestyle','none')
plot(ti,dataOut(1:N),'b.')
plot(ti(1:1:N), stats.bMean(1:1:N), 'y+', 'linewidth', 2)
plot(ti(1:1:N), stats.bMostProb(1:1:N), 'r.', 'linewidth', 2)
xlim([ti(1) ti(end)])
axis tight
ylim([min(pdfRanges) max(pdfRanges)])
xlabel('location id')
ylabel('depth (m)')
grid on
legend('p95','p90','p50', 'obs.', 'mean','most likely')
title(sprintf('%s, meane=%.2f, rmse=%.2f, sk=%.2f, LR=%.0f', testcasename,meaneM,rmseM,skM,sum(likelihoodRatio)),'interpreter','none')
figure(gcf)

% plot prob exceed
dthresh = 2.0;
threshid = find(pdfRanges<dthresh)
pthresh = sum(pdf(:,threshid),2);
id = find(dataOut<dthresh);
dbin = 0.2;
pbin = 0:dbin:1
pctbin = pbin*nan;
for i=1:length(pbin)
    idAll = find(abs(pthresh-pbin(i))/dbin<0.5);
    idCrit = find(abs(pthresh-pbin(i))/dbin<0.5 & dataOut<dthresh);
    pctbin(i) = length(idCrit)/(length(idAll)+eps);
end
figure
plot(pthresh,dataOut,'.', pthresh(id),dataOut(id),'r.', [0 , 1],[min(pdfRanges),max(pdfRanges)],'k--')
hold on
xlim([0 1])
hp=plotyy(nan,nan, pbin,pctbin)
ylabel('depth (m)')
xlabel('probability depth < threshold')
title(sprintf('Dcrit = %.1f m', dthresh))


% % % improved figure
% % % red is bad (mission failure requires depth>dcrit)
% % figure
% % [n0,x] = hist(dataOut, 0:1:20)
% % bar(x,n0/sum(n0),'facecolor','b')
% % hold on
% % id = find(pthresh>0.90);
% % [n]= hist(dataOut(id),x)
% % bar(x,n/sum(n0), 'facecolor','r')
% % id = find(pthresh<0.10);
% % [n]= hist(dataOut(id),x)
% % bar(x,n/sum(n0), 'facecolor','g')
% % xlabel('depth (m)')
% % xlim([-.5 10])
% % title(sprintf('Dcrit = %.1f m', dthresh))
% % ylabel('pct occurrence')



% improved figure
% red is bad (mission failure requires depth>dcrit)
figure
risk_level = 0.25; % level of failure permitted
[n0,x] = hist(dataOut, 0:1:20)
%bar(x,n0/sum(n0),'facecolor','b')
%hold on
idshallow = find(pthresh>(1-risk_level)); % % certain that depth<Dcrit
[nshallow]= hist(dataOut(idshallow),x)
bar(x,nshallow./n0, 'facecolor','r'); % fraction of test set in each class is 90% likely to be shallow.  Classify as failure (red)
xlabel('observed depth (m)')
ylabel('fraction of observations')
legend(sprintf('p(d<d_{crit})>%.0f%%',(1-risk_level)*100))
hold on
iddeep = find(pthresh<risk_level);
[ndeep]= hist(dataOut(iddeep),x)
hb= bar(x,ndeep./n0, 'facecolor','g');% fraction of test set in each class is 10% likely to be shallow.  Classify as success (green)
set(get(hb,'children'), 'FaceAlpha',1)
xlabel('depth (m)')
xlim([-.5 10])
legend(sprintf('p(d<d_{crit})>%.0f%%',(1-risk_level)*100), sprintf('p(d<d_{crit})<%.0f%%',(risk_level)*100))
title(sprintf('%s: Dcrit = %.1f m, %.0f%% risk level', testcasename, dthresh, risk_level*100), 'interpreter','none')
%ylabel('pct occurrence')
ylim([0 1])

% quantify
pctShallowTrue = sum(dataOut(idshallow)<dthresh)/sum(dataOut<dthresh) % predict shallow, is shallow at risk level
pctShallowFalse = sum(dataOut(iddeep)<dthresh)/sum(dataOut<dthresh) % predict deep, is shallow
pctTrueDeep = sum(dataOut(iddeep)>dthresh)/sum(dataOut>dthresh) % predict deep, is deep
pctFalseDeep = sum(dataOut(idshallow)>dthresh)/sum(dataOut>dthresh) % predict shallow, is deep

%% simplify
figure
risk_level = 0.25; % level of failure permitted
[n0,x] = hist(dataOut, 0:1:20)
bar(x,n0/sum(n0),'facecolor','b')
hold on
idshallow = find(pthresh>(1-risk_level)); % % certain that depth<Dcrit
[nshallow]= hist(dataOut(idshallow),x)
iddeep = find(pthresh<risk_level);
[ndeep]= hist(dataOut(iddeep),x)
hb= bar(x,(ndeep+nshallow)/sum(n0), 'facecolor','g');% fraction of test set in each class is 10% likely to be shallow.  Classify as success (green)
bar(x,nshallow/sum(n0), 'facecolor','r'); % fraction of test set in each class is 90% likely to be shallow.  Classify as failure (red)
xlabel('observed depth (m)')
ylabel('fraction of observations in each depth')
xlim([-.5 10.5])
legend('amgiguous', sprintf('predict deep: p(d<d_{crit})<%.0f%%',(risk_level)*100), sprintf('predict shallow: p(d<d_{crit})>%.0f%%',(1-risk_level)*100))
title(sprintf('%s: Dcrit = %.1f m, %.0f%% risk level', testcasename, dthresh, risk_level*100), 'interpreter','none')
%ylabel('pct occurrence')
%ylim([0 1])

%% just do above/below crit
% plot prob exceed
dthresh = 2.0;
threshid = find(pdfRanges<dthresh)
pthresh = sum(pdf(:,threshid),2);


figure
risk_level = 0.25; % level of failure permitted
n0=nan(2,1);
n0(1) = sum(dataOut<dthresh);
n0(2) = sum(dataOut>=dthresh);
x = dthresh+[-1,1]
bar(x,n0/sum(n0),'facecolor','b')
hold on
idshallow = find(pthresh>(1-risk_level)); % % certain that depth<Dcrit
nshallow = nan*n0;
nshallow(1) = sum(dataOut(idshallow)<dthresh);
nshallow(2) = sum(dataOut(idshallow)>=dthresh);
iddeep = find(pthresh<risk_level);
ndeep = nan*n0;
ndeep(1) = sum(dataOut(iddeep)<dthresh);
ndeep(2) = sum(dataOut(iddeep)>=dthresh);
hb= bar(x,(ndeep+nshallow)/sum(n0), 'facecolor','g');% fraction of test set in each class is 10% likely to be shallow.  Classify as success (green)
bar(x,nshallow/sum(n0), 'facecolor','r'); % fraction of test set in each class is 90% likely to be shallow.  Classify as failure (red)
xlabel('observed depth (m)')
ylabel('fraction of observations in each depth')
xlim([dthresh+[-1,1]])
legend('amgiguous', sprintf('predict deep: p(d<d_{crit})<%.0f%%',(risk_level)*100), sprintf('predict shallow: p(d<d_{crit})>%.0f%%',(1-risk_level)*100))
title(sprintf('%s: Dcrit = %.1f m, %.0f%% risk level', testcasename, dthresh, risk_level*100), 'interpreter','none')
%ylabel('pct occurrence')
%ylim([0 1])

