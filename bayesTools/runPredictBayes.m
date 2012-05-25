function [pred, dataIn, dataOut, bh]  = runPredictBayes(netName, nodeNamesIn, nodeNamesOut, data, dataNames, testcasename, dataErr)
% runPredictBayes(netName, nodeNamesIn, nodeNamesOut, data, dataNames)
% fix up inputs and drive them through the processor

% optional for test case
if(~exist('testcasename'))
    testcasename = '';
end

% get data
dataIn = [];
if(exist('dataErr'))
    % got to extract errors too
    dataInErr = [];
end
for j=1:length(nodeNamesIn)
    did = strmatch(nodeNamesIn{j},dataNames, 'exact');
    if(length(did)==1)
        dataIn(:,j) = data(:,did);
        if(exist('dataErr'))
            % package the errors
            dataInErr(:,j) = dataErr(:,did);
        end
    else
        error(sprintf('no data or multiple data matching %s', nodeNamesIn{j}))
    end
end
dataOut = [];
for j=1:length(nodeNamesOut)
    did = strmatch(nodeNamesOut{j},dataNames, 'exact');
    if(length(did)==1)
        dataOut(:,j) = data(:,did);
    else
        warning(sprintf('no data or multiple data matching %s', nodeNamesOut{j}))
        dataOut(:,j) = 0*data(:,1);
    end
end
[N,m] = size(dataIn);

% if(~exist('dataInErr') & ~isstruct(dataOut))
%     dataInErr = repmat(eps, N,m);
% end


%% make predictions in Netica
% pred = predictBayes(netName, nodeNamesIn, nodeNamesOut, dataIn, dataOut)
% using dataErr causes probabilistic input
if(exist('dataErr'))
    fprintf('predicting using data errors->pdfs\n')
    pred = predictBayes_v2dev(netName, nodeNamesIn, nodeNamesOut, dataIn, dataOut, dataInErr)
    pred.dataInErr = dataInErr;
else
    fprintf('predicting using data values\n')
    pred = predictBayes_v2dev(netName, nodeNamesIn, nodeNamesOut, dataIn, dataOut)
    pred.dataInErr = [];
end


% extract update and  prior
[fd,netRoot]=fileparts(netName);
for j = 1:length(nodeNamesOut)
    bh(j) = figure
    pdfRanges = getfield(pred.ranges, nodeNamesOut{j});
    pdf_prior = getfield(pred.pdf, ['prior_',nodeNamesOut{j}]);
    pdfp = repmat(pdf_prior,N,1);
    pdf = getfield(pred.pdf, [nodeNamesOut{j}]);
    if(length(pdfRanges)>length(pdf_prior(1,:)))
        % continuous
        pdfRangePlot = pdfRanges; %pdfRanges(2:end)-0.5*diff(pdfRanges) % for pcolor
    else
        % discrete
        pdfRangePlot = [pdfRanges; pdfRanges(end)+1]-0.5;
    end
    pdfData = getfield(pred.pdf, ['data_',nodeNamesOut{j}]);
    
    % get likelihoods
    probModelPrior = pred.probModelPrior(:,j);
    probModelUpdate = pred.probModelUpdate(:,j);
    likelihoodRatio = pred.logLikelihoodRatio(:,j);
    
    
    %% plot
    subplot 121
    pcolor(pdfRangePlot,1:N, double(repmat([pdf_prior,pdf_prior(end)],N,1))); shading flat; caxis([0 1])
    hold on
    plot(dataOut(:,j),1:N, 'k+')
    hold off
    title([netName,': prior'], 'interpreter', 'none')
    xlim([min(pdfRanges) max(pdfRanges)])
    xlabel(nodeNamesOut{j}, 'interpreter','none')
    colorbar
    
    subplot 122
    pcolor(pdfRangePlot,1:N, double([pdf,pdf(:,end)])); shading flat; caxis([0 1])
    hold on
    plot(dataOut(:,j),1:N, 'k+')
    hold off
    title(sprintf('%s update: LR= %.2f', testcasename, nansum(likelihoodRatio)), 'interpreter','none')
    xlim([min(pdfRanges) max(pdfRanges)])
    xlabel(nodeNamesOut{j}, 'interpreter','none')
    colorbar
    
    % save
    try
        % case where permisions are fouled up
        savename = [netRoot,'[',nodeNamesOut{j},']_LR=',num2str(round(nansum(likelihoodRatio))),'_',testcasename,'.',num2str(fix(now)),'.fig'];
        hgsave(savename)
        savename = [netRoot,'[',nodeNamesOut{j},']_LR=',num2str(round(nansum(likelihoodRatio))),'_',testcasename,'.',num2str(fix(now)),'.png'];
        print('-dpng', savename)
    catch
        warning('file save %s failed', savename)
    end
    
end

% save data
try
    % case where permisions are fouled up
    savename = [netRoot,'_LR=',num2str(round(nansum(likelihoodRatio))),'_',testcasename,'.',num2str(fix(now)),'.mat']
    save(savename)
catch
    warning('file save %s failed', savename)
end

