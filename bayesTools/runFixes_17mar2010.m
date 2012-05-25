% use the net to get what we need
netName = 'SB_Netica_test3.neta'
[netstuff] = netInfo(netName)
nodeNames = cellstr(char(netstuff.name))


% use "dataNames" match file to nodes--this lines up each column
dataNames = {'FID'	'X'	'Y'	'Cliff_Slope'	'HistRate'	'Retreat_98_05'	'Height'	'Impact_hours'	'geology'}


%% optimize the bin distributions
%nodeNames = {'HistRate'	'Cliff_Slope'	'Retreat_03_06' 'Impact_hours'	'geology' 'Height'}
% {'HistRate'	'Cliff_Slope'	'Retreat_03_06' 'Impact_hours'	'geology' 'Height'}
if(0) 
minbin    = [0.1,  10,   .05,   50,   1, 10]
M = 5; % shoot for this many bins
for j=1:length(nodeNames)
    % get the variable j
    %eval(sprintf('tmpdata=%s;', nodeNames{j}));
    did = strmatch(nodeNames{j},dataNames)
    tmpdata = data(:,did);
    tmpdata = tmpdata(find(~isnan(tmpdata)));
    unique(tmpdata)
    
    N = length(find(~isnan(tmpdata)))
    zmin = floor(min(tmpdata)/minbin(j))*minbin(j)
    zmax = ceil(max(tmpdata)/minbin(j))*minbin(j)
    %dz = (zmax-zmin)/M
    %z = zmin:dz:zmax
    [sortz,sortid] = sort(tmpdata);
    z = sortz(1:round(N/M):end); z(end)= zmax
    % adjust the bins
    %dz = max([minbin(j),median(diff(z))]);
    dz = ceil(median(diff(z))/minbin(j))*minbin(j)
    z = zmin:dz:zmax; z(end)=zmax;
    if(0)
        for i=2:(M-1)
            while((z(i)-z(i-1))< minbin(j) && z(i)<zmax)
                z(i:end-1) = z((i+1):end)
            end
        end
    end
    z = z([1;1+find(diff(z(:))>0)])';
    hist(tmpdata, z)
    title(nodeNames{j})
    xlabel(num2str(z(:)'))
    figure(gcf)
    pause
end
end


%% must choose from these
%nodenames for sb: 'Impact_hours'    'geology'    'Height'    'HistRate'    'Cliff_Slope'    'Retreat_98_05'
% list cases
allTestCases = {
    'updateInput_checkRetreat'            {'HistRate'	'Cliff_Slope'	'Impact_hours' 'geology' 'Height'}, {'Retreat_98_05'};
    'updateInputRetreat_checkRetreat'     {'Retreat_98_05' }, {'Retreat_98_05'};
    'updateLTR_checkRetreat'              {'HistRate'	}, {'Retreat_98_05'} ;
    'updateSlope_checkRetreat'            {'Cliff_Slope'	}, {'Retreat_98_05'};
    'updateImpact_checkRetreat'           {'Impact_hours'},  {'Retreat_98_05'};
    'updateGeology'                       {'geology'}, {'Retreat_98_05'};
    'updateHeight'                        {'Height'}, {'Retreat_98_05'} ;
    'updateLTRSlope'                      {'HistRate'	'Cliff_Slope'	}, {'Retreat_98_05'};
    'updateLTRSlopeImpact'                {'HistRate'	'Cliff_Slope'	'Impact_hours' }, {'Retreat_98_05'};
     'updateNotGeol_checkRetreat'         {'HistRate'	'Cliff_Slope'	'Impact_hours' 'Height'}, {'Retreat_98_05'};
    }


for i=1:length(allTestCases(:,1))
    % assign inputs
    testcase = allTestCases{i,1}; % string
    nodeNamesIn = allTestCases{i,2}; % cell
    nodeNamesOut = allTestCases{i,3};
    % get data
    dataIn = [];
    for j=1:length(nodeNamesIn)
        did = strmatch(nodeNamesIn{j},dataNames);
        if(length(did)==0)
            error(sprintf('no data matching %s', nodeNamesIn{j}))
        end
        dataIn(:,j) = data(:,did);
    end
    dataOut = [];
    for j=1:length(nodeNamesOut)
        did = strmatch(nodeNamesOut{j},dataNames);
        if(length(did)==0)
            error(sprintf('no data matching %s', nodeNamesOut{j}))
        end
        dataOut(:,j) = data(:,did);
    end
    N = length(dataOut);
    pred = predictBayes(netName, nodeNamesIn, nodeNamesOut, dataIn(1:N,:), dataOut(1:N,:))

    %% compare update to prior
    pdfRanges = getfield(pred.ranges, nodeNamesOut{1});
    pdf_prior = getfield(pred.pdf, ['prior_',nodeNamesOut{1}]);
    pdf = getfield(pred.pdf, [nodeNamesOut{1}]);
    pdfRangePlot = pdfRanges(2:end)-0.5*diff(pdfRanges)

    %% get likelihoods
    probModelPrior = pred.probModel0Data;
    probModelUpdate = pred.probModel2Data;
    likelihoodRatio = pred.likelihoodRatio;
    pdfp = repmat(pdf_prior,N,1);
 

    %% plot
    figure
    subplot 121
    pcolor(pdfRanges,ti(1:N), repmat([pdf_prior,pdf_prior(end)],N,1)); shading flat; caxis([0 1])
    hold on
    plot(dataOut(1:N),ti(1:N), 'm+')
    hold off
    title([netName,': prior'], 'interpreter', 'none')
    xlim([min(pdfRanges) max(pdfRanges)])
    xlabel(nodeNamesOut, 'interpreter','none')
    colorbar

    subplot 122
    pcolor(pdfRanges,ti(1:N), [pdf,pdf(:,end)]); shading flat; caxis([0 1])
    hold on
    plot(dataOut(1:N),ti(1:N), 'm+')
    hold off
    title(sprintf('%s update: LR= %.2f', testcase, sum(likelihoodRatio)), 'interpreter','none')
    xlim([min(pdfRanges) max(pdfRanges)])
    xlabel(nodeNamesOut, 'interpreter','none')
    colorbar
    
    %% save
    [fd,netRoot]=fileparts(netName);
    savename = [netRoot,'_LR=',num2str(round(sum(likelihoodRatio))),'_',testcase,'.',num2str(fix(now)),'.mat'];
    save(savename)
    savename = [netRoot,'_',testcase,'.',num2str(fix(now)),'.fig'];
    hgsave(savename)
    savename = [netRoot,'_',testcase,'.',num2str(fix(now)),'.eps'];
    print('-depsc2', savename)
    pause(1)
end

% compare to prior
figure
R = probModelUpdate./probModelPrior;
plot(ti, probModelPrior, 'r', ti,probModelUpdate,'k', ti,R,ti,0*ti+1,'g--');
ylabel('prob.')
xlabel('location')
xlim([min(ti) max(ti)])



%%
%% show priors
priornamesout=nodeNames;
priornamesout=[nodeNamesIn, nodeNamesOut];
predprior = predictBayes(netName, nodeNamesIn(1), priornamesout, 0,nanmean([dataIn,dataOut]))
idcomp = find(R==max(R)); idcomp = idcomp(1);
idcomp = find(R==min(R)); idcomp = idcomp(1);
% idcomp = find(R==nanmedian(R)); idcomp = idcomp(2);
idcomp = find(data(:,1)==201)
idcomp = find(data(:,1)==165)
idcomp = find(data(:,1)==79)

for k = 1:N
    idcomp=k;
    a = [dataIn,dataOut];
    a = a(idcomp,:)
    
    figure
    for i=1:length(priornamesout)
        pdfRanges = getfield(predprior.ranges, priornamesout{i});
        pdfp = getfield(predprior.pdf, ['prior_',priornamesout{i}]);
        if(i==1)
            title(num2str(k));
        end
        
        subplot(6,1,i)
        %patch([pdfRanges(:); pdfRanges(end);pdfRanges(1)], [pdf(1);pdf(:);0;0], 'r')
        for j=1:length(pdfp)
            try
                patch([pdfRanges(j:(j+1)); pdfRanges(j+1);pdfRanges(j)], [pdfp(j);pdfp(j);0;0], 'r')
                hold on
                pdf = getfield(pred.pdf, [priornamesout{i}]);
                pdf = pdf(idcomp,:);
                hp=patch([pdfRanges(j:(j+1)); pdfRanges(j+1);pdfRanges(j)], [pdf(j);pdf(j);0;0], 'k')
                set(hp, 'faceAlpha',0.5)
            catch
                % skip
                %         patch([pdfRanges(j); pdfRanges(j);pdfRanges(j)], [pdfp(j);pdfp(j);0;0], 'r')
                %         hold on
                %         pdf = getfield(pred.pdf, [priornamesout{i}]);
                %         pdf = pdf(idcomp,:);
                %         patch([pdfRanges(j); pdfRanges(j);pdfRanges(j)], [pdf(j);pdf(j);0;0], 'k')
                %
            end
            
        end
        plot(a(i),0.5,'m*')
        xlabel(priornamesout{i})
        xlim([min(pdfRanges) max(pdfRanges)])
        ylim([0 1])
    end
    
    pause
    clf
end


% make cleaner plots
pdf = pred.pdf.(nodeNamesOut{1});
pdfRanges = pred.ranges.(nodeNamesOut{1});
[stats] = PDF2Stats(pdfRanges, pdf, 0.5);
py975 = getPy(pdfRanges, pdf, .975);
py95 = getPy(pdfRanges, pdf, .95);
py75 = getPy(pdfRanges, pdf, .75);
py25 = getPy(pdfRanges, pdf, .25);
py05 = getPy(pdfRanges, pdf, .05);
py025 = getPy(pdfRanges, pdf, .025);
% get mostproblist for cheryl, including all most probs if there is a tie
% store, k,ti(k),mostprob
mostprobforch = [];
for k=1:N
    id = find(pdf(k,:)==max(pdf(k,:)));
    mostprobforch = [mostprobforch; [k+0*id(:),ti(k+0*id(:)),pdfRangePlot(id)]];
end
idgood = find(~isnan(dataOut) & ~isnan(stats.bMean(1:N))); M = length(idgood)
[bM,brmse1,skM] = regr_xzw([ones(M,1) stats.bMean(idgood) ],dataOut(idgood));%, 1./stats.bStd(idgood))
rmseM = sqrt(nanmean((stats.bMean(idgood)-dataOut(idgood)).^2))
meaneM = nanmean((stats.bMean(idgood)-dataOut(idgood)))
fprintf('LR=%.0f\n', sum(likelihoodRatio))


%% plot stats
figure(100)
clf
patch([ti;flipud(ti)],[py025;flipud(py975)],[0,0,0]+0.75)
hold on
patch([ti;flipud(ti)],[py05;flipud(py95)],[0,0,0]+0.5)
patch([ti;flipud(ti)],[py25;flipud(py75)],[0,0,0]+.25)
%plot(ti,dataOut(1:N),'b.')  
plot( ti,dataOut(1:N),'c') % for input 
plot(ti(1:1:N), stats.bMean(1:1:N), 'g+', 'linewidth', 2)
plot(ti(1:1:N), stats.bMostProb(1:1:N), 'mo', 'linewidth',2)
plot(mostprobforch(:,2),mostprobforch(:,3), 'mo', 'linewidth',2)
% plot(ti, interp2(xs,ti,Hrms_m,xi(end),ti), 'g.')
xlim([ti(1) ti(end)])
%datetick('x',6,'keeplimits')
axis tight
ylim([min(pdfRanges) max(pdfRanges)])
xlabel('id')
ylabel('rate (m/y)')
legend('p95','p90','p50', 'output obs.', 'mean', 'most prob.')
title(sprintf('%s update: LR= %.2f', testcase, sum(likelihoodRatio)), 'interpreter','none')
% save
savename = [netRoot,'_probseries_',testcase,'.',num2str(fix(now)),'.fig'];
hgsave(savename)
savename = [netRoot,'_probseries_',testcase,'.',num2str(fix(now)),'.eps'];
print('-depsc2', savename)
savename = [netRoot,'_probseries_',testcase,'.',num2str(fix(now)),'.eps'];
print('-dill', savename)
savename = [netRoot,'_probseries_',testcase,'.',num2str(fix(now)),'.ai'];
hgsave(savename)


%% colorful version
% plot stats
figure(101)
clf
patch([ti;flipud(ti)],[py025;flipud(py975)],[0,0,0]+0.75)
hold on
patch([ti;flipud(ti)],[py05;flipud(py95)],[0,0,0]+0.5)
patch([ti;flipud(ti)],[py25;flipud(py75)],[0,0,0]+.25)
plot( ti,dataOut(1:N),'c.','linewidth',4) % for input 
scatter(ti(1:1:N), stats.bMean(1:1:N),20,stats.bStd(1:1:N), 'o', 'linewidth', 2)
%scatter(ti(1:1:N), stats.bMostProb(1:1:N),20,stats.pMaxProb(1:1:N), 'x', 'linewidth',2)
plot(ti(1:1:N), stats.bMostProb(1:1:N), 'mx', 'linewidth',2)
% plot(ti, interp2(xs,ti,Hrms_m,xi(end),ti), 'g.')
xlim([ti(1) ti(end)])
%datetick('x',6,'keeplimits')
axis tight
ylim([-5 0])
xlabel('id')
ylabel('rate (m/y)')
legend('p95','p90','p50', 'output obs.', 'mean', 'most prob.')
title(sprintf('%s update: LR= %.2f', testcase, sum(likelihoodRatio)), 'interpreter','none')

savename = [netRoot,'_probseries_',testcase,'.',num2str(fix(now)),'.eps'];
print('-depsc2', savename)
savename = [netRoot,'_probseries_',testcase,'.',num2str(fix(now)),'.eps'];
print('-dill', savename)
savename = [netRoot,'_probseries_',testcase,'.',num2str(fix(now)),'.ai'];
hgsave(savename)

%%  morestatistics
figure
%plot(dataOut(1:N),stats.bMean(1:1:N), '.', dataOut(1:N),stats.bMostProb(1:1:N), '+', )
%plot(interp2(xs,ti,Hrms_m,xi(3),ti),dataOut(1:N),'b.', stats.bMean(1:1:N),dataOut(1:N),'ko', stats.bMostProb(1:1:N),dataOut(1:N),'r+',[0,2],[0,2],'k--')
plot( stats.bMean(1:1:N),dataOut(1:N),'ko', stats.bMostProb(1:1:N),dataOut(1:N),'r+',[-10 1],[-10 1],'k--')
axis equal
grid on
xlabel('predicted (m/y)')
ylabel('observed (m/y)')
%legend('TG83', 'Bayes mean', 'Bayes most prob.')
legend('Bayes mean', 'Bayes most prob.')


[bML,brmse1,skML] = regr_xzw([ones(N,1) stats.bMostProb(1:1:N) ],dataOut(1:N))
[bM,brmse1,skM] = regr_xzw([ones(N,1) stats.bMean(1:1:N) ],dataOut(1:N))
% [bTG,brmse1,skTG] = regr_xzw([ones(N,1) interp2(xs,ti,Hrms_m,xi(3),ti) ],dataOut(1:N))
rmseM = sqrt(nanmean((stats.bMean(1:1:N)-dataOut(1:N)).^2))
meaneM = nanmean((stats.bMean(1:1:N)-dataOut(1:N)))
rmseML = sqrt(nanmean((stats.bMostProb(1:1:N)-dataOut(1:N)).^2))
meaneML = nanmean((stats.bMostProb(1:1:N)-dataOut(1:N)))
% rmseTG = sqrt(nanmean((interp2(xs,ti,Hrms_m,xi(3),ti)-dataOut(1:N)).^2))
% meaneTG = nanmean((interp2(xs,ti,Hrms_m,xi(3),ti)-dataOut(1:N)))





% trick, save the stats from the real prediction as stats1 the, load the prediction where anser is provided as input. Compare stats1 to stats variables for a binned comparison
