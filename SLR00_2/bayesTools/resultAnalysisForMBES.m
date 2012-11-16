function outputstr = resultAnalysisForMBES(testName,jobName, sites, suData, burialData, compareType)

if(~exist('compareType'))
    compareType = 'UNIFORM'; % 'MIWDOC'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
try
    % impact
    a=load([testName,'\ShearStrength.mat'])
catch
    % scour
    a=load([testName,'\waveOrbitalVelocity.mat'])
end
try
    b=load([testName,'\waterDepth.mat'])
catch
    warning('no water depth data found')
    b.lon = nan;
    b.lat = nan;
end
try
    c=load([testName,'\mineType.mat'])
catch
    warning('no mine type data found')
end
d=load([testName,'\',jobName,'.mat'])

% process from matlab output
Lat = d.burialPrediction.lat;
Lon = d.burialPrediction.lon;
dlat = d.burialPrediction.dLat;
dlon = d.burialPrediction.dLon;
NLatLon = length(Lat(:));

% find analysis sites here
gridID = zeros(length(sites(:,1)),1);
for i=1:length(gridID)
    r = (d.Lat-sites{i,2}).^2 + (d.Lon-sites{i,3}).^2;
    gridID(i) = find(r==min(r(:)));
end
gridID

% plots
% locations
figure
%plot(d.Lon(:),d.Lat(:),'r+', a.lon(:),a.lat(:),'b.', b.lon(:),b.lat(:),'bo', burialData(:,2),burialData(:,1),'k*',d.Lon(gridID),d.Lat(gridID),'rs')
plot(d.Lon(:),d.Lat(:),'r+', a.lon(:),a.lat(:),'bo', b.lon(:),b.lat(:),'bx', burialData(:,2),burialData(:,1),'k*',d.Lon(gridID),d.Lat(gridID),'rs')
grid on ; axis equal; axis square; axis tight
xlabel('lat')
ylabel('lon')
text(d.Lon(gridID)+dlon,d.Lat(gridID),[char(sites(:,1)),num2str(gridID)],'fontSize',12,'FontWeight','bold','BackgroundColor','w')
title(testName, 'interpreter','none')
legend('grid loc.','ShearStrengths','depths','burials','analysis')
figure(gcf)
hgsave([testName,filesep,jobName,'map.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputFieldsToWrite = {
    'bMean'
    'bStd'
    'bMostProb'
    'b95'
    'medal'
    };
outputStats = computeStatsFromPDF(d.burialPrediction.pdfRanges, d.burialPrediction.pdf, .5, inputFieldsToWrite);

switch compareType
    case 'MIWDOC'
        % compare to NAVO
        [latNavo, lonNavo, pdfRangesNavo, pdfDataNavo, pdfClassNAVO, pdfCountNAVO] = ...
            NAVOSedshape2medalPDF('inputData\corpusTest\u_texas_regsedse_200209.shp', Lat,Lon,dlat,dlon);
        %sedshape2medalPDF('inputData\corpusTest\u_texas_regsedse_200209.shp', Lat,Lon,10*dlat,10*dlon);
        % NOTE: INTERPRET MIW DOCTRINE AS 95% confidence, this is implied here:

    case 'UNIFORM'
        % just plot the a uniform distribution for burial
        pdfRangesNavo = [0:5:100];
        pdfDataNavo = 0*pdfRangesNavo(2:end) + 1;
        pdfDataNavo = repmat(pdfDataNavo/sum(pdfDataNavo),NLatLon,1);
    case 'UNIFORMinit'
        % just plot the a uniform distribution for burial
        pdfRangesNavo = [0:5:100];
        pdfDataNavo = 0*pdfRangesNavo(2:end) + 1;
        pdfDataNavo = repmat(pdfDataNavo/sum(pdfDataNavo),NLatLon,1);
        
    case 'MIWBCAT'
        % get MIW prediction for this category
        pdfRangesNavo = [0:5:100];
        %pdfDataNavo = repmat(0*pdfRangesNavo(2:end),NLatLon,1) +0.05/length(pdfRangesNavo-1);
        pdfDataNavo = repmat(0*pdfRangesNavo(2:end),NLatLon,1);
        medalCats = outputStats.bMedalCat(gridID);
        minPctBurial = [0,0,10,20,75];
        maxPctBurial = [0,10,20,75,100];
        minPctBurial = minPctBurial(medalCats);
        maxPctBurial = maxPctBurial(medalCats);
        for i=1:length(gridID)
            idpdf = find(pdfRangesNavo(1:(end-1))>=minPctBurial(i) & pdfRangesNavo(2:(end))<=maxPctBurial(i));
            pdfDataNavo(gridID(i),idpdf) = 1;
            %pdfDataNavo(gridID(i),:) = 0.05/length(pdfRangesNavo-1) + 0.95*pdfDataNavo(gridID(i),:)/sum(pdfDataNavo(gridID(i),:));
            pdfDataNavo(gridID(i),:) = pdfDataNavo(gridID(i),:)/sum(pdfDataNavo(gridID(i),:));
        end
end

% get shear strength stats
inputFieldsToWrite = {
    'bMean'
    'bStd'
    'bMostProb'
    'b95'
    };

inputID = strmatch('shear',lower(char(d.burialPrediction.inputPDF.name)) );
if(isempty(inputID))
    % try grain size
    inputID = strmatch('waveorbitalvelocity',lower(char(d.burialPrediction.inputPDF.name)) );
end
inputStats = computeStatsFromPDF(d.burialPrediction.inputPDF(inputID).pdfRanges, d.burialPrediction.inputPDF(inputID).inputPDF, .5, inputFieldsToWrite);

% analysis
outputstr = {};
for i=1:length(gridID)
    % shear strength
    try
    figure
    subplot 211
    [xPlot,pPlot,hp] = plotPDF(d.burialPrediction.inputPDF(inputID).pdfRanges, d.burialPrediction.inputPDF(inputID).inputPDF(gridID(i),:));
    set(hp, 'faceColor','r')
    
    % show in situ data
    probSuData = nan;
    idData = find(abs(suData(:,1)-d.Lat(gridID(i)))<d.dlat & abs(suData(:,2)-d.Lon(gridID(i)))<d.dlon);
    if(length(idData)>0)
        if(length(idData)>1)
            xdata = 0.1*(max(suData(idData,3))-min(suData(idData,3)))*[0:10]+min(suData(idData,3));
            if(var(xdata)<eps)
                xdata = d.burialPrediction.inputPDF(inputID).pdfRanges'; % row vect
            end
        else
            xdata = suData(idData,3)*[0.9 1.1];
        end
        [pdata,xdata] = hist(suData(idData,3),xdata);
        hold on
        bar(xdata, max(pPlot)*pdata/max(pdata), 'b')
        % get prob. of data
        probSuData = getProbOfData(d.burialPrediction.inputPDF(inputID).pdfRanges, d.burialPrediction.inputPDF(inputID).inputPDF(gridID(i),:), suData(idData,3)+eps);
        %26jan2007: xunk = [min(d.burialPrediction.inputPDF(inputID).pdfRanges):1:max(d.burialPrediction.inputPDF(inputID).pdfRanges)];
        dunk = [max(d.burialPrediction.inputPDF(inputID).pdfRanges)-min(d.burialPrediction.inputPDF(inputID).pdfRanges)]/10;
        xunk = [min(d.burialPrediction.inputPDF(inputID).pdfRanges):dunk:max(d.burialPrediction.inputPDF(inputID).pdfRanges)];
        punk = 0*xunk(2:end)+1; punk = punk/sum(punk);
        probSuDataUNK = getProbOfData(xunk, punk, suData(idData,3));
        
        % get stats
        ylimNow = ylim;
        %text(0, 0.9*ylimNow(2), sprintf('   mean=%.2f, std=%.2f, P-ratio=%.3g', inputStats.bMean(gridID(i)),inputStats.bStd(gridID(i)),exp(log(prod(probSuData))-log(prod(probSuDataUNK)))))
        text(0, 0.9*ylimNow(2), sprintf('   mean=%.2f, std=%.2f, logP=%.3g', inputStats.bMean(gridID(i)),inputStats.bStd(gridID(i)),    (log(prod(probSuData))-log(prod(probSuDataUNK)))))
        outputstri = sprintf('   mean=%.2f, std=%.2f, logP=%.3g', inputStats.bMean(gridID(i)),inputStats.bStd(gridID(i)),    (log(prod(probSuData))-log(prod(probSuDataUNK))));
        text(0, 0.9*ylimNow(2), outputstri);
        outputstr{i,2} = outputstri;
        ylim(ylimNow)
        hold off
        
    end
    xlim([0 d.burialPrediction.inputPDF(inputID).pdfRanges(end-1)])
    grid on
    xlabel(d.burialPrediction.inputPDF(inputID).name)
    ylabel('prob.')
    legend('MBES','in situ obs.')
    titlestr = sprintf('[%s] site %s - grid location %d (Lat.=%.5f,Lon.=%.5f)', jobName,sites{i,1},gridID(i), d.Lat(gridID(i)),d.Lon(gridID(i)));
    title(titlestr, 'interpreter','none')
    outputstr{i,1} = titlestr;

    %%%%%%%%
    % burial
    subplot 212
    % add navo prediction
    [xnPlot,pnPlot,hp] = plotPDF(pdfRangesNavo,pdfDataNavo(gridID(i),:));
    hold on

    % add mbes
    [xPlot,pPlot,hp] = plotPDF(d.burialPrediction.pdfRanges,d.burialPrediction.pdf(gridID(i),:));
    set(hp, 'faceColor','r')
    %%alpha(hp,.75) % make transparent

    % add reality
    probBurialData=nan;
    probBurialDataNAVO=nan;
    idData = find(abs(burialData(:,1)-d.Lat(gridID(i)))<d.dlat & abs(burialData(:,2)-d.Lon(gridID(i)))<d.dlon);
    if(length(idData)>0)
        if(length(idData)>1 & max(diff(abs(burialData(idData,3))))>0)
            xdata = 0.1*(max(burialData(idData,3))-min(burialData(idData,3)))*[0:10]+min(burialData(idData,3));
            if(var(xdata)==0)
                xdata = d.burialPrediction.pdfRanges'; % row vect
            end
        else
            xdata = burialData(idData(1),3)*[0.99,1.1];
        end
        [pdata,xdata] = hist(burialData(idData,3),xdata);
        hold on
        hd = bar(xdata, max(pPlot)*pdata/max(pdata), 'b');
        % get prob. of data
        probBurialData = getProbOfData(d.burialPrediction.pdfRanges, d.burialPrediction.pdf(gridID(i),:), burialData(idData,3)+eps);
        probBurialDataNAVO = getProbOfData(pdfRangesNavo, pdfDataNavo(gridID(i),:), burialData(idData,3));
        xunk = [min(d.burialPrediction.pdfRanges):1:max(d.burialPrediction.pdfRanges)];
        punk = 0*xunk(2:end)+1; punk = punk/sum(punk);
        probBurialDataUNK = getProbOfData(xunk, punk, burialData(idData,3));
    else
        probBurialData = 0;
        probBurialDataNAVO = 0;
        probBurialDataUNK = 0;
    end

    % finish
    % add outline of the alternate prediction
    plot(xnPlot,pnPlot,'k--','linewidth',2)
    xlim([0 100])
    grid on
    xlabel(d.burialPrediction.name)
    ylabel('prob.')
    legend(compareType,'MBES', 'Obs.')

    % get stats from shapefile
    ylimNow = ylim;
    %text(0, 0.9*ylimNow(2), sprintf('   mean=%.2f, std=%.2f, b_{most prob.}=%.2f, b_{95}=%.2f, P-ratio_{MBES}=%.3f, P-ratio_{%s}=%.3f', ...
        %outputStats.bMean(gridID(i)),outputStats.bStd(gridID(i)),outputStats.bMostProb(gridID(i)),outputStats.b95(gridID(i)),...
        %exp(sum(log(probBurialData))-sum(log(probBurialDataUNK))),compareType, exp(sum(log(probBurialDataNAVO))-sum(log(probBurialDataUNK)))))
    %text(0, 0.9*ylimNow(2), sprintf('   mean=%.2f, std=%.2f, b_{most prob.}=%.2f, b_{95}=%.2f, logP-ratio_{MBES}=%.3f, logP-ratio_{%s}=%.3f', ...
    bmeani = outputStats.bMean(gridID(i));
    bstdi = outputStats.bStd(gridID(i));
    bmostprobi = outputStats.bMostProb(gridID(i));
    b95i = outputStats.b95(gridID(i));
    logPmbesi = (log(prod(probBurialData))-log(prod(probBurialDataUNK)));
    logPotheri =(log(prod(probBurialDataNAVO))-log(prod(probBurialDataUNK)));
    burialstri = sprintf('   mean=%.2f, std=%.2f, b_{most prob.}=%.2f, b_{95}=%.2f, logP{MBES}=%.3f, logP{%s}=%.3f', ...
        bmeani,bstdi,bmostprobi,b95i,...
        logPmbesi,compareType, logPotheri);
    title(burialstri)
    outputstr{i,3} = burialstri;
    ylim(ylimNow)
    hold off

    figure(gcf)
    orient tall
    figsavename = [testName,filesep,jobName,'_analysis_vs_',compareType,num2str(gridID(i)),'.fig'];
    if(exist(figsavename)==0)
        hgsave(figsavename)
        save([testName,filesep,jobName,'_analysis_vs_',compareType,num2str(gridID(i)),'.mat'], 'bmeani','bstdi','bmostprobi','b95i','logPmbesi','logPotheri')
    else
        figsavename = [testName,filesep,jobName,'_analysis_vs_',compareType,num2str(gridID(i)),'_v2.fig'];
        hgsave(figsavename)
        save([testName,filesep,jobName,'_analysis_vs_',compareType,num2str(gridID(i)),'_v2.mat'], 'bmeani','bstdi','bmostprobi','b95i','logPmbesi','logPotheri')
    end

    % finally, save output
    catch
        lasterr
        keyboard
    end
end
