function savedStats = BayesPostProc(netName,nodeNamesOut, pred, dataNames, ...
    data,  figdir,testcasename)

% import netica stuff
import norsys.netica.*;
% make the netica stuff global --> kinda kludgy but required...
global NETICA_ENV NETICA_NET
    VARNAMES = nodeNamesOut
    cnt = 1;
    N = length(data);
    did = 1:1:N;
    for j=1:length(VARNAMES)
        cnt = cnt +1;
        VARNAME = VARNAMES{j};
        varid = strmatch(VARNAME, dataNames);
        parid = strmatch(VARNAME,nodeNamesOut);
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
        LR = sum(pred.logLikelihoodRatio(:,parid));
        
        % regression
        [bML,brmse1,skML] = regr_xzw([ones(N,1) stats.bMostProb(1:1:N) ],data(did,varid)-nanmean(data(did,varid)), 1./(stats.bStd(1:1:N)+eps));
        [bM,brmse1,skM] = regr_xzw([ones(N,1) stats.bMean(1:1:N) ],data(did,varid)-nanmean(data(did,varid)), 1./(stats.bStd(1:1:N)+eps));
        rmseM = sqrt(nanmean((stats.bMean(1:1:N)-data(did,varid)).^2));
        meaneM = nanmean((stats.bMean(1:1:N)-data(did,varid)));
        rmseML = sqrt(nanmean((stats.bMostProb(1:1:N)-data(did,varid)).^2));
        meaneML = nanmean((stats.bMostProb(1:1:N)-data(did,varid)));
        
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