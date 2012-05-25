function [stats] = computeStatsFromPDF(pdfRanges, pdf, alpha, doFields)
% [stats] = computeStatsFromPDF(pdfRanges, pdf, alpha, doFields)
%
% compute statistics from pdf
% 
% Input
%   pdfRanges, the bins of the pdf
%   pdf, the pdf 
%   alpha, a probability level in order to compute prob(x<x_alpha) = alpha
%   optional doFields, list of desired outputs
%
% Output
%   stats, structure with fields
%       .palpha, the value whos prob(x<palpha)=alpha
%       .mean
%       .p025
%       .p975
%       .std
%       .maxProb, bin with most prob
%       etc.
%  NOTE: these  values are interpolated to an arbitrary discrete set of ranges to avoid too many values in shapefile output
% it appears that 12 char is limit to stats name

% to make all answers easy to categorize, interp results to these ranges
burialInterpValues = min([0,pdfRanges(1)]);
% changed the interp interval from 1/20 to 1/100 to resolve thresholds
burialInterpValues = burialInterpValues + [0:(1/100):1]*(pdfRanges(end) - burialInterpValues);
probabilityInterpValues = [0:0.05:1];

% make sure pdf is normalized
pdf=pdf./repmat(sum(pdf,2),1,size(pdf,2));

% do we do all possible statistics?
if(nargin==3)
    doFields = 'all';
end

% compute
% basically, interp all values to their class in the pdfRanges 
[Nlocs,Npdf] = size(pdf);
NpdfRanges = length(pdfRanges);
id = find(isnan(pdf(:,1)));
blank = ones(Nlocs,1);
blank(id) = nan;

% some statistics depend on whether we have discrete or continuous
if(Npdf==NpdfRanges)
    % discrete
    id1=1:Npdf;
    id2=id1;
elseif(Npdf==(NpdfRanges-1))
    id1 = 1:Npdf;
    id2 = id1+1;
else
    error('size of pdfRanges must equal size of pdf (discrete) or (size + 1) of pdf (continous)')
end

% get user-specified burial percentile
if(length(strmatch('bAlpha', doFields))>0 | strcmp('all',doFields))
    stats.bAlpha = pdfRanges(1+sum(cumsum(pdf')<=alpha)).*blank; % alpha prob. cumulative
    % for ALL OUTPUTS provide full range of ouput values in first 1:(number of ranges) bins
    stats.bAlpha(end+[1:NpdfRanges])=pdfRanges;
end

% get the burial tails of the user-specified percentile
if(length(strmatch('bAlphaPlus', doFields))>0 | strcmp('all',doFields))
    dalpha = (1-alpha)/2;
    % get the tails for the user-specified percentile
    stats.bAlphaPlus  = pdfRanges(1+sum(cumsum(pdf')<=(1-dalpha))).*blank; % +alpha/2 prob. cumulative
    stats.bAlphaPlus(end+[1:NpdfRanges])=pdfRanges;
    stats.bAlphaMinus = pdfRanges(1+sum(cumsum(pdf')<=  (dalpha))).*blank;; % -alpha/2 prob. cumulative
    stats.bAlphaMinus(end+[1:NpdfRanges])=pdfRanges;
end

% burial 2.5th percentiles
if(length(strmatch('b025', doFields))>0 | strcmp('all',doFields))
    stats.b025 = pdfRanges(1+sum(cumsum(pdf')<=0.025)).*blank;k
    stats.b025(end+[1:NpdfRanges])=pdfRanges;
    stats.b975 = pdfRanges(1+sum(cumsum(pdf')<=0.975)).*blank;
    stats.b975(end+[1:NpdfRanges])=pdfRanges;
end

% burial 95th percentiles
if(length(strmatch('b95', doFields))>0 | strcmp('all',doFields))
    stats.b95 = pdfRanges(1+sum(cumsum(pdf')<=0.95)).*blank;
    stats.b95(end+[1:NpdfRanges])=pdfRanges;
end

% burial 25th percentiles
if(length(strmatch('b25', doFields))>0 | strcmp('all',doFields))
    stats.b25 = pdfRanges(1+sum(cumsum(pdf')<=0.25)).*blank;
    stats.b25(end+[1:NpdfRanges])=pdfRanges;
    stats.b50 = pdfRanges(1+sum(cumsum(pdf')<=0.5)).*blank;
    stats.b50(end+[1:NpdfRanges])=pdfRanges;
    stats.b75 = pdfRanges(1+sum(cumsum(pdf')<=0.75)).*blank;
    stats.b75(end+[1:NpdfRanges])=pdfRanges;
end

% get burial class that is most probable
if(length(strmatch('bMostProb', doFields))>0 | strcmp('all',doFields))
    % bin with maxProb
    stats.bMostProb = repmat(nan,Nlocs,1);
    for i =1:Nlocs
        id = max(find(pdf(i,:)==max(pdf(i,:))));
        if(length(id)>0)
            id = id(end);
            stats.bMostProb(i) = 0.5*(pdfRanges(id1(id))+pdfRanges(id2(id))) ;
        end
    end
    stats.bMostProb(end+[1:(NpdfRanges-1)])= 0.5*(pdfRanges(id1)+pdfRanges(id2));
end

% do mean,std
if(length(strmatch('bMean', doFields))>0|length(strmatch('bStd', doFields))>0 | strcmp('all',doFields))
    stats.bMean = (pdf*[pdfRanges(id1)+pdfRanges(id2)]/2).*blank;
    stats.bStd = zeros(Nlocs,1);
    % get subbin scale variance
    binWidthVar = 0.5*pdf*[pdfRanges(id2)-pdfRanges(id1)]/sqrt(3);
    for i=1:Nlocs
        stats.bStd(i) = pdf(i,:)*(([pdfRanges(id1)+pdfRanges(id2)]/2-stats.bMean(i)).^2);
        % add subbin scale variance (only relevant if all in one bin)
        stats.bStd(i) =  sqrt(binWidthVar(i)+stats.bStd(i));
    end
    % force to known category
    stats.bMean = interp1(burialInterpValues,burialInterpValues,stats.bMean,'nearest');  
    stats.bMean(end+[1:length(burialInterpValues)])=burialInterpValues;
    % force to known category,  but force smallest value to be binwidth
    %stats.bStd = interp1(burialInterpValues,;burialInterpValues,stats.bStd,'nearest').*blank;
    stats.bStd = interp1(burialInterpValues,[burialInterpValues+sqrt(min(binWidthVar))],stats.bStd,'nearest').*blank;
    stats.bStd(end+[1:length(burialInterpValues)])=burialInterpValues;
end

% get this measure of prob. concentration, related to sampling
if(length(strmatch('pMostProb', doFields))>0 | strcmp('all',doFields))
    % get maxProb occuring in any bin
    stats.pMaxProb =interp1(probabilityInterpValues,probabilityInterpValues,(max(pdf')'), 'nearest').*blank;
    stats.pMaxProb(end+[1:length(probabilityInterpValues)])=probabilityInterpValues; 
end

% probability that burial exceeds 50%
if(length(strmatch('pBurialGT50', doFields))>0 | strcmp('all',doFields))
    bAlpha = 50;
    id = find(pdfRanges(id1)>=50);
    stats.pBurialGT50 = sum(pdf(:,id),2);
    % interp to static range
    stats.pBurialGT50 =interp1(probabilityInterpValues,probabilityInterpValues,stats.pBurialGT50, 'nearest').*blank;
    stats.pBurialGT50(end+[1:length(probabilityInterpValues)])=probabilityInterpValues; 
end

% get burial into medal categories
if(length(strmatch('medal', doFields))>0 | strcmp('all',doFields))
    [NAVO_MIW_REDUCED15, NAVO_MASTER_RECLASS, SHEPARD_ID_DEFINITION, NAVO_5CAT_BURIAL] = NAVOreclass;
    % NAVO_MASTER_RECLASS field 1 is the cat, field 5 is the 5cat burial
    navoBurialCat = [NAVO_5CAT_BURIAL{:,1}];
    minPctBurial = [NAVO_5CAT_BURIAL{:,3}];
    maxPctBurial = [NAVO_5CAT_BURIAL{:,4}];
    minPctBurial = minPctBurial(:);
    maxPctBurial = maxPctBurial(:);
    % get the navo burials
    %[minPctBurial, maxPctBurial, navoBurialCat, penColors] = NAVOburialCat;
    navoPDFRanges = [minPctBurial; maxPctBurial(end)];
    navoPDFRanges = [minPctBurial(~isnan(minPctBurial(:))); max(maxPctBurial)];
    navoPDFWidth = diff(navoPDFRanges);
    % get prob in each category
    medalCatsPDF = zeros( Nlocs, length(minPctBurial) );
    % put the burial prediction on a finer grid
    outputRanges = [0:10:100];
    outputPDF = pdfInterp(pdfRanges,pdf, outputRanges);
    outputRanges = 5:10:100;
    % scan each medal category
    docats = find(~isnan(minPctBurial));
    for i=docats(:)' % force row-vec
        % get prediction pdf categories that correspond to this medal category
        id = find(outputRanges>minPctBurial(i) & outputRanges<=maxPctBurial(i));
        medalCatsPDF(:,i) = sum(outputPDF(:,id),2); % here is prob in this range, at all locations
    end
    % get most likely medal category at each location
    stats.bMedalCat = repmat(nan,Nlocs,1);
    stats.pMedalProb = repmat(nan,Nlocs,1);
    stats.pMedalGTProb = repmat(nan,Nlocs,1);
    stats.pMedalLTProb = repmat(nan,Nlocs,1);
    %stats.pMedalRatio = repmat(nan,Nlocs,1);
    for i = 1:Nlocs
        % get category with maximum probability
        id = max(find(medalCatsPDF(i,:)==max(medalCatsPDF(i,:))));
        if(length(id)>0)
            id = id(end);
            % get largest burial if a tie
            stats.bMedalCat(i) = navoBurialCat(id); % this is a category (0-5)
            % get the likelihood here
            stats.pMedalProb(i) = medalCatsPDF(i,id);
            % now, operator wants to know what liklihood is that burial exeeds medal cat
            if(id<max(docats))
                stats.pMedalGTProb(i) = sum(medalCatsPDF(i,(id+1):end));
            else
                % already at largest bin, no prob of greater
                stats.pMedalGTProb(i) = 0;
            end
            % and what is prob that burial is less
            if(id>min(docats))
                stats.pMedalLTProb(i) = sum(medalCatsPDF(i,1:(id-1)));
            else
                % already at largest bin, no prob of greater
                stats.pMedalGTProb(i) = 0;
            end
            % get the power of this statistic
            if(0)%navoPDFWidth(id)>0)
                stats.pMedalRatio(i) = stats.pMedalProb(i)*100/navoPDFWidth(id);
            elseif(0)
                stats.pMedalRatio(i) = nan;
            end
        end
    end
    % force to some sensible categories
    stats.bMedalCat(end+[1:length(navoBurialCat)])=navoBurialCat; 
    stats.pMedalProb =interp1(probabilityInterpValues,probabilityInterpValues,stats.pMedalProb, 'nearest');
    stats.pMedalProb(end+[1:length(probabilityInterpValues)])=probabilityInterpValues; 
    stats.pMedalLTProb =interp1(probabilityInterpValues,probabilityInterpValues,stats.pMedalLTProb, 'nearest');
    stats.pMedalLTProb(end+[1:length(probabilityInterpValues)])=probabilityInterpValues; 
    stats.pMedalGTProb =interp1(probabilityInterpValues,probabilityInterpValues,stats.pMedalGTProb, 'nearest');
    stats.pMedalGTProb(end+[1:length(probabilityInterpValues)])=probabilityInterpValues; 
end

if(exist('stats')==0)
    fprintf('doFields did not match the possibilities'); doFields
end

