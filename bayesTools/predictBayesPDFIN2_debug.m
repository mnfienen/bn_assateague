function pred = predictBayesPDFIN2(netName, nodeNamesIn, nodeNamesOut, pdfIn, dataOut)

%%%%%%%%%%%%%%%%%%%%
% now, try the net
%netName = 'barrierIsland_RDHighLow_Dchange_v3.dne'
% start the network
netStatus = NetStartStop(netName)

% all is well if status is empty
if(isempty(netStatus))
    % get to the net
    global NETICA_ENV NETICA_NET NETICA_NETNAME NETICA_NODESTRUCT NETICA_STREAMER NETICA_ALL_NODES
else
    error(netStatus)
end

% get info
% display titles
nodetitles = {NETICA_NODESTRUCT.title}'
% display nodenames
nodenames = {NETICA_NODESTRUCT.name}'
% other info
%inputData(i).nodeID = strmatch(requiredInputNames{i,2}, char(NETICA_NODESTRUCT.title)); % id to talk to bayesian net
%inputData(i).pdfRanges = [NETICA_NODESTRUCT(inputData(i).nodeID).levels]; % range of values in the bayesian net
%inputData(i).pdfStates = {NETICA_NODESTRUCT(inputData(i).nodeID).state.name};
% get bins, either number of pdfRanges -1 (continuous) or number of pdfRanges (discrete)
%inputData(i).Nbins = NETICA_NODESTRUCT(inputData(i).nodeID).Nbeliefs; %length(inputData(i).pdfRanges)-1; % number of bins between these ranges
%inputData(i).pdf = zeros(1,inputData(i).Nbins);

% What will we predict?
predictID = zeros(length(nodeNamesOut),1);
for i=1:length(nodeNamesOut)
    j = strmatch(nodeNamesOut(i), nodenames);
    if(~isempty(j))
        predictID(i)=j;
    end
end
% what will we constrain
inputID = zeros(length(nodeNamesIn),1);
for i=1:length(nodeNamesIn)
    j = strmatch(nodeNamesIn(i), nodenames);
    if(~isempty(j))
        inputID(i)=j;
    end
end

% initialize the output 
flds = fieldnames(pdfIn)
m = length(flds)
ptmp = getfield(pdfIn,flds{1});
N = length(ptmp(:,1))
% [N,m] = size(dataIn);
if(m~=length(nodeNamesIn))
    error('dataIn must have same number of columns as nodeNameIn');
end
Nbins = zeros(length(nodenames),1);
pdfstruct = struct;
rangestruct = struct;
for j=1:length(nodenames)
    Nbins(j) = NETICA_NODESTRUCT(j).Nbeliefs;
    pdfstruct = setfield(pdfstruct, nodenames{j}, zeros(N,Nbins(j)));
    pdfLevels = [NETICA_NODESTRUCT(j).levels];
    rangestruct = setfield(rangestruct, nodenames{j}, pdfLevels)
 end

% get prior (before data - same for all data!)
% clear the net
NETICA_NET.retractFindings
for i=1:length(predictID)
    j = predictID(i); % j is index into net arrangement, i is index into data
    pdfstruct = setfield(pdfstruct, ['prior_',nodenames{j}], zeros(1,Nbins(j)));
    outputPDF = NETICA_ALL_NODES(j).getBeliefs;
    pdf = getfield(pdfstruct, ['prior_',nodenames{j}]);
    pdf(1,:) = outputPDF(:)';
    pdfstruct = setfield(pdfstruct, ['prior_',nodenames{j}], pdf);
end

% tack on some alternative pdfs
% prior prob. before data arrive
% update after data
% for j=1:length(predictID)
%     pdfstruct = setfield(pdfstruct, ['prior_',nodenames{predictID(j)}], zeros(N,Nbins(predictID(j))));
%     pdfstruct = setfield(pdfstruct, ['update_',nodenames{predictID(j)}], zeros(N,Nbins(predictID(j))));
% end

% loop through each input series and get netica prediction
for i=1:N
    fprintf('doing %d of %d\n', i,N)
    % assign inputs
    IDnum = i;

    % clear the net
    NETICA_NET.retractFindings

    % loop throught the nodenames, setting input
    dataid=0;
    for j = inputID(:)'
        % get data for this node
        dataid = dataid+1;
        newPDF = getfield(pdfIn, nodenames{j});
        nodenames{j}
        newPDF = newPDF(i,:)
        % update net if have data
        if(~(any(isnan(newPDF))))
            % enter data
            NETICA_ALL_NODES(j).enterLikelihood(newPDF);
        end
    end
    
    % store after all inputs!
    for j = inputID(:)'
        % get update for this variable
        nodenames{j}
        inputPDF = NETICA_ALL_NODES(j).getBeliefs
        % store
        pdf = getfield(pdfstruct, nodenames{j});
        pdf(i,:) = inputPDF(:)';
        pdfstruct = setfield(pdfstruct, nodenames{j}, pdf);
        pause
    end
    
    
    % update prediction
    for j = predictID';
        outputPDF = NETICA_ALL_NODES(j).getBeliefs;
        pdf = getfield(pdfstruct, nodenames{j});
        pdf(i,:) = outputPDF(:)';
        pdfstruct = setfield(pdfstruct, nodenames{j}, pdf);
    end
end
% done!
% should save output at this point
pred.pdf = pdfstruct;
pred.ranges=rangestruct;
%return


% below is post processing
%%%
% plot comparison of prediction to observation
% pdfRanges = rangestruct.DhighDt;
for j = predictID'
    pdfRanges = getfield(rangestruct, nodenames{j});
    pdfRangePlot = pdfRanges(2:end)-0.5*diff(pdfRanges)
    % show conturs at this Pvalue
    Pvalue = 0.90;

    %%%%%%%%%
    % prediction prior to updating (the mean for all situations)
    pdf0 = getfield(pdfstruct, ['prior_', nodenames{j}]);
    pdf0 = repmat(pdf0,N,1);
    pdf = pdf0;

    if(0)
        %    prior is same, repmat it
        figure(1); clf
        % remove floor
        pdf = pdf - min(pdf,[],2)*ones(size(pdf(1,:)))+1e-3; pdf =pdf./(sum(pdf,2)*ones(size(pdf(1,:))));
        % CORRECT FOR PCOLOR OFFSET!
        keyboard
        pcolor(1:N, pdfRanges, log10([pdf';zeros(1,N)])); shading flat
        caxis(log10([min(pdf(:)) max(pdf(:))])); colorbar;
        % add obs and percentiles
        hold on
        % do P%, P=Pvalue
        px = [[1:N]';flipud([1:N]')];
        id = find(isnan(px));id2 = find(~isnan(px));
        py = [getPy(pdfRanges, pdf, 1-Pvalue); flipud(getPy(pdfRanges, pdf, Pvalue))];
        px(id) = interp1(id2,px(id2), id);
        hp=patch(px, py, Pvalue*[1 1 1]); set(hp, 'FaceAlpha', 0.9)
        % observation
        plot(1:N,dataOut(:,j),  'linewidth',2)
        hold off
        xlabel('t')
        ylabel(nodetitles{j})
        title('prior prediction log prob.')
        hl=legend('log prob.',  sprintf('%.0f%% prob.',Pvalue*100), 'observed');
    end

    %%%%%%%%%
    % prediction prior to updating (the mean for all situations)
    pdf2 = getfield(pdfstruct, [nodenames{j}]);
    pdf = pdf2;
    if(0)
        figure(1); clf
        % remove floor
        pdf = pdf - min(pdf,[],2)*ones(size(pdf(1,:)))+1e-3; pdf =pdf./(sum(pdf,2)*ones(size(pdf(1,:))));
        % CORRECT FOR PCOLOR OFFSET!
        pcolor(1:N, pdfRanges, log10([pdf';zeros(1,N)])); shading flat
        caxis([min(pdfRanges) max(pdfRanges)]); colorbar;
        % add obs and percentiles
        hold on
        % do P%, P=Pvalue
        px = [[1:N]';flipud([1:N]')];
        id = find(isnan(px));id2 = find(~isnan(px));
        py = [getPy(pdfRanges, pdf, 1-Pvalue); flipud(getPy(pdfRanges, pdf, Pvalue))];
        px(id) = interp1(id2,px(id2), id);
        hp=patch(px, py, Pvalue*[1 1 1]); set(hp, 'FaceAlpha', 0.9)
        % observation
        plot(1:N,dataOut(:,j),  'linewidth',2)
        hold off
        xlabel('t')
        ylabel(nodetitles{j})
        title('updated prediction log prob.')
        hl=legend('log prob.',  sprintf('%.0f%% prob.',Pvalue*100), 'observed');
        %%%%%%%%%
    end
    
    % don't continue if no data for comparison
    if(exist('dataOut')~=1)
        return
    end


    % %%%%%%%%%
    % % updated prediction with forcing data
    % pdf2 = getfield(pdfstruct, nodenames{j});
    % pdf = pdf2;
    % figure(3); clf
    % pdf = pdf - min(pdf,[],2)*ones(size(pdf(1,:)))+1e-3; pdf =pdf./(sum(pdf,2)*ones(size(pdf(1,:))));
    % pcolor(x0-diff(x0(1:2))/2, pdfRanges, log10([pdf';zeros(1,N)])); shading flat
    % caxis([-3 0]); colorbar;hold on
    % px = [x0;flipud(x0)];
    % id = find(isnan(px));id2 = find(~isnan(px));
    % py = [getPy(pdfRanges, pdf, 1-Pvalue); flipud(getPy(pdfRanges, pdf, Pvalue))];
    % px(id) = interp1(id2,px(id2), id);
    % hp=patch(px, py, Pvalue*[1 1 1]); set(hp, 'FaceAlpha', 0.9)
    % % observation
    % plot(x0,dhighdt,'k', x0,d0,'b', x0,rlow,'g--',x0,rhigh,'g-',x0,dhighdt,'ko',  'linewidth',2)
    % hold off
    % xlabel('alongshore distance (km)')
    % ylabel('change in D_{high} (m)')
    % title('updated D_{high} and forcing log prob.')
    % hl=legend('log prob.',  sprintf('%.0f%% prob.',Pvalue*100), 'observed change', 'initial D_{high}','R_{low}','R_{high}');
    % % focus on sub region
    % xlim([29 31])
    %
    %%%%%%%%%%%%%%
    % scroll past the probability maps
    % for i=0:70
    %     figure(1)
    %     xlim(i + [-1 1]*1)
    %     figure(2)
    %     xlim(i + [-1 1]*1)
    %     pause
    % end

    %%%%%%%%%%%%
    %get prob data.
    jd = find(j==predictID)
    probModel0Data = zeros(N,1);
    probModel1Data = zeros(N,1);
    probModel2Data = zeros(N,1);
    likelihoodRatio = 0*pdfRangePlot;
    for k=1:length(pdfRangePlot)
        
        id = find(dataOut(:,jd)>pdfRanges(k) & dataOut(:,jd)<=pdfRanges(k+1));
        probModel0Data(id) = pdf0(id,k);
        %     probModel1Data(id) = pdf1(id,k);
        probModel2Data(id) = pdf2(id,k);
        likelihoodRatio(k) = sum(log10(probModel2Data(id)))-sum(log10(probModel0Data(id)));
    end
    likelihoodRatio
end
pred.probModel0Data  = probModel0Data;
pred.probModel2Data  = probModel2Data;
pred.likelihoodRatio = likelihoodRatio;
%return
return

%%%%%%%%%%%%%%
% compare pdfs at each location to each other and to data
figure
for i=find(probModel2Data./probModel0Data <1)'
    clf
    %patch(pdfRanges([1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]),...
    % pdf2(i,[1 1 2 2 3 3 4 4 5 5 6 6 7 7 1]),'r')
    patch(pdfRanges([1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]),...
        [0,pdf2(i,    [1 1 2 2 3 3 4 4 5 5 6 6 7 7]),0],'r')
    hold on
    %hp= patch(pdfRanges([1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]),...
    %pdf0(i,[1 1 2 2 3 3 4 4 5 5 6 6 7 7 1]),'k--')
    hp=patch(pdfRanges([1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]),...
        [0,pdf0(i,    [1 1 2 2 3 3 4 4 5 5 6 6 7 7]),0],'k--')
    set(hp, 'FaceAlpha', 0.5)
    plot(dataOut(i,jd),max([pdf0(i,:),pdf2(i,:)]),'bo')
    hold off
    ylim([0 1])
    xlabel('change in D_{high} (m)')
    ylabel('prob.')
    legend('update','prior','obs')
    % title(sprintf('i=%d, x=%.2f, Lratio=%.1f', i, x0(i), log(probModel2Data(i))-log(probModel0Data(i))))
    title(sprintf('i=%d, x=%.2f, Lratio=%.1f', i, x0(i), probModel2Data(i)/probModel0Data(i)))
    figure(gcf)
    pause
end
return


% include measurement error.
% there is prob(data=value)>0, all values
% measurement errors
sigma_obs = sqrt(2*(0.15)^2); % assume 15cm obs error. Error is double for difference estimates
%sigma_obs = 10
%for j=[eps, 1e-3, 1e-2,1e-1,2e-1,4e-1,8e-1,1,1e1, 1e2,1e3, 1e10] % find that 1e-1 is good
%sigma_obs =j
pdfModel0Obs = zeros(N,length(pdfRanges)-1);
pdfModel1Obs = zeros(N,length(pdfRanges)-1);
pdfModel2Obs = zeros(N,length(pdfRanges)-1);
for i=1:N
    pdfObs = diff(normcdf(pdfRanges,dhighdt(i),sigma_obs));
    pdfModel0Obs(i,:) = pdfObs'.*pdf0(i,:);
    pdfModel1Obs(i,:) = pdfObs'.*pdf1(i,:);
    pdfModel2Obs(i,:) = pdfObs'.*pdf2(i,:);
end
% LRvsLoc = nansum(log10(pdfModel2Obs)-log10(pdfModel1Obs),2);
LRvsLoc20 = log10(nansum(pdfModel2Obs,2))-log10(nansum(pdfModel0Obs,2));
nansum(LRvsLoc20)
LRvsLoc21 = log10(nansum(pdfModel2Obs,2))-log10(nansum(pdfModel1Obs,2));
nansum(LRvsLoc21)
LRvsLoc10 = log10(nansum(pdfModel1Obs,2))-log10(nansum(pdfModel0Obs,2));
nansum(LRvsLoc10)
%end

% compare improvements
plot(x0,LRvsLoc20,'r', x0,LRvsLoc10,'b', x0,LRvsLoc21,'g')

%%%%%%%%%%%%%%
% now, look at results
% compare pdfs at each location to each other and to data
figure
for i=find(LRvsLoc20>1)'
    clf
    hp=patch(pdfRanges([1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]),...
        [0,pdf2(i,    [1 1 2 2 3 3 4 4 5 5 6 6 7 7]),0],'r')
    set(hp,'LineStyle','none')
    hold on
    hp=patch(pdfRanges([1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]),...
        [-1,pdf0(i,    [1 1 2 2 3 3 4 4 5 5 6 6 7 7]),-1],'k--')
    set(hp, 'FaceAlpha', 1,'LineStyle','--', 'faceColor','none')

    plot(dhighdt(i),max([pdf0(i,:),pdf2(i,:)]),'bo')
    hold off
    ylim([0 1])
    xlabel('change in D_{high} (m)')
    ylabel('prob.')
    legend('update','prior','obs')
    title(sprintf('i=%d, x=%.2f, Lratio20=%.1f', i, x0(i),LRvsLoc20(i)))
    figure(gcf)
    pause
end

% compare
figure
patch(pdfRanges([1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]),...
    likelihoodRatio([1 1 2 2 3 3 4 4 5 5 6 6 7 7 1]), 'b')
xlabel('change D_{high}')
ylabel('likelihood ratio')

% explain
figure
subplot 211
plot(x0,d0,'b', x0,dhighdt,'r')
ylabel('elevation (m)')
ylim([-5 6])
legend('initial D_{high}', 'change  D_{high}')
xlim([min(x0) max(x0)])
xlim([29 31])
grid on
subplot 212
plot(x0,probData2,'b.',  x0,probData,'r')
ylabel('likelihood')
ylim([0 1])
legend('no forcing', 'use forcing')
grid on
xlim([29 31])

% where do we do well/poorly?
figure
plot(dhighdt,probData2,'b+', dhighdt,probData,'r.')
xlabel('change D_{high} (m)')
ylabel('prediction likelihood (skill)')
legend('prior', 'updated')


figure
plot(probData,probData2,'b+',[0 1],[0 1],'k')
grid on
xlabel('updated likelihood')
ylabel('prior likelihood')

% save ICP_BB_v8.mat



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save for Kara
% package the data
IvanData = struct;
%data0 = load('ivanData/fl02_fl04.0405r2.dhigh.mat')% (pre-storm, CHARTS, combined dates)
%data1 = load('ivanData/fl02_fl04.040919.dhigh.mat')% (post-storm, EAARL, Sep 19 2004)
IvanData.x0 = data0.dhigh(:,6); % easting (m)
IvanData.x1 = data1.dhigh(:,6);
IvanData.d0 = data0.dhigh(:,9); % use smoothed dhigh(9)
IvanData.d1 = data1.dhigh(:,9);
IvanData.dhighdt = d1-d0; % Dhigh difference - predict this!
% dataR = load('ivanData/rhigh.ivan.prof.mat');
IvanData.rhigh = max(dataR.Rdata(:,9,:), [],3); % 2% runup + surge + setup
IvanData.rlow = max(dataR.Rdata(:,7,:) + dataR.Rdata(:,8,:), [],3); % surge + setup

% package the prediction
IvanPrediction = struct;
IvanPrediction.pdfRanges = pdfRanges; % the dhighdt prediction bin boundaries
% 8x1 values (meters, negative is erosion):  [    -8    -6    -4    -2    -1     0     2     8]
IvanPrediction.pdf0 = pdf0; % the pdfs prior to Bayesian update
% Nx7 probabilities:     [0.0134    0.0207    0.0879    0.2606    0.4052    0.1899   0.0223]
IvanPrediction.pdf1 = pdf1; % the pdfs afterBayesian update of initial Dhigh
IvanPrediction.pdf2 = pdf2; % the pdfs afterBayesian update of initial Dhigh AND Rlow,Rhigh
IvanPrediction.LikelihoodRatio20 = LRvsLoc20; % log-likelihood of data given fully updated model compared to prior model
% Nx1 ratios: L>0 means updated model is better.
% L>1 means updated model is 10 times better than the prior model
IvanPrediction.LikelihoodRatio10 = LRvsLoc10; % log-likelihood of data given partially updated model (update Dhigh only) compared to prior model

% save to disk
%% save BayesianIvanPrediction.mat IvanData IvanPrediction

% check
load BayesianIvanPrediction.mat
% compare pdfs at each location to each other and to data
figure
% select a likelihood threshold to inspect
Lthresh = -0.5; % Lthresh<0 looks at bad, Lthresh>0 looks at good, Lthresh>>0 looks at really good,...
for i=find((sign(Lthresh)*IvanPrediction.LikelihoodRatio20)>(sign(Lthresh)*Lthresh))'
    clf
    % show updated histogram
    hp=patch(IvanPrediction.pdfRanges([1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]),...
        [0,IvanPrediction.pdf2(i,    [1 1 2 2 3 3 4 4 5 5 6 6 7 7]),0],'r')
    set(hp,'LineStyle','none')
    hold on
    hp=patch(IvanPrediction.pdfRanges([1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]),...
        [-1,IvanPrediction.pdf0(i,    [1 1 2 2 3 3 4 4 5 5 6 6 7 7]),-1],'k--')
    set(hp, 'FaceAlpha', 1,'LineStyle','--', 'faceColor','none')

    plot(IvanData.dhighdt(i),max([IvanPrediction.pdf0(i,:),IvanPrediction.pdf2(i,:)]),'bo')
    hold off
    ylim([0 1])
    xlabel('change in D_{high} (m)')
    ylabel('prob.')
    legend('update','prior','obs')
    title(sprintf('i=%d, x=%.2f, Lratio20=%.1f', i, IvanData.x0(i),IvanPrediction.LikelihoodRatio20(i)))
    figure(gcf)
    pause
end

