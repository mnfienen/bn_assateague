function pred = predictBayes(netName, nodeNamesIn, nodeNamesOut, dataIn, dataOut)
% pred = predictBayes(netName, nodeNamesIn, nodeNamesOut, dataIn, dataOut)
% input
%
% output
%

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
    j = strmatch(nodeNamesOut(i), nodenames,'exact');
    if(~isempty(j))
        predictID(i)=j;
    end
end
% what will we constrain
inputID = zeros(length(nodeNamesIn),1);
for i=1:length(nodeNamesIn)
    j = strmatch(nodeNamesIn(i), nodenames,'exact');
    if(~isempty(j))
        inputID(i)=j;
    end
end

% initialize the output
[N,m] = size(dataIn);
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
    rangestruct = setfield(rangestruct, nodenames{j}, pdfLevels);
    %    make sure data in bounds
    dataid = strmatch(nodenames{j},nodeNamesIn, 'exact');
    if(~isempty(dataid))
        id = find(dataIn(:,dataid)>max(pdfLevels));
        dataIn(id,dataid) = max(pdfLevels);
        id = find(dataIn(:,dataid)<min(pdfLevels));
        dataIn(id,dataid) = min(pdfLevels);
    end
end

% get prior (before data - same for all data!)
for i=1:length(predictID)
    j = predictID(i); % j is index into net arrangement, i is index into data
    pdfstruct = setfield(pdfstruct, ['prior_',nodenames{j}], zeros(1,Nbins(j)));
    outputPDF = NETICA_ALL_NODES(j).getBeliefs;
    pdf = getfield(pdfstruct, ['prior_',nodenames{j}]);
    pdf(1,:) = outputPDF(:)';
    pdfstruct = setfield(pdfstruct, ['prior_',nodenames{j}], pdf);
end

% loop through each input series and get netica prediction
fprintf('doing %d locations\n',N)
tic; 
for i=1:N
    if(toc>10)
        % provide feeback if running long
        fprintf('doing %d of %d\n', i,N)
        tic
    end
    % assign inputs
    IDnum = i;
    
    % clear the net
    NETICA_NET.retractFindings
    
    % loop throught the nodenames, setting input
    dataid=0;
    for j = inputID(:)'
        % get data for this node
        dataid = dataid+1;
        newData = dataIn(i,dataid);
        %fprintf('i=%d,j=%d,dataid=%d,%s=%.3f', i,j,dataid,nodenames{j},newData)
        % WARNING-treating data as perfect
        if(~isnan(newData))
            % enter data
            NETICA_ALL_NODES(j).enterValue(newData);
            % get update for this variable
            inputPDF = NETICA_ALL_NODES(j).getBeliefs;
            % store
            pdf = getfield(pdfstruct, nodenames{j});
            pdf(i,:) = inputPDF(:)';
            pdfstruct = setfield(pdfstruct, nodenames{j}, pdf);
            
            %             NETICA_ALL_NODES(j)
            %             keyboard
            %
        end
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

%% below is post processing
% plot comparison of prediction to observation
jcnt=0;
for j = predictID'
    jcnt=jcnt+1;
    pdfRanges = getfield(rangestruct, nodenames{j});
    pdfRangePlot = pdfRanges(2:end)-0.5*diff(pdfRanges);
    
    %%%%%%%%%
    % prediction prior to updating (the mean for all situations)
    pdf0 = getfield(pdfstruct, ['prior_', nodenames{j}]);
    pdf0 = repmat(pdf0,N,1);
    pdf2 = getfield(pdfstruct, [nodenames{j}]);
%    pdf = pdf2;
    %%%%%%%%%%%%
    %get prob data.
    %     jd = find(j==predictID);
    %     probModel0Data = zeros(N,1);
    %     probModel1Data = zeros(N,1);
    %     probModel2Data = zeros(N,1);
    %     binLikelihoodRatio = zeros(length(pdfRangePlot),1);
    %     for k=1:length(pdfRangePlot)
    %         id = find(dataOut(:,jd)>pdfRanges(k) & dataOut(:,jd)<=pdfRanges(k+1));
    %         probModel0Data(id) = pdf0(id,k);
    %         probModel2Data(id) = pdf2(id,k);
    %         binLikelihoodRatio(k) = sum(log10(probModel2Data(id)+eps))-sum(log10(probModel0Data(id)+eps));
    %         pdfData(id,k) = 1;
    %     end
    % end
    % avoid trajedy when the observation lies on the boundary between states:
    % distribute that point over two states.  This only affects the prob(data),
    % which is compared to prior and update pdfs identically
    pdfData  = makeInputPdf(pdfRanges, 'norm',[dataOut(:,jcnt),repmat(eps,N,1)]);
    probModelUpdate = nansum(pdfData.*pdf2,2);
    probModelPrior = nansum(pdfData.*pdf0,2);
    logLikelihoodRatio = log10(probModelUpdate+eps)-log10(probModelPrior+eps);
    % add to the output
    pred.probModelPrior(:,jcnt)  = probModelPrior;
    pred.probModelUpdate(:,jcnt)  = probModelUpdate;
    pred.logLikelihoodRatio(:,jcnt) = logLikelihoodRatio;
    pred.pdf.(['data_',nodenames{j}]) = pdfData;
end
return

