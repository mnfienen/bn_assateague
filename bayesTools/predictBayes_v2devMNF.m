function pred = predictBayes_v2devMNF(netName, nodeNamesIn, nodeNamesOut, dataIn, dataOut, dataInErr)
% pred = predictBayes(netName, nodeNamesIn, nodeNamesOut, dataIn, dataOut)
%% input
%   netName, a netica network (dne, neta)
%   nodeNamesIn, Mx1 cell array of strings listing node names to use as input
%   nodeNamesOut, Px1 cell array of strings listing node names to use as output
%   dataIn, NxM array of input values, collumns aranged in order specified by nodeNamesIn
%   (optional) dataOut, NxP array of verification data corresponding ot nodeNamesOut
%   (optional) dataInErr, NxM array of data value standard deviations (must use pdf if more complicated error model is used)
%
%% output
%   pred, structure with fields
%       .netName, the name of the net
%       .pdf, another structure with one NxR entry per node
%       .ranges, another structure with one Rx1 entry listing bins
%       . lots of other stuff
%
% 
%%% REPLACING NetStartStop for now with explicit starting and stopping 

% % start the network
% netStatus = NetStartStop(netName)
% % all is well if status is empty
% if(isempty(netStatus))
%     % get to the net
%     global NETICA_ENV NETICA_NET NETICA_NETNAME NETICA_NODESTRUCT NETICA_STREAMER NETICA_ALL_NODES
% else
%     error(netStatus)
% end
global NETICA_ENV NETICA_NET NETICA_NETNAME NETICA_NODESTRUCT NETICA_STREAMER NETICA_ALL_NODES

% import netica stuff
import norsys.netica.*;

% fire up the netica API with correct PWD
ifp = fopen('mikeppwd.txt','r');
PWD = fscanf(ifp,'%s');
NETICA_ENV = Environ(PWD);
NETICA_STREAMER = Streamer(netName);
NETICA_NET = Net(NETICA_STREAMER);
NETICA_NETNAME = char(NETICA_STREAMER.getFileName);
NETICA_ALL_NODES = NETICA_NET.getNodes.toArray;

% need to compile the net
NETICA_NET.compile;

% how many?
Nnodes = length(NETICA_ALL_NODES);

% suck info out
% store info in matlab struct
NETICA_NODESTRUCT = struct; 
for i=1:Nnodes
    % get its name (short) and title (human readable)
    NETICA_NODESTRUCT(i).name = char(NETICA_ALL_NODES(i).toString);
    NETICA_NODESTRUCT(i).title = char(NETICA_ALL_NODES(i).getTitle);
    % NETICA_ALL_NODES(i).getNumStates

    % get numeric values of the states
    NETICA_NODESTRUCT(i).levels = NETICA_ALL_NODES(i).getLevels;

    % get its beliefs
    NETICA_NODESTRUCT(i).beliefs = double(NETICA_ALL_NODES(i).getBeliefs);
    NETICA_NODESTRUCT(i).Nbeliefs = length(NETICA_NODESTRUCT(i).beliefs);

    % these beliefs correspond to discrete states
    NETICA_NODESTRUCT(i).state = struct;
    for j=1:NETICA_NODESTRUCT(i).Nbeliefs
        NETICA_NODESTRUCT(i).state(j).obj = NETICA_ALL_NODES(i).state(j-1); % increment from zero
        NETICA_NODESTRUCT(i).state(j).name = char(NETICA_NODESTRUCT(i).state(j).obj.getName);
        NETICA_NODESTRUCT(i).state(j).numeric = NETICA_NODESTRUCT(i).state(j).obj.getNumeric;
    end

    % get current liklihoods
    NETICA_NODESTRUCT(i).liklihood = double(NETICA_ALL_NODES(i).getLikelihood);
end

% what have we got 
char(NETICA_NODESTRUCT.title);

% status is ok
netStatus = [];



%%% back to Nathaniel's code now


% get info
% display titles
nodetitles = {NETICA_NODESTRUCT.title}'
% display nodenames
netNodeNames = {NETICA_NODESTRUCT.name}'
% other info
%inputData(i).nodeID = strmatch(requiredInputNames{i,2}, char(NETICA_NODESTRUCT.title)); % id to talk to bayesian net
%inputData(i).pdfRanges = [NETICA_NODESTRUCT(inputData(i).nodeID).levels]; % range of values in the bayesian net
%inputData(i).pdfStates = {NETICA_NODESTRUCT(inputData(i).nodeID).state.name};
% get bins, either number of pdfRanges -1 (continuous) or number of pdfRanges (discrete)
%inputData(i).Nbins = NETICA_NODESTRUCT(inputData(i).nodeID).Nbeliefs; %length(inputData(i).pdfRanges)-1; % number of bins between these ranges
%inputData(i).pdf = zeros(1,inputData(i).Nbins);

% check input 
if(isstruct(dataIn))
    % structure input means we have pdfs
    [N,m] = size(dataIn.(nodeNamesIn{1}));
    m = length(fieldnames(dataIn));
    fprintf('input is pdf with fields\n')
    fieldnames(dataIn)
else
    % we have a data matrix
    [N,m] = size(dataIn);
    if(m~=length(nodeNamesIn))
        error('dataIn must have same number of columns as nodeNameIn');
    end
end

% initialize output
pred = struct;
pred.netName = netName;
pred.pdf = struct;
pred.ranges = struct;
NETICA_NET.retractFindings
for j=1:length(netNodeNames)
    % intialize output pdf
    Nbins = NETICA_NODESTRUCT(j).Nbeliefs;
    pred.pdf.(netNodeNames{j}) =  zeros(N,Nbins);
    % initialize pdf ranges
    pred.ranges.(netNodeNames{j}) = [NETICA_NODESTRUCT(j).levels];
    % get the plotable version
    if(Nbins<length(pred.ranges.(netNodeNames{j})))
        % continuous, bin centers
        pred.rangesp.(netNodeNames{j}) = pred.ranges.(netNodeNames{j})(2:end)-0.5*diff(pred.ranges.(netNodeNames{j}));
    else
        % discrete, bin value
        pred.rangesp.(netNodeNames{j}) = pred.ranges.(netNodeNames{j});
    end
    % get priors
    % enforce as a row vector
    pdfp = NETICA_ALL_NODES(j).getBeliefs;
    pred.pdf.(['prior_',netNodeNames{j}]) = pdfp(:)';
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
    % clear the net
    NETICA_NET.retractFindings
    
    % loop through the INPUT nodes, setting input values or likelihoods
    for j = 1:length(nodeNamesIn)
        % and where is this in the net's node list?
        k = strmatch(nodeNamesIn(j), netNodeNames,'exact');
        if(isempty(k))
            fprintf('network does not contain a node with nodeNamesIn=%s\n', nodeNamesIn{j})
            fprintf('available netNodeNames are:\n')
            netNodeNames;
            error('nodeNamesIn dont match netNodeNames');
        end
        % need to know if continuous or discrete
        if(length(pred.ranges.(nodeNamesIn{j}))==length(pred.pdf.(nodeNamesIn{j})(i,:)));
            varType = 'discrete';
        else
            varType = 'continuous';
        end
        % get data for this node
        if(isstruct(dataIn))
            %% probabilities are input
            newPDF = dataIn.(nodeNamesIn{j})(i,:); % dynamic field name
            % save what we used as input (it will change with updating)
            pred.pdfIn.(nodeNamesIn{j})(i,:)=newPDF;
            % enter into the network
            if(~(any(isnan(newPDF))))
                % enter data
                NETICA_ALL_NODES(k).enterLikelihood(newPDF);
            end
        elseif(exist('dataInErr'))
            % data values and (optional) errors are input
            % convert data to likelihood using error model
            newPDF  = makeInputPdf(pred.ranges.(nodeNamesIn{j}), 'norm',[dataIn(i,j),eps+dataInErr(i,j)],varType);
            pred.pdfIn.(nodeNamesIn{j})(i,:)=newPDF;
            % enter into the network
            if(~(any(isnan(newPDF))) & max(newPDF)~=1)
                % enter data as a likelihood, since prob occurs in more than one bin
                try
                NETICA_ALL_NODES(k).enterLikelihood(newPDF);
                catch
                    lasterr
                    keyboard
                end
            elseif(~(any(isnan(newPDF))) & max(newPDF)==1)
                % this is same as perfect data, update valule (faster)
                %% BUT FOR NOW, DO IT THE SLOW WAY
                NETICA_ALL_NODES(k).enterLikelihood(newPDF);
%                 NETICA_ALL_NODES(k).enterValue(dataIn(i,j));
            else
                % nan, skip it
            end
        else
            % only data values are input, literally (notused)
            if(~isnan(dataIn(i,j)))
                % enter data
                % catch inconsistent
                try
                NETICA_ALL_NODES(k).enterValue(dataIn(i,j));
                catch
                    % skip this input but report error
                    fprintf('record %d: %s\n', i,lasterr)
                    %keyboard
                end
            end
        end
    end
    % end of updating
    
    % now, get all the beliefs 
    for j = 1:length(netNodeNames)
        try
        pred.pdf.(netNodeNames{j})(i,:) = NETICA_ALL_NODES(j).getBeliefs;
        catch
            lasterr
            keyboard
        end
    end
end
% done!

%% post processing
for j = 1:length(nodeNamesOut)
    pdfRanges = pred.ranges.(nodeNamesOut{j});
    pdfRangePlot = pdfRanges(2:end)-0.5*diff(pdfRanges);
    % prediction prior to updating (the mean for all situations)
    pdf0 = pred.pdf.(['prior_', nodeNamesOut{j}]);
    pdf0 = repmat(pdf0,N,1);
    pdf1 = pred.pdf.(nodeNamesOut{j});
    if(length(pdfRanges)==length(pdf0(1,:)));
        varType = 'discrete';
    else
        varType = 'continuous';
    end
    
    % parse the data
    if(isstruct(dataOut))
        % warning('cant handle pdf on output data yet')
        pdfData = dataOut.(nodeNamesOut{j});
    else
        % avoid trajedy when the observation lies on the boundary between states:
        % distribute that point over two states.  This only affects the prob(data),
        % which is compared to prior and update pdfs identically
        pdfData  = makeInputPdf(pdfRanges, 'norm',[dataOut(:,j),repmat(eps,N,1)], varType);
    end
    probModelUpdate = nansum(pdfData.*pdf1,2);
    probModelPrior = nansum(pdfData.*pdf0,2);
    logLikelihoodRatio = log10(probModelUpdate+eps)-log10(probModelPrior+eps);
    % add to the output
    pred.probModelPrior(:,j)  = probModelPrior;
    pred.probModelUpdate(:,j)  = probModelUpdate;
    pred.logLikelihoodRatio(:,j) = logLikelihoodRatio;
    pred.pdf.(['data_',nodeNamesOut{j}]) = pdfData;
end

% shut down the NETICA session here
NETICA_ENV.finalize

