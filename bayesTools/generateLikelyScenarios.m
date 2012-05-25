%% generate a bunch of boundary conditions for Delft3D at Cape Canaveral
addpath c:\svn\bayesTools

% use the net to get what we need
% netName = 'BayesWave_loaded_buoy_em_wis-wind-wavedir_in+r2p+ds+new.neta' % canaveral
netName = 'C:\SVN\BayesWaveTide_GulfCoast\gc040_007_Pilot_10years_buoyOnly.neta' % chandeleurs
[netstuff] = netInfo(netName)
nodeNames = cellstr(char(netstuff.name))

% use "dataNames" match file to nodes--this lines up each column
% all of them dataNames = {    'tide'    'windSpd'    'windDir'    'surge'    'dsurge'    'WaveDir'    'wavePer'    'waveHgt'    'slope'    'R2'
% order dataNames = { 'R2'  'slope' 'waveHgt'  'wavePer' 'WaveDir'  'windSpd'   'windDir'  'surge' 'dsurge' 'tide'
% use as controls
% dataNames = { 'R2'  'slope' 'waveHgt'  'wavePer' 'WaveDir'  'windSpd'   'windDir' }
%     'winddir1'
%     'windspd1'
%     'wavehgt1'
%     'waveper1'
%     'wavehgt2'
%     'waveper2'
%     'tide'
%     'wavedir1'
%     'wavedir2'

%% list the input and output nodes
defaultInputNodes = { 'wavehgt1'   }
defaultOutputNodes = {'wavehgt1','waveper1'  'wavedir1'  'windspd1' 'winddir1'  'wavehgt2' 'waveper2'  'wavedir2'  }
defaultTestCase = {defaultInputNodes, defaultOutputNodes}
defaultInputData = zeros(size(defaultInputNodes))
defaultOutputData = zeros(size(defaultOutputNodes))
% make a default prediction to use as a structure for populating
defaultPred = predictBayes(netName, defaultInputNodes, defaultOutputNodes, defaultInputData, defaultOutputData)


% %% populate with priors
% for j=1:length(nodeNames)
%     switch nodeNames{j}
%         case 'R2'
%             % update runup, want all cases with runup>threshhold
%             runupThresh = 1;
%             pdf = getfield(defaultPred.pdf, 'R2');
%             pdfRanges = getfield(defaultPred.ranges, 'R2');
%             id = find(pdfRanges(2:end)>runupThresh);
%             pdf = 0*pdf; % clear probs
%             pdf(id)=1;% set to uniform
%             defaultPred.pdf = setfield(defaultPred.pdf, 'R2',pdf); % replace
%
%         case 'slope'
%             % insert actual slope distributions
%             meanSlope = 0.05; % use real stats
%             stdSlope = 0.05;
%             pdf = getfield(defaultPred.pdf, 'slope');
%             pdfRanges = getfield(defaultPred.ranges, 'slope');
%             pdf = makeInputPdf(pdfRanges, pdf, 'norm',[meanSlope,stdSlope])
%             defaultPred.pdf = setfield(defaultPred.pdf, 'slope',pdf); % replace
%         otherwise
%             % otherwise, insert uncertain
%             nodeName = nodeNames{j};
%             pdf = getfield(defaultPred.pdf, nodeName);
%             pdf = 0*pdf +1;
%             defaultPred.pdf = setfield(defaultPred.pdf, nodeName,pdf); % replace
%     end
% end

%% store these for generating output
parameterNames = defaultOutputNodes; %{'waveHgt'  'wavePer' 'waveDir'  'windSpd'   'windDir'}
paramList = struct;
paramList.value=struct;
paramList.prob = struct;
for j=1:length(parameterNames)
    paramList.value =  setfield(paramList.value,parameterNames{j}, []);
    paramList.prob =  setfield(paramList.prob,parameterNames{j}, []);
end

%% generate the driving scenarios: wave height
waveHgtTests = 0.5:1:10;

% % initialize
% nodeNamesIn = defaultInputNodes;
% nodeNamesOut = defaultOutputNodes(1);
% % get needed input
% newPred = struct;
% for i=1:length(nodeNamesIn)
%     % use dynamic struct
%     newPred =  setfield(newPred,nodeNamesIn{i}, defaultPred.pdf.(nodeNamesIn{i}));
% end
%
% defaultOutputNodes = {'wavehgt2' 'waveper2'  'wavedir2'   'windspd1' 'winddir1' }
%


%% step through, updating, contstraining, updating...
cnt = 0;
dataOut = zeros(1,length(nodeNames)); % place holder
for i=1:length(waveHgtTests)
    % update wave height, get waveper
    nodeName1 = 'wavehgt1'; % initial constraint
    waveHgt = waveHgtTests(i);
    nodeNamesIn = {'wavehgt1'};
    nodeNamesOut = {'wavehgt1', 'waveper1'};% get the waveheight prediction too, so we get the prior
    % initialize
    %     pdf1 = getfield(newPred, nodeName1);
    rng1 = getfield(defaultPred.ranges, nodeName1);
    pdf1 = makeInputPdf(rng1, 'norm',[waveHgt,eps],'continuous');
    newPred = struct;
    newPred = setfield(newPred, nodeName1,pdf1); % replace
    
    % update
    pred1 = predictBayes_v2dev(netName, nodeNamesIn, nodeNamesOut, newPred, dataOut);
    % get the wave height prior prob
    pdf1 = getfield(pred1.pdf, nodeName1);
    pdf1p = getfield(pred1.pdf, ['prior_',nodeName1]);
    rng1 = getfield(pred1.ranges, nodeName1); rng1 = (rng1(1:end-1)+rng1(2:end))/2;
    % get the bin for the constrained variable
    [pdfsort,sortid1] = sort(-pdf1);
    %     % store
    %     cnt = cnt +1;
    %     paramList(cnt).value.(nodeNamesIn(1)) = rng1(sortid1(1)); % get most probable
    %     paramList(cnt).prob.(nodeNamesIn(1)) = pdf1p(sortid1(1)); % get most probable (prior)
    
    %% pass to wave period
    nodeName2 = 'waveper1';
    pdf2 = getfield(pred1.pdf, nodeName2);
    rng2 = getfield(pred1.ranges, nodeName2); rng2 = (rng2(1:end-1)+rng2(2:end))/2;
    % get top 3 values
    [pdfsort,sortid2] = sort(-pdf2);
    for j2=1:3
        % update period
        wavePer = rng2(sortid2(j2))
        pdf = 0*pdf2;
        pdf(sortid2(j2))=1;
        newPred = setfield(newPred, nodeName2,pdf); % replace
        % do it
        nodeNamesIn = {'wavehgt1','waveper1'};
        nodeNamesOut = {'wavedir1'};
        pred2 = predictBayes_v2dev(netName, nodeNamesIn, nodeNamesOut, newPred, dataOut);
        
        %% pass to wave direction
        nodeName3 = 'wavedir1';
        pdf3 = getfield(pred2.pdf, nodeName3);
        rng3 = getfield(pred2.ranges, nodeName3); rng3 = (rng3(1:end-1)+rng3(2:end))/2;
        [pdfsort,sortid3] = sort(-pdf3);
        for j3=1:3
            % update wave direction
            waveDir = rng3(sortid3(j3))
            pdf = 0*pdf3;
            pdf(sortid3(j3))=1;
            newPred = setfield(newPred, nodeName3,pdf); % replace
            nodeNamesIn = {'wavehgt1','waveper1','wavedir1'};
            nodeNamesOut = {'windspd1'};
            pred3 = predictBayes_v2dev(netName, nodeNamesIn, nodeNamesOut, newPred, dataOut);
            
            %% pass to wind speed
            nodeName4 = 'windspd1';
            pdf4 = getfield(pred3.pdf, nodeName4);
            rng4 = getfield(pred3.ranges, nodeName4); rng4 = (rng4(1:end-1)+rng4(2:end))/2;
            [pdfsort,sortid4] = sort(-pdf4);
            for j4=1:3
                % update wind speed
                windSpd = rng4(sortid4(j4))
                pdf = 0*pdf4;
                pdf(sortid4(j4))=1;
                newPred = setfield(newPred, nodeName4,pdf); % replace
                nodeNamesIn = {'wavehgt1','waveper1','wavedir1','windspd1'};
                nodeNamesOut = {'winddir1' };
                pred4 = predictBayes_v2dev(netName, nodeNamesIn, nodeNamesOut, newPred, dataOut);
                
                %% pass to wind dir
                nodeName5 = 'winddir1';
                pdf5 = getfield(pred4.pdf, nodeName5);
                rng5 = getfield(pred4.ranges, nodeName5); rng5 = (rng5(1:end-1)+rng5(2:end))/2;
                [pdfsort,sortid5] = sort(-pdf5);
                for j5=1:3
                    % update wind dir
                    windDir = rng5(sortid5(j5))
                    pdf = 0*pdf5;
                    pdf(sortid5(j5))=1;
                    newPred = setfield(newPred, nodeName5,pdf); % replace
                    nodeNamesIn = {'wavehgt1','waveper1','wavedir1','windspd1','winddir1'};
                    nodeNamesOut = {'wavehgt1','waveper1','wavedir1','windspd1','winddir1', 'wavehgt2', 'waveper2',  'wavedir2'};
                    pred5 = predictBayes_v2dev(netName, nodeNamesIn, nodeNamesOut, newPred, dataOut);
                    
                    %% offshore is fully specified, just collect them, we are done
                    cnt = cnt + 1;
                    paramList(cnt).rank = [1, j2, j3, j4, j5, 1]';
                    for k=1:length(nodeNamesOut)
                        pdf6 = getfield(pred5.pdf, nodeNamesOut{k});
                        rng6 = getfield(pred5.ranges, nodeNamesOut{k}); rng6 = (rng6(1:end-1)+rng6(2:end))/2;
                        [pdfsort,sortid6] = sort(-pdf6);
                        paramList(cnt).value.(nodeNamesOut{k}) = rng6(sortid6(1)); % get most probable
                        paramList(cnt).prob.(nodeNamesOut{k}) = pdf6(sortid6(1)); % get most probable
                    end
                    paramList(cnt).value
                    paramList(cnt).prob
                    paramList(cnt).rank % shows whether we pick most likely or not
                    pause(0.1);
                end
                % remove this one from the trial
                newPred = rmfield(newPred, nodeName4);
            end
            newPred = rmfield(newPred, nodeName3);
        end
        newPred = rmfield(newPred, nodeName2);
    end
end

% save BermWaveScenarios_v1.mat paramList netName nodeNames




