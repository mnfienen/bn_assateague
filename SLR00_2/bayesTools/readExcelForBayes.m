function [data,dataNames] = readExcelForBayes(netstuff, datafilename, textNames)
% [data,dataNames] = readExcelForBayes(netstuff, datafilename, textNames)
%
% Input
% netstuff, info returned by netInfo:
%   netstuff = 1x10 struct array with fields:
%     name
%     title
%     levels
%     beliefs
%     Nbeliefs
%     state
%     liklihood
% datafilename, and excel file
% textNames, cell of names that will have text that need to be processed to
% numeric states

% read the file
[a,b,c] = xlsread(datafilename);

% prepare to match nodenames to titles
N = length(a(:,1)); % number of data
dataTitles = b(1,:); % text data
% initialize and fill
dataNames = cell(1,length(netstuff));
data = nan*ones(N,1);
for i=1:length(netstuff)
    % deal with text input
    fprintf('getting data from %s\n', netstuff(i).name)
    switch netstuff(i).name
        % remap names as needed
        case  textNames
            fprintf('dealing with text input for %s\n', netstuff(i).name)
            id = strmatch(netstuff(i).name, char(dataTitles)); % which column in the data file
            stateNames = char(netstuff(i).state.name); % get text description
            stateNums = [(netstuff(i).state.numeric)]; % get numeric
            % convert to state value (index)
            tmp1 = c(2:N+1,id);
            % find these
            unique_data = unique(cellstr(tmp1))
            for j=1:length(stateNames(:,1))
                idd = strmatch( stateNames(j,:), tmp1);
                data(idd, i) = stateNums(j); % assign numeric
            end
            dataNames{i} = netstuff(i).name;
        otherwise
            % get from cell
            id = strmatch(netstuff(i).name, char(dataTitles));
            if(length(id)==1)
                dataNames{i} = netstuff(i).name;
                data(:,i) = [c{2:N+1,id}]';
            else
                error('no match to node name %s\n', netstuff(i).name)
            end
    end
end

