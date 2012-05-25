function netStatus = NetStartStop(netName, neticapw)
% netStatus = NetStartStop(netName)
%
% function to handle the creation and management of the java neticaJ package
%
% Input
%    netName, the name of the *.dne file to run
%       IF netName='stop', clear the net
%    (optional) neticapw, password for licensed version
%
% Output
%   netStatus, returns nothing if ok, else an error or warning message
%   globals: 
%       NETICA_ENV, object to kill the neticaJ session
%       NETICA_STREAMER, used to get access to the net 
%       NETICA_NETNAME, what a net calls itself 
%       NETICA_NET, net object, used to compile a particular net and get access to functions
%           - use this to retract finding, 
%       NETICA_NODESTRUCT, my structure to store info about each nodes, like its pdf ranges, its name, etc.
%       NETICA_ALL_NODES, list of node-objects, used to manipulate  each node (like getLikelihod)
%
% implements the java version of netica on the NETICA network
% can we do it all from matlab? YES
% classpath.txt must have, e.g., 
% C:\Users\Plant\Projects\...\NeticaJ_221\bin
% C:\Users\Plant\Projects\...\NeticaJ_221\bin\NeticaJ.jar
% YOU MUST START MATLAB IN THE WORKING DIRECTORY!!!
%
%

% net's java objects are available via this global 
global NETICA_ENV NETICA_NET NETICA_NETNAME NETICA_NODESTRUCT NETICA_STREAMER NETICA_ALL_NODES

% import netica stuff
import norsys.netica.*;
%import norsys.neticaEx.aliases.Node;

% if no arguments, just return status
if(nargin==0)
    if(isempty(NETICA_ENV))
        netStatus = sprintf('no net started yet');
    else
        netStatus = sprintf('using %s\n', char(NETICA_STREAMER.getFileName));
    end
    return;
end    

% check input for stop command
if(strcmp('stop', lower(netName)) | strcmp('finalize', lower(netName)))
    if(isempty(NETICA_ENV))
        netStatus = sprintf('no net started yet');
    else
        % stop it, let errors tell us what went wrong
        NETICA_ENV.finalize;
        netStatus = '';
        clear global NETICA_ENV NETICA_NET NETICA_NETNAME NETICA_NODESTRUCT NETICA_STREAMER NETICA_ALL_NODES
    end
    return
end

% process if we can find the net
if(exist(netName)==2)
    netStatus = sprintf('found %s, trying to start...\n', netName);
    % start session if not already started
    if(~isempty(NETICA_ENV))
        netStatus = sprintf('found existing net %s, trying to stop it...', NETICA_NETNAME);
        NETICA_ENV.finalize;
    end
    
    % lets init the environment with a method
    % see if one is out there
    %NETICA_ENV = Environ.getDefaultEnviron;
    % now set password
    NETICA_ENV = Environ('+FienenM/USGS/120,310-6/27866');
    %if(nargin==2)
    %    NETICA_ENV = Environ(neticapw);
    %else
    %    NETICA_ENV = Environ('');
    %end
    % kill by 
    %%  NETICA_ENV.finalize
    
    % get a net
    % NETICA_STREAMER = Streamer('impact.dne')
    NETICA_STREAMER = Streamer(netName);
    NETICA_NET = Net(NETICA_STREAMER);
    NETICA_NETNAME = char(NETICA_STREAMER.getFileName);
    
    % compile it so it can do some analysis
    NETICA_NET.compile;
    % NOTE: kill by 
    %%  NETICA_NET.finalize
    
    % ok, now, get the parts of each node that we need and orgnanize in a structure
    % get the nodes, as an array of objects
    NETICA_ALL_NODES = NETICA_NET.getNodes.toArray;
    
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
    
else
    netStatus = sprintf('could not find netName: %s\n', netName);
end

return;