function compile_net(NetName,newCaseFilename,voodooPar,outfilename)
% compile_net(NetName,newCaseFilename,voodoPar,outfilename)
% a m!ke@usgs joint <mnfienen@usgs.gov>
% function to build the CPT tables for a new CAS file on an existing NET
% (be existing, meaning that the nodes, edges, and bins are dialed)
% INPUT:
%       NetName --> a filename, including '.neta' extension
%       newCaseFilename --> new case file including '.cas' extension
%       voodoPar --> the voodoo tuning parameter for building CPTs
%       outfilename --> netica file for newly build net (including '.neta'
global NETICA_ENV NETICA_NET
% import netica stuff
import norsys.netica.*;

% fire up the netica API with correct PWD
ifp = fopen('mikeppwd.txt','r');
PWD = fscanf(ifp,'%s');
NETICA_ENV = Environ(PWD);

% open the net
NETICA_STREAMER = Streamer(NetName);
NETICA_NET = Net(NETICA_STREAMER);

% not sure if this compilation step is really necessary . . . 
NETICA_NET.compile;

NETICA_ALL_NODES = NETICA_NET.getNodes;
numNODES = length(NETICA_ALL_NODES);

for i = 1:numNODES
    cnode = NETICA_ALL_NODES.get(i);
    cnode.deleteTables;
end
% connect to the case file containing the new data
newCaseFile = Streamer(newCaseFilename);

% generate a new 
NETICA_NET.reviseCPTsByCaseFile(newCaseFile,NETICA_ALL_NODES,voodooPar);

% write out the new CPT info to a new .neta file
NETICA_NET.write(Streamer(outfilename));

% shut everything down
NETICA_NET.finalize
NETICA_ENV.finalize

disp('all done!')