function compile_net_EM(NetName,newCaseFilename,voodooPar,outfilename)
% compile_net(NetName,newCaseFilename,voodoPar,outfilename)
% a m!ke@usgs joint <mnfienen@usgs.gov>
% function to build the CPT tables for a new CAS file on an existing NET
% (be existing, meaning that the nodes, edges, and bins are dialed)
% INPUT:
%       NetName --> a filename, including '.neta' extension
%       newCaseFilename --> new case file including '.cas' extension
%       voodoPar --> the voodoo tuning parameter for building CPTs
%       outfilename --> netica file for newly build net (including '.neta'
global NETICA_ENV NETICA_NET NETICA_NETNAME NETICA_NODESTRUCT NETICA_STREAMER NETICA_ALL_NODES
global NETICA_CASEFILE_STREAMER NETICA_CS
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


% make a Learner, which is Neticaspeak for using a CPT learning algorithm
NETICA_LEARNER = Learner(Learner.EM_LEARNING);
NETICA_LEARNER.setMaxTolerance(1.0e-6);

NETICA_LEARNER.setMaxIterations(1000);

% make a caseset for the cases to use with the learner
NETICA_CASEFILE_STREAMER = Streamer(newCaseFilename);
NETICA_CS = Caseset;
NETICA_CS.addCases(NETICA_CASEFILE_STREAMER,1.0,[]);

NETICA_LEARNER.learnCPTs(NETICA_ALL_NODES,NETICA_CS,voodooPar);



% write out the new CPT info to a new .neta file
NETICA_NET.write(Streamer(outfilename));

% shut everything down
NETICA_CASEFILE_STREAMER.finalize
NETICA_LEARNER.finalize
NETICA_NET.finalize
NETICA_ENV.finalize


disp('New Net created by Learner.EM_LEARNING')