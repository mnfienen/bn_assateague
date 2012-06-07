

NetName = 'SLR00.neta';


% net's java objects are available via this global 
global NETICA_ENV NETICA_NET NETICA_NETNAME NETICA_NODESTRUCT NETICA_STREAMER NETICA_ALL_NODES

% import netica stuff
import norsys.netica.*;

% fire up the netica API with correct PWD
ifp = fopen('mikeppwd.txt','r');
PWD = fscanf(ifp,'%s');
NETICA_ENV = Environ(PWD);

% open the net
NETICA_STREAMER = Streamer(NetName);
NETICA_NET = Net(NETICA_STREAMER);
NETICA_NETNAME = char(NETICA_STREAMER.getFileName);

NETICA_NET.compile;

NETICA_ALL_NODES = NETICA_NET.getNodes;
numNODES = length(NETICA_ALL_NODES);

for i = 1:numNODES
    cnode = NETICA_ALL_NODES.get(i);
    cnode.deleteTables;
end

newCaseFile = Streamer('SEAWAT2NETICA_SLR_00_abbrev.cas');

NETICA_NET.reviseCPTsByCaseFile(newCaseFile,NETICA_ALL_NODES,10.0);

NETICA_NET.write(Streamer('testSLR.neta'));

NETICA_NET.finalize
NETICA_ENV.finalize