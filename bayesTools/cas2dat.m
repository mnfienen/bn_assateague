function cas2dat(root)
%% convert a cas file to a dat file
% only thing is to comment out the initial line
% N.B. --> this does not, yet, handle missing data
ifp = fopen([root,'.cas'],'r');
ofp= fopen([root,'.dat'],'w');
headerline = fgets(ifp);
fprintf(ofp,'%s',['% ',headerline]);
while ~feof(ifp)
    line = fgets(ifp);
    fprintf(ofp,'%s',line);
end
fclose('all')