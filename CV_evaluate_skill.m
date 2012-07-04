% Postprocessing to evaluate skill using the saved output from a
% cross-validation run. Required variables are:
% nfolds --> number of cross validation folds
% ntestcases --> number of test cases
% cas_dat_root --> root name of the cas/dat files used to create the net

nfolds = 10;
ntestcases = 2;
cas_dat_root = 'SEAWAT2NETICA_SLR_00';
calval = {'calibration','validation'};

%% first set up some text files for output
    for j = 1:ntestcases
      indatafile = [cas_dat_root, '_', num2str(1), '_fold-', ...
            num2str(1), '_ct-',num2str(j), '.mat'];
       load(indatafile)
       testcasename =  savedStats.valsavedStats(1).testcasename;
       numvars = length(savedStats.valsavedStats);
       
       for cv = 1:2
           for k = 1: numvars
               varnme = savedStats.valsavedStats(k).VARNAME;
               ofp = fopen([testcasename, '_CV_summary_',calval{cv},'_',varnme,'.dat'], 'w');
               fprintf(ofp,'Cross Validation Summary for %s: %s in %s\n',testcasename,calval{cv},varnme);
               fprintf(ofp,'%-16s %-16s %-16s %-16s\n', 'fold','skM','skML','LR');
               fclose(ofp);
           end; % k
       end %cv
    end % j in ntestcases

%% now, populate each file with the output data
for j = 1:ntestcases
       for cv = 1:2
           for k = 1: numvars
               varnme = savedStats.valsavedStats(k).VARNAME;
               testcasename =  savedStats.valsavedStats(1).testcasename;
               indatafile = [cas_dat_root, ...
                  '_ct-',num2str(j), '_NO_CV.mat'];
               load(indatafile);
              ofp = fopen([testcasename, '_CV_summary_',calval{cv},'_',varnme,'.dat'], 'a');
              fprintf(ofp,'%-16d %-16f %-16f %-16f\n', 0, ...
                  savedStats.ALLsavedStats(k).skM, ...
                  savedStats.ALLsavedStats(k).skML, ...
                  savedStats.ALLsavedStats(k).LR);
              fclose(ofp);
               for i = 1:nfolds
                clear savedStats;
                indatafile = [cas_dat_root, '_', num2str(i), '_fold-', ...
                    num2str(i), '_ct-',num2str(j), '.mat'];
                load(indatafile)
                       ofp = fopen([testcasename, '_CV_summary_',calval{cv},'_',varnme,'.dat'], 'a');
                      switch cv
                          case 1
                              fprintf(ofp,'%-16d %-16f %-16f %-16f\n', i, ...
                              savedStats.calsavedStats(k).skM, ...
                              savedStats.calsavedStats(k).skML, ...
                              savedStats.calsavedStats(k).LR);
                              fclose(ofp);  
                          case 2
                              fprintf(ofp,'%-16d %-16f %-16f %-16f\n', i, ...
                              savedStats.valsavedStats(k).skM, ...
                              savedStats.valsavedStats(k).skML, ...
                              savedStats.valsavedStats(k).LR);
                              fclose(ofp);
                      end % switch cv
                end % i in nfolds

               end % k
           end % cv
        
        
    
end % j in ntestcases