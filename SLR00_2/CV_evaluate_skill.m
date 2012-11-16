function CV_evaluate_skill(nfolds,ntestcases,cas_dat_root,voodoopar)
% Postprocessing to evaluate skill using the saved output from a
% cross-validation run. Required variables are:
% nfolds --> number of cross validation folds
% ntestcases --> number of test cases
% cas_dat_root --> root name of the cas/dat files used to create the net


calval = {'calibration','validation'};


global NETICA_ENV NETICA_NET;

cross_val_driver('SLR00_2.neta',voodoopar,[cas_dat_root,'.cas'],nfolds);

%% first set up some text files for output
    for j = 1:ntestcases
      indatafile = [cas_dat_root, '_', num2str(1), '_fold-', ...
            num2str(1), '_ct-',num2str(j), '.mat'];
       load(indatafile)
       testcasenameval =  savedStats.valsavedStats(j).testcasename;
       testcasenamecal =  savedStats.calsavedStats(j).testcasename;
       numvars = length(savedStats.valsavedStats);
       
           for k = 1: numvars
               varnme = savedStats.valsavedStats(k).VARNAME;
               ofp = fopen([testcasenamecal, '_CV_summary_',calval{1},'_',varnme,'_voodoo_', num2str(voodoopar),'.dat'], 'w');
               fprintf(ofp,'Cross Validation Summary for %s\n %s\nresponse variable:%s\nNfolds: %d\nvoodoopar: %d\n', ...
                   testcasenamecal,calval{1},varnme,nfolds,voodoopar);
               fprintf(ofp,'%-16s %-16s %-16s %-16s %-16s %-16s\n', 'fold','skM','skML','rmseM','rmseML','LR');
               fclose(ofp);
               ofp = fopen([testcasenameval, '_CV_summary_',calval{2},'_',varnme,'_voodoo_', num2str(voodoopar),'.dat'], 'w');
               fprintf(ofp,'Cross Validation Summary for %s\n %s\nresponse variable:%s\nNfolds: %d\nvoodoopar: %d\n', ...
                   testcasenameval,calval{2},varnme,nfolds,voodoopar);
               fprintf(ofp,'%-16s %-16s %-16s %-16s %-16s %-16s\n', 'fold','skM','skML','rmseM','rmseML','LR');
               fclose(ofp);
           end; % k
    end % j in ntestcases

%% now, populate each file with the output data
for j = 1:ntestcases
       for cv = 1:2
           for k = 1: numvars
               varnme = savedStats.valsavedStats(k).VARNAME;
               testcasenameval =  savedStats.valsavedStats(j).testcasename;
               testcasenamecal =  savedStats.calsavedStats(j).testcasename;
               indatafile = [cas_dat_root, ...
                  '_ct-',num2str(j), '_NO_CV.mat'];
               load(indatafile);
              switch cv
                  case 1
                      ofp = fopen([testcasenamecal, '_CV_summary_',calval{cv},'_',varnme,'_voodoo_', num2str(voodoopar),'.dat'], 'a');
                  case 2
                      ofp = fopen([testcasenameval, '_CV_summary_',calval{cv},'_',varnme,'_voodoo_', num2str(voodoopar),'.dat'], 'a');
              end %switch cv for opening file and writing all out
              fprintf(ofp,'%-16d %-16f %-16f %-16f %-16f %-16f\n', 0, ...
                  savedStats.ALLsavedStats(k).skM, ...
                  savedStats.ALLsavedStats(k).skML, ...
				  savedStats.ALLsavedStats(k).rmseM, ...
				  savedStats.ALLsavedStats(k).rmseML, ...
                  savedStats.ALLsavedStats(k).LR);
              fclose(ofp);
               for i = 1:nfolds
                clear savedStats;
                indatafile = [cas_dat_root, '_', num2str(i), '_fold-', ...
                    num2str(i), '_ct-',num2str(j), '.mat'];
                load(indatafile)
                      switch cv
                          case 1
                              ofp = fopen([testcasenamecal, '_CV_summary_',calval{cv},'_',varnme,'_voodoo_', num2str(voodoopar),'.dat'], 'a');

                              fprintf(ofp,'%-16d %-16f %-16f %-16f %-16f %-16f\n', i, ...
                              savedStats.calsavedStats(k).skM, ...
                              savedStats.calsavedStats(k).skML, ...
							  savedStats.calsavedStats(k).rmseM, ...
							  savedStats.calsavedStats(k).rmseML, ...
                              savedStats.calsavedStats(k).LR);
                              fclose(ofp);  
                          case 2
                              ofp = fopen([testcasenameval, '_CV_summary_',calval{cv},'_',varnme,'_voodoo_', num2str(voodoopar),'.dat'], 'a');

                              fprintf(ofp,'%-16d %-16f %-16f %-16f %-16f %-16f\n', i, ...
                              savedStats.valsavedStats(k).skM, ...
                              savedStats.valsavedStats(k).skML, ...
                              savedStats.valsavedStats(k).rmseM, ...
                              savedStats.valsavedStats(k).rmseML, ...
                              savedStats.valsavedStats(k).LR);
                              fclose(ofp);
                      end % switch cv
                end % i in nfolds

               end % k
           end % cv
        
        
    
end % j in ntestcases