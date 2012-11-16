nfolds = 10;
ntestcases = 2;
cas_dat_root = 'SEAWAT2NETICA_SLR_00';
voodoopar = [100];

for i = 1:length(voodoopar)
    CV_evaluate_skill(nfolds,ntestcases,cas_dat_root,voodoopar(i))
end