import numpy as np
import os, glob

folders_to_summarize = ['SLR00_2']

for cf in folders_to_summarize:
    ofp = open(cf + '_summary.dat','w')
    cpath = os.path.join(os.getcwd(),cf)
    for cfile in glob.glob(os.path.join(cpath,'*.dat')):
        indat = open(cfile,'r').readlines()
        tmp = indat.pop(0)
        if 'cross val' in tmp.lower():
            # this is a summary file
            scenario = tmp.strip().split()[-1]
            cal_val = indat.pop(0).strip()
            response_var = indat.pop(0).strip()
            nfolds = int(indat.pop(0).strip().split()[-1])
            voodoo = int(indat.pop(0).strip().split()[-1])
            junkus = indat.pop(0)
            tmp = indat.pop(0).strip().split()
            skM = float(tmp[1])
            skML = float(tmp[2])
            rmseM = float(tmp[3])
            rmseML = float(tmp[4])
            LR  = float(tmp[5])
            skM_cv = 0.0
            skML_cv = 0.0
            rmseM_cv = 0.0
            rmseML_cv = 0.0
            LR_cv = 0.0
            for i in np.arange(nfolds):
                tmp = indat[i].strip().split()
                skM_cv += float(tmp[1])
                skML_cv += float(tmp[2])
                rmseM_cv += float(tmp[3])
                rmseML_cv += float(tmp[4])
                LR_cv += float(tmp[5])                
            skM_cv /= nfolds
            skML_cv /= nfolds
            rmseM_cv /= nfolds
            rmseML_cv /= nfolds
            LR_cv /= nfolds
            # write to the summary file
            del indat
            ofp.write('\n' + '#'*36 + '\n')
            ofp.write('Scenario: %s\n' %(scenario))
            ofp.write('Response var: %16s %16s %d folds, voodoo = %d\n' %(response_var,cal_val,nfolds,voodoo))
            ofp.write('Overall metrics\n%16s%16s%16s%16s%16s\n' %('skM','skML','rmseM','rmseML','LR'))
            ofp.write('%16f%16f%16f%16f%16f\n' %(skM,skML,rmseM,rmseML,LR))       
            ofp.write('Cross_validation metrics\n%16s%16s%16s%16s%16s\n' %('skM','skML','rmseM','rmseML','LR'))
            ofp.write('%16f%16f%16f%16f%16f\n' %(skM_cv,skML_cv,rmseM_cv,rmseML_cv,LR_cv))
    ofp.close()