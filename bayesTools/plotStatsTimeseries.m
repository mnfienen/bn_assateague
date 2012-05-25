function hp = plotStatsTimeseries(t,pdf,pdfRanges)
% input
%  t, time or space increment
%  pdf, pdfranges from prediction

% make sure we can integrate
pdf = pdf + eps;

% get levels
py975 = getPy(pdfRanges, pdf, .975);
py95 = getPy(pdfRanges, pdf, .95);
py75 = getPy(pdfRanges, pdf, .75);
py25 = getPy(pdfRanges, pdf, .25);
py05 = getPy(pdfRanges, pdf, .05);
py025 = getPy(pdfRanges, pdf, .025);

% plot patches
hp=patch([t;flipud(t)],[py025;flipud(py975)],[0,0,0]+0.75); 
set(hp, 'linestyle','none')
hold on
hp(2)=patch([t;flipud(t)],[py05;flipud(py95)],[0,0,0]+0.5); 
set(hp(2), 'linestyle','none')
hp(3)=patch([t;flipud(t)],[py25;flipud(py75)],[0,0,0]+.25); 
set(hp(3), 'linestyle','none')
hp(4) = legend('95%', '90%', '50%');