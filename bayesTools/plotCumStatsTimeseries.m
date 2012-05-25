function hp = plotCumStatsTimeseries(t,pdf,pdfRanges, plevels)
% input
%  t, time or space increment
%  pdf, pdfranges from prediction


% get levels
if(~exist('plevels'))
    plevels = [98,95,50];
end
pstruct = struct;
for i=1:length(plevels)
    pstruct(i,1).level = plevels(i);
    pstruct(i,1).py = getPy(pdfRanges, pdf, pstruct(i).level/100);
    pstruct(i,1).label = sprintf('%.0f%%',pstruct(i).level);
    pshade = 1-i/(length(plevels)+1);
    hp=patch([t;flipud(t)],[pstruct(i).py;0*pstruct(i).py+pdfRanges(1)],[0,0,0]+pshade); 
    set(hp, 'linestyle','none')
    hold on
end
s = char(pstruct.label);
legend(s);

    
% % plot patches
% hp=patch([t;flipud(t)],[py99;0*py99+pdfRanges(1)],[0,0,0]+0.75); 
% set(hp, 'linestyle','none')
% hold on
% hp(2)=patch([t;flipud(t)],[py95;0*py95+pdfRanges(1)],[0,0,0]+0.5); 
% set(hp(2), 'linestyle','none')
% hp(3)=patch([t;flipud(t)],[py50;0*py50+pdfRanges(1)],[0,0,0]+.25); 
% set(hp(3), 'linestyle','none')
% 
% legend('99%', '95%', '50%')
