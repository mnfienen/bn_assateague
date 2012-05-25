function hp = plotPDF2(pdfRanges,pdf,facealpha,faceColor, lineStyle)

if(~exist('facealpha'))
    facealpha =1;
end
if(~exist('faceColor'))
    faceColor ='r';
end
if(~exist('lineStyle'))
    lineStyle ='-';
end

hp = zeros(size(pdf));
for j=1:length(pdf)
    %     patch([pdfRanges(j:(j+1)); pdfRanges(j+1);pdfRanges(j)], [pdfp(j);pdfp(j);0;0], 'r')
    %     hold on
    %     pdf = getfield(pred.pdf, [priornamesout{i}]);
    %     pdf = pdf(idcomp,:);
    if(length(pdf)~=length(pdfRanges))
        % continuous variable
        hp(j)=patch([pdfRanges(j:(j+1)); pdfRanges(j+1);pdfRanges(j)], [pdf(j);pdf(j);0;0], 'r');
    else
        % discrete variable
        hp(j)=patch([j;(j+1);(j+1);(j)]-1.5, [pdf(j);pdf(j);0;0], 'r');
    end
    set(hp(j), 'faceAlpha',facealpha, 'faceColor', faceColor, 'LineStyle', lineStyle);
end
xlim([min(pdfRanges) max(pdfRanges)])
ylim([0 1])