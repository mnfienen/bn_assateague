function pdf = makeInputPdf(pdfRanges, pdfType,pdfParam, varType)
% pdf = makeInputPdf(pdfRanges, pdfType,pdfParam, varType)
if(~exist('varType'))
    % assume continuous
    varType = 'continuous'; % or
end

[N,m] = size(pdfParam); % deal with multiple inputs
r = length(pdfRanges);
switch lower(varType)
    case 'continuous'
        pdf = zeros(N,r-1); % not discrete
        for i=1:N
            switch lower(pdfType)
                % assume continuous!
                case {'norm', 'gaussian', 'normal'}
                    cdf = normcdf(pdfRanges(:)', pdfParam(i,1), pdfParam(i,2));
                    cdf(end) = 1; % force cdf to 1
                    pdf(i,:) = diff(cdf);
            end
        end
        
    case 'discrete'
        pdf = zeros(N,r);
        for i=1:N
            % assume that error is to be interpreted as spread across adjacent bins
            pdf(i,:) = normpdf(pdfRanges(:)', pdfParam(i,1), pdfParam(i,2));
            pdf(i,:) = pdf(i,:)/sum(pdf(i,:));
        end
end
