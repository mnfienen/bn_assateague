function outputPDF = pdfInterp(inputRanges,inputPDF, outputRanges, interpMethod)
% pdf = pdfInterp(inputRanges,inputPDF, outputRanges)
%
% function to interpolate pdf on one set of bins to pdf on another
% always assumes that outputRanges include start/end limits of bins
% 
%
% Input
%   inputRanges, the bin locations of input
%   inputPDF, the pdf values
%   ouputRanges, the bin locations of ouput
%       outputPDF(i) = prob(outputRange(i)<b<outputRange(i+1)
%
% Output
%   pdf, the interpolated pdf

% interp method
if(~exist('interpMethod'))
    interpMethod = 'linear';
end

% check to see if inputRange==outputRange, and just return data
if(length(inputRanges)==length(outputRanges))
    if(all(inputRanges(:)==outputRanges(:)))
        outputPDF = inputPDF;
        return;
    end
end

% get dimension of input/ouput
[N,mi] = size(inputPDF);
mo = length(outputRanges);
inputRanges = inputRanges(:);
 
% deal with discrete case: which means a descretized continuous pdf!
if(length(inputRanges)==mi)
    % we have a discrete pdf and hope for the best
    % makes no sense to interplate truly discrete cases (i.e., what is halfway between an orange and an apple!)
    % so, make the ranges
    dr = diff(inputRanges(:));
    rstart = inputRanges(1) - 0.5*dr(1); % half-step back
    rend  = inputRanges(end) + 0.5*dr(end);
    inputRanges = [rstart; inputRanges(1:(end-1))+ 0.5*dr; rend];
end

% convert to cumulative distribution
% prob of value less than first range is zero, cause not allowed to be less!
inputPDF = [zeros(N,1), cumsum(inputPDF,2)];
% now, we have same number of inputs as ranges

% assume that the output pdf must capture all of the probability, so adjust boundaries
% fix cases where outputrange is inside input range
if(outputRanges(1)>inputRanges(1))
    % move output start to input start
    outputRanges(1) = inputRanges(1); 
elseif(outputRanges(1)<inputRanges(1))
   % force proper extraplation if the output range is outside inputRange
    inputRanges(1) = outputRanges(1);
end
if(outputRanges(end)<inputRanges(end))
    % move output end to input end
    outputRanges(end) = inputRanges(end);
elseif(outputRanges(end)>inputRanges(end))
    % move input end to output
    inputRanges(end) = outputRanges(end);
end

% interpolate
% transpose so interp recognizes thing correctly
outputPDF = interp1(inputRanges, inputPDF', outputRanges)';

% convert back from cumulative, retaining anything from the first output bin
% outputPDF = [outputPDF(:,1), diff(outputPDF,1,2)] ;
outputPDF = [diff(outputPDF,1,2)] ;
% this is one-element shorter than the outputRanges
%keyboard
