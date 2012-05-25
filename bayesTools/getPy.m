function py = getPy(pdfRanges, pdf, p)
% py is value of y with prob p
M = length(pdfRanges);
[N,Mb] = size(pdf);
f = cumsum(pdf,2);
f = [zeros(N,1),f];

if(Mb==length(pdfRanges))
    % use discrete
    
    % switch for upper/lower bounds
    if(p<=0.5)
        % lower bound
        f = f<p;
        f = sum(f,2);
    else
        % upper bound
        f = (-f)<(-p);
        f = M-sum(f,2)+1;
        f(f>M)=M; % catch overrun
    end
    py = pdfRanges(f);
else
    % use continous and interpolation
    py = nan*ones(N,1);
    for i=1:N
        try
            % deal with 0,1 at ends
            [funique, uniqueid] = unique(f(i,:));
            py(i) = interp1(funique,pdfRanges(uniqueid), p);
        catch
            lasterr
            keyboard
        end
    end
end