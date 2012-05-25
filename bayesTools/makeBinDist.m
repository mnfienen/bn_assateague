% get a robust bin distribution

M = 5; % shoot for this many bins
for j=1:length(nodeNames)
    % get the variable j
    eval(sprintf('tmpdata=%s;', nodeNames{j}));
    tmpdata = tmpdata(find(~isnan(tmpdata)));
    
    N = length(find(~isnan(tmpdata)))
    zmin = floor(min(tmpdata)/minbin(j))*minbin(j)
    zmax = ceil(max(tmpdata)/minbin(j))*minbin(j)
    %dz = (zmax-zmin)/M
    %z = zmin:dz:zmax
    [sortz,sortid] = sort(tmpdata);
    z = sortz(1:round(N/M):end); z(end)= zmax
    % adjust the bins
    %dz = max([minbin(j),median(diff(z))]);
    dz = ceil(median(diff(z))/minbin(j))*minbin(j)
    z = zmin:dz:zmax; z(end)=zmax;
    if(0)
        for i=2:(M-1)
            while((z(i)-z(i-1))< minbin(j) & z(i)<zmax)
                z(i:end-1) = z((i+1):end)
            end
        end
    end
    z = z([1;1+find(diff(z(:))>0)])';
    hist(tmpdata, z)
    title(nodeNames{j})
    xlabel(num2str(z(:)'))
    figure(gcf)
    pause
end
