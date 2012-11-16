function makeShapefileFromStructs(Lat,Lon,mbesmstats, shapename)
%makeShapefileFromStructs(matfilename, shapename)
%
%reads data from matfile and writes it to shapefile
%matfilename must contain variables Lat, Lon, and BurialStats.
%BurialStats must contain the following fields: palpha, mean, p025, p975,
%   var, maxProb.
%shapename should be only the name -- no extensions.
WAITBAR = 0; % turn off waitbar
BOGUSVALUE = -999;

%this function expects to read variables from a mat file named matfilename.
%these variables should be: Lat, Lon, and BurialStats. For now, BurialStats
%better have the following fields: palpha, mean, p025, p975, var, maxProb.

%load data
%load(matfilename)

% check lat,lon inputs
%if(size(Lat,2)==1 | size(Lat,1)==1)
    % expand to array
    %%[Lon,Lat]= meshgrid(Lon,Lat);
%end
dlat = max(Lat(:))-min(Lat(:)) + 1e-5; % only used for tacking on ordered set of output ranges
dlon = max(Lon(:))-min(Lon(:)) + 1e-5;

% get outputs
statnames = fieldnames(mbesmstats);

% only spit out cases where inputCount>0
if(isfield(mbesmstats,'inputCount'))
    inputCount = mbesmstats.inputCount;
    doID = find(inputCount>0);
    Lat = Lat(doID);
    Lon = Lon(doID);
end

% begin
% add one point object per Lat,Lon point
Nlatlon = length(Lat(:));
gs = struct;
gs(Nlatlon,1)=struct;
if(WAITBAR)
hw = waitbar(0,'preparing shapefile');
end
for i=1:Nlatlon
    if(WAITBAR)
    waitbar(i/Nlatlon,hw) ;
    end
    % assign the geo stuff
    gs(i).Geometry = 'Point';
    %gs(i).BoundingBox = [Lat(i) Lon(i);Lat(i) Lon(i)];
    gs(i).X = Lon(i);
    gs(i).Y = Lat(i);
end
% now do data
for k=1:length(statnames)
    if(WAITBAR)
    waitbar(0,hw, ['getting ',statnames{k}]) ;
    end
    data = getfield(mbesmstats,statnames{k});
    data = data(doID);
    id = find(isnan(data));
    data(id) = BOGUSVALUE;
    for i=1:Nlatlon
        if(WAITBAR)
        waitbar(i/Nlatlon,hw);
        end
        if(i==1)
            % intit the struct now
            gs = setfield(gs,{i}, statnames{k},data(i));
        else
            gs(i) = setfield(gs(i),statnames{k},data(i));
        end
    end
    % fill in extra values that have the ranges padded into them
    Nextra = max([length(data), length(gs)]);
    for i=(Nlatlon+1):Nextra
        if(length(gs)<i)
            gs = [gs;gs(end)];
        end
        if(i>length(data))
            % keep using last value
            gs(i) = setfield(gs(i),statnames{k},data(end));
        else
            % use padded values
            gs(i) = setfield(gs(i),statnames{k},data(i));
        end
        gs(i).Geometry = 'Point';
        gs(i).X = Lon(Nlatlon)+2*dlon;
        gs(i).Y = Lat(1)+(i-1-Nlatlon)*dlat/10;
        % there are some padded values used to set the scale
        % put them in first
    end
end

% write file
try
    shapewrite(gs,shapename);
    if(WAITBAR)
    delete(hw);
    end
catch
    if(WAITBAR)
    delete(hw);
    end
    warning(sprintf('error while trying to write the shapefile %s.  Do you have it open already?',shapename))
end
return