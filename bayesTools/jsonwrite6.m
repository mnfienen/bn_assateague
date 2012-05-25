function jsonwrite4(fname, S)
% make this format
%
% {
% predictions:[
%     [1,    -67.12310791,45.09410095,0,0,0,1,0,0,0],
%     [2,    -67.11375772,45.05,0,0,0,1,0,0,0],
%     [3,    -67.1,45.0351308,0,0,0,1,0,0,0],
%     [4,    -67.09468842,45.0348587,0,0,0,1,0,0,0],
%     [5,    -67.08534187,45.02030148,0,0,0,1,0,0,0],
%     [6,    -67.07081604,45,0,0,0,1,0,0,0],
%     [7,    -67.05,44.96456076,0,0,0,1,0,0,0],
%     ...
%     [11922,-68.79552626,44.54210292,0,0,0,1,0,0,0],
%     [11923,-68.80697738,44.50998369,0,0,0,1,0,0,0]
% ]
% }

flds = fieldnames(S);
fid = fopen(fname, 'wt');
% open the object
fprintf(fid, '{\n');
% loop through fields
cnt = length(S.X);
fprintf(fid, 'COUNT: %d,\n', cnt);
for j=1:length(flds)
            fprintf('processing %s...\n', flds{j})
    switch flds{j}
        case {'X'}
            % geo stuff
            StmpX = S.(flds{j});
            StmpY = S.(flds{j+1});
            [Npdf] = length(StmpX);
            fprintf(fid, 'XY:{\n');
            for i=1:Npdf
                fprintf(fid, '  %d:  [    ', (i));
                fprintf(fid, '%.9f, %.9f', StmpX(i), StmpY(i));
                fprintf(fid, ']'); % close the line
                if(i<Npdf)
                    fprintf(fid, ',\n');
                else
                    fprintf(fid,'\n'); % last line
                end
            end
            if flds{j} == 'X'
                fprintf(fid,'},\n');
            else
                fprintf(fid, '},\n');
            end
        case 'id'
            % do nothing
        case 'Y'
            % do nothing, done in X
            
        otherwise
            %% write a new object
            fprintf(fid, '%s:\n{', flds{j});
            Stmp = S.(flds{j});
            % get pdf
            pdf = Stmp.pdf; % get the data
            pdfRanges = Stmp.ranges;
            [Npdf,Mpdf] = size(pdf);
            %% write pdf
            fprintf(fid, 'pdf:{\n');
            for i=1:Npdf
                % all need an id
                fprintf(fid, ' %d:   [    ', (i));
                for k = 1:Mpdf
                    fprintf(fid, '%.3g', pdf(i,k));
                    if k<Mpdf
                        fprintf(fid, ',');
                    end
                end
                fprintf(fid, ']'); % close the line
                
                if(i<Npdf)
                    fprintf(fid, ',\n');
                else
                    fprintf(fid,'\n'); % last line
                end
            end
            fprintf(fid, '},\n');
            
            %% get cdf
            cpdf = Stmp.cdf; % get the data
            [Ncdf,Mcdf] = size(cpdf);
            % write pdf
            fprintf(fid, 'cdf:{\n');
            for i=1:Npdf
                % all need an id
                fprintf(fid, ' %d:   [    ', (i));
                for k = 1:Mcdf
                    fprintf(fid, '%.3g', cpdf(i,k));
                    if k<Mcdf
                        fprintf(fid, ',');
                    end
                end
                fprintf(fid, ']'); % close the line
                if(i<Ncdf)
                    fprintf(fid, ',\n');
                else
                    fprintf(fid,'\n'); % last line
                end
            end
            fprintf(fid, '},\n');
            
            % do pdf ranges
            m = length(pdfRanges);
            fprintf(fid, 'pdfRanges:[');
            for k = 1:m
                fprintf(fid, '%.3g ', pdfRanges(k));
                if k < m
                    fprintf(fid, ',');
                end
            end
            fprintf(fid,']\n},\n');
    end
end
fprintf(fid,'}\n');
fclose(fid);
return