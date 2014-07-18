function data = readLocalisations(filename,header,numCols)
        
    if ~isempty(header)
        try
            format = ['%f',' '];
            for ii = 1 : numCols-1
                format = [format,'%f',' '];
            end
            
            fid = fopen(filename, 'rt');
            data = textscan(fid, format, 'headerLines', header, 'CollectOutput', true);
            fclose(fid);
            data = cell2mat( data);
        catch err %#ok, unused arguments
            errordlg('File could not be found, or wrong number of headerlines or columns','modal');
        end
    else
        try
            format = ['%f',' '];
            for ii = 1 : numCols-1
                format = [format,'%f',' '];
            end
            fid = fopen(filename, 'rt');
            data = textscan(fid, format,'CollectOutput', true);
            fclose(fid);
            data = cell2mat( data);
        catch err %#ok, unused arguments
            errordlg('File could not be found, or wrong number of headerlines or columns','modal');
        end
    end
end