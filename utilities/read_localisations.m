function data = read_localisations(filename,header,numCols,delimiter)
        
    if ~isempty(header)
        try
            format = ['%f',' '];
            for ii = 1 : numCols-1
                format = [format,'%f',' '];
            end
            
            fid = fopen(filename, 'rt');
            data = textscan(fid, format, 'HeaderLines', header, 'CollectOutput', true, 'Delimiter',delimiter);
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