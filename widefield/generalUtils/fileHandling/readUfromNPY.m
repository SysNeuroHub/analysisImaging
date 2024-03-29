function U = readUfromNPY(uFilePath, varargin)
% readVfromNPY(uFilePath[, nSVtoRead])
% Reads a U matrix from an NPY file, optionally reads just a certain number
% of components
% 
% Inputs:
%   uFilepath: fullpath to U
% Outputs:
%   U: yPix x xPix x nTimePoints

if isempty(varargin)
    U = readNPY(uFilePath);
else
    nSVtoRead = varargin{1};
    [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(uFilePath);
    
    yPix = arrayShape(1);
    xPix = arrayShape(2);
    
    fid = fopen(uFilePath, 'r');
    try
        header = fread(fid, totalHeaderLength, '*char');
        U = fread(fid, yPix*xPix*nSVtoRead, ['*' dataType]);
        fclose(fid);
    catch me
        fclose(fid);
        rethrow(me);
    end
        
    U = reshape(U, yPix, xPix, nSVtoRead);
end
        