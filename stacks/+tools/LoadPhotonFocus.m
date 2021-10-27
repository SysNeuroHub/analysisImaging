function [nRows, nCols, TimeStamps, Stack] = LoadPhotonFocus(FileName)
% Loads a PCO file and returns either basic info or the whole stack
%
% [nRows, nCols, TimeStamps] = LoadPhotonFocus(FileName) 
% returns only the basic info (very quick).
%
% [nRows, nCols, TimeStamps, Stack] = LoadPhotonFocus(FileName)
% returns the whole stack of data (can be slow).
%
% 2013-10-13 Matteo Carandini
    
[~, ShortFileName, ~] = fileparts(FileName);

if nargout < 4
    GetInfoOnly = true;
    fprintf('Getting basic info on PhotonFocus file %s. ', ShortFileName)
else
    GetInfoOnly = false;
    fprintf('Loading stack from PhotonFocus file %s. ', ShortFileName)
end

fid = fopen(FileName, 'r');
if fid == -1, error(['Unable to open file: ' FileName]); end
pp = fread(fid,[6 1],'double');
[totFrames, extraMem, Offset , ~, nRows, nCols] = deal(pp(1),pp(2),pp(3),pp(4),pp(5),pp(6));
fclose(fid);

nAllocFrames = (totFrames + round(totFrames/(100/(extraMem-1))));

m = memmapfile(FileName);
m.offset = Offset;
m.format = {...
    'double' [1 nAllocFrames] 'tt';...
    'double' [1 nAllocFrames] 'ttAbs';...
    'uint16' [nRows nCols 1 nAllocFrames] 'Tensor'};

nFrames = nnz(m.Data(1).tt);

TimeStamps  = m.Data(1).tt(1:nFrames);

fprintf('\n');

if GetInfoOnly, return; end

Stack  = squeeze(m.Data(1).Tensor(:,:,1,1:nFrames)); % this can be very slow






