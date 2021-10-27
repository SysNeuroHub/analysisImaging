

function [ShutterMode, ScanMode, FrameRate, nRows_raw] = pcoedgeExtraInfo(Exps, rawDataDir)
% [ShutterMode, ScanMode, FrameRate, nRows_raw] = pcoedgeExtraInfo(Exps, rawDataDir)
% extracts information on pco.edge camera stored in Imager.Info

% 28/8/15 DS created
% 16/11/15 DS added 4th output, changed input from CamInfo to Exps
% 15/02/16 DS added option when Exps.Type == 'PCO_tif'. change 2nd input
% from FileString to rawDataDir for compatibility with the SVD scheme 


%if nargin < 2
    FileString = Exps.Cam.FileString;
%end

%CamInfo = 'ShutterMode:Rolling, ScanMode:slow, FrameRate:50';

CamInfo = Exps.Cam.Info;

C = textscan(CamInfo, '%s %s %s %*[^\n]', 'delimiter',',');

C1 = textscan(char(C{1}),'%s %s', 'delimiter', ':');
ShutterMode = char(C1{2});

C2 = textscan(char(C{2}),'%s %s', 'delimiter', ':');
ScanMode = char(C2{2});

C3 = textscan(char(C{3}),'%s %s', 'delimiter', ':');
FrameRate = str2num(char(C3{2}));

if nargin < 2
    %% extract nRows_raw
    rawDataDir = fullfile(...
        Exps.Cam.DataDir,[Exps.animal FileString],...
        num2str(Exps.iseries), num2str(Exps.iexp) );
end

switch Exps.Cam.Type
    case 'PCO'
        DirContents = dir(sprintf('%s/*.mat', rawDataDir));
    case 'PCO_tif'
        DirContents = dir(sprintf('%s/*.tif', rawDataDir));
end

nRows_cache = [];
for ifile = 1:length(DirContents) %2016/1/31 DS commented out
    try
        fname_raw = fullfile(rawDataDir, DirContents(ifile).name);
        switch Exps.Cam.Type
            case 'PCO'
                nRows_cache = tools.LoadPCO(fname_raw,1);
            case 'PCO_tif'
                %imstack = loadTiffStack(thisFile, 'tiffobj', statusDest);
                imstack = imread(fname_raw, 1);
                nRows_cache = size(imstack,1);
        end
    catch err
        continue
    end
    break;
end
nRows_imaged = nanmedian(nRows_cache);
nRows_raw = nRows_imaged * Exps.Cam.HardwareBinning; %number of row before binning
