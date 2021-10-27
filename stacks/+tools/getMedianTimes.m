function [mFrameRate, mnFrames, mDur] = getMedianTimes(DataDir, CamType)
% [mFrameRate, mnFrames, mDur] = getMedianTimes(DataDir, camType)
% returns median framerate, number of frames, duration from all raw data within
% directory DataDir. camType is either ‘PCO’ or ‘PhotonFocus’.
% 
% See also: StackSet.LoadOneStimulus (can be used for stackset with multiple repeats)

% 13-12-25 DS from loadOneStimulus

% sort the data files in ascending order
DirContents = dir(sprintf('%s/*.mat',DataDir));
Data.nFiles = length(DirContents);

date = zeros(Data.nFiles,1);
for ifile = 1:Data.nFiles
    date(ifile) = datenum(DirContents(ifile).date);
end
[~, Order] = sort(date);
Data.FileNames = DirContents(Order);

% figure out Data.nx, Data.ny, Data.nFrames, Data.Dur, Data.FrameRate, and Data.AllTimeStamps

Data.AllTimeStamps = cell(Data.nFiles,1);
NumTimeStamps = zeros(Data.nFiles,1);
for iFile = 1:Data.nFiles
    fname = fullfile(DataDir, Data.FileNames(iFile).name);
    switch CamType
        case 'PCO'
            [Data.nx,Data.ny,TimeStamps] = tools.LoadPCO(fname,true,true);%DS on 3/3/14 
        case 'PhotonFocus'
            [Data.nx,Data.ny,TimeStamps] = tools.LoadPhotonFocus(fname);
        otherwise
            error('We cannot read this camera type yet');
    end
    Data.AllTimeStamps{iFile} = TimeStamps;
    NumTimeStamps(iFile) = nnz(isfinite(Data.AllTimeStamps{iFile}));
end

mFrameRate = 1./nanmedian(diff(Data.AllTimeStamps{1}));
mnFrames = floor(median(NumTimeStamps));
mDur = (mnFrames-1)/mFrameRate;