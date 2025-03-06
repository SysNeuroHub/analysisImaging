function [result, mtxDir, tlDir] = uploadMtx(animal, session, expNumber, mightexDirectory)
% copty mightex experiment data (camera and DMD) to the directory where
% timeline data is saved
% then rename tifstack so the pipeline can run

p = dat.paths;

if nargin < 4
    mightexDirectory = 'C:\Users\experiment\Documents\Mightex\PolyScan3\User Data\Daisuke_Shimaoka\DMDcamera\20241029_Daisuke_Shimaoka';
end

if isnumeric(expNumber)
    expNumber = num2str(expNumber);
end

%% detect date & time of a specified experiment number from Timeline data
%timeline data:
%C:\LocalExpData\susanoo\2024-11-22_1\6\2024-11-22_1_6_susanoo_Timeline.mat
% mtx data: C:\Users\experiment\Documents\Mightex\PolyScan3\User Data\Daisuke_Shimaoka\DMDcamera\20241029_Daisuke_Shimaoka\20241122_191702
tlDir = fullfile(p.localRepository, animal, session, expNumber);
%tlInfo = dir(fullfile(tlDir, [session '_' expNumber '_' animal '_Timeline.mat']));
%tlTimeStamp = tlInfo.datenum;
dtl = System.IO.File.GetCreationTime(tlDir);
tlTimeStamp = datenum(double([dtl.Year, dtl.Month, dtl.Day, dtl.Hour, dtl.Minute, dtl.Second]));

mtxCandidateInfo = dir(mightexDirectory);
mtxCandidateInfo = mtxCandidateInfo([mtxCandidateInfo.isdir]);
mtxCandidateInfo = mtxCandidateInfo( ~ismember({mtxCandidateInfo.name}, {'.','..'}));
%mtxCandidateTimeStamp = [mtxCandidateInfo.datenum]; %NG last access date, not creation date
mtxCandidateTimeStamp = nan(numel(mtxCandidateInfo),1);
for idir = 1:numel(mtxCandidateInfo)
    d = System.IO.File.GetCreationTime(fullfile(mtxCandidateInfo(idir).folder, ...
        mtxCandidateInfo(idir).name, 'SessionEvents.dat'));
    mtxCandidateTimeStamp(idir) = datenum(double([d.Year, d.Month, d.Day, d.Hour, d.Minute, d.Second]));
end

[~, tgtDirIdx] = min(abs(mtxCandidateTimeStamp - tlTimeStamp)); %may not be robust
if isempty(tgtDirIdx)
    error('Could not find directory where mightex data has been saved');
end

%% identify the mightex directory corresponding to the experiment
mtxInfo = mtxCandidateInfo(tgtDirIdx);
mtxDir = fullfile(mtxInfo.folder, mtxInfo.name);
disp(['Mightex data @: ' mtxDir]);

%% check if a tifstack has been created (if not, create here?)
tifStackInfo = dir(fullfile(mtxDir, '*.tif'));
tifStackInfo = tifStackInfo([~tifStackInfo.isdir]);
if isempty(tifStackInfo)
    error('Could not find a tifstack');
else
    originalTif = fullfile(tifStackInfo.folder, tifStackInfo.name);
end

%% rename the tiffstack to (expnumber).tif
%renamedTif = fullfile(tlPath, [expNumber '.tif']);
renamedTif = fullfile(tifStackInfo.folder, [expNumber '.tif']);
if ~isequal(originalTif, renamedTif)
    result = movefile(originalTif, renamedTif);
    if ~result
        error(['Could not rename the tifstack to ' [expNumber '.tif']]);
    end
end

%% copy the tifstack and other mightx data 
result = copyfile(mtxDir, tlDir);
if result
    disp(['Mightex data successfully copied to ' tlDir]);
elseif ~result
    error(['Could not copy mightex data to ' tlDir]);
end

