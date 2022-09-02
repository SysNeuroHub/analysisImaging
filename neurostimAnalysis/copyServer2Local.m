function [nsPath_dest, oeDirPath_dest, vargout] = copyServer2Local(origin, destination, fullOEName, fullOIName)
% [nsPath_dest, oeDirPath_dest, vargout] = copyServer2Local(origin, destination, fullOEName, fullOIName)

%TODO:
%check if destination directory is same as the origin to determine copy or skip

if nargin < 4
    fullOIName = [];
end

pat = '^(?<subject>.+)\.(?<paradigm>.+)\.(?<hour>\d{2})(?<minute>\d{2})(?<second>\d{2})_(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})'; % <subject>.<paradigm>.<timestamp>
m = regexp(fullOEName,pat,'names');


%neurostim, openEphys and imaging data must be all stored under here
originDir = fullfile(origin, m.year, m.month, m.day);
destinationDir = fullfile(destination,m.year, m.month, m.day);

if ~exist(destinationDir,'dir')
    mkdir(destinationDir);
end

%% neurostim
fullNSName = sprintf('%s.%s.%s%s%s.mat',m.subject,m.paradigm,m.hour,m.minute,m.second);
nsPath_ori = fullfile(originDir, fullNSName);
nsPath_dest = fullfile(destinationDir, fullNSName);
if exist(nsPath_ori,'file')
    copyfile(nsPath_ori, destinationDir);
else
    disp(['NOT FOUND: ' nsPath_ori]);
end


%% openEphys
oeDirPath_ori = fullfile(originDir, fullOEName);
oeDirPath_dest = fullfile(destinationDir, fullOEName);
if exist(oeDirPath_ori,'dir')
    copyfile(oeDirPath_ori, oeDirPath_dest);
else
    disp(['NOT FOUND: ' oeDirPath_ori]);
end


if ~isempty(fullOIName)
    %% optical imaging
    oiDirPath_ori = fullfile(originDir, fullOIName);
    oiDirPath_dest = fullfile(destinationDir, fullOIName);
    if exist(oiDirPath_ori,'dir')
        copyfile(oiDirPath_ori, oiDirPath_dest);
    else
        disp(['NOT FOUND: ' oiDirPath_ori]);
    end
    vargout{1} = oiDirPath_dest;
end


