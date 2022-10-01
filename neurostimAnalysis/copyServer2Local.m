function [nsPath_dest, oeDirPath_dest, vargout] = copyServer2Local(origin, destination, fullOEName, fullOIName, overwrite)
% [nsPath_dest, oeDirPath_dest, vargout] = copyServer2Local(origin, destination, fullOEName, fullOIName)
% origin and destination assume yyyy/mm/dd

if nargin < 5
    overwrite = false;
end
if nargin < 4
    fullOIName = [];
end

pat = '^(?<subject>.+)\.(?<paradigm>.+)\.(?<hour>\d{2})(?<minute>\d{2})(?<second>\d{2})_(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})'; % <subject>.<paradigm>.<timestamp>
m = regexp(fullOEName,pat,'names');
if isempty(m) %for intan??
    pat = '^(?<subject>.+)\.(?<paradigm>.+)\.(?<hour>\d{2})(?<minute>\d{2})(?<second>\d{2})_(?<year>\d{2})(?<month>\d{2})(?<day>\d{2})'; % <subject>.<paradigm>.<timestamp>
    m = regexp(fullOEName,pat,'names');
    m.year = ['20' m.year];
end

%neurostim, openEphys and imaging data must be all stored under here
originDir = fullfile(origin, m.year, m.month, m.day);
destinationDir = fullfile(destination,m.year, m.month, m.day);

if ~exist(destinationDir,'dir')
    mkdir(destinationDir);
end

%% neurostim
NSName = sprintf('%s.%s.%s%s%s.mat',m.subject,m.paradigm,m.hour,m.minute,m.second);
nsPath_ori = fullfile(originDir, NSName);
nsPath_dest = fullfile(destinationDir, NSName);
if exist(nsPath_ori,'file')
    if exist(nsPath_dest,'file') && ~overwrite
        disp(['Skipped Copying Neurostim file']);
    else
        copyfile(nsPath_ori, destinationDir);
        disp(['Copied FROM ' nsPath_ori ' TO ' destinationDir]);
    end
else
    disp(['NOT FOUND: ' nsPath_ori]);
end


%% openEphys
oeDirPath_ori = fullfile(originDir, fullOEName);
oeDirPath_dest = fullfile(destinationDir, fullOEName);
if exist(oeDirPath_ori,'dir')
    if exist(oeDirPath_dest,'dir') && ~overwrite
        disp('Skipped copying OpenEphys directory');
    else
        copyfile(oeDirPath_ori, oeDirPath_dest);
        disp(['Copied FROM ' oeDirPath_ori ' TO ' oeDirPath_dest]);
    end
else
    disp(['NOT FOUND: ' oeDirPath_ori]);
end


if ~isempty(fullOIName)
    %% optical imaging
    oiDirPath_ori = fullfile(originDir, fullOIName);
    oiDirPath_dest = fullfile(destinationDir, fullOIName);
    if exist(oiDirPath_ori,'dir')
        if exist(oiDirPath_dest,'dir') && ~overwrite
            disp('Skipped copying imaging directory');
        else
            copyfile(oiDirPath_ori, oiDirPath_dest);
            disp(['Copied FROM ' oiDirPath_ori ' TO ' oiDirPath_dest]);
        end
    else
        disp(['NOT FOUND: ' oiDirPath_ori]);
    end
    vargout{1} = oiDirPath_dest;
end


