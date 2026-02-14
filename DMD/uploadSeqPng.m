function uploadSeqPng(seqFile, expt)
%uploadSeqPng(seqFile, expt)
% find png files specified in seqFile, then upload to the market server
% png files and seqFile are supposed to be saved in the same directory (seqFolder)

% expt.subject = 'Trajan';
% expt.expDate = '2025-11-14_2';
% expt.expNum = 3;

seqFile = char(seqFile);
if ~strcmp(seqFile(end-4:end), '.seq')
    seqFile = [seqFile '.seq'];
end
if ~isfile(seqFile)
    error('SEQ file not found: %s', seqFile);
end

info = read_polyscan_seq_pngs(seqFile);
%   info.pngNames        - PNG filenames in sequence order
%   info.found           - logical array (true if resolved)
%   info.fullPaths       - resolved full paths ('' if not found)
%   info.seqFolder       - folder containing the .seq file
%   info.missingNames    - filenames not found on disk

seqFolder = fileparts(seqFile);
if isempty(seqFolder)
    seqFolder = pwd;
end

thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
destination =fullfile(fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'Timeline', 'master')), 'DMD');
if ~exist(destination, 'dir')
    mkdir(destination);
end


success_seqFile = copyfile(fullfile(seqFolder, seqFile), fullfile(destination, seqFile));
if success_seqFile
    disp([seqFile ' was successfully copied']);
else
    disp([seqFile ' was Failed to copy']);
end

nPng = numel(info.pngNames);
for ip = 1:nPng
    copyfile(info.fullPaths{ip}, fullfile(destination, info.pngNames{ip}));
    if success_seqFile
        disp([info.pngNames{ip} ' was successfully copied']);
    else
        disp([info.pngNames{ip} ' was Failed to copy']);
    end
end