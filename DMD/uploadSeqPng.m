function uploadSeqPng(seqFileFullpath, expt)
%uploadSeqPng(seqFile, expt)
% find png files specified in seqFile, then upload to the market server
% png files and seqFile are supposed to be saved in the same directory (seqFolder)

% expt.subject = 'Trajan';
% expt.expDate = '2025-11-14_2';
% expt.expNum = 3;

seqFileFullpath = char(seqFileFullpath);
if ~strcmp(seqFileFullpath(end-3:end), '.seq')
    seqFileFullpath = [seqFileFullpath '.seq'];
end
if ~isfile(seqFileFullpath)
    error('SEQ file not found: %s', seqFileFullpath);
end

info = read_polyscan_seq_pngs(seqFileFullpath);
%   info.pngNames        - PNG filenames in sequence order
%   info.found           - logical array (true if resolved)
%   info.fullPaths       - resolved full paths ('' if not found)
%   info.seqFolder       - folder containing the .seq file
%   info.missingNames    - filenames not found on disk

[seqFolder, seqFileName, seqFileExt] = fileparts(seqFileFullpath);
disp(['Found ' num2str(numel(info.fullPaths)) ' files in ' seqFileName]);
if isempty(seqFolder)
    seqFolder = pwd;
end

thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
destination =fullfile(fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'Timeline', 'master')), 'DMD');
if ~exist(destination, 'dir')
    mkdir(destination);
end


success_seqFile = copyfile(seqFileFullpath, fullfile(destination, [seqFileName seqFileExt]));
if success_seqFile
    disp([seqFileFullpath ' was successfully uploaded']);
else
    disp([seqFileFullpath ' was failed to upload']);
end

nPng = numel(info.pngNames);
for ip = 1:nPng
    success_pngFile = copyfile(info.fullPaths{ip}, fullfile(destination, info.pngNames{ip}));
    if success_pngFile
        disp([info.pngNames{ip} ' was successfully uploaded']);
    else
        disp([info.pngNames{ip} ' was failed to upload']);
    end
end