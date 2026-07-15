function [seriesDescription, protocolName, dirName, dirSizeBytes] = getMRprotocol(rootDir)
% [seriesDescription, protocolName, dirName] = getMRprotocol(rootDir)

seriesDescription = cell(1);
protocolName = cell(1);
dirName = cell(1);
dirSizeBytes = 0;
folders = dir(rootDir);
for k = 1:length(folders)
    if folders(k).isdir && ~startsWith(folders(k).name,'.')

        files = dir(fullfile(folders(k).folder,folders(k).name,'pdata/1/dicom/*.dcm'));

        if ~isempty(files)
            info = dicominfo(fullfile(files(1).folder,files(1).name));

            fprintf('%s\n', folders(k).name);
            dirName = [dirName folders(k).name];
            dirSizeBytes = [dirSizeBytes folders(k).bytes];
            
            if isfield(info,'SeriesDescription')
                fprintf('  SeriesDescription: %s\n',info.SeriesDescription);
                seriesDescription = [seriesDescription info.SeriesDescription];
            end

            if isfield(info,'ProtocolName')
                fprintf('  ProtocolName: %s\n',info.ProtocolName);
                protocolName = [protocolName info.ProtocolName];
            end

            fprintf('\n');
        end
    end
end

theseIdx = ~cellfun(@isempty,dirName);
dirName = dirName(theseIdx);
protocolName = protocolName(theseIdx);
seriesDescription = seriesDescription(theseIdx);
dirSizeBytes = dirSizeBytes(theseIdx);
