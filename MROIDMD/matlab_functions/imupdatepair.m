function imupdatepair(referenceImage, imageFolder)
% How to use this function:
% load reference image data or file
% turn on amber LED SIG-IN mode
% in pco.camware File>Direct Record to File, Select a folder with arbitrary filename, as "16bit tiff file"
% this folder is 2nd input to the current function

currentDir = pwd;
% Load the reference image
if strcmp(referenceImage(end-3:end),'tif')
    referenceImage = imread(referenceImage);
end

cd(imageFolder);

% Create a figure window
fig = figure;

% Loop to continuously update the latest image
while isvalid(fig)
    % Get list of TIFF files in the folder
    imageFiles = dir('*.tif');

    if ~isempty(imageFiles)
        % Sort files by date
        [~, idx] = sort([imageFiles.datenum], 'descend');
        latestImageFile = imageFiles(idx(1)).name;

        if numel(idx)>1
            otherImageFiles = {imageFiles(idx(2:end)).name};
            delete(otherImageFiles{:});
        end

        % Read the latest image
        latestImage = imread(fullfile(imageFolder, latestImageFile));

        % Display reference and latest image side by side
        %imshowpair(referenceImage, latestImage, 'montage'); %side by side
        imshowpair(referenceImage, latestImage); %superimposed

        rho = corr(single(referenceImage(:)), single(latestImage(:)));
        title(['Comparing Reference (Green) with Latest Image: ', latestImageFile '(Magenda) rho: ' num2str(rho)]);
    else
        disp('No TIFF images found in the folder.');
    end

    pause(.1);  % Pause for 1 second before checking again
end
close all

disp('Figure closed. Function terminated.');
cd(currentDir);
end

