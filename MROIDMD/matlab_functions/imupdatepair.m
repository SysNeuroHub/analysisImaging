function imupdatepair(referenceImage, imageFolder)
% Define the folder containing the TIFF images
%imageFolder = 'C:\path\to\your\image\folder';  % <-- Change this to your folder path

% Load the reference image
%referenceImage = imread('C:\path\to\your\reference_image.tiff');  % <-- Change this too

% Create a figure window
figure;

% Loop to continuously update the latest image
while true
    % Get list of TIFF files in the folder
    imageFiles = dir(fullfile(imageFolder, '*.tiff'));
    
    if ~isempty(imageFiles)
        % Sort files by date
        [~, idx] = sort([imageFiles.datenum], 'descend');
        latestImageFile = imageFiles(idx(1)).name;
        
        % Read the latest image
        latestImage = imread(fullfile(imageFolder, latestImageFile));
        
        % Display reference and latest image side by side
        imshowpair(referenceImage, latestImage, 'montage');
        title(['Comparing Reference with Latest Image: ', latestImageFile]);
    else
        disp('No TIFF images found in the folder.');
    end
    
    pause(.1);  % Pause for 1 second before checking again
end

