
B = 10;  %Gaussian smooting factor (sd
C = 1; %erosion factor (disk & diamond)

subject = 'tmpB';
imgDir = '/home/daisuke/tmp/CCFindividual/';
mkdir(imgDir);

regDir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD/'; %will be uploaded to market
load(fullfile(regDir, subject, 'Atlas_reg_info.mat'), 'proj_brain','ROI_info');

image4DMD = [];
for roi_id = 1:60
    image4DMD(:,:,roi_id) = DMD_pattern_generation(roi_id,B,C,proj_brain,ROI_info);
end
image4DMD = uint8(round(double(intmax("uint8"))*image4DMD));

image4DMD_all = DMD_pattern_generation(1:60,B,C,proj_brain,ROI_info);
image4DMD_all = uint8(round(double(intmax("uint8"))*image4DMD_all));

save(fullfile(imgDir, subject), 'image4DMD','image4DMD_all'); %what else to save?


mkdir(fullfile(imgDir, subject));
saveEveryImages(image4DMD, fullfile(imgDir, subject)); %is this really needed?

