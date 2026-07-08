function create_pills_labels_manual(subjectName, MRIdir)
% create_pills_labels_manual(subjectName, MRIdir)
% manually locates pill position using UTE data
% saves result as pills_labels.nii

% load result of pattern_generation/align_UTE_to_Allen.sh
ute = niftiread(fullfile(MRIdir, subjectName,'ute_in_t2_m2k.nii.gz'));

pills = zeros(size(ute));

zidx = (52:71);
%% show images to label
ax = axes;
for iz = 1:numel(zidx)
    ax(iz)=subplot(4,5,iz);
    imagesc(squeeze(ute(:,zidx(iz),:)));axis xy equal tight
    hold on;
    title([num2str(zidx(iz))]);
end
linkaxes(ax);

while true
    nPills = input('Enter the number of pills an integer between 1 and 10: ');
    if isscalar(nPills) && nPills == floor(nPills) && nPills >= 1 && nPills <= 10
        break;
    end
    fprintf('Invalid input. Try again.\n');
end

%% draw pill regions manually
for iz = 1:numel(zidx)
    figure('position',[0 0 2600 1300]);
    imagesc(squeeze(ute(:,zidx(iz),:)));axis xy equal tight
    title(sprintf('%d / %d images from the bottom \n Draw upto %d ROIs. Press ESC to proceed.',iz, numel(zidx), nPills));
    pills_c = zeros(size(ute,1), size(ute,3));
    for ii = 1:nPills
        title(sprintf(['%d / %d images from the bottom \n ...' ...
            'Draw %d / %d ROIs. Press ESC to skip.'],iz, numel(zidx), ii, nPills));
        roiAhand = images.roi.Polygon('Color','r');
        draw(roiAhand);
        pills_c= pills_c + createMask(roiAhand);
    end

    pills(:,zidx(iz),:)  = pills_c;
    close;
end

%% confirmation
fig = figure('position',[0 0 1800 1000]);
for iz = 1:numel(zidx)
    subplot(4,5,iz);
    imagesc(squeeze(ute(:,zidx(iz),:)));axis xy equal tight
    hold on;
    title([num2str(iz) '/' num2str(numel(zidx))]);
    contour(squeeze(pills(:,zidx(iz),:)),'Color','r');
end
screen2png(fullfile(MRIdir,'create_pills_labels_manual_result'),fig);

ute_info = niftiinfo('ute_in_t2_m2k.nii.gz');
pills_info = ute_info;
pills_info.Filename = fullfile(MRIdir, subjectName,'pills_labels');
niftiwrite(int16(pills), fullfile(MRIdir, subjectName, 'pills_labels'), ute_info);

%% sanity check
% pills_loaded = niftiread('/home/daisuke/Documents/git/analysisImaging/MROIDMD/Confucious/pills_labels.nii');
% volumeViewer(1e5*pills_loaded+ute);