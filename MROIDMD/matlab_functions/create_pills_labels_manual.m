function create_pills_labels_manual(subjectName, MRIdir)
% subjectName = 'Confucious';
% MRIdir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD';


% load result of pattern_generation/align_UTE_to_Allen.sh
ute = niftiread(fullfile(MRIdir, subjectName,'ute_in_t2_m2k.nii.gz'));

pills = zeros(size(ute));

%% draw pill regions manually
zidx = 52:66;
for iz = 1:numel(zidx)
    figure('position',[0 0 2600 1300]);
    imagesc(squeeze(ute(:,zidx(iz),:)));axis equal tight
    title([num2str(iz) '/' num2str(numel(zidx))]);

    pills_c = zeros(size(ute,1), size(ute,3));
    for ii = 1:4
        roiAhand = images.roi.Polygon('Color','r');
        draw(roiAhand);
        pills_c= pills_c + createMask(roiAhand);
    end

    pills(:,zidx(iz),:)  = pills_c;
close;
end

%% confirmation
for iz = 1:numel(zidx)
    subplot(3,5,iz);
    imagesc(squeeze(ute(:,zidx(iz),:)));axis equal tight
    hold on;
    title([num2str(iz) '/' num2str(numel(zidx))]);
    contour(squeeze(pills(:,zidx(iz),:)),'Color','r');
end

ute_info = niftiinfo('ute_in_t2_m2k.nii.gz');
pills_info = ute_info;
pills_info.Filename = fullfile(MRIdir, subjectName,'pills_labels');
niftiwrite(int16(pills), fullfile(MRIdir, subjectName, 'pills_labels'), ute_info);

%% sanity check
 % pills_loaded = niftiread('/home/daisuke/Documents/git/analysisImaging/MROIDMD/Confucious/pills_labels.nii');
 % volumeViewer(1e5*pills_loaded+ute);