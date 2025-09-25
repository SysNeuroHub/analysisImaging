dir_DS = '/home/daisuke/Dropbox/GSlab/MRI_alignment/Kim2023_test';
dir_Chris = '/mnt/dshi0006_market/MRI/analysis/MRA056-selection';
subject = 'tmpB';

p=niftiread(fullfile(dir_DS, subject, 'pills_labels.nii'))>0;
b=niftiread(fullfile(dir_DS, subject, 'T2w_brain.nii'));
pb = b + max(b(:)).*p;
r=niftiread(fullfile(dir_DS, subject, 'T2w_resample.nii'));
%r=niftiread(fullfile(dir_DS, subject, 'average_swapped_m2k.nii.gz'));
%r=niftiread(fullfile(dir_DS, subject, 'input_swapped_m2k.nii.gz')); r= squeeze(r(:,:,:,1));

% %Allen space ... extracted pill and raw align
% p=niftiread(fullfile(dir_DS, subject, 'pills_labels_Allen.nii'))>0;
% b=niftiread(fullfile(dir_DS, subject, 'T2w_brain_mask_Allen.nii'));
% r=niftiread(fullfile(dir_DS, subject, 'input_swapped_m2k_Allen.nii.gz')); r= squeeze(r(:,:,:,1));

rpb = r .* int16(pb~=0);


p0=patch(isosurface(r));
p0.FaceColor = [0 1 0];
p0.EdgeColor = 'none';
p0.FaceAlpha = 0.4;

hold on;

p1=patch(isosurface(b));
p1.FaceColor = [0 0 1];
p1.EdgeColor = 'none';
p1.FaceAlpha = 0.4;

%% pills, detected by me
p2=patch(isosurface(p));
p2.FaceColor = [1 0 0];
p2.EdgeColor = 'none';
p2.FaceAlpha = 0.4;

%% pills, detected by Chris
p_ori=niftiread(fullfile(dir_Chris, subject, 'pills_labels.nii.gz'))>0;
p_ori=rot90(permute(rot90(p_ori,1),[2 3 1]),2);
p2_ori=patch(isosurface(p_ori));
p2_ori.FaceColor = [1 1 0];
p2_ori.EdgeColor = 'none';
p2_ori.FaceAlpha = 0.4;


camlight; lighting gouraud

axis equal tight off;
view(90,0);
legend('raw','brain','pill DS', 'pill Chris');


