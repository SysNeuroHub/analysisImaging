

dir_DS = '/home/daisuke/Dropbox/GSlab/MRI_alignment/Kim2023_test';
dir_Chris = '/mnt/dshi0006_market/MRI/analysis/MRA056-selection';

for ii = 1:3
    switch ii
        case 1
            subject = 'tmpB';
        case 2
            subject = 'tmpC';
        case 3
            subject = 'tmpD';
    end

    p_tmp(:,:,:,ii)=niftiread(fullfile(dir_DS, subject, 'pills_labels_Allen.nii'))>0;
end
p = sum(p_tmp,4);

b=niftiread(fullfile(dir_DS, subject, 'T2w_brain_mask_Allen.nii'));
p0=patch(isosurface(b));
p0.FaceColor = [0 1 0];
p0.EdgeColor = 'none';
p0.FaceAlpha = 0.4;

hold on;

p2=patch(isosurface(p));
p2.FaceColor = [1 0 0];
p2.EdgeColor = 'none';
p2.FaceAlpha = 0.4;


camlight; lighting gouraud

axis equal tight off;
view(90,0);