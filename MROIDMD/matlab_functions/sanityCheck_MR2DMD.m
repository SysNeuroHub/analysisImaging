subjectName = 'Gaius';
MRIdir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD';


cd(fullfile(MRIdir,subjectName));
volNames = {'input_swapped_m2k.nii.gz',...
    'average_swapped_m2k.nii.gz',...
    'T2w_refit.nii',...
    'T2w_resample.nii',...
    'T2w_brain.nii',...
    'T2w_brain_mask.nii.gz',...
    'T2w_brain_mask_Allen.nii',...
    'pills_labels_Allen.nii',...
    'Tem_to_T2Warped.nii.gz',...
    'ute_in_t2_m2k.nii.gz',... %result of align_UTE_to_Allen.sh
    'pills_labels.nii'};
nVolumes = numel(volNames);
for iv = 1:nVolumes
    vol{iv}=niftiread(volNames{iv});
end

%% automatically detected pills in Allen space
volumeViewer(vol{8}+uint8(vol{7}))


%% manually detected pills in T2 space
volumeViewer(max(vol{9}(:))*vol{11}+int16(vol{9}))


%NG
% f = figure;
% p1 = uipanel(f,'Position',[0,0,0.5,1]);
% p2 = uipanel(f,'Position',[0.5,0,0.5,1]);
% volshow(vol{1}, 'Parent',p1);
% volshow(vol{2}, 'Parent',p2);

%NG
% tiledlayout(3,3)
% 
% for iv = 2:nVolumes
%     ax = nexttile;
% 
%     p = patch(isosurface(vol{iv},0.5));
%     set(p,'EdgeColor','none')
% 
%     axis(ax,'equal')
%     view(ax,3)
%     camlight(ax);
%     title(volNames{iv})
% end
% 
% linkprop(findall(gcf,'Type','axes'),{'CameraPosition','CameraTarget'})