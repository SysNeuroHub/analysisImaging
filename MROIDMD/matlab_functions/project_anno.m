function [proj_anno, proj_anno_cortex, ROI_info] = project_anno(oriimg, angle)
%oriimg = niftiread('Atlas_anno_to_T2.nii');
% created from DMD_pattern_prep.m

if nargin < 2
    angle = [];
end

oriimg(oriimg==2000)=0; img=oriimg;
 
%% Remove layers (Right hemi)

for i = [107 114 100 37 51 72 93 79 65 325 332 346 360 150 122 353 346 171 206 136 143 185 164 199 178 86 192 199 213 238]
    for j = 1:6
        img(oriimg == i+j) = i;
    end
end

for i = [285 45 18 24 226 232 298 367 373]
    for j = 1:5
        img(oriimg == i+j) = i;
    end
end

for j=1:12
    img(oriimg==494+j)=494;
end

%% Remove layers (Left hemi)

for i = [107 114 100 37 51 72 93 79 65 325 332 346 360 150 122 353 346 171 206 136 143 185 164 199 178 86 192 199 213 238]+2000
    for j = 1:6
        img(oriimg == i+j) = i;
    end
end

for i = [285 45 18 24 226 232 298 367 373]+2000
    for j = 1:5
        img(oriimg == i+j) = i;
    end
end

for j=1:12
    img(oriimg==2494+j)=2494;
end

%% Projection to +z direction (Annotation)
anno = img;
if ~isempty(angle)
    anno = imrotate3(anno, angle(1), [1 0 0],'nearest','crop'); %roll
    anno = imrotate3(anno, angle(2), [0 1 0],'nearest','crop'); %pitch 
    anno = imrotate3(anno, angle(3), [0 0 1],'nearest','crop'); %yaw 
end

anno=flip(rot90(permute(anno,[3 1 2]),2),3);


proj_anno=zeros([size(anno,1) size(anno,2)]);
for i=1:size(anno,3)
    for j=1:size(anno,1)
        for k=1:size(anno,2)
            if anno(j,k,size(anno,3)+1-i)~=0
                proj_anno(j,k)=anno(j,k,size(anno,3)+1-i);
            end
        end
    end
end

%% ROI extraction

fprintf('Extracting ROIs...\n')
proj_anno_cortex=zeros(size(proj_anno));
surviveR=[18 24 44 51 65 72 79 86 93 100 122 136 143 150 164 171 178 185 192 199 206 213 226 238 298 325 332 346 353 360];
surviveL=surviveR+2000;
survive=[surviveR surviveL];
ROI_info=cell(length(survive),4);

for i=1:length(survive) % 60 regions 
proj_anno_cortex(find(proj_anno==survive(i)))=proj_anno(find(proj_anno==survive(i)));
A=Load_ATLAS_info(survive(i));
ROI_info{i,1}=A{1}; % Original value
ROI_info{i,2}=A{2}; % Name, abbr
ROI_info{i,3}=A{3}; % Name, full
ROI_info{i,4}=length(find(proj_anno_cortex==survive(i))); % # pixels
end