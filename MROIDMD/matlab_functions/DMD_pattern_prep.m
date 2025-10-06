function DMD_pattern_prep(mr_bead, mr_brain, oriimg, image2, image3, image4, angle)
% image1, mrimg (T2w_resample.nii): MRI image including reference capsuls, projected to x-y plane (will be
% supplied as .nii in future). Must be the same dimension to Atlas_anno_to_T2.nii,  CCF registered to individual MR (120x160 pixels)
% image2, (OI_bead.jpg): image taken by camera that captures reference capsuls.
% Must be from the same animal used for image 1
% image3, CCFBL_800x500_star.png (DMD_ref.jpg) : reference DMD image. MUST be
%800x500 pixels
% image4, CCFBL_800x500_star.tif' (OI_ref.jpg): image taken by camera during
%projection of image 3. Must be the same dimension to image 2

autoTform = 1; %MR-OI
autoTform2 = 1; %DMD

% if size(image2) ~= size(image4)
%     error(['image2 size ' numstr(size(image2)) ' does not match image 4 size ' num2str(size(image4))]);
% end

%% Data load

if(~exist('TF','var')) %for tform (tform2 = similarity)
    TF='affine';
    fprintf('TF not inserted; affine transformation will be applied.\n')
end

if(~exist('J','var'))
    J=70;
    fprintf('J not inserted; default J=70.\n')
end

 [~, proj_anno_cortex, ROI_info] = project_anno(oriimg, angle);

%% Transform to OI
if ~isempty(angle)
    mr_bead = imrotate3(mr_bead, angle(1), [1 0 0],'linear','crop'); %roll
    mr_bead = imrotate3(mr_bead, angle(2), [0 1 0],'linear','crop'); %pitch
    mr_bead = imrotate3(mr_bead, angle(3), [0 0 1],'linear','crop'); %yaw

    mr_brain = imrotate3(mr_brain, angle(1), [1 0 0],'linear','crop'); %roll
    mr_brain = imrotate3(mr_brain, angle(2), [0 1 0],'linear','crop'); %pitch
    mr_brain = imrotate3(mr_brain, angle(3), [0 0 1],'linear','crop'); %yaw
end

aa_bead=permute(mr_bead,[1 3 2]);  
% a1=max(aa_brain(:,:,J:end),[],3);
a1_bead=mean(aa_bead,3);

Vesselness = vesselness3D(mr_brain, 1, .5,.5,5000); %porus
[~,mrimg_brain] = getSurfaceData(Vesselness);
mrimg_brain = mrimg_brain/max(mrimg_brain(:));

%b1=rgb2gray(imread('OI_bead.jpg'));

borig=image2; %b1
mrimg_bead= double(fliplr(rot90(a1_bead)));
mrimg_bead =  mrimg_bead-prctile(mrimg_bead(:),5);
mrimg_bead = mrimg_bead./max(mrimg_bead(:));
mrimg_bead(mrimg_bead>1) = 1; mrimg_bead(mrimg_bead<0)=0;

figure;
subplot(121);imshow(mrimg_bead);title('Original MR bead image')
subplot(122);imagesc(log(borig)); axis equal tight; title('OI; draw brain boundary. Return to the first vertex then double click')
%BW=roipoly;
roiAhand = images.roi.AssistedFreehand;
draw(roiAhand);
BW = createMask(roiAhand);
mapconf=edge(double(BW));
close all;

disp('Computing transformation from MR to OI');
if autoTform %SLOW AT IMREGTFORM
    % tform = imregcorr(mrimg_bead, borig); %cannot be used for images with too different sizes
    borig_th = borig;
    borig_th(borig_th < prctile(borig_th(:), 70)) = 0;
    % mrimg_bead_th = mrimg_bead;
    % mrimg_bead_th(mrimg_bead_th > prctile(mrimg_bead_th(:),70)) = 1;

    [optimizer,metric] = imregconfig("multimodal");
    optimizer.MaximumIterations = 1e4;
    fixedRef  = imref2d(size(borig_th), 0.0104, 0.0104);  % example pixel sizes in mm
    movingRef = imref2d(size(mrimg_brain), 0.1, 0.1);
    tformCoarse = imregtform(mrimg_bead, movingRef, borig_th, fixedRef,"similarity",optimizer, metric);
    tform = imregtform(mrimg_bead, movingRef, borig_th, fixedRef,"affine",optimizer, metric, ...
        'InitialTransformation', tformCoarse);

    mrwarped = imwarp(proj_anno_cortex,movingRef,tform,'nearest','OutputView',fixedRef);
    beadwarped = imwarp(mrimg_bead,movingRef,tform,'cubic','OutputView',fixedRef);
    brainwarped = imwarp(mrimg_brain,movingRef,tform,'cubic','OutputView',fixedRef);
else
    fprintf('Mark movingPoints and fixedPoints...\n')
    [movingPoints,fixedPoints] = cpselect(mrimg_bead,borig,'Wait',true);
    tform = fitgeotrans(movingPoints,fixedPoints,TF);

    mrwarped = imwarp(proj_anno_cortex,tform,'nearest','OutputView',imref2d(size(borig)));
    beadwarped = imwarp(mrimg_bead,tform,'cubic','OutputView',imref2d(size(borig)));
    brainwarped = imwarp(mrimg_brain,tform,'cubic','OutputView',imref2d(size(borig)));
end
disp('...done');

% figure;imshow(borig);hold on;h=imshow(cat(3,ones(size(borig)),zeros(size(borig)), ...
% zeros(size(borig))));hold off;set(h,'AlphaData',beadwarped);title('Warped bead laid over OI')
figure; image(cat(3,borig,brainwarped, beadwarped) );axis equal tight off;


%% Transform to DMD

% a1=double(imread('CCFBL_800x500_star.tif'));%OI_ref.jpg;
% a1 = a1/max(a1(:));
a1 = image4; %DMD ref image captured by widefield camera
b1 =image3;%rgb2gray(imread('CCFBL_800x500_star.png'));%DMD_ref.jpg

disp('Computing transformation from OI to DMD');
if autoTform2
    tform2 = imregcorr(a1, b1);
else
        fprintf('Mark movingPoints and fixedPoints...\n')
        [movingPoints,fixedPoints] = cpselect(a1,b1,'Wait',true);
        tform2 = fitgeotrans(movingPoints,fixedPoints,TF);
end
disp('...done');

mrwarpedtoDMD = imwarp(mrwarped,tform2,'nearest','OutputView',imref2d(size(b1)));
% mrwarped2=[zeros(size(mrwarped,1),floor((size(borig,2)-size(mrwarped,2))/2)) mrwarped zeros(size(mrwarped,1),ceil((size(borig,2)-size(mrwarped,2))/2))];
% mrwarped3=[zeros(floor((size(borig,1)-size(mrwarped2,1))/2),size(mrwarped2,2));mrwarped2;zeros(ceil((size(borig,1)-size(mrwarped2,1))/2),size(mrwarped2,2))];

OIwarpedtoDMD = imwarp(a1,tform2,'cubic','OutputView',imref2d(size(b1)));
% beadwarped2=[zeros(size(beadwarped,1),floor((size(borig,2)-size(beadwarped,2))/2)) beadwarped zeros(size(beadwarped,1),ceil((size(borig,2)-size(beadwarped,2))/2))];
% beadwarped3=[zeros(floor((size(borig,1)-size(beadwarped2,1))/2),size(beadwarped2,2));beadwarped2;zeros(ceil((size(borig,1)-size(beadwarped2,1))/2),size(beadwarped2,2))];

mapconfwarpedtoDMD = 100*imwarp(mapconf,tform2,'cubic','OutputView',imref2d(size(b1)));
borigwarpedtoDMD = imwarp(borig, tform2,'cubic','OutputView',imref2d(size(b1)));

%%
proj_brain=mrwarpedtoDMD;
fprintf('Making brain images...\n')
Brainimage={};
for i=1:size(ROI_info,1)
    ROI_all=zeros(size(proj_brain));
    ROI_all(find(proj_brain==ROI_info{i,1}))=i;
    ROI_all2=imerode(ROI_all,strel('disk',1));
    ROI_all3=imerode(ROI_all2,strel('diamond',1));
    Brainimage{i}=ROI_all3;
end

ImageRow=[];
ROI_information=ROI_info;
for re=1:size(Brainimage,2)
    ImageRow=cat(3,Brainimage{re},ImageRow);
    ROI_information{re,1}=re;
end
ImageTotal=sum(ImageRow,3);

jetkey=jet(400);

figure;imagesc(ImageTotal+mapconfwarpedtoDMD);colormap([[1 1 1];jetkey]);title('Stimulation pattern')

fid=fopen('ROI_information.txt','w');
fprintf(fid,sprintf('Index  Abbr   Name   #Pixel (Image size:%d X %d) \n\n',size(proj_brain,1),size(proj_brain,2)));
for k=1:size(ROI_information,1)
fprintf(fid,'%s     %s     %s     %s     \n',string(ROI_information(k,1)),string(ROI_information(k,2)),string(ROI_information(k,3)),string(ROI_information(k,4)));
end
fclose(fid);

TotalBrainImage=ImageTotal;

save('Atlas_reg_info.mat','ROI_info','proj_brain','TotalBrainImage','tform','tform2','mapconfwarpedtoDMD',...
    'mrwarpedtoDMD',"OIwarpedtoDMD",'borigwarpedtoDMD','mrimg_bead','mrimg_brain','image2','image3','image4','proj_anno_cortex');
fprintf('Done! Execute DMD_pattern_generation([ROI indices])...\n')
