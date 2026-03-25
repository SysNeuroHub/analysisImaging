% created from stereo2DMD

resultServer = '/mnt/dshi0006_market';


%% subject info
subject = 'Confucious';
binary = 0;
extractSurface = 1;

regDir =fullfile('/home/daisuke/Documents/git/analysisImaging/MROIDMD', subject);
% regDir = fullfile('/mnt/dshi0006_market/Subjects',subject,'MR2DMDresult');%better to be somewhere as applyStereo2DMD needs to cd to this directory
load(fullfile(regDir, 'Atlas_reg_info.mat'), 'proj_brain','ROI_info',...
    'tform','tform2','mrwarpedtoDMD', 'mrwarped','image2','mrangle','autoTform','image2th');
OIsize = size(image2);



%% image in stereotaxic coordinates (common across subjects)
load('/home/daisuke/Documents/git/analysisImaging/DMD/references/camImg_20260214.mat','camImg');
MmPerPixel_oi = camImg.MmPerPixel;

modeIdx = [2 3 4 20 50 100];
modeFreq = [0.057 0.073 0.103 0.266 0.407 0.541];
for ii =1:numel(modeIdx)
    imgPrefix = ['whole_mode=' num2str(modeIdx(ii))];
    niifile = sprintf('eigenmode_allen100um_hemi=%s_frequency=%sHz.nii.gz',imgPrefix, num2str(modeFreq(ii)));
    Vmode = niftiread(fullfile('/home/daisuke/Documents/MATLAB/mouseEigenmode/hemisphere_whole/',niifile));
    Vmode = flip(permute(Vmode,[1 3 2]),3);

    %% convert stereo image into image for DMD
    [image4DMD_ori, image4OI_ori] = applyCCF2DMD(Vmode, tform, tform2, OIsize, MmPerPixel_oi, regDir, autoTform, extractSurface);

    %% separate into positive and negative nodes
    image4OI(:,:,1) = image4OI_ori.*(image4OI_ori>0); %positive nodes
    image4OI(:,:,2) = -image4OI_ori.*(image4OI_ori<0); %negative nodes

    image4DMD(:,:,1) = image4DMD_ori.*(image4DMD_ori>0);  %positive nodes
    image4DMD(:,:,2) = -image4DMD_ori.*(image4DMD_ori<0); %negative nodes

    %% scale to 1
    image4DMD = image4DMD./max(abs(image4DMD(:)));
    image4OI = image4OI./max(abs(image4OI(:)));

    %% save png images for DMD
    imgDir = fullfile('~/tmp',imgPrefix);
    mkdir(imgDir);

    save(fullfile(imgDir, [imgPrefix '_' subject]), ...
        'image4DMD_ori','image4OI_ori', 'image4DMD','image4OI', 'camImg','niifile'); %'image4OI_all_wCCF'
    imagesc(image4OI_ori);colormap(flipud(RedWhiteBlue));

    mkdir(fullfile(imgDir, subject));
    saveEveryImages(image4DMD, fullfile(imgDir, subject), binary, [imgPrefix '_' subject]);


    %% upload results to Market
    copyfile(imgDir, fullfile(resultServer, 'DMD images',imgPrefix));

    % mkdir(fullfile(resultServer, 'DMD images',imgPrefix,subject));
    % copyfile(fullfile(fullfile(imgDir, subject)), fullfile(resultServer, 'DMD images',imgPrefix,subject));

    clear image4OI image4DMD image4DMD_ori image4OI_ori
end