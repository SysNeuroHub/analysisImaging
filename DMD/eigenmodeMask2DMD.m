% created from stereo2DMD

resultServer = '/mnt/dshi0006_market';


%% subject info
subject = 'Confucious';
binary = 0;
extractSurface = 0;

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
nsites = 2;
for ii =4:numel(modeIdx)
    imgPrefix = ['mask_whole_mode=' num2str(modeIdx(ii)) '_numsites=' num2str(nsites)];
    %    niifile = sprintf('eigenmode_allen100um_hemi=%s_frequency=%sHz.nii.gz',imgPrefix, num2str(modeFreq(ii)));
    mage4DMD = []; image4OI = [];
    for pp = 1:2
        switch pp
            case 1
                polarity = 'positive';
            case 2
                polarity = 'negative';
        end
        niifile = sprintf('mask_allen100um_hemi=whole_mode=%d_frequency=%sHz_type=antinode_%s_numsites=%d_numneighbours=124.nii.gz',...
            modeIdx(ii), num2str(modeFreq(ii)),polarity, nsites);
        Vmode = niftiread(fullfile('/home/daisuke/Documents/MATLAB/mouseEigenmode/hemisphere_whole/',niifile));
        Vmode = flip(permute(Vmode,[1 3 2]),3);

        %% convert stereo image into image for DMD
        [image4DMD(:,:,pp), image4OI(:,:,pp)] = applyCCF2DMD(Vmode, tform, tform2, OIsize, MmPerPixel_oi, regDir, autoTform, extractSurface);
    end

    %% scale to 1
    image4DMD = image4DMD./max(abs(image4DMD(:)));
    image4OI = image4OI./max(abs(image4OI(:)));

    %% save png images for DMD
    imgDir = fullfile('~/tmp',imgPrefix);
    mkdir(imgDir);

    save(fullfile(imgDir, [imgPrefix '_' subject]), ...
       'image4DMD','image4OI', 'camImg','niifile'); %'image4OI_all_wCCF'

    mkdir(fullfile(imgDir, subject));
    saveEveryImages(image4DMD, fullfile(imgDir, subject), binary, [imgPrefix '_' subject]);


    %% upload results to Market
    copyfile(imgDir, fullfile(resultServer, 'DMD images',imgPrefix));

    % mkdir(fullfile(resultServer, 'DMD images',imgPrefix,subject));
    % copyfile(fullfile(fullfile(imgDir, subject)), fullfile(resultServer, 'DMD images',imgPrefix,subject));

    clear image4OI image4DMD image4DMD_ori image4OI_ori
end