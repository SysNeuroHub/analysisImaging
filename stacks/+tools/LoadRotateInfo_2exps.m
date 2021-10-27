function [t_concord, fliplr_input, fliplr_base, flipud_input, flipud_base, ...
    inputImage, inputImage_registered, baseImage] = ...
    LoadRotateInfo_2exps( ServerDir, p, ResizeFactors, method)
% [t_concord, fliplr_input, fliplr_base, flipud_input, flipud_base, ...
%  inputImage, inputImage_registered, baseImage] = ...
%    LoadRotateInfo_2exps(ServerDir, p) saves/loads rotation information between 2
% experiments, using projective method. 

%input:
%p: p-file for 2 exps (should provided as a 2D cell)
%p{1}.animal: input, p{2}.animal:base
%ResizeFactors: Resize factor for the outputs (ResizeFactor(1):input, ResizeFactor(2):base)
%if scalar, apply same resize factor for input and base images

%output:
%inputImage, baseImage: resized images, of which orientation is the same to the input data

% 2014-5-4 DS made from LoadRotateInfo
% 2015-2-25 DS for choice world data, save directory under series,   
% also make the file for base>input conversion
% 2015-10-8 DS deal with the case when resize factors are different between
%target and input

%TO DO: merge this into LoadRotateInfo.m
%TO DO: generalize for single-camera imaging (gcamp,intrinsic)
%TO DO: mark location of input/base points on the images
%TO DO: generalize suffix (not only _ratio)

useCam1 = 0;

if nargin < 4
    method = 'projective';
end

if nargin < 3
    ResizeFactors = [1 1];
end
if length(ResizeFactors) == 1
    ResizeFactors = [ResizeFactors ResizeFactors];
end

if nargin < 2
    buildRotateInfo = true;
    animal_input = [];
    animal_base = [];
else
    
    animal_input = p{1}.animal;
    animal_base = p{2}.animal;
    
    AnimalDir_base   = fullfile( ServerDir, [animal_base '_ratio']); % 13.10.18 DS
    if ischar(p{2}.iseries)
        SeriesDir_base   = fullfile(AnimalDir_base, sprintf('%s', p{2}.iseries));
        iseries_base = p{2}.iseries;
    else
        SeriesDir_base   = fullfile(AnimalDir_base, sprintf('%03d', p{2}.iseries));
        iseries_base = num2str(p{2}.iseries);
    end
    ExpDir_base      = fullfile(SeriesDir_base, sprintf('%03d', p{2}.iexp));
    if ~isdir(ExpDir_base)
        mkdir(ExpDir_base);
    end
    
    %     savename = fullfile(AnimalDir, [animal_base 'rotate_flip_cams_vs' animal_input '.mat']);%needs .mat extension for exist function
    %     savename = fullfile(ExpDir_base, [animal_base 'rotate_flip_cams_vs' animal_input '.mat']);%needs .mat extension for exist function
    sname_base = fullfile(ExpDir_base, [animal_base '_' iseries_base '_' num2str(p{2}.iexp) ...
        'rotate_flip_cams_vs' animal_input '_' num2str(p{1}.iseries) '_' num2str(p{1}.iexp)]);
    savename_ori = [sname_base '.mat'];

    savename = savename_ori;
    
    if exist( savename, 'file') 
        buildRotateInfo = false;
    else
        fprintf('%s\n was not found.\nSpecify file for image rotation/flip.\nIf Canceled, build from scratch\n',...
            sname_base);
        [savename_tra, savepath] = uigetfile(ExpDir_base);
        if savename_tra ~= 0
            buildRotateInfo = false;
            savename_tra  = fullfile(savepath, savename_tra);
            savename = savename_tra;
        else
            buildRotateInfo = true;
        end
    end
end



if ~buildRotateInfo
    fprintf('Found saved rotation info for %s\n', ...
        p{2}.animal);
    
    load(savename, 'fliplr_input', 'fliplr_base', 'flipud_input', 'flipud_base', 't_concord',...
        'input_points','base_points','inputImage','inputImage_registered','baseImage','method');
    
    if  ~exist( savename_ori, 'file')
        save(savename_ori, 'fliplr_input', 'fliplr_base', 'flipud_input', 'flipud_base', 't_concord',...
            'input_points','base_points','inputImage','inputImage_registered','baseImage','method');
        
        %should save counter part as well...
    end
    
    %     %backward compatiblity for VDAQ
    %     if ~exist('fliplrCam1','var')||~exist('flipudCam1','var')||~exist('fliplrCam2','var')||~exist('flipudCam2','var')
    %         fliplrCam1 = 'n'; %14/4/14 DS
    %         fliplrCam2 = 'n'; %14/4/14 DS
    %         flipudCam1 = 'n'; %15/4/14 DS
    %         flipudCam2 = 'n'; %15/4/14 DS
    %         display('Information on image flipping has been added.');
    %
    %         cam2 = flipud(cam2');  %15/4/14 DS
    %         cam1 = flipud(cam1'); %15/4/14 DS
    %         cam1_registered = flipud(cam1_registered');  %15/4/14 DS
    %     end
    
    
else
    fprintf('Making rotation from scratch.\n Need input at least 4 points.\n');
    
    
    AnimalDir_base   = fullfile( ServerDir, [animal_base '_ratio']); % 13.10.18 DS
    %SeriesDir_base   = fullfile(AnimalDir_base, sprintf('%03d', p{2}.iseries));
     if ischar(p{2}.iseries)
        SeriesDir_base   = fullfile(AnimalDir_base, sprintf('%s', p{2}.iseries));
        iseries_base = p{2}.iseries;
    else
        SeriesDir_base   = fullfile(AnimalDir_base, sprintf('%03d', p{2}.iseries));
        iseries_base = num2str(p{2}.iseries);
     end
     ExpDir_base      = fullfile(SeriesDir_base, sprintf('%03d', p{2}.iexp));

    if ~isdir(ExpDir_base)
        mkdir(ExpDir_base);
    end
    
    sname_base = fullfile(ExpDir_base, [animal_base '_' iseries_base '_' num2str(p{2}.iexp) ...
        'rotate_flip_cams_vs' animal_input '_' num2str(p{1}.iseries) '_' num2str(p{1}.iexp)]);
    
    savename_base = [sname_base '.mat'];%needs .mat extension for exist function
    savefigname = [sname_base '.png'];
    
    
    %% Load images...assuming the rotaion info between cam1-2 is already made for each exp.
    % input image
    ExpDir      = fullfile( ServerDir, [p{1}.animal '_ratio'], sprintf('%03d', p{1}.iseries), sprintf('%03d', p{1}.iexp));
    loadname = fullfile(ExpDir, [p{1}.animal '_' num2str(p{1}.iseries) '_' num2str(p{1}.iexp) 'rotate_flip_cams.mat']);
    %loadname = sprintf('%s/%s_ratio/%srotate_flip_cams',ServerDir,animal_input,animal_input);
    
    load(loadname, 'cam1_registered', 'cam2');
    %cam1_input = cam1_registered;
    inputImage = cam2;
    if useCam1 %23/7/2014
        inputImage = cam1_registered;
    end

    % base image
    if ischar(p{2}.iseries)
        ExpDir      = fullfile( ServerDir, [p{2}.animal '_ratio'], sprintf('%s', p{2}.iseries), sprintf('%03d', p{2}.iexp));
    else
        ExpDir      = fullfile( ServerDir, [p{2}.animal '_ratio'], sprintf('%03d', p{2}.iseries), sprintf('%03d', p{2}.iexp));
    end
    loadname = fullfile(ExpDir, [p{2}.animal '_' iseries_base '_' num2str(p{2}.iexp) 'rotate_flip_cams.mat']);
    %loadname = sprintf('%s/%s_ratio/%srotate_flip_cams',ServerDir,animal_base,animal_base);
    
    load(loadname, 'cam1_registered', 'cam2');
    %cam1_base = cam1_registered;
    baseImage = cam2;
    if useCam1 %23/7/2014   
        baseImage = cam1_registered;
    end
    
    [t_concord, input_points, base_points, fliplr_input, fliplr_base, flipud_input, flipud_base, inputImage_registered, fig] = ...
        tools.rotate_flip_cams(inputImage, baseImage, method);
    
    
    save(savename_base, 'fliplr_input', 'fliplr_base', 'flipud_input', 'flipud_base', 't_concord',...
        'input_points','base_points','inputImage','inputImage_registered','baseImage','method');
    
    print(fig,savefigname,'-dpng');
    %close all;

    %% save counter part .. 25/2/15
    fliplr_base_rev = fliplr_input;
    fliplr_input_rev = fliplr_base;
    flipud_base_rev = flipud_input;
    flipud_input_rev = flipud_base;
    base_points_rev = input_points;
    input_points_rev = base_points;
    baseImage_rev = inputImage;
    inputImage_rev = baseImage;
    
    t_concord_rev = cp2tform(input_points_rev,base_points_rev, method);
    
    if strcmp(fliplr_input_rev, 'y')
        inputImage_rev = fliplr(inputImage_rev);
    end
    if strcmp(flipud_input_rev, 'y')
        inputImage_rev = flipud(inputImage_rev);
    end
    
    inputImage_registered_rev = ...
        imtransform(inputImage_rev,t_concord_rev,'XData',[1 size(baseImage_rev,2)], 'YData',[1 size(baseImage_rev,1)]);

    fliplr_input = fliplr_input_rev;
    fliplr_base = fliplr_base_rev;
    flipud_input = flipud_input_rev;
    flipud_base = flipud_base_rev;
    input_points = input_points_rev;
    base_points = base_points_rev;
    inputImage = inputImage_rev;
    baseImage = baseImage_rev;
    
    t_concord = t_concord_rev;
    inputImage_registered = inputImage_registered_rev;
    
    AnimalDir_input   = fullfile( ServerDir, [animal_input '_ratio']); % 13.10.18 DS
    SeriesDir_input   = fullfile(AnimalDir_input, sprintf('%03d', p{1}.iseries));
    ExpDir_input      = fullfile(SeriesDir_input, sprintf('%03d', p{1}.iexp));

    if ~isdir(ExpDir_input)
        mkdir(ExpDir_input);
    end

    sname_input = fullfile(ExpDir_input, [animal_input '_' num2str(p{1}.iseries) '_' num2str(p{1}.iexp) ...
        'rotate_flip_cams_vs' animal_base '_' iseries_base '_' num2str(p{2}.iexp)]);
    
    savename_input = [sname_input '.mat'];%needs .mat extension for exist function

    save(savename_input, 'fliplr_input', 'fliplr_base', 'flipud_input', 'flipud_base', 't_concord',...
        'input_points','base_points','inputImage','inputImage_registered','baseImage','method');

end


%adjust images according to ResizeFactor
if min(ResizeFactors) < 1
    %recalculate rotation info
    t_concord = cp2tform(ResizeFactors(1)*input_points, ResizeFactors(2)*base_points, ...
        method);
    
    inputImagex = ceil( size(inputImage,2) * ResizeFactors(1));
    inputImagey = ceil( size(inputImage,1) * ResizeFactors(1));
    inputImage = imresize(inputImage, [inputImagey, inputImagex]);
    
    inputImagerx = ceil( size(inputImage_registered,2) * ResizeFactors(1));
    inputImagery = ceil( size(inputImage_registered,1) * ResizeFactors(1));
    %inputImage_registered = imresize(inputImage_registered, [inputImagery,
    %inputImagerx]); %this is wrong!
    
    
    baseImagex = ceil( size(baseImage,2) * ResizeFactors(2));
    baseImagey = ceil( size(baseImage,1) * ResizeFactors(2));
    baseImage = imresize(baseImage, [baseImagey, baseImagex]);

    %for debugging
    inputImage_cache = inputImage;
    if strcmp(fliplr_input, 'y')
        inputImage_cache = fliplr(inputImage_cache);
    end
    if strcmp(flipud_input, 'y')
        inputImage_cache = flipud(inputImage_cache);
    end
    
    inputImage_registered = ...
        imtransform(inputImage_cache, t_concord,'XData',[1 size(baseImage,2)], ...
        'YData',[1 size(baseImage,1)]);

end



