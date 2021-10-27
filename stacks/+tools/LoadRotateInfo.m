function [t_concord, fliplrCam1, fliplrCam2, flipudCam1, flipudCam2, cam1, cam1_registered, cam2] = ...
    LoadRotateInfo( ServerDir, p, ResizeFactor)
% [t_concord, fliplrCam1, fliplrCam2, flipudCam1, flipudCam2, cam1, cam1_registered, cam2] = ...
%     LoadRotateInfo( ServerDir, p, ResizeFactor)
%loads rotation info for stackset
%if the info is not found, make a new one using rotate_flip_cams.m and save
%'rotate_flip_cams.mat' file under ExpDir.
%currently just for ratiometric signal
% [t_concord, fliplrCam1, fliplrCam2, flipudCam1, flipudCam2, cam1, cam1_registered, cam2] =
% LoadRotateInfo( ServerDir, p, ResizeFactor) loads previously stored
% rotation info
% [...] = LoadRotateInfo( ServerDir) makes the rotation info from scratch
%
% input:
%       ServerDir: directory where rotationInfo is saved, usually where the
%       stackset is saved as well
% output:
%       t_concord
%       fliplrCam1(2) : 'y' > image was flipped left-right
%       flipudCam1(2) : 'y' > image was flipped upside down
%       cam1(2) : images before flipping and rotation
%       cam1_registered : image of cam1 after flipping and rotation

% TO DO:
% bundle outputs to one structure
% merge this with LoadRotateInfo_2exps
% use trials with no locomotion

if nargin < 3
    ResizeFactor = 1;
end

if nargin < 2
    buildRotateInfo = true;
else
    
    %AnimalDir   = fullfile(ServerDir, [p.animal Cam.FileString]); % 13.10.18 DS
    AnimalDir   = fullfile( ServerDir, [p.animal '_ratio']); % 13.10.18 DS
    if isstr(p.iseries)
        SeriesDir   = fullfile(AnimalDir, sprintf('%03s', p.iseries));
    else
        SeriesDir   = fullfile(AnimalDir, sprintf('%03d', p.iseries));
    end
    ExpDir      = fullfile(SeriesDir, sprintf('%03d', p.iexp));
    
    %ExpDir   = fullfile( ServerDir, [p.animal '_ratio']); % 13.10.18 DS
    if ~isdir(AnimalDir)
        mkdir(AnimalDir);
    end
    if ~isdir(SeriesDir)
        mkdir(SeriesDir);
    end
    if ~isdir(ExpDir)
        mkdir(ExpDir);
    end
    
    savename_a = fullfile(AnimalDir, [p.animal 'rotate_flip_cams.mat']);%needs .mat extension for exist function
   if isstr(p.iseries)
    savename_s = fullfile(SeriesDir, [p.animal '_' (p.iseries) 'rotate_flip_cams.mat']);%needs .mat extension for exist function
   savename_e = fullfile(ExpDir, [p.animal '_' (p.iseries) '_' num2str(p.iexp) 'rotate_flip_cams.mat']);
   else
       savename_s = fullfile(SeriesDir, [p.animal '_' num2str(p.iseries) 'rotate_flip_cams.mat']);%needs .mat extension for exist function
       savename_e = fullfile(ExpDir, [p.animal '_' num2str(p.iseries) '_' num2str(p.iexp) 'rotate_flip_cams.mat']);
   end
    
    if exist(savename_e, 'file')
        buildRotateInfo = false;
        savename = savename_e;
        %     elseif exist(savename_s, 'file')
        %         buildRotateInfo = false;
        %         savename = savename_s;
        %     elseif exist(savename_a, 'file')
        %         buildRotateInfo = false;
        %         savename = savename_a;
    else
        fprintf('Rotation/flip info was not found.\nSpecify file for image rotation/flip.\nIf Canceled, build from scratch\n');
        [savename, savepath] = uigetfile(ExpDir);
        if savename ~= 0
            buildRotateInfo = false;
            savename  = fullfile(savepath, savename);
        else
            buildRotateInfo = true;
        end
    end
end

if ~buildRotateInfo
    fprintf('Found saved rotation/flip info %s\n', savename);
    
    load(savename, 't_concord', 'fliplrCam1', 'fliplrCam2', 'flipudCam1', 'flipudCam2',...
        'input_points','base_points','cam1','cam1_registered','cam2');
    
    if ~exist(savename_e,'file')
        save(savename_e, 't_concord', 'fliplrCam1', 'fliplrCam2', 'flipudCam1', 'flipudCam2',...
            'input_points','base_points','cam1','cam1_registered','cam2');
    end
    %copy the loaded data to the new directory??
    
    %backward compatiblity for VDAQ
    if ~exist('fliplrCam1','var')||~exist('flipudCam1','var')||~exist('fliplrCam2','var')||~exist('flipudCam2','var')
        fliplrCam1 = 'n'; %14/4/14 DS
        fliplrCam2 = 'y'; %14/4/14 DS
        flipudCam1 = 'y'; %15/4/14 DS
        flipudCam2 = 'y'; %15/4/14 DS
        display('Information on image rotation/flipping has been added.');
        
        cam2 = flipud(cam2');  %15/4/14 DS
        cam1 = flipud(cam1'); %15/4/14 DS
        cam1_registered = flipud(cam1_registered');  %15/4/14 DS
    end
    
    
else
    fprintf('Making rotation/flip info from scratch...\n');
    
    prompt={'Animal name:',...
        'Serie num',...
        'Experiment num',...
        'Stimulus num ',...
        'Trial num '...
        'firstFrameBase for image avg',...
        'lastFrameBase for image avg',...
        'camera magnification factor',...
        'camera binning factor',...
        };
    name='data to build rotation/flip info';
    numlines=1;
    %     defaultanswer={sprintf('M%6.0f_SD',str2num(datestr(date,'yymmdd'))),...
    %         '1','1','1','1','100','200','1.6','1','n'};
    answer={sprintf('M%6.0f_SD',str2num(datestr(date,'yymmdd'))),...
        '1','1','1','1','100','200','1.6','1'};
    
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    
    validanswer = 0;
    while ~validanswer
        try
            
            answer = inputdlg(prompt,name,numlines,answer,options);
            
            Exps.animal = answer{1};
            %Exps.iseries = str2num(answer{2});
            Exps.iseries = (answer{2});
            Exps.iexp = str2num(answer{3});
            
            iStim_input = str2num(answer{4});
            iRepeat_input = str2num(answer{5});
            firstFrameBase = str2num(answer{6});
            lastFrameBase = str2num(answer{7});
            Magnification = str2num(answer{8});
            HardwareBinning = str2num(answer{9});
           
            try
                protocolName = fullfile('\\zserver\Data\trodes\', Exps.animal, (Exps.iseries), num2str(Exps.iexp), 'Protocol.mat');
                load(protocolName);
                p = Protocol;
            catch err
                p = ProtocolLoadDS(Exps); %2015-09-04
            end
            rawDir = tools.getServerDirectory(p.animal);
            Cam = tools.Imager('PCO', rawDir, '',Magnification,HardwareBinning);
            
            AnimalDir   = fullfile( ServerDir, [p.animal '_ratio']); % 13.10.18 DS
            if isstr(p.iseries)
                SeriesDir   = fullfile(AnimalDir, sprintf('%s', p.iseries));
            else
                SeriesDir   = fullfile(AnimalDir, sprintf('%d', p.iseries));
            end
            ExpDir      = fullfile(SeriesDir, sprintf('%03d', p.iexp));
            if ~isdir(AnimalDir)
                mkdir(AnimalDir);
            end
            if ~isdir(SeriesDir)
                mkdir(SeriesDir);
            end
            if ~isdir(ExpDir)
                mkdir(ExpDir);
            end
            
            %     savename = fullfile(AnimalDir, [p.animal 'rotate_flip_cams.mat']);%needs .mat extension for exist function
            %     savefigname = fullfile(AnimalDir, [p.animal 'rotate_flip_cams.png']);
           if isstr(p.iseries)
               savename = fullfile(ExpDir, [p.animal '_' (p.iseries) '_' num2str(p.iexp) 'rotate_flip_cams.mat']);
               savefigname = fullfile(ExpDir, [p.animal '_' (p.iseries) '_' num2str(p.iexp) 'rotate_flip_cams.png']);
           else
               savename = fullfile(ExpDir, [p.animal '_' num2str(p.iseries) '_' num2str(p.iexp) 'rotate_flip_cams.mat']);
               savefigname = fullfile(ExpDir, [p.animal '_' num2str(p.iseries) '_' num2str(p.iexp) 'rotate_flip_cams.png']);
           end
            
            %% load cam1
            Cam.FileString = '_cam1';
            MyStack_cam1 = StackSet.LoadStacks( ServerDir, p, 1, iStim_input, iRepeat_input, Cam);
            cam1 = mean(MyStack_cam1.Values(:,:,firstFrameBase:lastFrameBase),3);
            cam1 = cam1/mean(cam1(:));
            MyStack_cam1 = [];
            
            %% load cam2
            Cam.FileString = '_cam2';
            MyStack_cam2 = StackSet.LoadStacks( ServerDir, p, 1, iStim_input, iRepeat_input, Cam);
            cam2 = mean(MyStack_cam2.Values(:,:,firstFrameBase:lastFrameBase),3);
            cam2 = cam2/mean(cam2(:));
            MyStack_cam2 = [];
            
            [t_concord, input_points, base_points, fliplrCam1, fliplrCam2, flipudCam1, flipudCam2, cam1_registered, fig] = ...
                tools.rotate_flip_cams(cam1, cam2);
            
            
            %TO DO: for choice world data, saving directory should be under series
            save(savename, 'fliplrCam1', 'fliplrCam2', 'flipudCam1', 'flipudCam2', 't_concord','p',...
                'input_points','base_points','cam1','cam1_registered','cam2',...
                'Exps','iStim_input','iRepeat_input','firstFrameBase','lastFrameBase');
            
            print(fig,savefigname,'-dpng');
            close all;
            validanswer = 1;
        catch err
            disp(err)
            validanswer = 0;
        end
    end
end

%adjust images according to ResizeFactor
if ResizeFactor < 1
    %recalculate rotation info
    t_concord = cp2tform(ResizeFactor*input_points, ResizeFactor*base_points, ...
        'nonreflective similarity');
    
    %use rotation info of the original images ... does not work properly
    % cam1_resize_registered = ...
    %     imtransform(cam1_resize,t_concord,'XData',[1 size(cam2_resize,2)], 'YData',[1 size(cam2_resize,1)]);
    cam1x = ceil( size(cam1,2) * ResizeFactor);
    cam1y = ceil( size(cam1,1) * ResizeFactor);
    cam1 = imresize(cam1, [cam1y, cam1x]);
    cam1rx = ceil( size(cam1_registered,2) * ResizeFactor);
    cam1ry = ceil( size(cam1_registered,1) * ResizeFactor);
    cam1_registered = imresize(cam1_registered, [cam1ry, cam1rx]);
    cam2x = ceil( size(cam2,2) * ResizeFactor);
    cam2y = ceil( size(cam2,1) * ResizeFactor);
    cam2 = imresize(cam2, [cam2y, cam2x]);
end