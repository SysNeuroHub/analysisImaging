movDir = '/mnt/dshi0006_market/natural/nishimoto2023';
% saveFormat = 'Motion JPEG 2000'; %'MPEG-4'; NG for linux
saveEachFig = 0;

nTotFrames = 10*60*30;%original movie frame number
durPerSnippet = 4;%[s]
nTotSnippets = 10;%9000/durPerSnippet; %=2.5hours
frameRate = 15; %[Hz] original data is 30Hz
frameLength = frameRate*durPerSnippet;
dropFrames = 1; %#frames to drop per 1frame
%0: effective 15Hz
%1: effective 30Hz
%CompressRatio = 10;
spaceSampleFac = 2;

%% parameter about stereotaxic coords
scale = 1/3;
width =  400;
height = 300;
bregma = [scale*(380-20) width/2+1]+0.5; %[y x]
lambda = [scale*(825-20) width/2+1]+0.5;

brainImage = zeros(height,width);
MmPerPixel_img = 0.0104 / scale; %measured w scale 27/1/25 from getMmPerPix.m

ioriMov = 2;%0;

saveDir = '/home/daisuke/tmp/';
saveName = sprintf('natural_%dx%dpix_%d', width, height, ioriMov);
mkdir(fullfile(saveDir, saveName, 'stereo'));


tic;
for imov = 1%:nTotSnippets
imageStereo = [];

    %saveName = fullfile(movDir, saveSuffix, num2str(imov));

    % if nTotFrames <= imov*(1+dropFrames)*frameLength - (ioriMov-1)*nTotFrames
    %     ioriMov = ioriMov + 1;
        loadName = sprintf('mtrn0%02d_10min_re10sec.avi', ioriMov);
        loadData = VideoReader(fullfile(movDir, loadName));
    % end

    % w = VideoWriter(saveName, saveFormat);
    % w.LosslessCompression = false;
    % w.CompressionRatio = CompressRatio;
    % %w.Quality = 100; not available for motion jpeg 2000
    % w.FrameRate = frameRate;
    %
    % disp([num2str(ioriMov) '-' num2str(imov)])

    %open(w);
    for i = 1:frameLength *(1+dropFrames)
        thisFrame = rgb2gray(readFrame(loadData));
        if mod(i,dropFrames+1)==0
            %writeVideo(w, thisFrame);

           
            ypix = (size(thisFrame,1)-spaceSampleFac*height)/2+1:spaceSampleFac:(size(thisFrame,1)+spaceSampleFac*height)/2;
            xpix = (size(thisFrame,2)-spaceSampleFac*width)/2+1:spaceSampleFac:(size(thisFrame,2)+spaceSampleFac*width)/2;

            trimmedFrame = thisFrame(ypix, xpix);
            normFrame = double(trimmedFrame)./double(intmax('uint8'));

            imageStereo = cat(3, imageStereo, normFrame);

            if saveEachFig
                f=figure;
                f.InnerPosition = [1 1 width height];

                image(trimmedFrame);colormap(gray);

                axis ij image off
                xlim([1 size(brainImage,2)]);
                ylim([1 size(brainImage,1)]);

                ax = gca;
                ax.Position = [0 0 1 1];

                thisName = [num2str(i)];
                screen2png(fullfile(saveDir, saveName, 'stereo', thisName),f);
                close(f);
            end

        end
    end

    %play the created video
    % v = VideoReader(fullfile([saveName '.mj2']));
    % figure('position',[0 0 1920 1080]);
    %currAxes = axes;
    % while hasFrame(v)
    %     vidFrame = readFrame(v);
    %     image(vidFrame,"Parent",currAxes)
    %     currAxes.Visible = "off";
    %     pause(1/v.FrameRate)
    % end

    %     clear frames frames_m
    t = toc

    save(fullfile(saveDir, saveName, [saveName '_stereo']), 'imageStereo','bregma', 'MmPerPixel_img');
end