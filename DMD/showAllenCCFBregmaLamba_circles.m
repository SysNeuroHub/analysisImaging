addpath(genpath('C:\Users\dshi0006\allenCCF'));
addpath(genpath('C:\Users\dshi0006\git\analysisImaging'));

%load('M:\Subjects\himiko\2025-01-23_1\dataSummary_amber.mat', 'dataSummary');
%brainImage = dataSummary.meanImage;

binarise = 1;
th_coverage = 0.5; % minimu ratio of pixels that are actually on the cortex compared to the desired

%% camera image
binning = 2; % scale = 1/binning
width = 1168/binning;
height = 900/binning;
bregma = [360/binning width/binning+1];
lambda = [805/binning width/binning+1];
MmPerPixel = 0.0104 * binning; %measured w scale 27/1/25 from getMmPerPix.m

camImg = cameraImageInfo([height width], bregma, lambda, MmPerPixel);


%% stimulus position
xfrombregma = -4.5:1:-0.5; %[mm]
yfrombregma = -4:1:3; %A>0, P<0
radiusmm = MmPerPixel*5; %[mm];


%% save CCF image with bregma, lambda and scale
saveName = ['CCFBL_' num2str(camImg.imageSize(2)) 'x' num2str(camImg.imageSize(1)) 'pix_' num2str(numel(xfrombregma)) 'x' num2str(numel(yfrombregma)) 'circle_bk'];
saveName_s = [saveName '_stereo'];
saveServer = '~/tmp';
saveDir = fullfile(saveServer, saveName);
mkdir(saveDir);
% mkdir(fullfile(saveDir, 'stereo'));
% screen2png(fullfile(saveDir, [saveName '_stereo']), f);
% close(f);

f=figure;
image(zeros(camImg.imageSize));
addAllenCtxOutlines(camImg.bregmapix, camImg.lambdapix, 'w', MmPerPixel);%this looks at lambda and shrinks the CCF
scatter(camImg.bregmapix(2), camImg.bregmapix(1), 600/binning,'markerEdgecolor','w','MarkerFaceColor','w','marker','x');
scatter(camImg.lambdapix(2), camImg.lambdapix(1), 600/binning,'markerEdgecolor','w','MarkerFaceColor','w','marker','x');
line([camImg.imageSize(2)/2-50/binning camImg.imageSize(2)/2-50/binning+1/MmPerPixel], [800 800]/binning,'linewidth',2,'color','w');
exportPng4DMD(fullfile(saveDir, [saveName_s '_ref']), f, binarise);
close(f);

%% superimpose grid on CCF
nPatches = numel(xfrombregma)*numel(yfrombregma);
disp([num2str(nPatches) ' patches'])

xfrombregmapix = 1/MmPerPixel * xfrombregma + camImg.bregmapix(2);
yfrombregmapix = -1/MmPerPixel * yfrombregma + camImg.bregmapix(1);
radiuspix = 1/MmPerPixel * radiusmm;

% grids on CCF ... sometimes it vanishes
% f=figure;
% image(zeros(camImg.imageSize));
% addAllenCtxOutlines(camImg.bregmapix, camImg.lambdapix, 'w', MmPerPixel);%this looks at lambda and shrinks the CCF
% hline(yfrombregmapix, gca,'-','w');
% vline(xfrombregmapix, gca,'-','w');
% exportPng4DMD([saveName '_grid'], f, 1);


%% check if stimulation is within the cortex
f=figure;
image(zeros(camImg.imageSize));
addAllenCtxOutlines(camImg.bregmapix, camImg.lambdapix, 'w', MmPerPixel);%this looks at lambda and shrinks the CCF
exportPng4DMD(fullfile(saveDir, 'test'), f, 1);close(f);
l=imread(fullfile(saveDir, 'test.png'));
delete(fullfile(saveDir, 'test.png'));
lf = imfill(l,'holes');


%% show only patches along the grid
set(0, 'DefaultFigureVisible', 'off');
mkdir(fullfile(saveDir,'stereo'));

patchNumber = 1;
imageName = cell(1);
position = [];

%% initial image = blank
fpatch=figure;
image(zeros(camImg.imageSize));%colormap(gray);

thisName = [num2str(patchNumber) ];
exportPng4DMD(fullfile(saveDir, 'stereo', [thisName '_' saveName_s]), fpatch, binarise);
close(fpatch);

imageStereo = imread(fullfile(saveDir, 'stereo', [thisName '_' saveName_s '.png']));
imageName{patchNumber} = thisName;
position(1,1:4) = nan;

for yy = 1:numel(yfrombregma)
    for xx = 1:numel(xfrombregma)
        for rr = 1:numel(radiusmm)

            disp([num2str(patchNumber) '/' num2str(nPatches+1)]);


            xval = xfrombregmapix(xx);
            yval = yfrombregmapix(yy);
            rval = radiuspix(rr);
            thisPosition = [xval-rval yval-rval 2*rval 2*rval];

            fpatch=figure;
            image(zeros(camImg.imageSize));%colormap(gray);

            hold on;
            %         viscircles([xfrombregmapix(xx) yfrombregmapix(yy)], radiuspix(rr), 'color','w');
            rectangle('position',thisPosition ,'curvature',[1 1], 'facecolor','w');



            exportPng4DMD(fullfile(saveDir, 'stereo', 'tmp'), fpatch, binarise);
            ltmp=imread(fullfile(saveDir, 'stereo','tmp.png'));

            %% check if the projected image is
            if sum(ltmp.*lf) / sum(ltmp) < th_coverage continue; end


            patchNumber = patchNumber + 1;
            thisName = [num2str(patchNumber) ];
            exportPng4DMD(fullfile(saveDir, 'stereo', [thisName '_' saveName_s]), fpatch, binarise);
            close(fpatch);

            %% load
            imageLoaded = imread(fullfile(saveDir, 'stereo', [thisName '_' saveName_s '.png']));

            imageStereo = cat(3, imageStereo, imageLoaded);
            imageName{patchNumber} = thisName;
            position(patchNumber, :) = thisPosition;
        end
    end
end

%% save everything in one file
save(fullfile(saveDir, saveName_s), 'imageStereo','camImg','position');

delete(fullfile(saveDir, 'stereo','tmp.png'));
set(0, 'DefaultFigureVisible', 'on');

%% stimulus positions w CCF
load(fullfile(saveDir, saveName_s), 'imageStereo');
fpatch=figure;
imagesc(sum(imageStereo,3));
hold on;
addAllenCtxOutlines(camImg.bregmapix, camImg.lambdapix, 'w', MmPerPixel);%this looks at lambda and shrinks the CCF
exportPng4DMD(fullfile(saveDir, [saveName_s '_all_wCCF' ]), fpatch, binarise);close(fpatch);
                

