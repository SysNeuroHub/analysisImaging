

function pixelRFViewer_test(mimg, RF, aziAxis, elevAxis, thisPixel)
 %pixelRFViewer_test(mimg, RF, aziAxis, elevAxis, thisPixel)
 % created from pixelTuningCurveViewer

 % 2020-11-17 added 4th input
 
% function pixelTuningCurveViewer(allFrames, cLabels, timePoints)



Ypix = size(mimg,1);
Xpix = size(mimg,2);

if nargin < 4 || isempty(thisPixel)
    thisPixel = [round(Ypix/2) round(Xpix/2)];
end

% cax = autoCax(allFrames, thisPixel);

f = figure; 

set(f, 'UserData', thisPixel);
set(f, 'KeyPressFcn', @(f,k)tcViewerCallback(f, k, mimg, RF, aziAxis, elevAxis));
%set(f, 'WindowScrollWheelFcn', @(f,k)tcViewerCallbackWheel(f,k, allFrames, cLabels, timePoints));

redrawFig(mimg, RF, aziAxis, elevAxis, thisPixel);


function tcViewerCallbackClick(f, keydata, mimg, RF, aziAxis, elevAxis)


figHand = get(f, 'Parent');
ud = get(figHand, 'UserData');
% thisTimePoint = ud(1); thisPixel = ud(2:3); thisCond = ud(4); cax = ud(5:6);
thisPixel = ud(1:2);

clickX = keydata.IntersectionPoint(1);
clickY = keydata.IntersectionPoint(2);

switch get(f, 'Tag')
    case 'brainImage'
        thisPixel = round([clickY clickX]);
%     case 'traces'
%         thisTimePoint = min(length(timePoints), find(timePoints>clickX,1));
%     case 'tuningCurves'
%         dists = abs(clickX-cLabels);
%         thisCond = find(dists==min(dists),1);
        
end

set(figHand, 'UserData', thisPixel);%[thisTimePoint thisPixel thisCond cax]);
redrawFig(mimg, RF, aziAxis, elevAxis, thisPixel);



function redrawFig(mimg, RF, aziAxis, elevAxis, thisPixel)
%(allFrames, cLabels, timePoints, thisPixel, thisTimePoint, thisCond, cax)

%nConditions = size(allFrames,4);

% plot the brain image with a marker where the selected pixel is
thisAx = subplot(1,2,1); 
q = imagesc(mimg); set(q, 'HitTest', 'off');
hold on;
q = plot(thisPixel(2), thisPixel(1), 'ro'); set(q, 'HitTest', 'off');
hold off;
%caxis(cax);
caxis(prctile(mimg(:),[1 99]));
axis equal tight;
title(sprintf('pixel %d, %d selected', thisPixel(1), thisPixel(2)));
set(thisAx, 'ButtonDownFcn', @(f,k)tcViewerCallbackClick(f, k, mimg, RF, aziAxis, elevAxis));
set(thisAx, 'Tag', 'brainImage');
colormap(thisAx, 'gray');

% RF in visual field
thisAx = subplot(1,2,2); 
thisRF = squeeze(RF(thisPixel(1),thisPixel(2),:,:));
cax = prctile(abs(thisRF(:)),99);
q = imagesc(aziAxis, elevAxis, thisRF); set(q, 'HitTest', 'off');
caxis([-cax cax]);
colormap(thisAx, 'jet');colorbar;
xlabel('azimuth [deg]'); ylabel('elevation [deg]');
axis equal tight;

% function cax = autoCax(allFrames, thisPixel)
% 
% minVal = min(min(allFrames(thisPixel(1), thisPixel(2), :,:)));
% maxVal = max(max(allFrames(thisPixel(1), thisPixel(2), :,:)));
% 
% midVal = (maxVal+minVal)/2; range = (maxVal-minVal)*1.1;
% cax = [midVal-range/2 midVal+range/2];

