function [Vout, ScaleFactor2D] = HemoCorrectLocal_simple(U, V, Vaux, fS, FreqRange, pixSpace, regressFiltered)
% function [Vout, ScaleFactor] = HemoCorrectLocal(U, V, Vaux, fS, FreqRange, pixSpace)
% 
% Does local hemodynamic correction for widefield imaging, in SVD space.
%
% U and V is an SVD representation of is the neural signal you are trying to correct
%
% Vaux is the other non-neural signals that you are using to measure hemodynmaics 
% (e.g. green channel reflectance). It should be compressed by the same U.
%
% V and Vaux: nComps x time
%
% Note that if these were recorded using alternating illumination, you will need to 
% resample them using SubSampleShift
%
% fS is sampling frequency. 
%
% FreqRange is the heartbeat range, in which the correction factors are estimated. 
% Default 9 to 13.
%
% pixSpace is the spacing between the pixel subgrid on which to compute the gain
% factors. Larger means faster but possibly less accurate. (Default 3)
%
% Outputs: Vout is corrected signal
%
% 2016-09-07 ScaleV is called inside the function

% Transpose to allow for conventional nSVDs x nTimes input
V = V';
Vaux = Vaux';

if nargin<6
    regressFiltered = false;
end

if nargin<5
    FreqRange = [9 13];
end

if nargin<6
    pixSpace = 3;
end

% first subtract out means so the filters don't go nuts
zV = bsxfun(@minus, V, mean(V));
zVaux = bsxfun(@minus, Vaux, mean(Vaux));

% now filter for heart-frequency
[b, a] = butter(2,FreqRange/(fS/2));
% fV = filter(b,a,zV); %27/11/17
% fVaux = filter(b,a,zVaux); %27/11/17
fV = single(filtfilt(b,a,double(zV)));
fVaux = single(filtfilt(b,a,double(zVaux)));

U = imresize(U, 1/pixSpace);

% make the pixel subgrid and compute submatrices etc.
ySize = size(U,1); xSize = size(U,2);
Usub = reshape(U, ySize*xSize,[]); % P x S
%[pixY, pixX] = meshgrid(1:pixSpace:ySize, 1:pixSpace:xSize); 
% pixInd = sub2ind([ySize, xSize], pixY, pixX);
% Usub = Uflat(pixInd,:);


% compute single pixel time traces:
pixTrace = fV*Usub';
pixAux = fVaux*Usub';

%temporal filtering here??

% now compute regression coefficient for each pixel. Since they have mean 0, this is
% Cov(v1, v2)/var(v2), for each column

ScaleFactor = sum(pixTrace.*pixAux) ./ sum(pixAux.*pixAux);
ScaleFactor2D = reshape(ScaleFactor, [ySize xSize]);

ScaleFactor2D(isinf(ScaleFactor2D(:)))=0;
ScaleFactor2D(isnan(ScaleFactor2D(:)))=0;

% plot it
% figure;
% imagesc(1:pixSpace:xSize, 1:pixSpace:ySize,reshape(ScaleFactor, size(pixY)));
% caxis([-1 1]*max(abs(ScaleFactor(:))));
% colormap(colormap_blueblackred);
% colorbar
% title('Correction Scale factor');

% now compute the corresponding V-space transformation matrix
if regressFiltered
    newV = ScaleV(U, fVaux', ScaleFactor2D, 1);
else
    newV = ScaleV(U, zVaux', ScaleFactor2D, 1);
end

% now make the prediction
Vout = V - newV';

% Transpose to return conventional nSVs x nTimes output
Vout = Vout';

