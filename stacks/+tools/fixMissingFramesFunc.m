function [Tensor, tt, tinfo] = fixMissingFramesFunc(Tensor, tt, method, ResizeFac, HardwareBinning)
% [Tensor] = fixMissingFramesFunc(Tensor, tt, method)
% returns Tensor, in which missed or skipped frames are filled with nans (method = 'nans') 
% or interpolated with adjacent frames(method = 'interp3).
% tt is the PCO camera time stamps from the memory map file (3rd output of
% LoadPCO/LoadPhotonFocus).
% 
% [...] = fixMissingFramesFunc(Tensor, tt, method, ResizeFac,
% HardwareBinning) usess timestamps (instead of tt) from image itself, only when
% ResizeFac=1 & HardwareBinning=1
% 
% [Tensor, tt] = fixMissingFramesFunc(...)
% also returns time stamps after the fixation
% 
% [Tensor, tt, tinfo] = fixMissingFramesFunc(...)
% also returns tinfo, in which
% 
%    tinfo.elim indicates whether the first(last) frame is eliminated or not
%            tinfo.elim(1)=1 > eliminate 1st frame
%            tinfo.elim(2)=1 > eliminate last frame
% 
%   tinfo.interp indicates time (in steps) when frame(s) is interpolated
% 
% See also: tools.LoadPCO, tools.LoadPhotonFocus, 
% StackSet.TrimDarkEnds (simpler method just to eliminate the last frames)
% 
% Note, this function was tested only with PCO.edge camera

%2013-12-05 DS created from fixMisingFrames.m
%2014-03-02 DS modefied. more precise estimate of time stamp
%               added tinfo output
%2014-08-10 DS added  HardwareBinning input to judge if to use tstampImage



if nargin < 5
    tstampImage = false;
end
if nargin < 4
    ResizeFac = 0;
end
if ResizeFac == 1 && HardwareBinning == 1;
    tstampImage = true;
else 
    tstampImage = false;
end

if size(tt,1) < size(tt,2)
    tt = tt';
end

FrameInterval = median(diff(tt));
[nr nc nt] = size(Tensor);


%% Substituting erraneous time stamp (M140527_SD)
%read time-stamp from image (only when binning factor and resize factor = 1)
tIdx = [];
if tstampImage 
    %time stamp should be saved as 8 bits
    timeRecord = Tensor(1:14,1,:); %this is valid only when tensor was not flipped/rotated
    [pixelIdx tIdx] = find(timeRecord > 255);
end

rr_error_center = unique(tIdx);
rr_error_center = setdiff(rr_error_center, [1 2 length(tt)- 1 length(tt)]);%eliminate the initial and end frames


if ~isempty(rr_error_center)
    fprintf('Fixing erraneous time stamp and tensor (mmap)\n');
   for ii = 1:length(rr_error_center)
       tt(rr_error_center(ii)) = tt(rr_error_center(ii)-1) + FrameInterval;
   end   
   
   [xi yi zi] = meshgrid(1:nc,1:nr,linspace(1,2,3));
   for ridx = 1:length(rr_error_center)
       
       disp([num2str(ridx) ' frame out of ' num2str(length(rr_error_center))]);
       iT = interp3(single(Tensor(:,:,[rr_error_center(ridx)-2 rr_error_center(ridx)+2])),xi,yi,zi);
       Tensor = cat(3, cat(3, Tensor(:,:,1:rr_error_center(ridx)-2), iT), Tensor(:,:,rr_error_center(ridx)+2:end));
   end
   
   tinfo.replace = sort([rr_error_center-1 rr_error_center rr_error_center+1]);%replace 1step before and after
end

%% Eliminating 1st/last frame (from LoadPCO)
tinfo.elim = zeros(1,2);
%if diff(tt(1:2)) > 2.0 * FrameInterval %this is not robust, when sampling rate is other than 50Hz
if diff(tt(1:2)) > FrameInterval + 0.01 %for M140808_SD_1_6 33Hz recording
    tt = tt(2:end) - tt(2) + 1.5*FrameInterval; %rough estimate. DS on 3/2/14
    Tensor = Tensor(:,:,2:end);
    tinfo.elim(1) = 1;
    fprintf('Eliminating the first frame (mmap)\n');
end

% if there is zero gap between last two frames, drop the last
if diff(tt(end-1:end)) < 0.5 *FrameInterval
    tt = tt(1:end-1);
    Tensor = Tensor(:,:,1:end-1);
    tinfo.elim(2) = 1;
    fprintf('Eliminating the last frame (mmap)\n');
end


%% interpolation
tinfo.interp = [];

new_tt = tt; % figure; clf; plot(diff(new_tt),'ko-')
%th = median(diff(new_tt)) + single(1*std(diff(new_tt)));%DS on 19/3/14. 3>1.
th = 1.5 * median(diff(new_tt)); %DS on 8/4/14...assuming outliers are double of 1/framerate
if (th-median(diff(new_tt)))<median(diff(new_tt))/100,
    th = median(diff(new_tt))+median(diff(new_tt))/100; %AP to avoid infinite loop
end
[r c v] = find(diff(new_tt)>th);

% NOTE: linspace does not work!
% uncomment below if interpolation is done later
% [x y z] = meshgrid(1:S.nCols,1:S.nRows,1:2);

if ~isempty(r), fprintf('Fixing missing frames (mmap)\n'); end

while ~isempty(r)
    % find how many frames are missing in that gap
    gapSize  = round( (new_tt(r(1)+1) - new_tt(r(1)))/median(diff(new_tt)) )-1;
    % create a vector of missing times
    missVect = new_tt(r(1)) + [1:gapSize].*median(diff(new_tt));
    % add the time vector
    new_tt   = [new_tt(1:r(1)); missVect(:); new_tt(r(1)+1:end)];
    tinfo.interp = [tinfo.interp missVect];
    
    if strcmp(method, 'nan')
        % add NaN-frames in the gap.
        missTensor = ones(nr,nc,gapSize).*NaN;%uint16(ones(nr,nc,gapSize).*NaN);
        Tensor = cat(3,cat(3,Tensor(:,:,1:r(1)), missTensor),Tensor(:,:,r(1)+1:end));
        display('NANs were interpolated to the missed Frames (mmap)');
        
    elseif strcmp(method, 'interp3')
        %Use this instead for tensor interpolation
        [xi yi zi] = meshgrid(1:nc,1:nr,linspace(1,2,length(missVect)+2));
        iT = interp3(single(Tensor(:,:,r(1):r(1)+1)),xi,yi,zi);%28/5/14 DS changed to single
        iT = iT(:,:,2:end-1); % take only the interpolated ones
        Tensor = cat(3,cat(3,Tensor(:,:,1:r(1)), iT), Tensor(:,:,r(1)+1:end));
        display('Missed Frames were interpolated with interp3  (mmap)');
    end
    
    % check if there are still gaps with missing frames
    [r c v]  = find(diff(new_tt) > th);
end

tt = new_tt;

%back to index
[~, tinfo.interp] = intersect(tt, tinfo.interp); %1/4/14

end



