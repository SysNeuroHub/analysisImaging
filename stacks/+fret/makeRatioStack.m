function [S, ratioInfo, registerInfo, SUM] = makeRatioStack(SMaster, SSlave, ratioInfo, varargin)
% make stackset of ratiometric using noise-gain equalization for each
% pixel. Images from two cameras are aligned in this function, using
% stackset.align.
%
% S = makeRatioStack(SMaster, SSlave, ratioInfo)
% returns stackset of ratiometric from stacksets of master and slave.
% ratioInfo is a structure containing parameters needed to process ratiometric signal:
%   properties
%   .rotationInfo: rotation information (cp2tform) of master
%   image to match the image from slave. If empty, no rotation
%   .flipSlave: true>flip the slave image
%   .masterFactor: master correction factor for noise-gain equalization. If
%   this field does not exist, calculate correction factor from the data
%   .master(slave)BaseImage: avg master and slave for noise-gain
%   equalization. If this field does not exist, calculate baseline images from the data
%
%   .correctFreq: minimum and maximum temporal frequency [Hz] to correct the
%   hemodynamic noise
%   correctSpace: spatial filter size to calculate the pixel-wise correction
%   factor
%   .firstFrameBase,lastFrameBase: initial and last time steps to calculate
%   the correction factor. If these are left empty, first and last time
%   points will be automatically substituted.
%
% S = makeRatioStack(SMaster, SSlave, ratioInfo, registerInfo)
% does image registration before taking the ratio
% registerInfo is a structure containing information needed for image registration
% This needs to be 2x1 vector, each specifies cam1 and cam2, respectively.
% properties:
%   .registerImage: image used as target. Size should be the same as SMaster
%   and SSlave. If one of the two registerImage is empty, do registration
%   using a single channel. If both are empty, do not apply registration.
%   .roix, roiy: specifies ROI in the registerImage (after image rotation)
%
% [S, ratioInfo] = makeRatioStack(...)
% also returns ratioInfo, to which following properties are added:
% ratioInfo.hbFreq: heart beat frequency at which gain-equalization is calculated
% ratioInfo.masterFactor: master correction factor
%
% [S, ratioInfo, registerInfo] = makeRatioStack(...)
% also returns registerInfo, to which following properties are added:
% registerInfo.Dx
% registerInfo.Dy
%
% [S, ratioInfo, registerInfo] = makeRatioStack(...,'suffix', 'detrend')
% apply linear detrend using [ratioInfo.firstFrameBase ratioInfo.lastFrameBase]
% see also. tools.DetrendStack
%
% See also. StackSet.Align, fret.makeratio, tools.RapidReg
%
% 25/3/14 - DS added varargin input for image registration
% 23/10/14  DS added doDifferential input (default:true)
% 07/06/15  DS modified when ratioInfo.masterFactor along with ratioInfo.master(slave)BaseImage
% are given in input, use it for calculating ratioSequence
% 15/7/15   DS moved contents in makeratio.m into this function. makeratio.m will not
% be used anymore.

% if nargin < 4 %tis not valid
%     registerInfo = [];
% end

if ~isempty(varargin)
    vidx = [];
    for vv = 1 : length(varargin)
        if isstr(varargin{vv})
            vidx = [vidx vv];
        end
    end
    
    %information for image registration in varargin
    for vv = vidx
        if any(strfind(varargin{vv}, 'registerInfo'))
            registerInfo = varargin{vv+1};
            break
        end
    end
    
    for vv = vidx
        if any(strfind(varargin{vv}, 'doDifferential'))
            doDifferential = varargin{vv+1};
            break
        end
    end
    
    for vv = vidx
        if any(strfind(varargin{vv}, 'suffix'))
            suffix = varargin{vv+1};
            break
        end
    end
    
end

if ~exist('suffix','var')
    suffix = '';
end

if ~isempty(regexp(suffix, regexptranslate('wildcard','*detrend*'), 'once'))%2015/2/23
    doDetrend = true;
else
    doDetrend = false;
end


if exist('registerInfo','var')
    doRegister = true;
    %registerInfo = varargin{1,2}{1};
    marginSize = 0.2;
    
    if ~isfield(registerInfo, 'registerImage') || ~isfield(registerInfo, 'roix') || ~isfield(registerInfo, 'roiy')
        error('registerInfo needs registerImage, roix and roix fields')
    else
        if ~isempty(registerInfo(1).registerImage) && ~isempty(registerInfo(2).registerImage)
            %     if ~isempty(registerInfo(1).roix) && ~isempty(registerInfo(1).roiy) && ~isempty(registerInfo(2).roix) && ~isempty(registerInfo(2).roiy)
            display('register using two cams');
            registerOption = 3;
        elseif ~isempty(registerInfo(1).registerImage) && isempty(registerInfo(2).registerImage)
            display('register using cam1');
            registerOption = 1;
        elseif isempty(registerInfo(1).registerImage) && ~isempty(registerInfo(2).registerImage)
            display('register using cam2');
            registerOption = 2;
        else
            error('cannot recognize registerInfo');
        end
    end
else
    doRegister = false;
end

if ~exist('doDifferential','var')
    doDifferential = true;
end


S = StackSet; % an empty object for ratiometric
S.TimeUnits = 's';
S.SpaceUnit = 'mm';
S.nRows = SSlave.nRows;
S.nCols = SSlave.nCols;
S.PixelSize = SSlave.PixelSize;
S.Info = SSlave.Info;
S.nConds = 1; %23/10/14

%added 19/5/14 DS
S.nFrames = SSlave.nFrames;
S.TimeVec = SSlave.TimeVec;


if ~isequal(SMaster.ResizeFactor, SSlave.ResizeFactor)
    error('Resize factor for master and slave should be equal.');
else
    S.ResizeFactor = SSlave.ResizeFactor;
end

masterSize = size(SMaster.Values);
slaveSize = size(SSlave.Values);

if ~isfield(ratioInfo, 't_concord')
    ratioInfo.t_concord = [];
end

if ~isfield(ratioInfo, 'correctFreq');
    ratioInfo.correctFreq = [6 15];
end
if isempty(ratioInfo.correctFreq)
    ratioInfo.correctFreq = [6 15];
end

if ~isfield(ratioInfo, 'correctSpace')
    ratioInfo.correctSpace = 15; %pixels
end
if isempty(ratioInfo.correctSpace)
    ratioInfo.correctSpace = 15;
end

if isfield(ratioInfo, 'masterFactor')
    globalEqFactor = true; %% use gain-equalization factors from outside
else
    globalEqFactor = false;
end

if isfield(ratioInfo, 'masterBaseImage') + isfield(ratioInfo, 'slaveBaseImage') == 2
    globalBaseImage = true; %%use baseline images from outside
elseif isfield(ratioInfo, 'masterBaseImage') + isfield(ratioInfo, 'slaveBaseImage') == 1
    warning('Either ratioInfo.masterBaseImage or ratioInfo.slaveBaseImage is empty. obtain baseline image from this trial.');
    globalBaseImage = false;
else
    globalBaseImage = false;
end

if ~globalBaseImage && (~isfield(ratioInfo, 'firstFrameBase') || ~isfield(ratioInfo, 'lastFrameBase')...
        || isempty(ratioInfo.firstFrameBase) || isempty(ratioInfo.lastFrameBase))
    h1=figure;
    plot(squeeze(mean(mean(SMaster.Values(:,:,:,1)))),'.');
    title('Average of master across pixels');
    xlabel('Steps')
    
    message = sprintf('input firstFrameBase idx: ');
    ratioInfo.firstFrameBase = str2num(input(message,'s'));
    if isempty(ratioInfo.firstFrameBase)
        ratioInfo.firstFrameBase = 1;
    end
    message = sprintf('input lastFrameBase idx: ');
    ratioInfo.lastFrameBase = str2num(input(message,'s'));
    if isempty(ratioInfo.lastFrameBase)
        ratioInfo.lastFrameBase = masterSize(3);
    end
    close(h1);
end

disp('Adjusting master images for ratio...');
SMaster = Align(SMaster, ratioInfo.fliplrCam1, ratioInfo.flipudCam1, ...
    ratioInfo.t_concord, SSlave.nRows, SSlave.nCols);%22.7.2014 DS SMaster.nRows>SSlave.nRows
masterRAW = SMaster.Values;
SMaster.Values = [];

disp('Adjusting slave images for ratio...');
SSlave = Align(SSlave, ratioInfo.fliplrCam2, ratioInfo.flipudCam2);
slaveRAW = SSlave.Values;
SSlave.Values = [];


%% hack to adjust the length of the two stacks
%Is this needless or should be replaced by tools.alignTimeStamp??
nframes_m = size(masterRAW,3);
nframes_s = size(slaveRAW,3);
if ~isequal(nframes_m, nframes_s);
    disp('temporal size not consistent between master and slave')
    if nframes_s > nframes_m
        slaveRAW = slaveRAW(:,:,1:nframes_m);
        S.FrameRate = SMaster.FrameRate;
        S.nFrames = SMaster.nFrames;
        S.TimeVec = SMaster.TimeVec;
        disp('adjusted slaveRAW temporal size');
    elseif nframes_m > nframes_s
        masterRAW = masterRAW(:,:,1:nframes_s);
        S.FrameRate = SSlave.FrameRate;
        S.nFrames = SSlave.nFrames;
        S.TimeVec = SSlave.TimeVec;
        disp('adjusted masterRAW temporal size');
    end
else
    S.FrameRate = SSlave.FrameRate;
    S.nFrames = SSlave.nFrames;
    S.TimeVec = SSlave.TimeVec;
end

%% registration
if doRegister
    
    
    switch registerOption
        case 1
            registerInfo(1).roix = registerInfo(1).roix(1):min([S.nCols registerInfo(1).roix(end)]);%2014/12/1
            registerInfo(1).roiy = registerInfo(1).roiy(1):min([S.nRows registerInfo(1).roiy(end)]);%2014/12/1
            
            rec = masterRAW(registerInfo(1).roiy, registerInfo(1).roix,:);
            registerImage_cam1 = registerInfo(1).registerImage(registerInfo(1).roiy, registerInfo(1).roix);
            [Dx, Dy] = tools.RapidReg(rec, registerImage_cam1, marginSize, 100, 'nopar');
        case 2
            registerInfo(2).roix = registerInfo(2).roix(1):min([S.nCols registerInfo(2).roix(end)]);%2014/12/1
            registerInfo(2).roiy = registerInfo(2).roiy(1):min([S.nRows registerInfo(2).roiy(end)]);%2014/12/1
            
            rec = slaveRAW(registerInfo(2).roiy, registerInfo(2).roix,:);
            registerImage_cam2 = registerInfo(2).registerImage(registerInfo(2).roiy, registerInfo(2).roix);
            [Dx, Dy] = tools.RapidReg(rec, registerImage_cam2, marginSize, 100, 'nopar');
        case 3
            %        mergeStack = [masterRAW(registerInfo.roiy,registerInfo.roix,:) slaveRAW(registerInfo.roiy,registerInfo.roix,:)];
            rec1 = masterRAW(registerInfo(1).roiy, registerInfo(1).roix,:);
            rec2 = slaveRAW(registerInfo(2).roiy, registerInfo(2).roix,:);
            for tt = 1:size(rec1,3)
                rec1t = rec1(:,:,tt);
                rec1(:,:,tt) = (rec1t - mean(rec1t(:))) ./ std(rec1t(:));
                rec2t = rec2(:,:,tt);
                rec2(:,:,tt) = (rec2t - mean(rec2t(:))) ./ std(rec2t(:));
            end
            %mergeStack = masterRAW(registerInfo.roiy,registerInfo.roix,:) + slaveRAW(registerInfo.roiy,registerInfo.roix,:);
            mergeStack = rec1 + rec2;
            registerImage_merge = registerInfo(1).registerImage(registerInfo(1).roiy,registerInfo(1).roix) ...
                + registerInfo(2).registerImage(registerInfo(2).roiy,registerInfo(2).roix);
            [Dx, Dy] = tools.RapidReg(mergeStack, registerImage_merge, marginSize, 100, 'nopar');
    end
    
    if median(abs(Dx)) > 10
        warning('Registration is NG (|Dx| > 10)');
    end
    if median(abs(Dy)) > 10
        warning('Registration is NG (|Dy| > 10)');
    end
    
    
    masterRAW = tools.imageReg(masterRAW, Dx, Dy, marginSize);
    slaveRAW = tools.imageReg(slaveRAW, Dx, Dy, marginSize);
    registerInfo(1).Dx = Dx;
    registerInfo(1).Dy = Dy;
    registerInfo(2).Dx = Dx;
    registerInfo(2).Dy = Dy;
    %registerInfo.registerImage = registerImage;
else
    registerInfo = [];
end

if globalEqFactor * globalBaseImage == 0
    if ~isempty(ratioInfo.firstFrameBase)
        firstFrameBase = ratioInfo.firstFrameBase;
    else
        firstFrameBase = 1;
    end
    
    if ~isempty(ratioInfo.lastFrameBase)
        lastFrameBase = ratioInfo.lastFrameBase;
    else
        lastFrameBase = S.nFrames;
    end
else
    ratioInfo.firstFrameBase = [];
    ratioInfo.lastFrameBase = [];
end

if globalEqFactor
    masterCorrection = ratioInfo.masterFactor;
    slaveCorrection = 0.5*(1+1./(2*masterCorrection - 1));
else
    
    slaveRAWMasked = reshape(slaveRAW(:,:,firstFrameBase:lastFrameBase), ...
        size(slaveRAW,1)*size(slaveRAW,2), lastFrameBase-firstFrameBase+1);
    
    slaveRAWMasked = slaveRAWMasked - (meshgrid(mean(slaveRAWMasked,2),ones(size(slaveRAWMasked,2),1)))';
    
    frameT = 1000/S.FrameRate;
    
    L = lastFrameBase - firstFrameBase;
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    faxis = 1000/frameT/2*linspace(0,1,NFFT/2+1);
    Y = fft(slaveRAWMasked, NFFT, 2)/L;
    amp = abs(Y(:,1:NFFT/2+1));
    meanLogAmp = mean(log10(amp),1);
    idx = intersect(find(faxis >= ratioInfo.correctFreq(1)), ...
        find(faxis <= ratioInfo.correctFreq(2)));
    [dum, maxidx] = max(meanLogAmp(idx));
    hbFreq = faxis(maxidx + min(idx) - 1);
    if isempty(hbFreq) %2015/2/23
        ratioInfo.hbFreq = (ratioInfo.correctFreq(1)+ratioInfo.correctFreq(2))/2;
    else
        ratioInfo.hbFreq = hbFreq;
    end
    disp(['Target freq = ', num2str(ratioInfo.hbFreq)]);
    clear amp slaveRAWMasked Y
    
    ratioInfo.correctFreq(1) = ratioInfo.hbFreq-1;
    ratioInfo.correctFreq(2) = ratioInfo.hbFreq+1;
    
    [ slaveCorrection, masterCorrection ] = ...
        fret.slaveMasterAdjustFactor2D( slaveRAW, masterRAW, ...
        [ratioInfo.correctFreq(1), ratioInfo.correctFreq(2)], 1/(frameT/1000), ratioInfo.correctSpace, ...
        firstFrameBase, lastFrameBase);
    
    ratioInfo.masterFactor = masterCorrection;
end

if globalBaseImage
    masterBaseImage = ratioInfo.masterBaseImage;
    slaveBaseImage = ratioInfo.slaveBaseImage;
else
    masterBaseImage = squeeze(mean(masterRAW(:,:,firstFrameBase:lastFrameBase), 3));
    slaveBaseImage = squeeze(mean(slaveRAW(:,:,firstFrameBase:lastFrameBase), 3));
    ratioInfo.masterBaseImage = masterBaseImage;
    ratioInfo.slaveBaseImage = slaveBaseImage;
end

masterRAW = fret.sequenceAdjust( masterRAW, masterCorrection, [], [], 1, size(masterRAW,3), masterBaseImage);
slaveRAW = fret.sequenceAdjust( slaveRAW, slaveCorrection, [], [], 1, size(slaveRAW, 3), slaveBaseImage);


ratioSequence = fret.ratioOfSequences(slaveRAW, masterRAW);
%should be 1, right?
%     ratioSequence = tools.differentialSequence(ratioSequence, firstFrameBase, lastFrameBase);

%2015/7/17
masterEqBaseImage = masterBaseImage;
slaveEqBaseImage = slaveBaseImage;
ratioEqBaseImage = masterEqBaseImage./slaveEqBaseImage;
ratioSequence = tools.differentialSequence(ratioSequence, [], [], 1, size(ratioSequence,3), ratioEqBaseImage);

masterSequence = tools.differentialSequence(masterRAW, [], [], 1, size(masterRAW,3), masterBaseImage);
slaveSequence = tools.differentialSequence(slaveRAW, [], [], 1, size(slaveRAW,3), slaveBaseImage);

sumSequence = 1/2 * (masterSequence + slaveSequence);


S.Values = ratioSequence;

SUM = S;
SUM.Values = sumSequence;

if doDetrend %added on 27/1/15
    %remove linear trend
    SUM = tools.DetrendStacks(SUM, 1, [SUM.TimeVec(ratioInfo.firstFrameBase) SUM.TimeVec(ratioInfo.lastFrameBase)]);
    %remove linear trend
    S = tools.DetrendStacks(S, 1, [S.TimeVec(ratioInfo.firstFrameBase) S.TimeVec(ratioInfo.lastFrameBase)]);
    ratioInfo.doDetrend = true;
end



end
