function [ratioSequence, heartBeatFreq, slaveCorrection, masterCorrection,mask, sumSequence] = ...
    makeratio(masterRAW,slaveRAW,HBmin,HBmax,frameT, firstFrame,lastFrame,...
    firstFrameBase,lastFrameBase,mask,suffix,filterSizeHB,doDifferential)
% returns sequence of ratiometric, from master and slave sequences (F
% instead of dF/F).
%
% [ratioSequence] = makeratio(masterRAW, slaveRAW, HBmin, HBmax, frameT)
% returns sequence of ratiometric, corrected for each pixel.
% INPUTS:
%   masterRAW, slaveRAW; x-y-t tensor after rotation/flip, after subtracting baseimage,
%   size of these matrix must be identical
%   HBmin, HBmax: frequency range at which the "noise" is calculated
%   frameT: interval between frames in ms.
%
% [...] = makeratio(masterRAW, slaveRAW, HBmin, HBmax, frameT, firstFrame,...
% lastFrame,firstFrameBase,lastFrameBase)
% lets you specify period of sequence
% INPUTS:
%   firstFrame,lastFrame: first and last frame number to calculate the
%   ratio (DEFAULT: 1, size(slaveRAW,3)).
%   firstFrameBase,lastFrameBase: first and last frame number to calculate
%   the noise-gain correction factor (DEFAULT: 1, size(slaveRAW,3)).
%
% [...] = makeratio(masterRAW, slaveRAW, HBmin, HBmax, frameT, firstFrame,...
% lastFrame,firstFrameBase,lastFrameBase, mask)
% lets you specify ROI to calculate the noise-gain equalization factor 
%   empty > define from the image via GUI
%   scalar > calculate the factor for each pixel
%
% [ratioSequence, heartBeatFreq, slaveCorrection, masterCorrection,mask] = makeratio(...)
% also returns relating information on the noise-gain equalization
% OUTPUTS:
%   heartBeatFrequency: estimated frequency of "noise", at which the
%   correction factor is calculated [Hz]
%   slaveCorrection, masterCorrection: noise-gain equalization factors. See
%   Akemann 2012 (JNP) for details.
%   mask: roi used to calculate the noise-gain equalization factor
% 
% See also: fret.makeRatioStack, fret.slavemasterAdjustFactor(2D),
% fret.SequenceAdjust, tools.differentialSequence

% 2014-10-23 DS added 13th input
% 2014-10-24 DS added 6th output

if nargin < 13
    doDifferential = true;
end
if nargin < 12
    filterSizeHB = [];
end
if nargin < 11
    suffix = [];
end
if nargin < 10
    mask = 1;
end
if nargin < 9
    lastFrameBase = size(slaveRAW,3);
end
if nargin < 8
    firstFrameBase = 1;
end
if nargin < 7
    lastFrame = size(slaveRAW,3);
end
if nargin < 6
    firstFrame = 1;
end

if isempty(filterSizeHB)
    filterSizeHB = 5;      %2014-1-7 DS
end

if isempty(mask)
    disp('*** Define area to be used for evaluating the correction factor ***');
    mask = tools.makeMask(masterRAW(:,:, firstFrameBase));
    pause(0.05);
    close all;
end

slaveRAWMasked = reshape(slaveRAW(:,:,firstFrameBase:lastFrameBase), ...
    size(slaveRAW,1)*size(slaveRAW,2), lastFrameBase-firstFrameBase+1);

if length(mask) > 1
    mask2 = find(mask);
    slaveRAWMasked = slaveRAWMasked(mask2,:);
end

slaveRAWMasked = slaveRAWMasked - (meshgrid(mean(slaveRAWMasked,2),ones(size(slaveRAWMasked,2),1)))';

L = lastFrameBase - firstFrameBase;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
faxis = 1000/frameT/2*linspace(0,1,NFFT/2+1);
Y = fft(slaveRAWMasked, NFFT, 2)/L;
amp = abs(Y(:,1:NFFT/2+1));
meanLogAmp = mean(log10(amp),1);
idx = intersect(find(faxis >= HBmin), find(faxis <= HBmax));
[dum, maxidx] = max(meanLogAmp(idx));
heartBeatFreq = faxis(maxidx + min(idx) - 1);
if isempty(heartBeatFreq) %2015/2/23
    heartBeatFreq = (HBmin+HBmax)/2;
end
disp(['Target freq = ', num2str(heartBeatFreq)]);
clear amp slaveRAWMasked Y

if length(mask)==1
    disp('Calculating correction factors for each pixel...');
 
      
    [ slaveCorrection, masterCorrection ] = ...
        fret.slaveMasterAdjustFactor2D( slaveRAW, masterRAW, ...
        [heartBeatFreq-1, heartBeatFreq+1], 1/(frameT/1000), filterSizeHB, ...
        firstFrameBase, lastFrameBase);

      
    %gain equalized sequences
    disp('Adjusting sequence for each pixel...')
    masterRAW = fret.sequenceAdjust( masterRAW, masterCorrection, firstFrameBase, lastFrameBase, firstFrame, lastFrame);
    slaveRAW = fret.sequenceAdjust( slaveRAW, slaveCorrection, firstFrameBase, lastFrameBase, firstFrame, lastFrame);

    

else
    if strcmp(suffix, '_slowGain')%want to combine these two into one
        [ slaveCorrection, masterCorrection ] = ...
            fret.slaveMasterAdjustFactor_slowGain( slaveRAW, masterRAW, heartBeatFreq, 1/(frameT/1000), mask, firstFrameBase, lastFrameBase, firstFrame, lastFrame);
    elseif isempty(suffix)
        [ slaveCorrection, masterCorrection ] = ...
            fret.slaveMasterAdjustFactor( slaveRAW, masterRAW, [heartBeatFreq-1, heartBeatFreq+1], 1/(frameT/1000), mask, firstFrameBase, lastFrameBase, firstFrame, lastFrame);
    else
        disp('Unknown suffix')
    end
    
    disp(['Proposed master correction factor = ', num2str(masterCorrection)]);
    disp(['Proposed slave correction factor = ', num2str(slaveCorrection)]);
    
    %gain equalized
    masterRAW = fret.sequenceAdjust( masterRAW, masterCorrection, firstFrameBase, lastFrameBase, firstFrame, lastFrame);
    slaveRAW = fret.sequenceAdjust( slaveRAW, slaveCorrection, firstFrameBase, lastFrameBase, firstFrame, lastFrame);
    
%     mmaster = squeeze(mean(mean(masterRAW(110:130, 130:150,:))));
%     mslave = squeeze(mean(mean(slaveRAW(110:130, 130:150,:))));
%     plot(mmaster, 'r');hold on
%     plot(mslave + mean(mmaster-mslave),'g');


end

%spatial filtering before the differentiation...slighty smaller vessel artifact than
%filtering after the differentiation


% if doSF
%     filterSize = 3;
%     disp(['(Master) Doing spatial filtering']);
%     masterRAW = averageFilter2(masterRAW, filterSize, 2);
%     
%     disp(['(Slave) Doing spatial filtering']);
%     slaveRAW = averageFilter2(slaveRAW, filterSize, 2);
% end

ratioSequence = fret.ratioOfSequences(slaveRAW, masterRAW, firstFrame, lastFrame);

%spatial filtering before the differentiation...slighty smaller vessel artifact than
%filtering after the differentiation
% if doSF
%     filterSize = 3;
%     disp(['(Ratio) Doing spatial filtering']);
%     ratioSequence = averageFilter2(ratioSequence, filterSize, 2);
% end

%dR/R
% this should come before taking ratio?? 15/2/16 
if doDifferential
    ratioSequence = tools.differentialSequence(ratioSequence, firstFrameBase, lastFrameBase, firstFrame, lastFrame);
end


if nargout > 5
    masterSequence = tools.differentialSequence(masterRAW, firstFrameBase, lastFrameBase, firstFrame, lastFrame);
    slaveSequence = tools.differentialSequence(slaveRAW, firstFrameBase, lastFrameBase, firstFrame, lastFrame);

    sumSequence = 1/2 * (masterSequence + slaveSequence);
end

% %% test taking differential before taking the ratio
% masterSequence = tools.differentialSequence(masterRAW, firstFrameBase, lastFrameBase, firstFrame, lastFrame);
% slaveSequence = tools.differentialSequence(slaveRAW, firstFrameBase, lastFrameBase, firstFrame, lastFrame);
% % 
%  mratio = squeeze(mean(mean(ratioSequence(50:80, 50:70, :))));
%  mmaster = squeeze(mean(mean(masterSequence(50:80, 50:70, :))));
%  mslave = squeeze(mean(mean(slaveSequence(50:80, 50:70, :))));
% % 
   
%    fit1 = fit(mmaster(101:200),mslave(101:200),'poly1')
%   fit2 = fit(mmaster(201:450),mslave(201:450),'poly1')
%   fit3 = fit(mmaster(451:650),mslave(451:650),'poly1')
% 
% plot(mmaster(100:200),mslave(100:200),'.');hold on
% plot(fit1,'y');hold on
% plot(mmaster(201:450),mslave(201:450),'r.');
% plot(fit2,'r');
% plot(mmaster(451:end-2),mslave(451:end-2),'c.');
% plot(fit3,'c');
% axis equal 
% grid on
% xlabel('dM/M');ylabel('dS/S');
% %legend('prestim','during stim','poststim');
% title('correction factors calculated during poststim period')

%   
%    subplot(131);
%    plot(mratio(100:200),mslave(100:200),'y.');hold on
%    plot(mratio(100:200),mmaster(100:200),'r.');
%    axis equal tight; grid on;
% 
%    
%    subplot(132);
%    plot(mratio(201:450),mslave(201:450),'y.');hold on
%    plot(mratio(201:450),mmaster(201:450),'r.');
%    axis equal tight; grid on;
%    
%    subplot(133);
%    plot(mratio(451:650),mslave(451:650),'y.');hold on
%    plot(mratio(451:650),mmaster(451:650),'r.');
%    axis equal tight; grid on;
