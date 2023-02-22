classdef StackSet
    % STACKSET a set of image stacks
    %
    % Methods:
    % AssignCondition(S,Values, conditions)
    % AverageConditions( S, weight )
    % Duration( S )
    % GetOneCondition( S, iCond )
    % PlayCondition( S, iCond, clim, overlayImage, cmap, roi)
    % Resize( S, scale )
    % Crop( S, RowFromTo, ColFromTo )
    % Trim( S, TimeClips )
    % SubtractCondition( S, jCond )
    % Differential( S, firstBaseTime, lastBaseTime)
    % FilterInTime( S, LoCutFreq, HiCutFreq )
    % FilterInSpace( S, LoCutFreq, HiCutFreq, DoGraphics )
    % FilterInSpace_gauss(S,sigma)
    % ExtractEvent( S, EventTime, tlims, marginFrames)
    % MeanEvent( S, EventTimes, tlims)
    % SpaceAverages(S, Conditions, TimeFrom, TimeTo, roi)
    % TimeAverages(S, Flag)
    % CrossCorr(S, Trace, iCond, nlags)
    % FourierMaps(S, ff)
    % ZScore(S, firstBaseTime, lastBaseTime, meanValue, stdValue)
    % Shift_Regress(S,predictor,t1,t2,down_sample,nodetrend_flag)
    % Align(S, doFliplr,doFlipud,t_concord, nRows, nCols)
    % InterpolateNan(S)
    % Register(S, Dy, Dx, RowFromTo, ColFromTo, targetImage, marginSize, precision)
    % WavletCoef(S, freq, nCycles)
    % SubtractBase(S, firstBaseTime, lastBaseTime)
    % Resample(S,NewTimeVec)
    % PixelTuningCurve(S, cLabels)
    % FillDims(S)
    % DerRec(S)
    
    %
    % 2013-04 Matteo Carandini
    % 2013-09 Xudong, Daisuke and Matteo
    
    %TODO:
    %loads of operations can be done with bsxfun/arrayfun
    
    properties
        nRows
        nCols
        nFrames
        TimeVec
        nConds
        Values      % the values in single precision, an array nRows x nCols x nFrames x nConds
        FrameRate
        PixelSize
        Info        % add Cam.FileString?
        ResizeFactor
        Description = ''
        TimeUnits = 's'
        SpaceUnit = 'mm'
        RowVec = [] % a vector of values for the rows
        ColVec = [] % a vector of values for the columns
    end
    
    methods (Static)
        
        S = LoadAVI(FileName, RowFromTo, ColFromTo, FrameFromTo);
        
        S = LoadOneStimulus(Cam, protocol, iStim, iTr, ResizeFactor, defineIndividual);
        
        S = LoadAllStimuliAverageRepeats(Cam, p, ResizeFactor);
        
        S = LoadOneStimulusAverageRepeats(Cam, p, ResizeFactor, iStim, registerImage);
        
        S = LoadAllStimuliSelectedRepeat(Cam, p, ResizeFactor, iStim);
        
        [S, ratioInfo, registerInfo] = LoadStacks( ServerDir, p, ResizeFactor, iStim, RepeatList, Cam, varargin);
        
        [S, ratioInfo, registerInfo, h, trace] = LoadOneStimulus_ratio(ServerDir, Cam, p, iStim, RepeatList, ResizeFactor, registerImage, ratioInfo, varargin);%DS 25/3/14
        
        [S, ratioInfo] = LoadOneStimulusAverageRepeats_ratio(ServerDir, Cam, p, ResizeFactor, iStim, ratioInfo);
        
        % to read structs that have not been saved as StackSets
        % Made by Matteo 2015-10
        function S = loadobj(S)
            if isstruct(S)
                
                newS = StackSet; % Call default constructor
                
                % Assign property values from struct
                newS.nRows = S.nRows;
                newS.nCols = S.nCols;
                newS.nFrames = S.nFrames;
                newS.TimeVec = S.TimeVec;
                newS.nConds = S.nConds;
                newS.Values = S.Values;
                newS.FrameRate = S.FrameRate;
                newS.PixelSize = S.PixelSize;
                newS.Info = S.Info;
                newS.ResizeFactor = S.ResizeFactor;
                newS.Description = S.Description;
                newS.TimeUnits = S.TimeUnits;
                newS.SpaceUnit = S.SpaceUnit;
                
                S = newS;
            end
        end
        
        
    end
    
    methods
        
        
        function S = StackSet()
        end
        
        function S = AssignCondition(S,Values, conditions)
            % S = S.AssignCondition(Values, conditions)
            % Keeps or Resizes Values according to the dimensions of the stack S,
            % and assign it to iCond of S.
            % If the incoming stack has more frames than needed, the last
            % ones are dropped. If it has fewer frames than needed, the
            % the last ones are set to NaN
            
            % convert uint16 to double
            % TO DO:
            % convert to single
            % assign in time?
            
            % 2014-12-22 DS allowed 2nd input to be a vector
            
            [nr,nc,nf] = size(Values);
            nf = min(S.nFrames,nf);
            
            for iii = 1:length(conditions)
                
                iCond = conditions(iii);
                
                if nr==S.nRows && nc==S.nCols
                    fprintf('Assigning condition %d\n', iCond);
                else
                    fprintf('Resizing from %d,%d to %d, %d and assigning condition %d\n', ...
                        nr, nc, S.nRows, S.nCols, iCond);
                end
                
                for iFrame = 1:nf
                    S.Values(:,:,iFrame,iCond) = ...
                        imresize(Values(:,:,iFrame,iii), [S.nRows, S.nCols]);
                end
                
            end
            S.nConds = size(S.Values,4);
        end
        
        function AvgStack = AverageConditions( S, weight )
            % Averages across the conditions
            %
            % AvgStack = S.AverageConditions()
            % returns the average across conditions of the StackSet
            %
            % AvgStack = S.AverageConditions(weight)
            % lets you specify weights for the averaging.
            % sum(weight) should be 1.
            
            % TO DO: allow weight in time as well
            
            AvgStack = S.GetOneCondition( 1 ); %initialize all but Values
            
            % S = S.InterpolateNan;
            
            if nargin < 2
                AvgStack.Values = mean( S.Values, 4 );
            else
                if (sum(weight) < 0.99) || (sum(weight) > 1.01)
                    error('Sum of weight needs to be 1.');
                end
                
                values_cache = zeros(S.nRows, S.nCols, S.nFrames);
                for iCond = 1:S.nConds
                    values_cache = values_cache + weight(iCond) * S.Values(:,:,:,iCond);
                end
                AvgStack.Values = values_cache;
            end
        end
        
        
        function Dur = Duration( S )
            % Dur = S.Duration()
            % calculates duration of StackSet in seconds, from StackSet.TimeVec
            
            % Dur =  S.TimeVec(end)-S.TimeVec(1); underestimates
            
            ActualFramerate = 1/median(diff(S.TimeVec));
            Dur = S.nFrames/ActualFramerate;
            
            % Dur = S.nFrames/S.FrameRate; less robust!!
        end
        
        
        
        function OneStack = GetOneCondition( S, iCond )
            % Extracts one stack from the stackset
            %
            % OneStack = S.GetOneCondition( iCond )
            
            OneStack = StackSet; % an empty object
            
            % boring: fill all properties one by one
            OneStack.nRows = S.nRows;
            OneStack.nCols = S.nCols;
            OneStack.nFrames = S.nFrames;
            OneStack.TimeVec = S.TimeVec;
            OneStack.FrameRate = S.FrameRate;
            OneStack.PixelSize = S.PixelSize;
            OneStack.Info = S.Info;
            OneStack.ResizeFactor = S.ResizeFactor;
            OneStack.Description  = S.Description;
            OneStack.TimeUnits  = S.TimeUnits;
            OneStack.SpaceUnit  = S.SpaceUnit;
            
            OneStack.nConds = 1;
            OneStack.Values = S.Values(:,:,:,iCond);
            
        end
        
            
        function PlayCondition( S, iCond, clim, overlayImage, cmap, roi, TitleString)
            % view and record one of the conditions as a movie.
            % Left click on the movie specifies ROI for time traces above the movie.
            % Right click on the movie allows to export ROI
            %
            % S.PlayCondition(iCond )
            % plays a movie of iCond condition
            %
            % S.PlayCondition(iCond, clim )
            % specify color range of the movie ([min max])
            %
            % PlayCondition( S, iCond, clim, overlayImage)
            % superimpose image on top of the movie (using transparency)
            % size(overlayImage) should be [S.nRows, S.nCols]
            %
            % S.PlayCondition(iCond,[],[],[],roi)
            % imports roi ( column of 4 numbers, specifying [xbegin ybegin xlength ylength])
            %
            % see also: \\zserver\Code\MovieInspect
            
            % 2014-08-28 DS deleted output
            % 2015-04-10 DS added roi input
            
            %addpath('\\zserver.cortexlab.net\Code\MovieInspect');
            
            % close all
            
%             if iCond<0 || iCond>S.nConds
%                 disp('This condition does not exist');
%                 return
%             end
            if nargin < 7
                TitleString = sprintf('Condition %d',iCond);
            end
            if nargin < 6
                roi = [];
            end
            if nargin < 5
                cmap = [];
            end
            
            if nargin < 4
                overlayImage = [];
            end
            
            Stmp = S;
            Stmp.nConds = numel(iCond);
            Stmp.nCols = Stmp.nCols*Stmp.nConds;
            Stmp.Values = permute(Stmp.Values(:,:,:,iCond), [1 2 4 3]);
            Stmp.Values = reshape(Stmp.Values,Stmp.nRows, ...
                Stmp.nCols, Stmp.nFrames);

            roi_tmp = roi;
            
            if nargin < 3
                % establish a good range
                MedianImg = squeeze(nanmedian(Stmp.Values,3));
                m = nanmedian(MedianImg(:));
                StdImg = squeeze(nanstd(Stmp.Values(:,:,2:end-1),[],3));
                StdImg(isinf(StdImg)) = NaN;
                s = nanmean(StdImg(:));
                clim = [ m-2*s m+2*s ];
                if diff(clim)<=0, error('Cannot establish a range'); end
            end
            
            MovieInspector( Stmp.Values, roi_tmp, Stmp.FrameRate,...
                clim, TitleString, Stmp.TimeVec, overlayImage, cmap);
        end
        
        function S2 = Resize( S, scale )
            % resizes in space (useful for downsampling)
            %
            % S2 = S.Resize(SCALE)
            % returns a stackset that is SCALE times the size of S
            
            % TODO change implementation of nr and nc so they are consistent
            % with assigncondition
            
            [nr, nc ] = size( imresize( S.Values(:,:,1,1), scale ) );
            
            S2 = S;
            S2.nRows = nr;
            S2.nCols = nc;
            S2.Values = zeros(nr,nc,S.nFrames,S.nConds);
            S2.ResizeFactor = S.ResizeFactor*scale;
            S2.PixelSize = S.PixelSize/scale;
            
            for iCond = 1:S.nConds
                for iFrame = 1:S.nFrames
                    S2.Values(:,:,iFrame,iCond) = ...
                        imresize( S.Values(:,:,iFrame,iCond), scale );
                end
            end
        end
        
        function S = Crop( S, RowFromTo, ColFromTo )
            % crops a StackSet in space
            %   S = S.Crop([RowFrom RowTo], [ColFrom ColTo] )
            %   S = S.Crop([], [ColFrom ColTo]) does not crop in rows
            %   S = S.Crop([RowFrom RowTo],[]) does not crop in columns
            % RowFromTo and ColFromTo in pixels
            %
            % S = S.Crop() specifies ROI to crop via GUI.
            %
            % see also: tools.makeMask
            
            if nargin < 2
                mask = tools.makeMask(S.Values(:,:,1,1));
                pause(0.05);
                close all;
                
                Rows = (find(sum(mask,2)));
                Cols = (find(sum(mask,1)));
                
                RowFromTo = [Rows(1) Rows(end)];
                ColFromTo = [Cols(1) Cols(end)];
            else
                if isempty(RowFromTo)
                    RowFromTo = [0 S.nRows+1];
                end
                if isempty(ColFromTo)
                    ColFromTo = [0 S.nCols+1];
                end
            end
            
            %2014/12/23
            RowFromTo = round(RowFromTo);
            ColFromTo = round(ColFromTo);
            
            fprintf('Cropping in space row:%d-%d column:%d-%d ...\n', ...
                RowFromTo(1), RowFromTo(2),ColFromTo(1),ColFromTo(2));
            RowFrom = max(0,RowFromTo(1));
            RowTo   = min( S.nRows+1, RowFromTo(2));
            
            ColFrom = max(0,ColFromTo(1));
            ColTo   = min( S.nCols+1, ColFromTo(2));
            
            S.Values( [1:RowFrom RowTo:end],: , :, : ) = [];
            S.Values( :, [1:ColFrom ColTo:end], :, : ) = [];
            
            S.nRows = RowTo - RowFrom - 1; %DS 31/3/14
            S.nCols = ColTo - ColFrom - 1; %DS 31/3/14
        end
        
        function S = Trim( S, TimeClips )
            % Trim trims a StackSet in time in seconds
            %
            % S = S.Trim([from to] )
            
            iFrom = find(diff(S.TimeVec >  TimeClips(1)));
            iTo   = find(diff(S.TimeVec >= TimeClips(2)));
            if ~isempty(iFrom) && ~isempty(iTo) && iFrom<=iTo
                fprintf('Trimming from time %f to time %f\n',TimeClips);
                S.TimeVec = S.TimeVec(iFrom:iTo);
                S.nFrames = length(S.TimeVec);
                S.Values = S.Values(:,:,iFrom:iTo,:);
            else
                disp('Strange parameters , will not trim...');
            end
        end
        
        function S = SubtractCondition( S, jCond )
            % S = S.SubtractCondition(iCond)
            % Subtracts Stackset.Values by iCond condition for each time
            % step
            
            % Fixed by MC on 2014-09-04: it is essential that it is jCond
            % and not iCond, or the output becomes all zeros!!!
            
            % TO DO: subtract by an image (make another method?)
            
            if jCond>0 && jCond<=S.nConds
                thisImage = S.Values(:,:,:,jCond); %31/10/17
                for iCond = 1:S.nConds
                    S.Values(:,:,:,iCond) = S.Values(:,:,:,iCond)-thisImage;
                end
            end
        end
        
        
        function S = Differential(S,firstBaseTime, lastBaseTime)
            % S = S.Differential()
            % divide Values by average values across
            % whole period for each condition.
            %
            % dF/F = S.Differential - 1
            %
            % S = S.Differential(firsstBaseTime, lastBaseTime)
            % takes differential, relative to average values between
            % firstBaseTime and lastBaseTime (in seconds).
            %
            % See also: tools.differentialSequence
            
            if nargin < 2
                firstBaseFrame = 1;
            else
                [~,firstBaseFrame] = min(abs(S.TimeVec - firstBaseTime));
            end
            
            if nargin < 3
                lastBaseFrame = S.nFrames;
            else
                [~,lastBaseFrame] = min(abs(S.TimeVec - lastBaseTime));
            end
            
            for iCond = 1:S.nConds
                S.Values(:,:,:,iCond) = tools.differentialSequence(S.Values(:,:,:,iCond), ...
                    firstBaseFrame, lastBaseFrame);
            end
        end
        
        % added by Matteo 2015-11
        function C = MeanSubtractAndDivide(S)
            % Subtract and divide each condition by its temporal mean.
            % Useful to compute the contrast of an image or dF/F in
            % imaging.
            C = S;
            TimeMeans = S.TimeAverages;
            for iCond = 1:S.nConds
                MeanStack = repmat(TimeMeans(:,:,iCond),[1 1 S.nFrames]);
                C.Values(:,:,:,iCond) = (S.Values(:,:,:,iCond) - MeanStack)./MeanStack;
            end
        end
        
        % added by Matteo 2015-11
        function [Res, Sep] = RemoveSeparable(S, showgraphics)
            % RemoveSeparable removes the spacetime separable component of the responses
            %
            % [Res, Sep] = RemoveSeparable(S, showgraphics)
            %
            
            if nargin<2, showgraphics = 0; end
            resp = reshape(S.Values,S.nRows*S.nCols,S.nFrames*S.nConds);
            [~,~,~,Model,Residual] = tools.MakeSeparable(resp,showgraphics);
            Res = S; Res.Values = reshape(Residual,Res.nRows,Res.nCols,Res.nFrames,Res.nConds);
            Sep = S; Sep.Values = reshape(   Model,Res.nRows,Res.nCols,Res.nFrames,Res.nConds);
        end
        
        function S = FilterInTime( S, LoCutFreq, HiCutFreq )
            % Bandpass filters in time
            %
            % S = FilterInTime( S, LoCutFreq, HiCutFreq )
            % Filters between LoCutFreq-HiCutFreq (in Hz).
            %
            % If you set LoCutFreq to 0 it will remove the mean.
            % Set it to -Inf if you don't want that to happen.
            
            % 2016-1-8 DS added mirror frames to minimize boundary artefact
            
            % TO DO: check unit of amplitude
            
            if LoCutFreq>HiCutFreq
                disp('You got the low cutoff and the high cutoff switched around');
                return
            end
            
            if any(isnan(S.Values))
                error('Cannot filter StackSets that contain NaNs');
            end
            
            msize = 0.5;
            mFrames = round(S.nFrames*msize);
            nFrames_c = S.nFrames + 2*mFrames;
            
            Values_c = cat(3, flipdim(S.Values(:,:,1:mFrames,:),3), S.Values);
            Values_c = cat(3, Values_c, flipdim(S.Values(:,:,end-mFrames+1:end,:),3));
            
            
            %frequencies corresponding to the elements of an fft
            %taken from \\zserver.cortexlab.net\Code\Matteobox\freq.m
            nsamples = nFrames_c;
            duration = nFrames_c/S.FrameRate;
            if mod(nsamples,2)
                %odd (the original freq by MC did always the odd case)
                ff = [0:floor(nsamples/2) , -floor((nsamples-1)/2):1:-1]/duration;
            else
                %even
                ff = [0:(floor(nsamples/2)-1) , -floor((nsamples-1)/2)-1:1:-1]/duration;
            end
            
            
            fprintf(1, 'Filtering in time using FFT...');
            ftTensor = fft( Values_c, [], 3 );
            KillFreqs = (abs(ff)<=LoCutFreq) | (abs(ff)>=HiCutFreq);
            ftTensor(:,:,KillFreqs,:) = 0;
            Values_c = real(ifft( ftTensor, [], 3 ));
            S.Values = Values_c(:,:,mFrames+1:mFrames+S.nFrames,:);
            fprintf(1, 'done\n');
            
        end
        
        function S = FilterInSpace( S, LoCutFreq, HiCutFreq, DoGraphics )
            % Bandpass filters the stack set in space
            %
            % S.FilterInSpace( LoCutFreq, HiCutFreq ) lets you specify the
            % low and high cutoffs LoCutFreq, HiCutFreq in cycles/mm.
            % Set LoCutFreq to 0   if you don't want to low-pass filter
            % Set HiCutFreq to Inf if you don't want to high-pass filter
            
            % FilterInSpace( LoCutFreq, HiCutFreq, DoGraphics ) lets you
            % specify whether to look at the filter shape (DEFAULT: 0)
            
            %2014-11-30 DS add mirror images to minimize boundary effect
            
            %% parse the parameters
            if isnan(S.PixelSize)
                error('PixelSize = Nan. Cannot apply spatial filter');
            end
            
            if nargin < 4
                DoGraphics = 0;
            end
            
            if nargin < 2, error('Need to specify at least two inputs'); end
            
            if LoCutFreq==0 && isinf(HiCutFreq)
                return
            end
            
            if LoCutFreq>HiCutFreq,
                warning('You got the low cutoff and the high cutoff switched around');
            end
            
            if LoCutFreq == 0
                fLo = eps;
            else
                fLo = (2*S.PixelSize)*LoCutFreq; %convert from cycles/mm to cycles/pixel
            end
            
            if isinf(HiCutFreq)
                fHi = 1-eps;
            else
                fHi = (2*S.PixelSize)*HiCutFreq; %convert from cycles/mm to cycles/pixel
            end
            
            if fHi<=fLo
                error('Highpass and lowpass specs do not make sense - not doing highpass');
            end
            if fHi>1 || fLo>1
                error('Cutoff frequencies must be below Nyquist frequency');
            end
            
            % add mirror images to avoid fft artifact ripple - 2014.1.1 DS
            msize = 0.5;
            xmsize = ceil(msize*S.nCols)-1; %margin size in x-axis
            ymsize = ceil(msize*S.nRows)-1; %margin size in y-axis
            
            
            %% Create the filter
            
            ntap = min(S.nRows + 2*(ymsize+1), S.nCols + 2*(xmsize+1))-1;
            if rem(ntap,2)
                ntap = ntap - 1;
            end
            if ntap == 0
                error('Something funny is going on');
            end
            
            BPfilter1D = fir1(ntap,[fLo fHi]);
            TheFilter = ftrans2(BPfilter1D);
            
            % show the filter and its frequency response
            if DoGraphics
                [hBP,w] = freqz(BPfilter1D,1,512);
                figure; clf
                subplot(2,1,1);
                plot(w/pi /2 /S.PixelSize,abs(hBP),'k'); hold on
                plot(fLo /2 /S.PixelSize*[1 1], [0 1], 'b:');
                plot(fHi /2 /S.PixelSize*[1 1], [0 1], 'b:');
                xlabel('Frequency (cycles/mm)');
                subplot(2,1,2);
                plot(linspace(-ntap/2,ntap/2,ntap+1)*S.PixelSize,BPfilter1D,'ko-', 'markerfacecolor','k')
                xlabel('Space (mm)');
                axis tight
                drawnow
            end
            
            %% Get the Fourier Transform of the filter
            
            n = size(TheFilter,1);
            n = (n-1)/2;
            PadFilter = zeros(S.nRows + 2*(ymsize+1), S.nCols + 2*(xmsize+1));
            PadFilter( round((S.nRows + 2*(ymsize+1))/2)+[-n:n], round((S.nCols + 2*(xmsize+1))/2)+[-n:n] ) = TheFilter;
            fftFilter = fft( fft( PadFilter,[], 1 ), [], 2 );
            fftFilterTensor = repmat( fftFilter, [ 1 1 S.nFrames] );
            
            %% Do the filtering (multiplication in the frequency domain)
            
            fprintf(1, 'Filtering in space ');
            
            
            
            for iCond = 1: S.nConds
                fprintf(1, '.');
                
                movie = S.Values(:,:,:,iCond);
                
                %large movie with mirror images
                lmovie = zeros(S.nRows + 2*(ymsize+1), S.nCols + 2*(xmsize+1), S.nFrames);
                for t = 1:S.nFrames
                    lmovie_c = [flipdim(movie(:,1:xmsize+1,t),2) movie(:,:,t) flipdim(movie(:,S.nCols-xmsize:end,t),2)];
                    lmovie(:,:,t) = [flipdim(lmovie_c(1:ymsize+1,:),1); lmovie_c; flipdim(lmovie_c(S.nRows-ymsize:end,:),1)];
                end
                
                fftDataTensor = fft(fft(lmovie,[],1),[],2);
                foo2 = real(ifft( ifft( fftFilterTensor .* fftDataTensor, [], 1 ), [], 2));
                lfiltered = fftshift(fftshift(foo2,1),2);
                S.Values(:,:,:,iCond) = lfiltered(ymsize+2:ymsize+S.nRows+1, xmsize+2:xmsize+S.nCols+1, :);
            end
            
            fprintf(1, 'done\n');
            
        end
        
        function [S, idx_new] = ExtractEvent(S, EventTime, tlims, marginFrames)
            % ExtractEvent(S, EventTime, tlims) extracts tensor triggered
            % by a single event at EventTime. Note that if the time points outside tlims are
            % interpolated with nans (this is not the case in StackSet.MeanEvent)
            %
            % ExtractEvent(S, EventTime, tlims, marginFrames) eliminates
            % initial and last frames
            %
            % inputs:
            % EventTime: time to trigger tensor (scalar)
            % tlims: begin and end time to extract tensor, from EventTime (in second)
            %
            % outputs:
            % idx_new: time indexes in which Values are from S. S.Values of other indexes
            % are filled with nans
            %
            % See also. S.MeanEvent
            %
            % 2015/4/2 DS created
            % 2015/4/16 DS added 4th input
            
            % alignment in time is not precise due to the use of "round"
            
            if nargin < 4
                marginFrames = 0;
            end
            
            dur =  S.TimeVec(end)-S.TimeVec(1);
            
            Frm_Tlims  = round((tlims/dur*S.nFrames)); % convert tlims to frames
            
            TimeVec_event = ( Frm_Tlims(1):Frm_Tlims(2) ) * dur/S.nFrames; % a vector of times
            
            nFrmEvent  = length(TimeVec_event);
            
            EventValues   = nan( S.nRows, S.nCols, nFrmEvent, S.nConds, 'single' );
            
            frmT     = round((EventTime/dur*S.nFrames)); % time of event in cam frames
            frmLeft  = frmT + Frm_Tlims(1);
            frmRight = frmT + Frm_Tlims(2);
            
            
            idx_original = frmLeft:frmRight;
            
            %time index outside the original tensor is trimmed
            idx_trimmed = idx_original(intersect(find(idx_original >= 1 + marginFrames), find(idx_original <= S.nFrames - marginFrames)));
            idx_new = idx_trimmed - frmLeft + 1;
            EventValues(:,:,idx_new,:) = S.Values(:,:,idx_trimmed,:);
            
            S.Values = EventValues;
            S.TimeVec = TimeVec_event;
            S.nFrames = nFrmEvent;
            S.Description = [S.Description 'triggered at time ' num2str(EventTime)];
        end
        
        
        
        function [Sev, validEventTimes] = ExtractAllEvents(S, eventTimes, tlims, marginFrames, alignMethod)
            % [Sev, validEventTimes] = ExtractAllEvents(S, eventTimes, tlims)
            % returns event triggered tensor of every specified events.
            %
            % 2017-1-30 created rom StackEvAll
            
            if nargin < 4
                marginFrames = 0;
            end
            if nargin < 5
                alignMethod = 'interp';
            end
            
            validEvent = ones(1,length(eventTimes));%
            validEventTimes = [];
            
            Sev = StackSet;
            for ev = find(validEvent > 0)
                
                %ExtractEvent extracts epoch, without aligning in time
                [EventStack, validFrameIdx] = S.ExtractEvent(eventTimes(ev), ...
                    tlims, marginFrames);
                
                if validFrameIdx(end) == S.nFrames
                    EventStack.nFrames = EventStack.nFrames - 1;
                    validFrameIdx = validFrameIdx(1:end-1);
                    EventStack.Values = EventStack.Values(:,:,1:EventStack.nFrames,:);
                    EventStack.TimeVec = EventStack.TimeVec(1:EventStack.nFrames);
                end
                invalidFrameIdx = setxor(1:EventStack.nFrames, validFrameIdx);
                
                %initialization
                if isempty(Sev.nConds)
                    Sev.TimeVec = EventStack.TimeVec;
                    Sev.FrameRate = EventStack.FrameRate;
                    Sev.nFrames = EventStack.nFrames;
                    Sev.nRows = EventStack.nRows;
                    Sev.nCols = EventStack.nCols;
                    Sev.PixelSize = EventStack.PixelSize;
                    Sev.ResizeFactor = EventStack.ResizeFactor;
                    Sev.nConds = 0;
                    
                    stackinfo_single = whos('EventStack');
                    estimatedBytes = length(find(validEvent > 0)) * stackinfo_single.bytes;
                    
                    disp(['Stackset.ExtractAllEvents: estimated stack size is ' num2str(estimatedBytes*1e-6) 'M bytes']);
                    %fuse
                    if estimatedBytes > 1e10 %10GB
                        disp(sprintf('ABORTED since stack size will be huge \n Consider using reResizing (5th input) or StackEvRepeats.m'));
                        return;
                    end
                end
                
                % missed frames is not yet implemented
                
                
                EventStack = EventStack.InterpolateNan; %20/7/2015
                
                if any(isnan(EventStack.Values(:))) %4/2/16
                    continue;
                else
                    
                    %this stuck if Stackset.Values include nans
                    [EventStack.Values, idx_align] = ...
                        tools.alignTimeStamp(EventStack.Values, EventStack.TimeVec, Sev.TimeVec, alignMethod);
                    
                    
                    EventStack.nFrames = length(Sev.TimeVec);
                    
                    if strcmp(alignMethod, 'nearest')
                        [~, invalidFrameIdx] = intersect(idx_align, invalidFrameIdx);
                    elseif strcmp(alignMethod, 'interp')
                        invalidFrames_cache = [];
                        for mmm = 1:length(invalidFrameIdx)
                            [~, timeDifIdx] = sort(abs(Sev.TimeVec - EventStack.TimeVec(invalidFrameIdx(mmm))));
                            invalidFrames_cache = [invalidFrames_cache timeDifIdx(1:2)];
                        end
                        invalidFrameIdx = unique(invalidFrames_cache);
                    end
                    validFrameIdx = setxor(1:EventStack.nFrames, invalidFrameIdx);
                    EventStack.Values(:,:,invalidFrameIdx,:) = nan; %make sure invalid frame has NAN values
                    
                    Sev.nConds = Sev.nConds + 1;
                    Sev = Sev.AssignCondition(EventStack.Values, Sev.nConds);
                    
                    validFrames_cache = ones(1, EventStack.nFrames); %modefied on 31/8/2015
                    validFrames_cache(invalidFrameIdx) =  0;  %modefied on 31/8/2015
                    
                    if ~exist('validFrames', 'var')
                        validFrames = validFrames_cache;
                    else
                        validFrames = [validFrames; validFrames_cache];
                    end
                    
                    
                    % S_validEvents = S_validEvents.AssignCondition(MyStack_event.Values, S_validEvents.nConds+1);
                    if isempty(validEventTimes)
                        validEventTimes = [eventTimes(ev)];
                    else
                        validEventTimes = [validEventTimes; eventTimes(ev)];
                    end
                    
                end
            end
        end
        
        function [Smean, Sall] = MeanEvent( S, EventTimes, tlims)
            % Computes the average event
            %
            % Smean = S.MeanEvent( EventTimes, tlims )
            % returns event-avg stackset
            %
            % [Sman, Sall] = S.MeanEvent( EventTimes, tlims )
            % also returns event-triggered stackset of every event
            %
            % inputs:
            % EventTimes: times of events (s). You can specify the times
            % individually across conditions as EventTimes{iCond}
            % tlims: beginning and end of the event-avg stackset (s).
            % Note that if The tensor are
            % not within [tt(i)-dt tt(i)+dt] for all i, the output S.Values is all 0
            %
            % returns the avg event for each condition, if EventTimes is a
            % cell.
            %
            % see also: TensorEvent
            
            % TO DO: also return SE across events
            
            %addpath('\\zserver.cortexlab.net\Code\BrainPics'); %TensorEvent
            
            
            dur =  S.TimeVec(end)-S.TimeVec(1);
            
            mEventValues = cell(S.nConds,1);
            EventValues = cell(S.nConds,1);
            for iCond = 1:S.nConds
                
                if iscell(EventTimes)
                    et_cache = EventTimes{iCond};
                else
                    et_cache = EventTimes;
                end
                [mEventValues{iCond}, tt, EventValues{iCond}] = TensorEvent( ...
                    S.Values(:,:,:,iCond), et_cache - S.TimeVec(1), tlims, dur );
            end
            Smean = S;
            Smean.Values = cat(4,mEventValues{:});
            Smean.nConds = size(Smean.Values,4);
            Smean.TimeVec = tt;
            Smean.nFrames = length(tt);
            
            Sall = S;
            Sall.Values = shiftdim(cat(4,EventValues{:}),1);
            Sall.nConds = size(Sall.Values,4);
            Sall.TimeVec = tt;
            Sall.nFrames = length(tt);
            
            function validEvent = validWindow(eventTimes,tlims,dur,nframe)
                % validEvent = validWindow(eventTimes,tlims,dur,nframe)
                % returns train of 1/0 of whether event window is within the record
                
                % eventTimes: event time(s) from onset of the recording [s]
                % tlims: time window before and after the event [begin-end] [s]
                % dur: duration of the record [s]
                % nframes: number of frames of the record
                
                %2014-12-19 DS from TensorEvent
                
                Frm_Tlims  = round((tlims/dur*nframe)); % convert tlims to frames
                
                for it = 1:length(eventTimes)
                    
                    frmT     = round((eventTimes(it)/dur*nframe)); % in cam frames
                    frmLeft  = frmT + Frm_Tlims(1);
                    frmRight = frmT + Frm_Tlims(2);
                    
                    if( frmLeft > 0 && frmRight < nframe )
                        validEvent(it)=1;
                    else
                        validEvent(it)=0;
                    end
                end
            end
            function [ EventTensor, EventTimes, allEvents] = TensorEvent( TheTensor, tt, tlims, dur, meansub )
                % TensorEvent makes an event-related tensor
                %
                % EventTensor = TensorEvent( TheTensor, tt, dt, dur, meansub ) TensorEvent cuts
                % frame sequences in tensor 'TheTensor', with dimensions (nx,ny,nframes).
                % It loops through time points (events) specified in 'tt' (secs) and
                % selects the sequence of frames between 'tt(i)-dt' and 'tt(i)+dt' (dt in
                % secs). It needs to know the duration dur (in seconds) of the tensor.
                % Returns the event-tensor 'EventTensor'. Note that if The tensor are
                % not within [tt(i)-dt tt(i)+dt] for all i, the output EventTensor is all 0
                %
                % If dt is a vector [t1 t2], it selects the sequence of frames between
                % 'tt(i)+t1' and 'tt(i)+tt2'. A classic is [ -0.1 0.4 ].
                %
                % EventTensor = TensorEvent( ...,  'meansub' )
                % subtracts the mean from each event
                %
                % [ EventTensor EventTimes] = TensorEvent( ... )
                % returns a vector of times
                %
                % [ EventTensor EventTimes allEvents] = TensorEvent( ... )
                % returns all the events
                %
                % SEE ALSO TensorMakeMovie
                
                % 2004-05 Andrea Benucci (TensorEventMake)
                % 2005-10 Matteo Carandini (TensorEvent) made it independent of protocol
                % 2006-04 MC added output EventTimes
                % 2006-04 MC added possibility of dt = [t1 t2], corrected normalization
                % 2006-07 MC added check that dur is a number
                % 2012-07 AB added allEvents as output
                % 2013-02 AB added keepAll nargout flag
                % 2014-12 DS added display output
                % 2014-12 DS split the function to TensorEvent and validWindow
                %%
                
                keepAll = 1; if nargout<3, keepAll = 0; end
                
                if nargin<4
                    error('Need to specify at least the first 4 parameters');
                end
                
                if ~isfinite(dur), error('Gotta give a number for dur'); end
                
                if nargin<5
                    meansub = 0;
                else
                    meansub = 1;
                end
                
                if isempty(tt)
                    EventTensor = 0;
                    EventTimes = [];
                    return
                end
                
                if length(tlims) == 1
                    tlims = [-tlims tlims];
                end
                
                %% Parameters
                
                nt    = length(tt);
                [ny,nx,nframe] = size(TheTensor);
                
                Frm_Tlims  = round((tlims/dur*nframe)); % convert tlims to frames
                
                EventTimes = ( Frm_Tlims(1):Frm_Tlims(2) ) * dur/nframe; % a vector of times
                
                nFrmEvent  = length(EventTimes);
                
                %% Main loop
                EventTensor   = zeros( ny,nx,nFrmEvent,'single' );
                
                cnt   = 0; allEvents = [];
                for it = 1:nt
                    frmT     = round((tt(it)/dur*nframe)); % in cam frames
                    frmLeft  = frmT + Frm_Tlims(1);
                    frmRight = frmT + Frm_Tlims(2);
                    
                    validEvent = validWindow(tt(it),tlims,dur,nframe);%2014/12/19 DS
                    
                    if validEvent
                        tensorInt = TheTensor(:,:,frmLeft:frmRight);
                        if meansub
                            tensorInt = tensorInt - mean(tensorInt(:));
                        end
                        EventTensor = EventTensor + tensorInt;
                        if keepAll
                            allEvents(end+1,:,:,:) = tensorInt;
                        end
                        cnt = cnt + 1;
                        %%%%%%%%%%%%%%%%%% DEBUG
                        % bpviewframes(TheTensor(:,:,frmLeft:frmRight));
                        % Mfr(it) = getframe(gcf);
                        %%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                end
                
                if cnt > 0
                    EventTensor = EventTensor./cnt;
                end
                disp([num2str(cnt) '/' num2str(nt) 'events were averaged.'])
                return
                
                %------------------ end of function ------------------
            end
            
        end
        
        function [SpaceAvg, timevec] = SpaceAverages(S, Conditions, TimeFrom, TimeTo, roi)
            % Returns time vectors of averages across space (time x condition x roi).
            %
            % SpaceAvg = S.SpaceAverages
            % SpaceAvg = S.SpaceAverages(Conditions )
            % SpaceAvg = S.SpaceAverages(Conditions, TimeFrom )
            % SpaceAvg = S.SpaceAverages(Conditions, TimeFrom, TimeTo )
            % SpaceAvg = S.SpaceAverages(Conditions, TimeFrom, TimeTo, roi)
            % returns avg time-course within specified roi
            % roi: 3D matrix of {0, 1}, size should be (S.nRows, S.nCols, number of ROI)
            % OR row of numROI x column of 4 numbers, specifying [xbegin ybegin xlength ylength]
            
            % 9/5/2014 DS added 5th input
            % 12/11/2014 DS added another roi format
            
            if nargin<5 roi = ones(S.nRows, S.nCols); end
            if nargin<4 || isempty(TimeTo    ), TimeTo   = S.TimeVec(end); end %modified 25/8/15
            if nargin<3 || isempty(TimeFrom  ), TimeFrom = S.TimeVec(1);   end %modified 25/8/15
            if nargin<2 || isempty(Conditions), Conditions = 1:S.nConds;   end
            
            %added on 25/8/15
            [~, FrameFrom] = min(abs(TimeFrom - S.TimeVec));
            [~, FrameTo] = min(abs(TimeTo - S.TimeVec));
            
            FrameFrom = max(1, FrameFrom);
            FrameTo   = min(S.nFrames, FrameTo);
            
            %if (size(roi,1) ~= size(S.Values,1)) || (size(roi,2) ~= size(S.Values,2))
            %    error('ROI size is not consistent with that of S.Valaues.');
            %end
            
            if size(roi,2) == 4
                numroi = size(roi,1);
                roi = round(roi);
                roi_cache = zeros(S.nRows, S.nCols, numroi);
                
                
                for rr = 1:numroi
                    if roi(rr,4) > S.nRows - roi(rr,2)
                        warning('roi y size excees S.nRows. Modify');
                        roi(rr,4) = S.nRows - roi(rr,2);
                    end
                    if roi(rr,3) > S.nCols - roi(rr,1)
                        warning('roi x size excees S.nCols. Modify');
                        roi(rr,3) = S.nRows - roi(rr,1);
                    end
                    roi_cache(roi(rr,2):roi(rr,2)+roi(rr,4), roi(rr,1):roi(rr,1)+roi(rr,3),rr) = 1;
                end
                roi = roi_cache;
                clear roi_cache
            else
                roi = single(roi);
                numroi = size(roi, 3);
            end
            SpaceAvg = zeros(FrameTo-FrameFrom+1, length(Conditions), numroi,'single');
            for rr = 1:numroi
                [roicache] = meshgrid(roi(:,:,rr), FrameFrom:FrameTo);
                roicache = reshape(roicache', S.nRows, S.nCols, FrameTo-FrameFrom+1);
                roisize = squeeze(sum(sum(roi(:,:,rr))));
                for cc = 1:length(Conditions)
                    SpaceAvg(:,cc,rr) = squeeze(nansum(nansum(roicache .* ...
                        squeeze(S.Values(:,:,FrameFrom:FrameTo,Conditions(cc)))))) / roisize;%DS 17/04/14. mean>nanmean
                end
            end
            
            if nargout > 1
                timevec = S.TimeVec(FrameFrom:FrameTo);
            end
        end
        
        function TimeAvg = TimeAverages(S, Flag)
            % TimeAvg = S.TimeAverages
            % returns images of mean across time
            %
            % TimeAvg = S.TimeAverages('median')
            % returns images of median across time
            
            if nargin == 1 || strcmp(Flag, 'mean')
                TimeAvg = squeeze(nanmean(S.Values,3));
            elseif strcmp(Flag, 'median')
                TimeAvg = squeeze(nanmedian(S.Values,3));
            end
        end
        
        
        function  [CrossCorrMaps, pvalMaps] = CrossCorr(S, Trace, nlags)
            %[CrossCorrMaps, pvalMaps] = CrossCorr(S, Trace)
            %calculates map of the cross-correlation between S.Values and
            %Trace with 0 lag
            %
            
            %[CrossCorrMaps, pvalMaps] = CrossCorr(S, Trace, nlags)
            %lets you specify lags between the two
            %
            % The output CrossCorrMaps has the length 2*nlags+1
            % note that the sampling rate should be the same for S.Values and Trace
            
            % 2015-12-31 DS modified when nlags=1, now much more efficient
            
            % TODO: use svd to make this faster?
            
            
            if nargin < 3
                nlags = 0;
            end
            
            if size(Trace,1)>size(Trace,2)
                Trace = Trace'; %16/10/18
            end
            CrossCorrMaps = zeros(size(S.Values,1),size(S.Values,2), 2*nlags + 1, S.nConds);
            
            for iCond = 1:S.nConds
                
                if nlags > 1
                    for iRow = 1: size(S.Values,1)
                        for iCol = 1: size(S.Values,2)
                            CrossCorrMaps(iRow,iCol,:,iCond) = ...
                                xcorr(S.Values(iRow,iCol,1:nlags+1,iCond),Trace(1:nlags+1), nlags,'unbiased');%'coeff'
                        end
                    end
                    pvalMaps = nan;%not yet implemented
                else
                    
                    Values_cat = reshape(S.Values(:,:,:,iCond), S.nRows*S.nCols, S.nFrames);
                    Values_cat = [Trace; Values_cat]';
                    [Rvector, Pvector] = corrcoef(Values_cat); %P is sometimes 0..
                    
                    %correlation to trace, in 1D
                    corr_cat = Rvector(1,2:end);
                    p_cat = Pvector(1,2:end);
                    
                    %correlation to trace, in 2D
                    CrossCorrMaps(:,:,iCond) = reshape(corr_cat, S.nRows, S.nCols);
                    pvalMaps_c = reshape(p_cat, S.nRows, S.nCols);
                    
                    %bonferroni correction
                    pvalMaps(:,:,iCond) = pvalMaps_c * S.nRows * S.nCols;
                end
            end
        end
        
        function  ParCorrMaps = ParCorrMaps(S, Traces)
            %CorrMaps = ParCorrMaps(S, Traces)
            %computes correlation of traces to stackset.values (0 lag only)
            %
            % Traces: variables x time
            % 2017-6-16 DS created from CrossCorr
            
            nVariables = size(Traces, 1);
            
            ParCorrMaps = zeros(size(S.Values,1),size(S.Values,2), nVariables, S.nConds);
            
            for iCond = 1:S.nConds
                Values_cat = reshape(S.Values(:,:,:,iCond), S.nRows*S.nCols, S.nFrames);
                %Values_cat = [Trace; Values_cat]';
                for tgtIdx = 1:nVariables
                    ctrlIdx = setxor(tgtIdx, 1:nVariables);
                    [Rvector] = partialcorr(Values_cat', Traces(tgtIdx,:)', Traces(ctrlIdx,:)');
                    
                    %correlation to trace, in 1D
                    %corr_cat = Rvector(1,2:end);
                    
                    %correlation to trace, in 2D
                    ParCorrMaps(:,:,tgtIdx,iCond) = reshape(Rvector, S.nRows, S.nCols);
                end
            end
        end
        
        function [ComplexMaps, AbsMaps, AngleMaps] = FourierMaps(S, ff)
            % Get a frequency component of a stackset
            %
            % [ComplexMaps, AbsMaps, AngleMaps] = FourierMaps(S, ff) lets
            % you specify the frequency. ff (in Hz) can be a number or a vector
            % with one value per condition.
            %
            % NOTE: The outputs are arrays nRows x nCols x nConds
            % (MC made this happen 4 Sept 2014, sorry if it breaks older
            % code, but it was a pain to have them as cell arrays)
            
            % TO DO: check unit of AbsMaps
            % right now ff needs to a scalar if nConds == 1
            
            % if the ff are all the same, make them become one number
            if length(unique(ff))==1, ff = ff(1); end
            
            yy = zeros(1,1,S.nFrames);
            aaa = zeros(S.nRows,S.nCols,S.nFrames);
            ComplexMaps = zeros(S.nRows,S.nCols,S.nConds);
            
            for iCond = 1:S.nConds
                
                if iCond==1 || length(ff)>1
                    if ff(iCond)==0
                        yy(:) =  ones(1,S.nFrames);
                    else
                        yy(:) = 2*exp(- (1:S.nFrames)/S.nFrames * ...
                            (S.TimeVec(end) - S.TimeVec(1)) *2*pi*1i*ff(iCond));
                    end
                    aaa = repmat(yy,[S.nRows,S.nCols,1]);
                end
                ComplexMaps(:,:,iCond) = double(mean(S.Values(:,:,:,iCond).*aaa, 3));
            end
            
            AbsMaps  = abs(  ComplexMaps );
            AngleMaps= angle(ComplexMaps );
            
        end
        
        function S = ZScore(S, firstBaseTime, lastBaseTime, meanValue, stdValue)
            % S = S.ZScore()
            % returns normalized tensor, by subtracting mean and dividing by
            % std over time for each pixel and condition
            %
            % S = S.ZScore(firsstBaseTime, lastBaseTime)
            % calculates mean and std from between firstBaseTime and lastBaseTime
            % (in seconds)  for each pixel and condition.
            %
            % S = S.ZScore([], lastBaseTime)
            % calculates mean and std from the begining of tensor to lastBaseTime
            %
            % S = S.ZScore(firstBaseTime, [])
            % calcualtes meand std from firstBaseTime to the end of the tensor
            %
            % S = S.ZScore([],[], meanValue, stdValue)
            %   = S.ZScore(firstBaseTime,lastBaseTime, meanValue, stdValue)
            %   = S.ZScore(firstBaseTime,[], meanValue, stdValue)
            %   = S.ZScore([],lastBaseTime, meanValue, stdValue)
            % uses meanValue and stdValue for calculation of zscore (ignore first and lastBaseTime inputs).
            % dimension: nRows x nCols x nConds.
            %
            % S = S.ZScore(firstBaseTime,[], meanValue, [])
            % calculates std from [firstBaseTime (end of the tensor)].
            %
            % S = S.ZScore(firstBaseTime,[], [], stdValue)
            % calculates mean from [firstBaseTime (end of the tensor)]
            %
            % S = S.ZScore([], lastBaseTime, meanValue, [])
            % calculates std from [(begin of the tensor) lastBaseTime]
            %
            % S = S.ZScore([], lastBaseTime, [], stdValue)
            % calculates mean from [(begin of the tensor) lastBaseTime]
            %
            % S = S.ZScore([], [], meanValue, [])
            % calculates std from whole period of the tensor
            %
            % S = S.ZScore([], [], [], stdValue)
            % calculates mean from whole period of the tensor
            
            
            % See also: stackset.SubtractBase
            
            % 2014-12-23 DS added 2nd and 3rd inputs
            % 2015-01-28 DS added "meanValue" and "stdValue" inputs
            
            if nargin < 2 || isempty(firstBaseTime)
                firstBaseTime = S.TimeVec(1);
            end
            
            if nargin < 3 || isempty(lastBaseTime)
                lastBaseTime = S.TimeVec(end);
            end
            
            if nargin < 4
                meanValue = [];
            end
            
            if nargin < 5
                stdValue = [];
            end
            
            baseStack = S.Trim([firstBaseTime lastBaseTime]);
            
            for iCond = 1:S.nConds
                
                if isempty(meanValue)
                    meanValue_cache = mean(baseStack.Values(:,:,:,iCond), 3);
                else
                    meanValue_cache = meanValue(:,:,iCond);
                end
                if isempty(stdValue)
                    stdValue_cache = std(baseStack.Values(:,:,:,iCond), 0, 3);
                else
                    stdValue_cache = stdValue(:,:,iCond);
                end
                
                MeanAcrossTime = repmat(meanValue_cache, [1 1 S.nFrames]);
                StdAcrossTime = repmat(stdValue_cache, [1 1 S.nFrames]);
                S.Values(:,:,:,iCond) = (S.Values(:,:,:,iCond) - MeanAcrossTime)./StdAcrossTime;
            end
        end
        
        
        function [Kernel, TimeVecKernel] = Shift_Regress(S,predictor,t1,t2,down_sample,nodetrend_flag)
            % [Kernel, TimeVecKernel] = Shift_Regress(S,predictor,t1,t2,down_sample,nodetrend_flag)
            % Compute the kernel between S.Values and the predictor for each condition. to Allows to down_sample and detrend.
            % Returns Kernel and TimeVecKernel
            %
            % See also: \\zserver\code\tools\shift_regress
            
            %AP created Jan '14
            
            if nargin<3
                tt=1:size(S.Values,3);
            end
            if nargin<4
                t1 = -5; %s
                t2 = 5; %s
            end
            if nargin<6
                down_sample = 1;
            end
            if nargin<7
                nodetrend_flag = 1;
            end
            
            % Kernel      % the kernel, an array nTimesCorr x nRows x nCols x nConds AP Jan '14
            % TimeVecKernel % nTimesCorr x nConds
            for icond=1:S.nConds,
                [Kernel(:,:,:,icond),TimeVecKernel(:,icond)] = shift_regress(predictor(:,icond),permute(S.Values(:,:,:,icond),[3 1 2]),S.TimeVec,t1,t2,down_sample,nodetrend_flag);
            end
        end
        
        function S = Align(S, doFliplr,doFlipud,t_concord, nRows, nCols)
            % rotate and/or flip Stackset
            %
            % S = A.Align(doFliplr, doFlipud) flips flips stackset image for each time.
            % doFliplr = 'y' flips left-right
            % doFlipud = 'y' flips up-down
            %
            % S = A.Align([],[],t_concord, nRows, nCols) rotates
            % stackset, where
            % t_concord is a spatial transformation structure returned by
            % maketform or cp2tform.
            % nRows and nCols specifies matrix size afte the conversion.
            %
            % S = A.Align(doFliplr, doFlipud, t_concord, nRows, nCols)
            % flips and then rotates stackset.
            %
            % see also: flipdim, imrotate
            % cf. tools.rotate_flip_cams, tools.LoadRotateInfo
            %
            % DS 16/5/14
            
            % TO DO: when the resize factors of 2 exps are different,
            % imresize one so that the resize factors will be same, and
            % then apply rotation
            %
            
            if isempty(doFliplr)
                doFliplr = 'n';
            end
            if isempty(doFlipud)
                doFlipud = 'n';
            end
            
            if nargin == 6 && ~isempty(t_concord) && ~isempty(nCols)
                doRotate = true;
            else
                doRotate = false;
            end
            
            if nargin == 4 || nargin == 5
                error('Missing information for rotation');
            end
            
            
            if doRotate
                S.nRows = nRows;
                S.nCols = nCols;
            end
            
            Values_aligned = zeros(S.nRows, S.nCols, S.nFrames, S.nConds,'single');
            
            
            disp('Aligning stackset...');
            if strcmp(doFliplr, 'y')
                %    S.Values = fliplr(S.Values(:,:,tt));
                S.Values = flipdim(S.Values, 2);
            end
            if strcmp(doFlipud, 'y')
                %    S.Values = flipud(S.Values(:,:,tt));
                S.Values = flipdim(S.Values, 1);
            end
            
            if doRotate
                % R = makeresampler({'cubic', 'nearest'}, 'replicate');
                R = makeresampler({'cubic', 'nearest'}, 'symmetric');
                % R = makeresampler({'cubic', 'nearest'}, 'circular');
                for icond = 1:S.nConds
                    for tt = 1:S.nFrames
                        %tformarray will make this faster but there are
                        %no XData/YData options...
                        % probably this is the reason: http://www.mathworks.com/matlabcentral/newsreader/view_thread/63635
                        
                        Values_aligned(:,:,tt,icond) = imtransform(S.Values(:,:,tt,icond), ...
                            t_concord, R, ...
                            'XData',[1 S.nCols], 'YData',[1 S.nRows]);
                    end
                end
                S.Values = Values_aligned;
            end
            
        end
        
        function [S, idxRow, idxCol, idxFrame] = InterpolateNan(S)
            % S  = InterpolateNan()
            % finds pixels of 'NAN' and 'INF', and interpolates with median
            % this tool is supposedly useful to apply before filtering
            %
            % [S, idxRow, idxCol, idxFrame] = InterpolateNan()
            % also returns index of interpolated pixels
            
            % 2015-1-26 DS added 2nd-4th output
            
            disp('Interpolating NAN pixels...');
            
            
            for icond = 1:S.nConds
                Values_cache = S.Values(:,:,:,icond);
                
                nanidx = [find(isnan(Values_cache(:))); find(isinf(Values_cache(:)))];
                Values_cache(nanidx) = nanmedian(Values_cache(:));
                S.Values(:,:,:,icond) = Values_cache;
                
                [idxRow{icond}, idxCol{icond}, idxFrame{icond}] = ind2sub(size(S.Values), nanidx);
            end
        end
        
        function [S, Dx, Dy] = Register(S, Dy, Dx, RowFromTo, ColFromTo, targetImage, marginSize, precision)
            % S = S.Register(Dy,Dx) register images according to the input (Dx, Dy)
            %
            % [S, Dx, Dy] = Register(S, [], [], RowFromTo, ColFromTo, targetImage)
            % finds (Dx, Dy) against targetImage, then register Images
            %
            % [...] = Register(S, [], [], RowFromTo, ColFromTo, targetImage, marginSize)
            % allows to specify size of margin around the targetImage, to
            % reduce boundary artifact using tools.RapidReg (max:1,
            % default:0.2)
            %
            % [...] = Register(S, [], [], RowFromTo, ColFromTo, targetImage, marginSize, precision)
            % allows to specify precision of image registration using
            % tools.RapidReg (default: 100)
            %
            % See also: tools.RapidReg, tools.ImageReg
            
            
            if size(Dy,1) < size(Dy,2) %DS added on 7/9/15
                Dy = Dy';
            end
            
            if size(Dx,1) < size(Dx,2) %DS added on 7/9/15
                Dx = Dx';
            end
            
            if nargin < 8
                precision = 100;
            end
            
            if nargin < 7
                marginSize = 0.2;
            end
            
            if nargin > 3 %find Dx and Dy
                
                S_crop = S.Crop(RowFromTo, ColFromTo);
                Values_crop = S_crop.Values;
                targetImage_crop = targetImage(RowFromTo(1)+1:RowFromTo(end)-1, ColFromTo(1)+1:ColFromTo(end)-1);
                
                for icond = 1:S.nConds
                    [Dx(:,icond), Dy(:,icond)] = ...
                        tools.RapidReg(Values_crop(:,:,:,icond), targetImage_crop, marginSize, precision, 'nopar');
                end
            end
            
            for icond = 1:S.nConds
                %                 addpath('\\zserver\Code\Register\development\');
                %                 [movie, validX, validY] = img.translate(S.Values(:,:,:,icond), Dx(:,icond), Dy(:,icond));
                %S.Values(:,:,:,icond) = tools.imageReg(S.Values(:,:,:,icond), Dx(:,icond), Dy(:,icond), marginSize);
                
                %slow, but without artifact near border of window
                Dx(isnan(Dx)) = 0;
                Dy(isnan(Dy)) = 0;
                Values_c = [];
                for tt = 1:S.nFrames
                    Values_c(:,:,tt) = imtranslate(S.Values(:,:,tt,icond),[Dx(tt,icond) Dy(tt,icond)]);
                end
                S.Values(:,:,:,icond) = Values_c;
            end
        end
        
        function S = WavletCoef(S, freq, nCycles)
            % S. WavletCoef(S, freq, nCycles)
            % returns wavelet coefficient (complete Morlet wavelet) at
            % given frequncy and number of cycles
            
            % TO DO: allow freq to be a vector?
            
            Fb = 2; %fixed. other values not allowed for zero-sum
            delta = 1/S.FrameRate; %sampling interval(s)
            Fc=1/2/pi * nCycles; %[cycles]. modefied on 9.20.2011
            scale = Fc ./ delta ./ freq; %this ensures linear relation bet. frequency and cycles(Fc)
            
            
            %% wavlet transform
            matlabPath = '\\ZOMBIE\Users\daisuke\Documents\MATLAB\dsbox\mywavelet';
            path(path, matlabPath);
            
            wavename = sprintf('cmoc%f-%f',Fb, Fc);
            
            coef = zeros(S.nRows, S.nCols, S.nFrames,'single');
            for icond = 1:S.nConds
                
                disp(['Calculating wavelet coefficient ' num2str(icond) '/' num2str(S.nConds)]);
                for rr = 1:S.nRows
                    for cc = 1:S.nCols
                        signal = squeeze(S.Values(rr,cc,:,icond))';
                        coef(rr,cc,:,icond) = (mycwt(signal, scale, wavename));
                    end
                end
            end
            
            S.Values = coef;
        end
        
        % added by Matteo 2015-10
        function T = Resample(S,NewTimeVec)
            % S = S.Resample(NewTimeVec)
            
            T = S;
            vv = shiftdim( S.Values, 2);
            vv_interp = interp1( S.TimeVec, vv, NewTimeVec );
            if S.nConds == 1 %added on 5/11/15 by DS
                T.Values = shiftdim(vv_interp,1);
            else
                T.Values = shiftdim(vv_interp,2);
            end
            T.TimeVec = NewTimeVec;
            T.nFrames = length(NewTimeVec);
            T.FrameRate = 1./median(diff(NewTimeVec));
        end
        
        function S = SubtractBase(S, firstBaseTime, lastBaseTime)
            % S = S.SubtractBase()
            % subtract average values across
            % whole or specified period for each condition.
            %
            % If 1st input is an image, subtract S.Values by
            % this image.
            %
            % S = S.SubtractBase(firstBaseTime, lastBaseTime)
            % subtract average values between
            % firstBaseTime and lastBaseTime (in seconds).
            %
            % S = S.SubtractBase([], lastBaseTime)
            % subtract from the begining of tensor to lastBaseTime
            %
            % S = S.SubtractBase(firstBaseTime, [])
            % subtract from firstBaseTime to the end of the tensor
            %
            % See also: StackSet.Differential
            
            % 2014-12-02 DS created
            % 2014-12-21 DS allowed inputs to be empty
            % 2015-9-25 DS if 1st input is an image, subtract S.Values by
            % this image
            
            if nargin < 2 || isempty(firstBaseTime)
                firstBaseTime = S.TimeVec(1);
            end
            
            if nargin < 3 || isempty(lastBaseTime)
                lastBaseTime = S.TimeVec(end);
            end
            
            
            if (size(firstBaseTime,1) == S.nRows) && (size(firstBaseTime,2) == S.nCols)
                baseImage = firstBaseTime;
            else
                baseStack = S.Trim([firstBaseTime lastBaseTime]);
                baseImage = baseStack.TimeAverages;
            end
            
            for icond = 1:S.nConds
                baseCube = meshgrid(baseImage(:,:,icond),1:S.nFrames);
                baseCube = reshape(baseCube, S.nFrames, S.nRows, S.nCols);
                baseCube = shiftdim(baseCube,1);
                S.Values(:,:,:,icond) = S.Values(:,:,:,icond) - baseCube;
            end
        end
        
        function f = PixelTuingCurve(S, cLabels, alphadata)
            % pixelTuningCurveViewer(allFrames, cLabels, alphadata)
            % Displays: A) the image of the brain for a certain time point and
            % condition; B) the traces across time for a certain pixel and all
            % conditions; C) the value at that pixel and time point across conditions
            %
            % input:
            % - cLabels: 1 x nConditions, the "value" for each condition, e.g. the
            % contrast of the stimulus. Default is 1:nConditions. needs
            % to be numerical array
            % See also: tools.pixelTuningCurveViewer
            %
            % 2016-1-4 DS
            
            if nargin<3
                alphadata = [];
            end
            if nargin<2
                cLabels = 1:S.nConds;
            end
            f = tools.pixelTuningCurveViewer(S.Values, cLabels, S.TimeVec, alphadata);
        end
        
        
        function S = FilterInSpace_gauss(S,sigma)
            %S = FilterInSpace_gauss(sigma)
            for icond = 1:S.nConds
                for tt = 1:S.nFrames
                    S.Values(:,:,tt,icond) = imgaussfilt(S.Values(:,:,tt,icond),sigma);
                end
            end
        end
        
        function S = FillDims(S)
            % S = FillDims
            % updates S.nRows, S.nCols, S.nFrames, S.nConds using dimensions S.Values
            tensorSize = size(S.Values);
            S.nRows = tensorSize(1);
            S.nCols = tensorSize(2);
            if length(tensorSize)>2
                S.nFrames = tensorSize(3);
            else
                S.nFrames = 1;
            end
            if length(tensorSize)>3
                S.nConds = tensorSize(4);
            else
                S.nConds = 1;
            end
        end
        
        function S = DerRec(S)
            %take derivative then rectify in each condition
            
            S.Values = cat(3, zeros(S.nRows, S.nCols, 1, S.nConds), max(diff(S.Values,[], 3),0));
        end
    end
end
