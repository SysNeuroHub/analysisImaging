function [S, coef] = DetrendStacks(S,TimeClips, mask, fitDeg)
% detrend S.Values with polyfit, using data between TimeClips
%
% S = DetrendStacks(S) detrends stackset using linear regression coefficient (leaving time-invariant component)
%
% S = DetrendStacks(S,TimeClips) detrends stackset using the
% data with specified period TimeClips (in seconds)
%
% S = DetrendStacks(S,TimeClips, mask, fitDeg) calculates the coefficient
% for an averaged time-course within mask, and subtract the fit from each
% pixel.
% mask size should be (S.nRow, S.nCols)
%
% [S, coef] = DetrendStacks(...) also returns estimated coefficients.
% if mask is specified > (order of coefficients) x (condition in stackset)
% if mask is not specified > (nRows) x (nCols) x (order of coefficients) x
% (condition in stackset)

% 30-05-2014 DS created
% 03-06-2014 DS added 2nd output
% 2015-01-28 DS changed algorithm for linear case (use linear regression)
% 2015-08-21 DS overhaul. replaced 2nd output with coefficients. leave time-invariant component.
% changed order of inputs

%TO DO: exponential fit
% incopolate this to Stackset.Differential??

useInput = false;

if nargin < 4
    fitDeg = 1;
end
if nargin < 3
    mask = [];
end
if nargin < 2 || isempty(TimeClips)
    TimeClips = [S.TimeVec(1) S.TimeVec(end)];
end

%hack...why does this happen?
if size(S.TimeVec,1) > size(S.TimeVec,2)
    S.TimeVec = S.TimeVec';
end

firstFrameBase = find(diff(S.TimeVec >  TimeClips(1)));
lastFrameBase   = find(diff(S.TimeVec >= TimeClips(2)));
% [~, firstFrameBase] = min(abs(S.TimeVec-TimeClips(1)));
% [~, lastFrameBase] = min(abs(S.TimeVec-TimeClips(2)));

if ~isempty(firstFrameBase) && ~isempty(lastFrameBase) && firstFrameBase<lastFrameBase
    
    if isempty(mask)
        
        %% linear regression for each pixel
        coef = [];
        for iCond = 1:S.nConds
            disp(['Detrending stack ' num2str(iCond) '/' num2str(S.nConds) '...']);
            Values_cache = reshape(S.Values(:,:,:,iCond), ...
                S.nRows*S.nCols, S.nFrames);
            Values_cache = Values_cache'; %time x pixels
            
            taxis = [S.TimeVec' ones(S.nFrames,1)];
            
            %estimate linear-regression coefficient using [firstFrameBase:lastFrameBase]
            %note this will also subtract offset to be 0
            coef_part = (taxis(firstFrameBase:lastFrameBase,:)\Values_cache(firstFrameBase:lastFrameBase,:));
            regressed_part = taxis(:,1) * coef_part(1,:);%leave baseline
            Values_det = (Values_cache - regressed_part)';
            
                     
            %Values_det = detrend_regress(Values_cache, 1)';
            
            S.Values(:,:,:,iCond) = reshape(Values_det, S.nRows, ...
                S.nCols, S.nFrames);
            
            coef(:,:,:,iCond) = reshape(coef_part', S.nRows, S.nCols, 2);
            %coef(:,:,1): 1st-order coefficient
            %coef(:,:,2): time-invariant component or "baseline"
        end
        
        
        %% fit just once using pixels 
    else
        Values_cache = S.SpaceAverages([],[],[],mask);
        
        coef = zeros(fitDeg+1, S.nConds);
        for iCond = 1:S.nConds
            if ~useInput
                mp = polyfit(S.TimeVec(firstFrameBase:lastFrameBase), squeeze(Values_cache(firstFrameBase:lastFrameBase,iCond))', fitDeg);
            end
            
            for jjj = 1:length(mp)
                coef(jjj,iCond) = mp(jjj);
            end
            
            mp(end) = 0;% do not subtract time-invariant component
            
            mfit = polyval(mp, S.TimeVec);
            [~,~,MFIT] = ndgrid(1:S.nRows, 1:S.nCols, mfit);
            
            disp(['Detrending stack ' num2str(iCond) '/' num2str(S.nConds) '...']);
            S.Values(:,:,:,iCond) = S.Values(:,:,:,iCond) - MFIT + mean(mfit(firstFrameBase:lastFrameBase));%mean(S.Values) is not changed
        end
    end
    
else
    disp('Strange parameters , will not detrend...');
end