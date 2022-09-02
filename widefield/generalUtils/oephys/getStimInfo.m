function stimInfo = getStimInfo(c)
switch c.paradigm
    case {'oriXYZ','oriXYZc'}
        
        %     pos = get(c.gabor.prms.X, 'atTrialTime',inf);
        posX = get(c.gabor.prms.X, 'atTrialTime',inf);
        posY = get(c.gabor.prms.Y, 'atTrialTime',inf);
        
        posXparam = unique(posX);
        posYparam = unique(posY);
        if numel(posXparam)==1 && numel(posYparam)>1
            stimInfo.labelDescription = 'vertical position [deg]';
            pos = posY;
        elseif numel(posXparam) >= 1 && numel(posYparam)==1
            stimInfo.labelDescription = 'horizontal position [deg]';
            pos = posX;
        else
            error('both X&Y positions are variable?');
        end
            
%         
%         posXYparam = zeros(length(posXparam)*length(posYparam),1);
%         istim = 1;
%         for ix = 1:length(posXparam)
%             for iy = 1:length(posYparam)
%                 posXYparam(istim) = posXparam(ix) +1i*posYparam(iy);
%                 istim = istim+1;
%             end
%         end
        
        duration = get(c.gabor.prms.duration,'atTrialTime',inf);
        stimInfo.duration = duration(1)/1e3;
        %calcWin = [-0.2*duration duration*1.2];
        stimInfo.stimLabels = pos;
        stimInfo.condLabels = unique(pos);
        stimInfo.tgtFreq = [];
    case 'kalatsky'
        contDir = get(c.gabor.prms.contDir, 'atTrialTime', inf);
        duration = get(c.gabor.prms.duration,'atTrialTime',inf);
        stimInfo.duration = duration(1)/1e3;
        %calcWin = [-0.2*duration duration*1.2];
        stimInfo.stimLabels = contDir;
        stimInfo.condLabels = unique(contDir);
        stimInfo.labelDescripiton = 'stimulus direction';
        blankDur = get(c.gabor.prms.bdur,'atTrialTime',inf,'trial',1);
        dur = get(c.gabor.prms.dur,'atTrialTime',inf,'trial',1);
        stimInfo.tgtFreq = 1/get(c.gabor.prms.dur); %[hz] for frequency analysis. Skip if empty
        %tgtFreq = 1/(dur + blankDur);
        stimInfo.stimDur = dur - blankDur;
        stimInfo.vfRange = [45 105];%TODO determine from cic
    case 'runPassiveMovies'
        stimInfo.tgtFreq = [];
        stimInfo.stimLabels = [];
        duration = get(c.movie.prms.duration,'atTrialTime',inf);
        stimInfo.duration = duration(1)/1e3;
        %calcWin = [-0.2*duration duration*1.2];
end
