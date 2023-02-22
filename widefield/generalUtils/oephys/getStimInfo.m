function stimInfo = getStimInfo(c)
%stimInfo = getStimInfo(c)
% .labelDescription
% .duration
% .stimLabels
% .condLabels
% .tgtFreq

switch c.paradigm
    case {'oriXYZ','oriXYZc'}
        
        %     pos = get(c.gabor.prms.X, 'atTrialTime',inf);
        posX = get(c.gabor.prms.X, 'atTrialTime',inf);
        posY = get(c.gabor.prms.Y, 'atTrialTime',inf);
        orientation = get(c.gabor.prms.orientation, 'atTrialTime',inf);
       
        posXparam = unique(posX);
        posYparam = unique(posY);
        oriParam = unique(orientation);
        if numel(posXparam)==1 && numel(posYparam)>1
            stimInfo.labelDescription = 'vertical position [deg]';
            pos = posY;
        elseif numel(posXparam) > 1 && numel(posYparam)==1
            stimInfo.labelDescription = 'horizontal position [deg]';
            pos = posX;
        elseif numel(posXparam) == 1 && numel(posYparam)==1 && numel(oriParam)>1
            pos = orientation;
            stimInfo.labelDescription = 'orientation [deg]';
        else
            pos = posX;
            %error('both X&Y positions are variable?');
        end
            
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
        orientation = unique(get(c.gabor.prms.orientation,'atTrialTime',inf));
        if orientation==0
            stimInfo.labelDescription = 'kalatsky Y (-1: downward, +1:upward)';
            stimInfo.vfRange = sort(get(c.gabor.prms.Y,'atTrialTime',inf));
        elseif orientation==90
            stimInfo.labelDescription = 'kalatsky X (-1: leftward, +1:rightward)';
            stimInfo.vfRange = sort(get(c.gabor.prms.X,'atTrialTime',inf));
        end
        stimInfo.blankDur = get(c.gabor.prms.bdur,'atTrialTime',inf,'trial',1);
        dur = get(c.gabor.prms.dur,'atTrialTime',inf,'trial',1);
        stimInfo.tgtFreq = 1/get(c.gabor.prms.dur); %[hz] for frequency analysis. Skip if empty
        %tgtFreq = 1/(dur + blankDur);
        stimInfo.stimDur = dur - stimInfo.blankDur;

    case 'kalatsky_rotation'
        contDir = get(c.contour.prms.contDir, 'atTrialTime', inf);
        duration = c.trialDuration; %get(c.contour.prms.duration,'atTrialTime',inf);
        stimInfo.duration = duration(1)/1e3;
        stimInfo.stimLabels = contDir;
        stimInfo.condLabels = unique(contDir);
        stimInfo.stimDur = get(c.contour.prms.dur,'atTrialTime',inf,'trial',1);
        stimInfo.blankDur = 0;
        stimInfo.tgtFreq = 2/stimInfo.stimDur; %2F to detect orientation column
        stimInfo.labelDescription = '-1: CW, +1: CCW';
        stimInfo.vfRange = [0 360];%
    case 'kalatsky_wedge'
        contDir = get(c.contour.prms.contDir, 'atTrialTime', inf);
        duration = c.trialDuration; %get(c.contour.prms.duration,'atTrialTime',inf);
        stimInfo.duration = duration(1)/1e3;
        stimInfo.stimLabels = contDir;
        stimInfo.condLabels = unique(contDir);
        stimInfo.stimDur = get(c.contour.prms.dur,'atTrialTime',inf,'trial',1);
        stimInfo.blankDur = 0;
        stimInfo.tgtFreq = 1/stimInfo.stimDur; %2F to detect orientation column
        stimInfo.labelDescription = '-1: CW, +1: CCW';
        stimInfo.vfRange = [0 360];%
    case 'freqTag'
        fields = fieldnames(c);
        jj=1;
        for ii = 1:numel(fields)
           if regexp(fields{ii}, regexptranslate('wildcard','gabor*'));
            tgtfields{jj} = fields{ii};
           jj=jj+1;
           end
        end
        tgtfields = sort(tgtfields);
        
        dur = []; xpos = []; ypos = [];
        for ff = 1:numel(tgtfields)
            dur(ff) = get(c.(tgtfields{ff}).prms.dur,'atTrialTime',inf);
            xpos(ff) = get(c.(tgtfields{ff}).prms.X,'atTrialTime',inf);
            ypos(ff) = get(c.(tgtfields{ff}).prms.Y,'atTrialTime',inf);
        end
        
        stimInfo.tgtFreq = 1./dur; 
        stimInfo.stimLabels=1;
        stimInfo.blankDur = 0;
        
        if numel(unique(xpos))==1 && numel(unique(ypos))>1
            stimInfo.labelDescription = 'y [deg]';
            stimInfo.condLabels = ypos;
        elseif numel(unique(xpos))>1 && numel(unique(ypos))==1
            stimInfo.labelDescription = 'x [deg]';
            stimInfo.condLabels = xpos;
        end
        duration = c.trialDuration; %get(c.contour.prms.duration,'atTrialTime',inf);
        stimInfo.duration = duration(1)/1e3; %nealy 20s longer than specified

    case 'runPassiveMovies'
        stimInfo.tgtFreq = [];
        stimInfo.stimLabels = [];
        duration = get(c.movie.prms.duration,'atTrialTime',inf);
        stimInfo.duration = duration(1)/1e3;
        
        fileNames_full = get(c.movie.prms.filename,'atTrialTime',inf);
        [fileDirectory,fileNames_idx] = cellfun(@(x)(fileparts(x)), fileNames_full, 'UniformOutput', false);
        stimInfo.stimLabels = cellfun(@(x)(str2num(x)), fileNames_idx);
        stimInfo.condLabels = unique(stimInfo.stimLabels); %sorted incremental order
        stimInfo.labelDescription = ['movie name in' fileDirectory{1}];
        
        %calcWin = [-0.2*duration duration*1.2];
        
        stimInfo.screenPix = [c.movie.prms.yPixels.value c.movie.prms.xPixels.value]; %#pixels
        stimInfo.width = c.movie.prms.width.value; % width of stimulus window in [deg]
        stimInfo.height = c.movie.prms.height.value; % height of stimulus window in [deg]
end
