
%specifications
% needs to be fast so can be used during experiment
% work on SVD data
% analyze response to sparse noise 
% compatible with 2p (map and trace) and wf 
% options for post processing on imaging data such as dF/F, spikes

addpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis'));
rmpath(genpath('\\ad.monash.edu\home\User006\dshi0006\Documents\MATLAB\master\Analysis\wfAnalysis\obs'));

%% experiment info
expt.subject = 'test';%'dummy_wf';
expt.expDate = '2020-06-25_1';%'2020-06-16_1';
expt.expNum = 3;%9;

%% analysis info
nSV=1000;
params.movieSuffix = 'blue';% 'purple''corr_dFF'; %'blue'  %U for corr_dFF can have NANs..
params.useCorrected = 0;
respWin = [0.1 0.5];%for 6s AP [0.05 0.2]%for 6f AP %[0 0.18];

%% data location
thisDate = expt.expDate(1:10);
thisSeries = str2num(expt.expDate(12:end));
expPath = fileparts(dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'widefield','master'));
saveVpath = expPath;

%% get stimulus info
expt = grabStimTimesWF(expt,1);
expt.timelineFile = dat.expFilePath(expt.subject, thisDate, thisSeries, expt.expNum, 'timeline','master');
expt = grabSparseFrames(expt); 
%< upto this point, every experiment in 2p/wf should be able to pass 

stimInd = find([1 diff(expt.stimFrames.ss.ImageSequence) > 0]);
stimTimes = stimTimes(stimInd); %SF
stimImg = nan(size(expt.stimFrames.ss.ImageTextures{1},1),...
    size(expt.stimFrames.ss.ImageTextures{1},2));
for t = 1:length(stimInd)
    stimImg(:,:,t) = expt.stimFrames.ss.ImageTextures{expt.stimFrames.ss.ImageSequence(stimInd(t))};
end


%% check case for mismatch between photodiode and stimulus in Protocol
frameTimes_c = expt.stimTimes.frameTimes{1};


        stim_screen = cat(3, expt.stimFrames.ss.ImageTextures{:});
        
        % odd number of stimuli, but one extra photodiode flip to come back down
        if mod(size(stim_screen,3),2) == 1 && ...
                length(frameTimes_c) == size(stim_screen,3) + 1
            frameTimes_c(end) = [];
            warning('Odd number of stimuli, removed last photodiode');
        end
        
        % If there's still a mismatch, break
        if size(stim_screen,3) ~= length(frameTimes_c)
            warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                num2str(length(frameTimes_c)) ' photodiode pulses']);
            
            % Try to estimate which stim were missed by time difference
            photodiode_diff = diff(frameTimes_c);
            max_regular_diff_time = prctile(diff(frameTimes_c),99);
            skip_cutoff = max_regular_diff_time*1.5;%2;
            photodiode_skip = find(photodiode_diff > skip_cutoff);
            est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
            stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                1:length(photodiode_skip),'uni',false));
            
            if isempty(est_n_pulse_skip) || length(frameTimes_c) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                error('Can''t match photodiode events to stimuli')
            end
        end
        
        stim_times = frameTimes_c;
        

%% retrieve imaging data (SVD format)
%[U,V,t] = loadSVD(expt);%not yet
[U, V, expt.frameTimes, mimg] = quickLoadUVt(expPath, nSV, saveVpath, params);

%% some post processing of imaging data (optional)

%% compute retinotopy & VFS map, visualize them
%stimPositions = ??


[signMap, xMap, yMap] = sparseRetinotopy(U, V, expt.frameTimes, stimPositions, stimTimes, mimg, respWin); %function by NS. different to AP/SF?

