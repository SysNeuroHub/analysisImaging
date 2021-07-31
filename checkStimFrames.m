function expt = checkStimFrames(expt, nframes_protocol)
%expt = checkStimFrames(expt, nframes_protocol)
%returns expt.stimTimes.frameTimes after checking
%case for mismatch between photodiode and stimulus in Protocol

%10/7/20 created from AP_sparseRetinotopy

%TODO: make the script compatible with multiple repeats

nrepeats = size(expt.stimTimes.frameTimes,2);
nstims = size(expt.stimTimes.frameTimes,1);

for irepeat = 1:nrepeats
    for istim = 1:nstims
        frameTimes_pd = expt.stimTimes.frameTimes{istim, irepeat};
        nframes = nframes_protocol(istim);
        
        % odd number of stimuli, but one extra photodiode flip to come back down
        if mod(nframes,2) == 1 && ...
                length(frameTimes_pd) == nframes + 1
            frameTimes_pd(end) = [];
            warning('Odd number of stimuli, removed last photodiode');
        end
        
        % If there's still a mismatch, break
        if nframes ~= length(frameTimes_pd)
            warning([num2str(nframes) ' stimuli, ', ...
                num2str(length(frameTimes_pd)) ' photodiode pulses']);
            
            % Try to estimate which stim were missed by time difference
            photodiode_diff = diff(frameTimes_pd);
            max_regular_diff_time = prctile(diff(frameTimes_pd),99);
            skip_cutoff = max_regular_diff_time*1.5;%2;
            photodiode_skip = find(photodiode_diff > skip_cutoff);
            est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
            stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                1:length(photodiode_skip),'uni',false));
            
            if isempty(est_n_pulse_skip) || length(frameTimes_pd) + sum(est_n_pulse_skip) ~= nframes_protocol
                error('Can''t match photodiode events to stimuli')
            end
        end
        
        expt.stimTimes.frameTimes{istim,irepeat} = frameTimes_pd;
    end
end
