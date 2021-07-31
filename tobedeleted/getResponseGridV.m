function [response_grid, surround_time, response_t] = getResponseGridV(expt, stim_screen, V, respWin)
% INPUT:
% expt.frameTimes: time of recording
% expt.stimTimes.frameTimes{iStim, iRepeat}: time of stimulus
% stim_screen: [yScreenPix x xScreenPix x time x cond]
% V: [nSV x time]
% respWin: [begin end] (s)
%
% OUTPUT:
% response_grid
% surround_time: avg within respWin {yScreenPix, xScreenPix}
% response_t: avg across events [delay x nSV? x yScreenPix*xScreenPix]

framerate = 1./nanmedian(diff(expt.frameTimes));

ny = size(stim_screen,1);
nx = size(stim_screen,2);

% NS:
% [stimTimeInds, stimPositions, stimArray] = computeSparseNoiseSignals(block)
% this function assumes stimlus is presented at one location at one time
% [signMap, xMap, yMap] = sparseRetinotopy(U, V, expt.frameTimes, stimPositions, stimTimes, mimg, respWin); %function by NS. different to AP/SF

%AP:
% Get average response to each stimulus
%INPUT: stim_screen, V, respWin
%OUTPUT: response_grid{px_y,px_x}(nSV, #presentation)
surround_samplerate = 1/(framerate*1);
surround_time = respWin(1):surround_samplerate:respWin(2);
response_n = nan(ny,nx);
response_grid = cell(ny,nx);
response_t = [];
nstims = size(expt.stimTimes.frameTimes,1);
nrepeats = size(expt.stimTimes.frameTimes,2);
for iStim = 1:nstims
    for iRepeat = 1:nrepeats
        for px_x = 1:nx
           for px_y = 1:ny
                 
                % Use first frame of dark or light stim
                align_stims = (stim_screen(px_y,px_x,2:end, iStim)~= 0) & ...
                    (diff(stim_screen(px_y,px_x,:, iStim),[],3) ~= 0);
                align_times = expt.stimTimes.frameTimes{iStim, iRepeat}(find(align_stims)+1);
                %         keyboard
                align_times = align_times(round(length(align_times)/2):end);
                if size(align_times,1) < size(align_times,2)
                    align_times = align_times';
                end
                
                response_n(px_y,px_x) = length(align_times);
                %         keyboard
                % Don't use times that fall outside of imaging
                align_times(align_times + surround_time(1) < expt.frameTimes(2) | ...
                    align_times + surround_time(2) > expt.frameTimes(end)) = [];
                %         keyboard
                % Get stim-aligned responses, 2 choices:
                
                % 1) Interpolate times (slow - but supersamples so better)
                %         align_surround_times = bsxfun(@plus, align_times, surround_time);
                %         peri_stim_v = permute(mean(interp1(expt.frameTimes,V',align_surround_times),1),[3,2,1]);
                
                % 2) Use closest frames to times (much faster - not different)
                align_surround_times = bsxfun(@plus, align_times, surround_time);
                frame_edges = [expt.frameTimes,expt.frameTimes(end)+1/framerate];
                align_frames = discretize(align_surround_times,frame_edges);
                
                align_frames(any(isnan(align_frames),2),:) = [];
                %         keyboard
                % Define the peri-stim V's as subtracting first frame (baseline)
                peri_stim_v = bsxfun(@minus, ...
                    reshape(V(:,align_frames)',size(align_frames,1),size(align_frames,2),[]), ...
                    reshape(V(:,align_frames(:,1))',size(align_frames(:,1),1),size(align_frames(:,1),2),[]));
                %         keyboard
                mean_peri_stim_v = permute(mean(peri_stim_v,2),[3,1,2]); %[nV, #presentation?] temporal mean within surround_time
                %         keyboard
                % Save V's
                response_grid{px_y,px_x} = mean_peri_stim_v; %temporal avg of every events
                
                response_t = cat(3,response_t, squeeze(mean(peri_stim_v,1))');%avg across events
            end
        end
    end
end
