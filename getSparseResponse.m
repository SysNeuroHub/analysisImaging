function [response_t, response_grid] = getSparseResponse(V, frameTimes, ...
    stim_screen, stimFrameTimes, respWin, polarity)
% Get average response to each stimulus of sparse noise
%INPUT:
% V: [nSV x frames]
% frameTimes: [frames]
% stim_screen: [px_y x px_x x stimFrames] 
% stimFrameTimes: time of stimulus presentation [s] (no longer cell)
% respWin: [min max] (T)
% polarity: stimulus color of interest. if empty, both (white/black/'')
% (default:white)
%
%OUTPUT:
% response_grid{px_y,px_x}(nSV, #presentation)
% response_t: avg across events [delay x nSV? x yScreenPix*xScreenPix]

% 16/7/20 created from AP_sparsenoise...
% TODO separate times for baseline and response. cf SF's code
% this code cannot distintguish white/black stim. hence cannot distinguish
% on/off subfields

if nargin<6
    polarity = 'white';
end
subtractFirstFrame = 0;

framerate = 1./nanmedian(diff(frameTimes));
ny = size(stim_screen,1);
nx = size(stim_screen,2);

surround_samplerate = 1/framerate;
% surround_time = respWin(1):surround_samplerate:respWin(2);
surround_time = linspace(respWin(1),respWin(2),round(diff(respWin)/surround_samplerate+1));
response_n = nan(ny,nx);
response_grid = cell(ny,nx);
response_t = [];
for px_x = 1:nx
    for px_y = 1:ny
        
        switch polarity
            case 'white'
                align_stims = (stim_screen(px_y,px_x,2:end)== 1) & ...
                    (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
            case 'black'
                align_stims = (stim_screen(px_y,px_x,2:end)== 0) & ...
                    (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
            otherwise
                % Use first frame of dark or light stim
                %                 align_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
                %                     (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
                 align_stims = ((stim_screen(px_y,px_x,2:end)== 1) | ...
                     (stim_screen(px_y,px_x,2:end)== 0)) & ...
                    (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
        end
        
        align_times = stimFrameTimes(find(align_stims)+1);
        %align_times = stimFrameTimes(find(align_stims==0)+1);
        
        align_times = align_times(round(length(align_times)/2):end);
        %<trim the first half of the events ...why?
        
        if size(align_times,1) < size(align_times,2)
            align_times = align_times';
        end
        
        response_n(px_y,px_x) = length(align_times);
        %         keyboard
        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < frameTimes(2) | ...
            align_times + surround_time(2) > frameTimes(end)) = [];
        %         keyboard
        % Get stim-aligned responses, 2 choices:
        
        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(expt.frameTimes,V',align_surround_times),1),[3,2,1]);
        
        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = bsxfun(@plus, align_times, surround_time);
        frame_edges = [frameTimes,frameTimes(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);
        
        align_frames(any(isnan(align_frames),2),:) = [];
        %         keyboard
        % Define the peri-stim V's as subtracting first frame (baseline)
        peri_stim_v = reshape(V(:,align_frames)',size(align_frames,1),size(align_frames,2),[]);
        if subtractFirstFrame
            peri_stim_v = bsxfun(@minus, peri_stim_v, ...
                reshape(V(:,align_frames(:,1))',size(align_frames(:,1),1),size(align_frames(:,1),2),[]));
        end
        %         keyboard
        mean_peri_stim_v = permute(mean(peri_stim_v,2),[3,1,2]); %[nV, #presentation?] temporal mean?
        %         keyboard
        % Save V's
        response_grid{px_y,px_x} = mean_peri_stim_v;
        
        response_t = cat(3,response_t, squeeze(mean(peri_stim_v,1))');
    end
end
