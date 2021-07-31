function [stim_screen, ny, nx] = getStimScreen(stimFrames)
%stim_screen = getStimScreen(stimFrames)

nstim = size(stimFrames,1);

for istim = 1:nstim
    if strcmp(stimFrames(istim).ss.Type, 'stimSparseNoiseUncorrAsync') %AP
        stim_screen(:,:,:,istim) = cat(3, stimFrames(istim).ss.ImageTextures{:});
    elseif strcmp(stimFrames(istim).ss.Type, 'stimSparseNoiseUncorr2') %SF
        stim_screen(:,:,:,istim) = nan(size(stimFrames(istim).ss.ImageTextures{1},1),...
            size(stimFrames(istim).ss.ImageTextures{1},2),stimFrames(istim).ss.nFrames);
        for t = 1:stimFrames.ss.nFrames
            stim_screen(:,:,t,istim) = stimFrames(istim).ss.ImageTextures{stimFrames(istim).ss.ImageSequence(t)};
        end
    else
        error([stimFrames(istim).ss.Type ' is currently not compatible.']);
    end
end

ny = size(stim_screen,1);
nx = size(stim_screen,2);
