function [S_dB, F, T] = spectrogramSVD(U,V,t, mask)

nRows = size(U,1);
nCols = size(U,2);
nFrames = size(V,2);
fs = 1/median(diff(t)); %[Hz]

window_length = round(4 * fs);   % 50 ms window (adjust as needed)
overlap = round(0.9 * window_length);
nfft = max(256, 2^nextpow2(window_length));

Wave = reshape(svdFrameReconstruct(U(1,1,:), V), [nFrames 1]);
[~, F, T] = spectrogram(Wave, window_length, overlap, nfft, fs);
     
S_dB = nan(nRows, nCols, numel(F),numel(T));
for irow = 1:nRows
    for icol = 1:nCols
        if mask(irow,icol) == 0; continue; end
        Wave = reshape(svdFrameReconstruct(U(irow,icol,:), V), [nFrames 1]);
        
        % --- Compute spectrogram ---
        [S, F, T] = spectrogram(Wave, window_length, overlap, nfft, fs);
        
        % --- Convert to power (dB) ---
        S_dB(irow, icol,:,:) = 20*log10(abs(S) + eps);
    end
end