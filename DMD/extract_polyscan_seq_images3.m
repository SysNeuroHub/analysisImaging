function frames = extract_polyscan_seq_images3(seqFile)
% Extract raw 1-bit image buffers from PolyScan .seq
% Assumes:
%   Resolution: 800 x 500
%   Bit depth: 1-bit
%   Number of frames: 3

    WIDTH  = 800;
    HEIGHT = 500;
    BYTES_PER_FRAME = WIDTH * HEIGHT / 8;  % 50,000

    fid = fopen(seqFile, 'rb');
    assert(fid > 0, 'Cannot open seq file');
    raw = fread(fid, inf, 'uint8');
    fclose(fid);

    N = numel(raw);
    fprintf('SEQ size: %.2f MB\n', N/1e6);

    candidates = [];

    % Slide window with safe bounds (MATLAB 1-based indexing)
    step = 256;
    for i = 1 : step : (N - BYTES_PER_FRAME + 1)
        block = raw(i : i + BYTES_PER_FRAME - 1);

        % Heuristic: 1-bit images have very few unique byte values
        if numel(unique(block)) <= 8
            candidates(end+1) = i; %#ok<AGROW>
        end
    end

    % Cluster nearby hits (same frame detected multiple times)
    candidates = unique(round(candidates / BYTES_PER_FRAME) * BYTES_PER_FRAME);
    candidates(candidates < 1) = [];

    fprintf('Candidate frame offsets found: %d\n', numel(candidates));

    % We expect exactly 3 frames
   % assert(numel(candidates) >= 3, 'Could not find enough frames');

    % candidates = candidates(1:3);

    % Extract frames
    frames = cell(1,numel(candidates));
    for k = 1:numel(candidates)-1
        startIdx = candidates(k);
        frames{k} = raw(startIdx : startIdx + BYTES_PER_FRAME - 1);
    end
end
