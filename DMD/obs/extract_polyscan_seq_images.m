function frames = extract_polyscan_seq_images(seqFile)
% Extract embedded 1-bit images from a PolyScan .seq file
% Assumes:
%   Resolution: 800 x 500
%   Bit depth: 1-bit
%   Number of frames: 4

    WIDTH  = 800;
    HEIGHT = 500;
    BYTES_PER_FRAME = WIDTH * HEIGHT / 8;  % 50,000

    fid = fopen(seqFile, 'rb');
    raw = fread(fid, inf, 'uint8');
    fclose(fid);

    fprintf('SEQ size: %.2f MB\n', numel(raw)/1e6);

    candidates = [];

    % sliding window search
    for i = 1 : 512 : (numel(raw) - BYTES_PER_FRAME)
        block = raw(i : i + BYTES_PER_FRAME - 1);

        % Heuristic: 1-bit images have low entropy
        uniqueVals = numel(unique(block));

        if uniqueVals < 16   % very conservative
            candidates(end+1) = i; %#ok<AGROW>
        end
    end

    % Deduplicate nearby hits
    candidates = unique(round(candidates / BYTES_PER_FRAME) * BYTES_PER_FRAME);

    fprintf('Found %d candidate blocks\n', numel(candidates));

    % Keep first 4 plausible frames
    candidates = candidates(1:4)+1;

    % Extract frames
    frames = cell(1,4);
    for k = 1:4
        frames{k} = raw(candidates(k) : candidates(k) + BYTES_PER_FRAME - 1);
    end
end
