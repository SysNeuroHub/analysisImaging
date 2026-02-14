function frames = extract_seq_blocks(seqFile, nFrames)

    W = 800;
    H = 500;
    BYTES = W*H/8;   % 50000

    fid = fopen(seqFile,'rb');
    raw = fread(fid,inf,'uint8');
    fclose(fid);

    entropy = zeros(1, numel(raw)-BYTES);

    for i = 1:512:(numel(raw)-BYTES)
        entropy(i) = numel(unique(raw(i:i+BYTES-1)));
    end

    [~,idx] = sort(entropy);

    starts = sort(idx(1:nFrames));

    frames = cell(1,nFrames);

    for k = 1:nFrames
        s = starts(k);
        frames{k} = raw(s:s+BYTES-1);
    end
end
