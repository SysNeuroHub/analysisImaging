function frames = extract_polyscan_seq_images2(seqFile)

    WIDTH  = 800;
    HEIGHT = 500;
    BYTES_PER_FRAME = WIDTH * HEIGHT / 8;

    fid = fopen(seqFile,'rb');
    raw = fread(fid,inf,'uint8');
    fclose(fid);

    candidates = [];

    for i = 1:512:(numel(raw)-BYTES_PER_FRAME)
        block = raw(i:i+BYTES_PER_FRAME-1);
        if numel(unique(block)) < 32
            candidates(end+1) = i; %#ok<AGROW>
        end
    end

    candidates = unique(round(candidates/BYTES_PER_FRAME)*BYTES_PER_FRAME)+1;

    frames = cell(1,min(4,numel(candidates)));

    for k = 1:numel(frames)
        frames{k} = raw(candidates(k):candidates(k)+BYTES_PER_FRAME-1);
    end
end
