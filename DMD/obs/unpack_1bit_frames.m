function imgs = unpack_1bit_frames(frames)
% Convert raw 1-bit buffers into images

    WIDTH  = 800;
    HEIGHT = 500;

    imgs = cell(size(frames));

    for k = 1:numel(frames)
        bits = reshape(dec2bin(frames{k}, 8)' - '0', [], 1);
        img  = reshape(bits(1:WIDTH*HEIGHT), WIDTH, HEIGHT)';
        imgs{k} = logical(img);

        figure;
        imshow(imgs{k});
        title(sprintf('Recovered frame %d', k));
    end
end
