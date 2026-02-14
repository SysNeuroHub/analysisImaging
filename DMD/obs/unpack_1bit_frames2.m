function imgs = unpack_1bit_frames2(frames)
% Convert raw 1-bit buffers into logical images (800x500)

    WIDTH  = 800;
    HEIGHT = 500;

    imgs = cell(size(frames));

    for k = 1:numel(frames)
        % Convert bytes -> bits (MSB first)
        bits = reshape(de2bi(frames{k}, 8, 'left-msb')', [], 1);

        % Keep only needed bits
        bits = bits(1 : WIDTH * HEIGHT);

        % Reshape: [row, col] = [y, x]
        img = reshape(bits, WIDTH, HEIGHT)';
        imgs{k} = logical(img);

        % Display
        figure;
        imshow(imgs{k});
        title(sprintf('Recovered frame %d', k));
    end
end
