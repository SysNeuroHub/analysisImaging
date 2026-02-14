function imgs = unpack_polyscan_1bit(frames)

    WIDTH  = 800;
    HEIGHT = 500;

    imgs = cell(size(frames));

    for k = 1:numel(frames)

        buf = frames{k};

        img = false(HEIGHT, WIDTH);

        idx = 1;

        % Column-major, 8-row blocks
        for x = 1:WIDTH
            for yb = 1:(HEIGHT/8)

                b = buf(idx);
                idx = idx + 1;

                % LSB-first vertical bits
                bits = bitget(b,1:8);

                y0 = (yb-1)*8;

                for j = 1:8
                    img(y0+j, x) = bits(j);
                end
            end
        end

        % Flip vertically (common for DMD)
        img = flipud(img);

        imgs{k} = img;

        figure;
        imshow(img);
        title(['Recovered frame ' num2str(k)]);
    end
end
