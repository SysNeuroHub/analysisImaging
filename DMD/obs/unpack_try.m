function img = unpack_try(buf, mode)

    W = 800;
    H = 500;

    img = false(H,W);
    idx = 1;

    switch mode

        % Row-wise, MSB
        case 1
            for y=1:H
                for x=1:W/8
                    b = buf(idx); idx=idx+1;
                    bits = bitget(b,8:-1:1);
                    img(y,(x-1)*8+(1:8)) = bits;
                end
            end

        % Row-wise, LSB
        case 2
            for y=1:H
                for x=1:W/8
                    b = buf(idx); idx=idx+1;
                    bits = bitget(b,1:8);
                    img(y,(x-1)*8+(1:8)) = bits;
                end
            end

        % Column-wise, MSB
        case 3
            for x=1:W
                for y=1:H/8
                    b = buf(idx); idx=idx+1;
                    bits = bitget(b,8:-1:1);
                    img((y-1)*8+(1:8),x) = bits;
                end
            end

        % Column-wise, LSB
        case 4
            for x=1:W
                for y=1:H/8
                    b = buf(idx); idx=idx+1;
                    bits = bitget(b,1:8);
                    img((y-1)*8+(1:8),x) = bits;
                end
            end
    end
end
