function [mov, ops, t_subsample] = readSubsampleTensor(ops)
%ops with the following fields
% NavgFramesSVD
% Nframes
% RegFile
% mimg
% verbose

% created from get_svdcomps

[Ly, Lx] = size(ops.mimg);

if ~isfield(ops, 'verbose')
    ops.verbose = false;
end

rawDType = 'uint16';

ntotframes          = ceil(sum(ops.Nframes));
ops.NavgFramesSVD   = min(ops.NavgFramesSVD, ntotframes);
nt0 = ceil(ntotframes / ops.NavgFramesSVD);


ops.NavgFramesSVD = floor(ntotframes/nt0);
nimgbatch = nt0 * floor(1000/nt0);

ix = 0;
t_subsample = [];
mov = zeros(Ly, Lx, ops.NavgFramesSVD, 'single');

if ops.verbose
    fprintf(1, 'loading data\n');
end

try 
    fid = fopen(ops.RegFile, 'r');
    while 1
        if ops.verbose
            fprintf(1, '   frame %d out of %d\n', ix*nt0, ops.Nframes);
        end
        
        data = fread(fid,  Ly*Lx*nimgbatch, ['*' rawDType]);
        if isempty(data)
            break;
        end
        data = single(data);
        data = reshape(data, Ly, Lx, []);      
    
        irange = 1:nt0*floor(size(data,3)/nt0);
        data = data(:,:, irange);

        data = reshape(data, Ly, Lx, nt0, []);
        davg = single(squeeze(mean(data,3)));

        mov(:,:,ix + (1:size(davg,3))) = davg;

        ix = ix + size(davg,3);
        if isempty(t_subsample)
            t_subsample = -nt0/2+nt0*(1:size(data,4));
        else
            t_subsample = [t_subsample t_subsample(end)-nt0/2+nt0*(1:size(data,4))];
        end
    end
catch me
    fclose(fid);
    rethrow(me);
end
fclose(fid);
toc
mov(:, :, (ix+1):end) = [];