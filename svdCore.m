function [U, Sv, V] = svdCore(mov, ops)
%[U, Sv, V] = svdCore(mov, ops)
%
%ops.nSVD
%ops.useGPU
%ops.roi [0/1]
%mov: row x column x time (single)

%created from get_svdcomps
% number of SVD components kept
%ops.nSVD = min(ops.nSVD, size(mov,3));
%
ops.nSVD = min(ops.nSVD, size(mov,3));
mov             = reshape(mov, [], size(mov,3));

% If an ROI for the brain was selected, zero all outside pixels
% (AP 160804)
if isfield(ops,'roi') && ~isempty(ops.roi)
    mov(~ops.roi(:),:) = 0;
end

% mov             = mov./repmat(mean(mov.^2,2).^.5, 1, size(mov,2));
COV             = mov' * mov/size(mov,1);

% total variance of data. If you ask for all Svs back then you will see
% this is equal to sum(Sv). In this case Sv are the singular values *of the
% covariance matrix* not of the original data - they are equal to the Sv of
% the original data squared (the variances per dimension). 
totalVar = sum(diag(COV)); 
                            

ops.nSVD = min(size(COV,1)-2, ops.nSVD);
%toc
if ops.nSVD<1000 || size(COV,1)>1e4
    [V, Sv]          = eigs(double(COV), ops.nSVD);
else
    if ops.useGPU
        [V, Sv]         = svd(gpuArray(double(COV)));
        V = gather(V);
        Sv = gather(Sv);
    else
         [V, Sv]         = svd(COV);
    end
    V               = V(:, 1:ops.nSVD);
    Sv              = Sv(1:ops.nSVD, 1:ops.nSVD);
end

clear COV
U               = normc(mov * V);
clear mov
U               = single(U);
Sv              = single(diag(Sv));
