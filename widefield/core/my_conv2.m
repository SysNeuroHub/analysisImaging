function S1 = my_conv2(S1, sig, varargin)
% S1 = my_conv2(S1, sig)
% applies gaussian filter(s) on each dimension sequentially
% S1 = my_conv2(S1, sig, idims)
% takes an extra argument (idims) which specifies which dimension to filter on
% extra argument can be a vector with all dimensions that need to be
% smoothed, in which case sig can also be a vector of different smoothing
% constants

%called in smooth_movie

if sig>.25
    idims = 2;
    if ~isempty(varargin)
        idims = varargin{1};
    end
    if numel(idims)>1 && numel(sig)>1
        sigall = sig;
    else
        sigall = repmat(sig, numel(idims), 1);
    end
    
    for i = 1:length(idims)
        sig = sigall(i);
        
        idim = idims(i);
        Nd = ndims(S1);
        
        S1 = permute(S1, [idim 1:idim-1 idim+1:Nd]);

        S1 = my_conv(S1, sig);
%         dsnew = size(S1);
%         
%         S1 = reshape(S1, size(S1,1), []);
%         dsnew2 = size(S1);
%                 
%         tmax = ceil(4*sig);
%         dt = -tmax:1:tmax;
%         gaus = exp( - dt.^2/(2*sig^2));
%         gaus = gaus'/sum(gaus);
%                 
%         cNorm = filter(gaus, 1, cat(1, ones(dsnew2(1), 1), zeros(tmax,1)));
%         cNorm = cNorm(1+tmax:end, :);
%         
%         S1 = filter(gaus, 1, cat(1, S1, zeros([tmax, dsnew2(2)])));
%         S1(1:tmax, :) = [];
%         S1 = reshape(S1, dsnew);
%         
%         S1 = bsxfun(@rdivide, S1, cNorm);
        
        S1 = permute(S1, [2:idim 1 idim+1:Nd]);
    end
end