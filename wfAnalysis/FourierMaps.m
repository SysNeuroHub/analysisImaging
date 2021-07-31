function [ComplexMaps, AbsMaps, AngleMaps] = FourierMaps(Tensor, srate, ff)
% Get a frequency component of a tensor
%
% [ComplexMaps, AbsMaps, AngleMaps] = FourierMaps(S, ff) lets
% you specify the frequency. ff (in Hz) can be a number or a vector
% with one value per condition.
%
% NOTE: The outputs are arrays nRows x nCols x nConds
% (MC made this happen 4 Sept 2014, sorry if it breaks older
% code, but it was a pain to have them as cell arrays)

% 12/10/20 DS created from stackset.FourierMaps

% TO DO: check unit of AbsMaps
% search for faster way for SVD format


nRows = size(Tensor,1);
nCols = size(Tensor,2);
nFrames = size(Tensor, 3);
duration = 1/srate * nFrames; 

yy = zeros(1,1,nFrames);
%aaa = zeros(nRows,nCols,nFrames);
%ComplexMaps = zeros(nRows,nCols);


yy(:) = 2*exp(- (1:nFrames)/nFrames * ...
    duration *2*pi*1i*ff);
aaa = repmat(yy,[nRows,nCols,1]);
ComplexMaps = double(mean(Tensor.*aaa, 3));

AbsMaps  = abs(  ComplexMaps );
AngleMaps= angle(ComplexMaps );

end