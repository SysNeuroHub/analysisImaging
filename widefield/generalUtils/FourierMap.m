function [ComplexMap, AbsMap, AngleMap] = FourierMap(Tensor, TimeVec, ff)
% Get a frequency component of a tensor
%
% [ComplexMaps, AbsMaps, AngleMaps] = FourierMaps(Tensor, TimeVec, ff) lets
% you specify the frequency. ff (in Hz) can be a number
%
% NOTE: The outputs are arrays nRows x nCols
% code, but it was a pain to have them as cell arrays)

% 2026-02-13 created from Stackset.FourierMaps

[nRows, nCols, nFrames] = size(Tensor);

yy = zeros(1,1,nFrames);
%aaa = zeros(nRows,nCols,nFrames);
%ComplexMap = zeros(nRows,nCols,S.nConds);

yy(:) = 2*exp(- (1:nFrames)/nFrames * ...
    (TimeVec(end) - TimeVec(1)) *2*pi*1i*ff);
aaa = repmat(yy,[nRows,nCols,1]);
ComplexMap = double(mean(Tensor.*aaa, 3));

AbsMap  = abs(  ComplexMap );
AngleMap= angle(ComplexMap );

end