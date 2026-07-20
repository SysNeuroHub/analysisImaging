function Icorr = applyStretchedExp(I,param, frameNumbers)
%Icorr = applyStretchedExp(I,param)

k = param.k;
beta  = param.beta;

[nTrace,nT] = size(I);

if nargin < 3
    frameNumbers = (1:nT)';
else
    if numel(frameNumbers) ~= nT
        error('frameNumbers must match size(I,3)');
    end
end

B = zeros(nTrace,nT);

for ii=1:nT

    B(:,ii)=exp(-k.*(frameNumbers(ii).^beta));

end


B = max(B,0.2);


Icorr = I./B;

if isfield(param,'scaleFac') && ~isempty(param.scaleFac)
    Icorr = Icorr.*repmat(param.scaleFac,[1 size(Icorr,2)]);
end