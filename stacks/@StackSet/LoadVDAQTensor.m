function S = LoadVDAQTensor( S, animal, iseries, iexp, ScaleFac )
% Loads a VDAQ tensor, for backward compatibility

global VDAQ

if nargin == 0 % a constructor must handle this case
    return
end

if nargin<4
    ScaleFac = 1;
end

%  TensorBuild(animal,iseries,iexp,[],[],1);
% run this the first time you analyze imaging data,
% then the tensor is saved and you can just load it

TensorLoad(animal,iseries,iexp);

% now import the data from the global VDAQ,
% scaling it by the factor ScaleFac

if isempty(VDAQ), error('Could not load VDAQ'); end

S.nRows = round(VDAQ.nx * ScaleFac);
S.nCols = round(VDAQ.ny * ScaleFac);
S.nConds = VDAQ.nstim;

nn = zeros(S.nConds,1);
for iCond = 1:S.nConds
    nn(iCond) = size(VDAQ.tensor{iCond},3);
end

S.nFrames = max(nn);
S.FrameRate     = VDAQ.FrameRate;
S.TimeVec       = (1:S.nFrames)/S.FrameRate;
S.PixelSize     = VDAQ.MmPerCameraPix / ScaleFac;
S.Description   = sprintf('%s-%d-%d',VDAQ.animal,VDAQ.iseries,VDAQ.iexp);

S.Values = nan(S.nRows,S.nCols,S.nFrames,S.nConds);

fprintf('Importing VDAQ tensor');
if ScaleFac ~= 1, fprintf(' and scaling by %2.2f', ScaleFac); end
for iCond = 1:S.nConds
    fprintf('.');
    for iFrame = 1:nn(iCond)
        S.Values(:,:,iFrame,iCond) = ...
            imresize(VDAQ.tensor{iCond}(:,:,iFrame), [S.nRows, S.nCols]);
    end
end
fprintf('done\n');

VDAQ = []; % to save memory

end
