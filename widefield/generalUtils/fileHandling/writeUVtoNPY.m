function writeUVtoNPY(U, V, uFilePath, vFilePath)
% writeUVtoNPY(U, V, uFilePath, vFilePath)
% Writes the U and V from the SVD analysis to an npy file. 
% Inputs
%   U: yPix x xPix x nSV
%   V: nSV x nFrames
% For U this is uncomplicated but for V you want to transpose first so that
% each component's time course will be written in order. This facilitates
% loading just a limited number of components later. 

if ~isempty(uFilePath)
    writeNPY(U, [uFilePath '.npy']);
end
if ~isempty(vFilePath)
    writeNPY(V', [vFilePath '.npy']);
end