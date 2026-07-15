function nE_UTE = searchUTEdirectory(rootDir)
%returns directory names whose protocol include T2star
[seriesDescription, protocolName, dirName, dirSizeBytes] = getMRprotocol(rootDir);

nE_UTE_idx = find(contains(seriesDescription, 'UTE3D'));
nE_UTE = dirName(nE_UTE_idx);

if numel(nE_UTE)> 1
    [~, biggerDirectoryIdx] = max(dirSizeBytes(nE_UTE_idx));
    nE_UTE = nE_UTE(biggerDirectoryIdx);
end

nE_UTE = str2num(nE_UTE{1});