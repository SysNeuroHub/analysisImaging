function nE_T2 = searchT2directory(rootDir)
%returns directory names whose protocol include T2star
[seriesDescription, protocolName, dirName, dirSizeBytes] = getMRprotocol(rootDir);

nE_T2_idx = find(contains(protocolName, 'T2star'));
nE_T2 = dirName(nE_T2_idx);

if numel(nE_T2)> 1
    [~, biggerDirectoryIdx] = max(dirSizeBytes(nE_T2_idx));
    nE_T2 = nE_T2(biggerDirectoryIdx);
end
nE_T2 = str2num(nE_T2{1});