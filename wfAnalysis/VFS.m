function signMap = VFS(xMap, yMap)
% 1) get gradient direction
[~,Vdir] = imgradient(imgaussfilt(yMap,1));
[~,Hdir] = imgradient(imgaussfilt(xMap,1));

% 3) get sin(difference in direction) if retinotopic, H/V should be
% orthogonal, so the closer the orthogonal the better (and get sign)
angle_diff = sind(Vdir-Hdir);
angle_diff(isnan(angle_diff)) = 0;

signMap = angle_diff;