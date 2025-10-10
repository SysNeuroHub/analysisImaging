function TexImg = getSufraceData2(TexVol_T2, brainMask)

volumeMask = nan(size(brainMask));
volumeMask(brainMask>0)=1;
TexImg = rot90(squeeze(nansum(TexVol_T2.*volumeMask, 2)));

% % cf. suface mask
%  [~,~,~,S] = vol2Surf(T2w_brain_us>0, 50*scaleFactor);
% surfaceMask = nan(size(S));
% surfaceMask(S>0)=1;
% TexImg = fliplr(rot90(squeeze(nansum(TexVol_T2.*surfaceMask, 2))));




