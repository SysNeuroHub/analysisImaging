function C = Load_ATLAS_info(ROI)

%ROI: Integer; 18 24 45 51 65 72 79 86 93 100 122 136 143 150 164 171 178 185 192 199 206 213 226 238 298 325 332 346 353 360;
%2018 2024 2045 2051 2065 2072 2079 2086 2093 2100 2122 2136 2143 2150 2164 2171 2178 2185 2192 2199 2206 2213 2298 2346 2353 2360;

C{1,1} = ROI;

switch ROI
    
        % MOTOR
        case 18
        C{1,2}="Mop"; C{1,3}="Right Primary motor area";
        case 24
        C{1,2}="Mos"; C{1,3}="Right Secondary motor area";
        
        % SOMATOSENSORY
        case 44
        C{1,2}="SSP-n"; C{1,3}="Right Primary somatosensory area, nose";
        case 51
        C{1,2}="SSP-bfd"; C{1,3}="Right Primary somatosensory area, barrel field";
        case 65
        C{1,2}="SSP-ll"; C{1,3}="Right Primary somatosensory area, lower limb";
        case 72
        C{1,2}="SSP-m"; C{1,3}="Right Primary somatosensory area, mouth";
        case 79
        C{1,2}="SSP-ul"; C{1,3}="Right Primary somatosensory area, upper limb";
        case 86
        C{1,2}="SSP-tr"; C{1,3}="Right Primary somatosensory area, trunk";
        case 93
        C{1,2}="SSP-un"; C{1,3}="Right Primary somatosensory area, unassigned";
        case 100
        C{1,2}="SSs"; C{1,3}="Right Supplemental somatosensory area";
        
        % AUDITORY
        case 122
        C{1,2}="AUDd"; C{1,3}="Right Dorsal auditory area";
        case 136
        C{1,2}="AUDp"; C{1,3}="Right Primary auditory area";
        case 143
        C{1,2}="AUDpo"; C{1,3}="Right Posterior auditory area";
        case 150
        C{1,2}="AUDv"; C{1,3}="Right Ventral auditory area";
        
        % VISUAL
        case 164
        C{1,2}="VISal"; C{1,3}="Right Anterolateral visual area";
        case 171
        C{1,2}="VISam"; C{1,3}="Right Anteromedial visual area";
        case 178
        C{1,2}="VISl"; C{1,3}="Right Lateral visual area";
        case 185
        C{1,2}="VISp"; C{1,3}="Right Primary visual area";
        case 192
        C{1,2}="VISpl"; C{1,3}="Right Posterolateral visual area";
        case 199
        C{1,2}="VISpm"; C{1,3}="Right Posteromedial visual area";
        
        case 206
        C{1,2}="VISli"; C{1,3}="Right Laterointermediate area";
        case 213
        C{1,2}="VISpor"; C{1,3}="Right Postrhinal area";
        case 226
        C{1,2}="ACAd"; C{1,3}="Right Anterior cingulate area, dorsal part";
        case 238
        C{1,2}="PL"; C{1,3}="Right Prelimbic area";
        case 298
        C{1,2}="RSPagl"; C{1,3}="Right Retrosplenial area, lateral agranular part";
        case 325
        C{1,2}="RSPd"; C{1,3}="Right Retrosplenial area, dorsal part";
        case 332
        C{1,2}="RSPv"; C{1,3}="Right Retrosplenial area, ventral part";
        case 346
        C{1,2}="VISa"; C{1,3}="Right Anterior area";
        case 353
        C{1,2}="VISrl"; C{1,3}="Right Rostrolateral visual area";
        case 360
        C{1,2}="TEa"; C{1,3}="Right Temporal association areas";
        
         % MOTOR
        case 2018
        C{1,2}="Mop"; C{1,3}="Left Primary motor area";
        case 2024
        C{1,2}="Mos"; C{1,3}="Left Secondary motor area";
        
        % SOMATOSENSORY
        case 2044
        C{1,2}="SSP-n"; C{1,3}="Left Primary somatosensory area, nose";
        case 2051
        C{1,2}="SSP-bfd"; C{1,3}="Left Primary somatosensory area, barrel field";
        case 2065
        C{1,2}="SSP-ll"; C{1,3}="Left Primary somatosensory area, lower limb";
        case 2072
        C{1,2}="SSP-m"; C{1,3}="Left Primary somatosensory area, mouth";
        case 2079
        C{1,2}="SSP-ul"; C{1,3}="Left Primary somatosensory area, upper limb";
        case 2086
        C{1,2}="SSP-tr"; C{1,3}="Left Primary somatosensory area, trunk";
        case 2093
        C{1,2}="SSP-un"; C{1,3}="Left Primary somatosensory area, unassigned";
        case 2100
        C{1,2}="SSs"; C{1,3}="Left Supplemental somatosensory area";
        
        % AUDITORY
        case 2122
        C{1,2}="AUDd"; C{1,3}="Left Dorsal auditory area";
        case 2136
        C{1,2}="AUDp"; C{1,3}="Left Primary auditory area";
        case 2143
        C{1,2}="AUDpo"; C{1,3}="Left Posterior auditory area";
        case 2150
        C{1,2}="AUDv"; C{1,3}="Left Ventral auditory area";
        
        % VISUAL
        case 2164
        C{1,2}="VISal"; C{1,3}="Left Anterolateral visual area";
        case 2171
        C{1,2}="VISam"; C{1,3}="Left Anteromedial visual area";
        case 2178
        C{1,2}="VISl"; C{1,3}="Left Lateral visual area";
        case 2185
        C{1,2}="VISp"; C{1,3}="Left Primary visual area";
        case 2192
        C{1,2}="VISpl"; C{1,3}="Left Posterolateral visual area";
        case 2199
        C{1,2}="VISpm"; C{1,3}="Left Posteromedial visual area";
        
        case 2206
        C{1,2}="VISli"; C{1,3}="Left Laterointermediate area";
        case 2213
        C{1,2}="VISpor"; C{1,3}="Left Postrhinal area";
        case 2226
        C{1,2}="ACAd"; C{1,3}="Left Anterior cingulate area, dorsal part";
        case 2238
        C{1,2}="PL"; C{1,3}="Left Prelimbic area";
        case 2298
        C{1,2}="RSPagl"; C{1,3}="Left Retrosplenial area, lateral agranular part";
        case 2325
        C{1,2}="RSPd"; C{1,3}="Left Retrosplenial area, dorsal part";
        case 2332
        C{1,2}="RSPv"; C{1,3}="Left Retrosplenial area, ventral part";
        case 2346
        C{1,2}="VISa"; C{1,3}="Left Anterior area";
        case 2353
        C{1,2}="VISrl"; C{1,3}="Left Rostrolateral visual area";
        case 2360
        C{1,2}="TEa"; C{1,3}="Left Temporal association areas";
end
 

