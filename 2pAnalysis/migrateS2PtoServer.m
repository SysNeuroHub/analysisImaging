function migrateS2PtoServer(expt, TSDir, S2PDir, marketDir, vaultDir, deleteOriginals)
% migrateS2PtoServer(expt, TSDir, S2PDir, marketDir, vaultDir)
% copies raw data and S2P data to the server

% migrateS2PtoServer(expt, TSDir, S2PDir, marketDir, vaultDir, 1)
% moves raw data and S2P data

% 15/7/20

if nargin < 6
    deleteOriginals = 0;
end

theseExps = expt.expNum;

for eee = 1:length(theseExps)
    S2PDataDir = sprintf('%s/%s/%s/%d/', S2PDir, expt.subject, expt.expDate, theseExps(eee));
    TSDataDir = sprintf('%s/%s/%s/%d/', TSDir, expt.subject, expt.expDate, theseExps(eee));
    marketDataDir = sprintf('%s/%s/%s/%d/', marketDir, expt.subject, expt.expDate, theseExps(eee));
    vaultDataDir = sprintf('%s/%s/%s/%d/', vaultDir, expt.subject, expt.expDate, theseExps(eee));
    
    tsname = fullfile(TSDataDir, 'Episode001.h5'); %fullfile(TSDir, subject, expDate, num2str(expr), 'Episode001.h5');
    tsname2 = fullfile(TSDataDir, 'ThorRealTimeDataSettings.xml');
    tsname3 = fullfile(TSDataDir, 'Experiment.xml');
    
    if deleteOriginals
        %S2P data to Market
        movefile(S2PDataDir, marketDataDir);
    else
        copyfile(S2PDataDir, marketDataDir);
    end
    
    %raw TS data to Market (because vault is slow to read)
    copyfile(tsname, marketDataDir);
    copyfile(tsname2, marketDataDir);
    copyfile(tsname3, marketDataDir);
    
    if deleteOriginals
        %raw TI data to Vault
        movefile(TSDataDir, vaultDataDir);
    else
        copyfile(TSDataDir, vaultDataDir);
    end
end